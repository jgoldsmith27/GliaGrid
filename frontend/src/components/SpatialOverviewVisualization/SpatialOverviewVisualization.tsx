import React, { useState, useMemo, useCallback, useEffect, useRef } from 'react';
import { ScatterplotLayer, PolygonLayer, GeoJsonLayer, PathLayer } from '@deck.gl/layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters, Color, PickingInfo, MapViewState, Layer } from '@deck.gl/core';
import { scaleOrdinal } from 'd3-scale';
import { schemeCategory10 } from 'd3-scale-chromatic'; // Example color scheme
import IconButton from '@mui/material/IconButton'; // Keep this if other IconButtons are used, otherwise remove
import Icon from '@mui/material/Icon'; // Use base Icon component
import styles from './SpatialOverviewVisualization.module.css';
import { SharedDataStore, useSharedData } from '../../services/data/SharedDataStore'; // Import SharedDataStore
import LoadingSpinner from '../LoadingSpinner/LoadingSpinner'; // Import LoadingSpinner
import useJobStatus from '../../hooks/useJobStatus'; // Correct: Default import
// Import turf/boolean-point-in-polygon or similar for intersection checks
import booleanPointInPolygon from '@turf/boolean-point-in-polygon'; // Correct: Default import
import { polygon as turfPolygon, point as turfPoint, feature, polygon, lineString as turfLineString } from '@turf/helpers'; // Correct: Named imports + types
import { Feature as GeoJsonFeature, Point as GeoJsonPoint, GeoJsonProperties, Polygon } from 'geojson'; // Import GeoJSON types
import ScaleBar from '../VisualizationHelpers/ScaleBar'; // Import the new ScaleBar component
import Tooltip from '@mui/material/Tooltip'; // Import Tooltip

// Define types locally or import from a central types file
// Removed SpatialPoint definition as PointFeature is imported

// Type for layer boundary data (array of [x, y] coordinates)
type LayerBoundary = [number, number][];

// Helper to get initial view state based on boundaries or points
const getInitialViewState = (validBoundaries: Record<string, LayerBoundary>) => {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    let hasData = false;

    if (validBoundaries && Object.keys(validBoundaries).length > 0) {
        Object.values(validBoundaries).forEach(boundary => {
            if (boundary) { // Check if boundary is not null/undefined
                boundary.forEach(([x, y]) => {
                    minX = Math.min(minX, x);
                    minY = Math.min(minY, y);
                    maxX = Math.max(maxX, x);
                    maxY = Math.max(maxY, y);
                    hasData = true;
                });
            }
        });
    }

    if (!hasData) {
        return { target: [0, 0, 0] as [number, number, number], zoom: 1, minZoom: -10, maxZoom: 10 };
    }

    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;
    // Add padding to width/height calculation
    const padding = 50; // Adjust as needed
    const width = Math.max(1, maxX - minX) + padding * 2;
    const height = Math.max(1, maxY - minY) + padding * 2;
    
    const safeWindowWidth = window.innerWidth > 0 ? window.innerWidth : 1000; // Default if window size is 0
    const safeWindowHeight = window.innerHeight > 0 ? window.innerHeight : 800;
    
    const zoomX = Math.log2(safeWindowWidth / width);
    const zoomY = Math.log2(safeWindowHeight / height);
    let zoom = Math.min(zoomX, zoomY) - 0.5; // Adjust zoom level
    zoom = isNaN(zoom) ? 1 : zoom; // Handle potential NaN

    return {
        target: [centerX, centerY, 0] as [number, number, number],
        zoom: Math.max(-10, Math.min(10, zoom)), // Clamp initial zoom too
        minZoom: -10,
        // Set maxZoom limit. Rationale:
        // - Coordinate units are typically ~2 units/µm.
        // - Orthographic zoom means pixels = worldUnits * 2^zoom.
        // - At zoom=6, 1 world unit (0.5 µm) = 2^6 = 64 pixels.
        // - Zooming further overly magnifies the discrete coordinate space 
        //   without revealing more biological detail (points become huge blocks).
        // - This value provides a balance between close inspection and meaningful resolution.
        maxZoom: 6 
    };
};

// Define component props
interface SpatialOverviewVisualizationProps {
  jobId: string;
  onLassoSelect?: (selectedCoords: [number, number][] | null) => void; // MODIFIED: Expect coords or null
  onAnalyzeSelection?: () => void; // ADDED: Callback to trigger analysis
}

const SpatialOverviewVisualization: React.FC<SpatialOverviewVisualizationProps> = ({ 
  jobId, 
  onLassoSelect, 
  onAnalyzeSelection // Destructure the prop
}) => {
  // ADDED: Log received prop
  console.log("[SpatialOverviewVisualization] Received onAnalyzeSelection:", typeof onAnalyzeSelection);

  // Use useJobStatus to get the final analysis results
  const { jobStatus, isLoading: isJobStatusLoading, error: jobStatusError } = useJobStatus(jobId);

  // --- State --- 
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [visibleLayers, setVisibleLayers] = useState<Set<string>>(new Set());
  const [layerColors, setLayerColors] = useState<Record<string, Color>>({});
  const [uniqueLayerNames, setUniqueLayerNames] = useState<string[]>([]);
  const [rawLayerBoundaries, setRawLayerBoundaries] = useState<Record<string, LayerBoundary | null>>({});
  const [isDrawingLasso, setIsDrawingLasso] = useState(false);
  const [lassoPoints, setLassoPoints] = useState<[number, number][]>([]);

  // --- Ruler State ---
  const [isRulerActive, setIsRulerActive] = useState(false);
  const [rulerPoints, setRulerPoints] = useState<[number, number][]>([]);
  const [rulerDistanceMicrons, setRulerDistanceMicrons] = useState<number | null>(null);
  const [rulerHoverCoord, setRulerHoverCoord] = useState<[number, number] | null>(null); // State for hover coordinate

  // --- Refs --- 
  const isDragging = useRef(false);
  const deckContainerRef = useRef<HTMLDivElement>(null);

  // Process JOB RESULTS to find unique layers, assign colors, and get boundaries
  useEffect(() => {
    const results = jobStatus?.results;
    const outputs = results?.outputs;
    const boundaries = outputs?.layer_boundaries as Record<string, LayerBoundary | null> | undefined; // Type assertion with undefined possibility

    if (outputs && boundaries) {
        setRawLayerBoundaries(boundaries);
        
        // Extract unique layer names *from boundaries* where boundary is not null
        const uniqueLayers = Object.keys(boundaries).filter(layer => boundaries[layer] !== null);
        setUniqueLayerNames(uniqueLayers);
        
        // Assign colors using d3 scale and convert to RGBA array inline
        const colorScale = scaleOrdinal<string, string>(schemeCategory10);
        const colors: Record<string, Color> = {};
        uniqueLayers.forEach((layer) => {
            const hexColor = colorScale(layer);
            // Simple hex to RGBA conversion (assuming #RRGGBB)
            const r = parseInt(hexColor.slice(1, 3), 16);
            const g = parseInt(hexColor.slice(3, 5), 16);
            const b = parseInt(hexColor.slice(5, 7), 16);
            colors[layer] = [r, g, b, 180]; // RGBA color with some transparency
        });
        setLayerColors(colors);
        
        // Initialize visibility - make all layers with non-null boundaries visible initially
        setVisibleLayers(new Set(uniqueLayers));
        console.log("[SpatialOverviewVisualization] Layers found (from boundaries):", uniqueLayers);
    } else {
        // Reset if no boundaries found in results
        setUniqueLayerNames([]);
        setLayerColors({});
        setVisibleLayers(new Set());
        setRawLayerBoundaries({});
    }
  }, [jobStatus?.results]); // Depend specifically on the results part of jobStatus

  // Calculate initial view state based ONLY on VALID layer boundaries
  const initialViewState = useMemo(() => {
      console.log(`[SpatialOverviewVisualization] Calculating initial view state from boundaries`);
      // Filter out null boundaries before passing to calculation
      const validBoundaries: Record<string, LayerBoundary> = {};
      Object.entries(rawLayerBoundaries).forEach(([layer, boundary]) => {
          if (boundary) {
              validBoundaries[layer] = boundary;
          }
      });
      return getInitialViewState(validBoundaries);
  }, [rawLayerBoundaries]); // Depend on the raw boundaries state

  // Initialize state with computed initial values and bounds
  const [viewState, setViewState] = useState<OrthographicViewState>({
      ...initialViewState, // Spread calculated target and initial zoom
      minZoom: initialViewState.minZoom ?? -10,
      maxZoom: initialViewState.maxZoom ?? 6,
  });

  // Reset view state when initial view state calculation changes (e.g., data loads)
  useEffect(() => {
       // Check if viewState and target exist before comparing
       if (viewState && viewState.target && 
           (initialViewState.target[0] !== viewState.target[0] || 
            initialViewState.target[1] !== viewState.target[1] || 
            initialViewState.zoom !== viewState.zoom)) {
           setViewState(initialViewState);
       }
  }, [initialViewState, viewState]); // Add viewState dependency

  const onViewStateChange = useCallback(({ viewState: vs }: ViewStateChangeParameters<OrthographicViewState>) => {
    if (!isDrawingLasso || !isDragging.current) {
        let newZoom = vs.zoom;
        // Clamp zoom only if it's a number, respecting the defined min/max bounds.
        if (typeof newZoom === 'number') {
            newZoom = Math.max(initialViewState.minZoom ?? -10, Math.min(initialViewState.maxZoom ?? 6, newZoom));
        }
        // Update state, ensuring bounds are preserved
        setViewState({...vs, zoom: newZoom, minZoom: initialViewState.minZoom ?? -10, maxZoom: initialViewState.maxZoom ?? 6 });
    }
  }, [isDrawingLasso, initialViewState.minZoom, initialViewState.maxZoom]);

  const toggleLayerVisibility = useCallback((layerName: string) => {
    setVisibleLayers(prev => {
        const next = new Set(prev);
        if (next.has(layerName)) {
            next.delete(layerName);
        } else {
            next.add(layerName);
        }
        return next;
    });
  }, []);

  // --- Lasso Event Handlers --- 
  const handleDragStart = useCallback((info: PickingInfo, event: any) => {
      if (!isDrawingLasso) return;
      if (event.leftButton) {
          isDragging.current = true;
          setLassoPoints([info.coordinate as [number, number]]); 
      }
  }, [isDrawingLasso]);

  const handleDrag = useCallback((info: PickingInfo, event: any) => {
      if (!isDrawingLasso || !isDragging.current) return;
      if (event.leftButton) {
          setLassoPoints(prev => {
              const next = [...prev, info.coordinate as [number, number]];
              console.log('[handleDrag] Adding point, new count:', next.length, 'Coords:', info.coordinate); // LOGGING
              return next;
          });
      }
  }, [isDrawingLasso]);

  const handleDragEnd = useCallback((info: PickingInfo, event: any) => {
      if (!isDrawingLasso || !isDragging.current) return;
      
      isDragging.current = false;
      
      const finalPoints = [...lassoPoints, lassoPoints[0]]; // Close the polygon
      setLassoPoints(finalPoints); // MODIFIED: Save the final closed points to state

      if (finalPoints.length >= 4) { 
          // MODIFIED: Call callback with lasso coordinates
          console.log("[SpatialOverview] Lasso completed with points:", finalPoints);
          if (onLassoSelect) {
              onLassoSelect(finalPoints);
          }
      } else {
          // If not enough points, clear selection via callback
          if (onLassoSelect) {
              onLassoSelect(null);
          }
      }
      setIsDrawingLasso(false); // Still deactivate drawing MODE
  }, [isDrawingLasso, lassoPoints, onLassoSelect]); 

  // --- Lasso UI Controls --- 
  const toggleLasso = () => {
      if (isDrawingLasso) {
          clearLasso();
      }
      setIsDrawingLasso(!isDrawingLasso);
  };

  const clearLasso = () => {
      setIsDrawingLasso(false);
      setLassoPoints([]);
      // MODIFIED: Call callback with null to clear selection
      if (onLassoSelect) {
        onLassoSelect(null);
      }
  };

  // --- Ruler Controls ---
  const toggleRuler = useCallback(() => {
    setIsRulerActive(prev => {
        const nextIsActive = !prev;
        // Clear points and distance when toggling ruler off
        if (!nextIsActive) { 
          setRulerPoints([]);
          setRulerDistanceMicrons(null);
        }
        // Ensure lasso is off when ruler is activated
        if (nextIsActive) {
            setIsDrawingLasso(false); 
            setLassoPoints([]); 
            if (onLassoSelect) onLassoSelect(null); // Clear external lasso state too
        }
        return nextIsActive;
    });
  }, [onLassoSelect]);

  const clearRuler = useCallback(() => {
    setRulerPoints([]);
    setRulerDistanceMicrons(null);
    // Keep ruler active after clearing, user can click again or toggle off
  }, []);

  // --- DeckGL Click Handler for Ruler ---
  const handleDeckClick = useCallback((info: PickingInfo, event: any) => {
    if (isRulerActive && info.coordinate) { 
      const clickedCoord: [number, number] = [info.coordinate[0], info.coordinate[1]];
      setRulerHoverCoord(null); // Reset hover coord on click
      setRulerPoints(prev => {
        let nextPoints = [...prev];
        if (nextPoints.length >= 2) {
          nextPoints = [clickedCoord]; 
          setRulerDistanceMicrons(null);
        } else {
          nextPoints.push(clickedCoord);
        }
        if (nextPoints.length === 2) {
          const p1 = nextPoints[0];
          const p2 = nextPoints[1];
          const dx = p2[0] - p1[0];
          const dy = p2[1] - p1[1];
          const distanceUnits = Math.sqrt(dx * dx + dy * dy);
          const distanceMicrons = distanceUnits / 2;
          setRulerDistanceMicrons(distanceMicrons);
          // Maybe deactivate ruler here? setIsRulerActive(false);
        }
        return nextPoints;
      });
      return; 
    }
  }, [isRulerActive]);

  const handleHover = useCallback((info: PickingInfo, event: any) => {
    if (isRulerActive && rulerPoints.length === 1 && info.coordinate) {
      setRulerHoverCoord([info.coordinate[0], info.coordinate[1]]);
    } else {
        // Clear hover coord if ruler not active or 0/2 points selected
        if (rulerHoverCoord !== null) { // Avoid unnecessary state updates
             setRulerHoverCoord(null);
        }
    }
  }, [isRulerActive, rulerPoints.length, rulerHoverCoord]); // Add rulerHoverCoord to dependencies

  // Define DeckGL layers based on BOUNDARIES and visibility state
  const layers = useMemo(() => {
      if (!rawLayerBoundaries || Object.keys(rawLayerBoundaries).length === 0) return [];

      const boundaryLayers = Object.entries(rawLayerBoundaries)
          // Filter for layers that are visible AND have non-null boundaries
          .filter(([layerName, boundary]) => boundary && visibleLayers.has(layerName))
          .map(([layerName, boundary]) => {
              // Type guard to ensure boundary is not null here
              if (!boundary) return null; 
              
              // Define type for polygon data item
              interface PolygonDataItem { polygon: LayerBoundary; layerName: string; }
              
              return new PolygonLayer<PolygonDataItem>({
                  id: `boundary-${layerName}`,
                  data: [{ polygon: boundary, layerName: layerName }], // Data format for PolygonLayer
                  getPolygon: (d: PolygonDataItem) => d.polygon,
                  getFillColor: layerColors[layerName] || [128, 128, 128, 50],
                  getLineColor: [255, 255, 255, 150], // White outline
                  getLineWidth: 1,
                  lineWidthMinPixels: 1,
                  pickable: !isDrawingLasso && !isRulerActive, // Disable picking if ruler is active too
                  autoHighlight: !isDrawingLasso && !isRulerActive, // Disable auto-highlight while drawing
                  highlightColor: [255, 255, 255, 100],
                  updateTriggers: { // Trigger update when selection changes
                    // REMOVED: selectedLayerNames triggers
                    // getFillColor: [selectedLayerNames],
                    // getLineWidth: [selectedLayerNames]
                  }
              });
          })
          .filter(layer => layer !== null); // Filter out any null layers from map
          
       // Add Lasso Layer
       const lassoLayer = lassoPoints.length >= 2 ?
            new GeoJsonLayer({
                id: 'lasso-layer',
                // MODIFIED: Use LineString while dragging OR if < 4 points.
                // Use Polygon only when drag is finished AND there are enough points.
                data: isDragging.current || lassoPoints.length < 4
                    ? turfLineString(lassoPoints, { layerType: 'lasso' })
                    : turfPolygon([[...lassoPoints]], { layerType: 'lasso' }), // Only called when !isDragging and length >= 4
                stroked: true,
                // MODIFIED: Fill only when it's a finished polygon
                filled: !isDragging.current && lassoPoints.length >= 4,
                getLineColor: [255, 0, 0, 200],
                getFillColor: [255, 0, 0, 50],
                getLineWidth: 2,
                lineWidthMinPixels: 2,
                pickable: false
            }) : null;

       const finalLayers = [
           ...boundaryLayers,
           lassoLayer
       ].filter(Boolean); // Filter out null lasso layer

       console.log(`[SpatialOverviewVisualization] Rendering ${boundaryLayers.length} boundary layers. Lasso active: ${isDrawingLasso}, Dragging: ${isDragging.current}`); // Updated log
       // LOGGING added below
       if (lassoLayer) {
           console.log('[Layers Memo] Lasso Layer Data:', lassoLayer.props.data);
       } else {
           console.log('[Layers Memo] No Lasso Layer. Points:', lassoPoints.length);
       }

       // --- Add Ruler Layers (Path and Points) ---
       let rulerLayers: Layer[] = [];
       let rulerPathData: [number, number][] | null = null;

       if (rulerPoints.length === 1 && rulerHoverCoord) {
           // Line follows cursor
           rulerPathData = [rulerPoints[0], rulerHoverCoord];
       } else if (rulerPoints.length === 2) {
           // Fixed line after second click
           rulerPathData = rulerPoints;
       }

       if (rulerPathData) {
           rulerLayers.push(new PathLayer({
               id: 'ruler-path',
               data: [{ path: rulerPathData }],
               getPath: d => d.path,
               getColor: [0, 200, 200, 200], // Cyan with some transparency
               getWidth: 2,
               widthMinPixels: 2,
               billboard: true, 
               pickable: false,
           }));
       }
       
       // Layer for the points (only show clicked points, not hover coord)
       if (rulerPoints.length > 0) {
           rulerLayers.push(new ScatterplotLayer({
               id: 'ruler-points',
               data: rulerPoints.map(p => ({ coordinates: p })), 
               getPosition: d => d.coordinates,
               getFillColor: [0, 200, 200, 255], // Solid Cyan color for points
               getRadius: 5, 
               radiusScale: 1, 
               radiusMinPixels: 3, 
               radiusMaxPixels: 10, 
               pickable: false,
           }));
       }

       const finalLayersWithRuler = [
           ...finalLayers,
           ...rulerLayers
       ].filter(Boolean); // Filter out null layers

       return finalLayersWithRuler as Layer[];

  }, [rawLayerBoundaries, visibleLayers, layerColors, lassoPoints, isDrawingLasso, rulerPoints, isRulerActive, rulerHoverCoord]); // Added rulerHoverCoord dependency

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  const onDeckError = useCallback((error: Error, layer: any) => {
    console.error('[DeckGL Error]', error, 'Layer:', layer?.id);
    // Optionally, you could add state to display a user-friendly error message
  }, []);

  // Tooltip function - update to show layer name when hovering over polygon
  const getTooltipText = useCallback((info: PickingInfo | null | undefined): string | null => { // Correct type: PickingInfo
      if (isDrawingLasso) return null; // No tooltips while drawing
      if (info?.layer?.id.startsWith('boundary-') && info?.object) {
           // Access layerName from the data item we structured for PolygonLayer
           const layerName = (info.object as { layerName: string }).layerName;
           return `Layer: ${layerName}`;
       }
      // Add tooltip for lasso layer if needed (using properties we added)
      // else if (info?.layer?.id === 'lasso-layer' && info?.object?.properties?.layerType === 'lasso') {
      //    return `Lasso Selection Area`;
      // }
      return null;
  }, [isDrawingLasso]);

  // --- Loading and Error States --- 
  // Check if job status is loading AND results haven't arrived yet
  const isLoading = isJobStatusLoading && !jobStatus?.results;
  // Use job status error first, then potentially add checks for missing boundaries later
  const error = jobStatusError;

  if (isLoading) {
    return (
        <div className={`${styles.spatialOverviewContainer} ${styles.loadingOrError}`}>
            <LoadingSpinner message={`Loading analysis results...`} />
        </div>
    );
  }

  if (error) {
    return (
        <div className={`${styles.spatialOverviewContainer} ${styles.loadingOrError}`}>
            <p className={styles.errorMessage}>Error loading analysis results: {error}</p>
        </div>
    );
  }
  
  // Check if results are present but boundaries are missing (potential backend issue)
  if (jobStatus?.status === 'success' && jobStatus?.results && !jobStatus?.results?.outputs?.layer_boundaries) {
      return (
        <div className={`${styles.spatialOverviewContainer} ${styles.loadingOrError}`}>
            <p className={styles.errorMessage}>Analysis complete, but layer boundary data is missing in the results.</p>
        </div>
    );
  }

  // Main Render
  return (
    <div ref={deckContainerRef} className={`${styles.spatialOverviewContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
      <div className={styles.deckGlWrapper}> 
        <DeckGL
          views={new OrthographicView({id: 'ortho-overview-view'})} // Unique view id
          viewState={viewState} // Use the controlled viewState for zoom and target
          onViewStateChange={onViewStateChange}
          controller={true} // Enable default DeckGL controller
          onError={onDeckError}
          layers={layers as any[]}
          getTooltip={getTooltipText} 
          onDragStart={handleDragStart}
          onDrag={handleDrag}
          onDragEnd={handleDragEnd}
          onClick={handleDeckClick} // Use updated handler
          onHover={handleHover} // Add the hover handler
          getCursor={({isDragging: deckIsDragging}) => 
              isRulerActive ? 'crosshair' : 
              (isDrawingLasso ? 'crosshair' : 
              (deckIsDragging ? 'grabbing' : 'grab'))
          }
          style={{ width: '100%', height: '100%', position: 'relative' }}
        />
        {/* Add ScaleBar here, ensure viewState and zoom are valid */}
        {viewState && typeof viewState.zoom === 'number' && (
          <ScaleBar 
            currentZoom={viewState.zoom} 
            unitsPerMicron={2} // 2 coordinate units = 1 micron
            targetPixelWidth={100} // Optional: aim for a 100px wide bar
          />
        )}

        {/* Ruler Distance Display - Improved Styling */} 
        {rulerDistanceMicrons !== null && (
            <div className={styles.measurementDisplay}> {/* Use a dedicated class */} 
                <Icon sx={{ fontSize: 16, mr: 0.5, color: '#00c8c8' }}>straighten</Icon> {/* Optional Icon */} 
                <span>{`Distance: ${rulerDistanceMicrons.toFixed(1)} µm`}</span>
                <Tooltip title="Clear Measurement">
                    <IconButton size="small" onClick={clearRuler} sx={{ ml: 0.5, p: '2px' }}>
                        <Icon fontSize="inherit">close</Icon>
                    </IconButton>
                </Tooltip>
            </div>
        )}
      </div>
      {/* Controls Overlay */} 
      <div className={styles.controlsOverlay}>
          <div className={styles.controlSection}>
              <h4>Layers</h4>
              <div className={styles.layerLegendContainer}>
                {uniqueLayerNames.length > 0 ? (
                  uniqueLayerNames.map((layerName) => (
                    <div key={layerName} className={styles.layerControlItem} onClick={() => toggleLayerVisibility(layerName)} title={`Toggle ${layerName}`}>
                      <input
                        type="checkbox"
                        readOnly // Control checked state via parent div click
                        className={styles.layerToggleCheckbox}
                        checked={visibleLayers.has(layerName)}
                        // Remove onChange handler, handled by div click
                        disabled={isDrawingLasso} 
                      />
                      {/* Ensure color swatch uses the correct class */}
                      <span 
                        className={styles.legendColorSwatch} 
                        style={{ backgroundColor: `rgba(${layerColors[layerName]?.join(',')})` }}
                      ></span>
                      <span className={styles.layerNameText}>{layerName}</span>
                    </div>
                  ))
                ) : (
                  <p>No layers found.</p>
                )}
              </div>
          </div>
          {/* Tools Section */} 
          <div className={styles.toolsSection}> 
                {/* Lasso Controls */} 
                <div className={`${styles.controlSection} ${styles.toolGroup}`}> 
                    <h4>Selection</h4>
                    <div className={styles.toolButtons}> 
                        <Tooltip title={isDrawingLasso ? "Cancel Lasso" : "Start Lasso"}>
                            <span> {/* Span needed for tooltip on disabled button */} 
                            <IconButton onClick={toggleLasso} disabled={isRulerActive} color={isDrawingLasso ? 'secondary' : 'default'}>
                                <Icon>gesture</Icon> {/* Example lasso icon */} 
                            </IconButton>
                            </span>
                        </Tooltip>
                        <Tooltip title="Clear Selection">
                             <span> {/* Span needed for tooltip on disabled button */} 
                             <IconButton onClick={clearLasso} disabled={lassoPoints.length === 0 || isRulerActive} size="small">
                                 <Icon fontSize="inherit">clear_all</Icon> {/* Example clear icon */} 
                             </IconButton>
                             </span>
                        </Tooltip>
                        <Tooltip title="Analyze Selection">
                             <span> {/* Span needed for tooltip on disabled button */} 
                             <IconButton onClick={onAnalyzeSelection} disabled={lassoPoints.length < 4 || isDrawingLasso || isRulerActive} size="small">
                                 <Icon fontSize="inherit">analytics</Icon> {/* Example analyze icon */} 
                             </IconButton>
                             </span>
                        </Tooltip>
                    </div>
                </div>
                {/* Ruler Controls */} 
                <div className={`${styles.controlSection} ${styles.toolGroup}`}> 
                    <h4>Measure</h4>
                    <div className={styles.toolButtons}> 
                        <Tooltip title={isRulerActive ? "Deactivate Ruler" : "Activate Ruler"}>
                             <span> {/* Span needed for tooltip on disabled button */} 
                            <IconButton onClick={toggleRuler} disabled={isDrawingLasso} color={isRulerActive ? "primary" : "default"}>
                                <Icon>straighten</Icon>
                            </IconButton>
                             </span>
                        </Tooltip>
                        <Tooltip title="Clear Measurement">
                             <span> {/* Span needed for tooltip on disabled button */} 
                             <IconButton onClick={clearRuler} disabled={rulerPoints.length === 0} size="small">
                                 <Icon fontSize="inherit">close</Icon>
                             </IconButton>
                             </span>
                        </Tooltip>
                    </div>
                </div>
          </div>
      </div>
      {/* Fullscreen Toggle Button - Use standard button with unicode characters */}
      <button 
          onClick={toggleFullscreen} 
          className={styles.fullscreenButton} // Apply the style class
          title={isFullscreen ? 'Exit Fullscreen' : 'Enter Fullscreen'}
      >
        {isFullscreen ? '✕' : '⛶'}
      </button>
    </div>
  );
};

export default SpatialOverviewVisualization; 