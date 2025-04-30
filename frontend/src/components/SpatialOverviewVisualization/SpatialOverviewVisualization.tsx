import React, { useState, useMemo, useCallback, useEffect, useRef } from 'react';
import { ScatterplotLayer, PolygonLayer, GeoJsonLayer } from '@deck.gl/layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters, Color, PickingInfo, MapViewState } from '@deck.gl/core';
import { scaleOrdinal } from 'd3-scale';
import { schemeCategory10 } from 'd3-scale-chromatic'; // Example color scheme
import IconButton from '@mui/material/IconButton'; // Keep this if other IconButtons are used, otherwise remove
import styles from './SpatialOverviewVisualization.module.css';
import { SharedDataStore, useSharedData } from '../../services/data/SharedDataStore'; // Import SharedDataStore
import LoadingSpinner from '../LoadingSpinner/LoadingSpinner'; // Import LoadingSpinner
import useJobStatus from '../../hooks/useJobStatus'; // Correct: Default import
// Import turf/boolean-point-in-polygon or similar for intersection checks
import booleanPointInPolygon from '@turf/boolean-point-in-polygon'; // Correct: Default import
import { polygon as turfPolygon, point as turfPoint, feature, polygon } from '@turf/helpers'; // Correct: Named imports + types
import { Feature as GeoJsonFeature, Point as GeoJsonPoint, GeoJsonProperties, Polygon } from 'geojson'; // Import GeoJSON types

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
        zoom: Math.max(-10, Math.min(10, zoom)), // Clamp zoom within bounds
        minZoom: -10,
        maxZoom: 10
    };
};

// Define component props
interface SpatialOverviewVisualizationProps {
  jobId: string;
  onLassoSelect?: (selectedLayerNames: string[]) => void; // Optional callback
}

const SpatialOverviewVisualization: React.FC<SpatialOverviewVisualizationProps> = ({ jobId, onLassoSelect }) => {
  // Use useJobStatus to get the final analysis results
  const { jobStatus, isLoading: isJobStatusLoading, error: jobStatusError } = useJobStatus(jobId);

  const [isFullscreen, setIsFullscreen] = useState(false);
  const [visibleLayers, setVisibleLayers] = useState<Set<string>>(new Set());
  const [layerColors, setLayerColors] = useState<Record<string, Color>>({});
  const [uniqueLayerNames, setUniqueLayerNames] = useState<string[]>([]);
  const [rawLayerBoundaries, setRawLayerBoundaries] = useState<Record<string, LayerBoundary | null>>({});

  // --- Lasso State --- 
  const [isDrawingLasso, setIsDrawingLasso] = useState(false);
  const [lassoPoints, setLassoPoints] = useState<[number, number][]>([]);
  const [selectedLayerNames, setSelectedLayerNames] = useState<string[]>([]);
  const isDragging = useRef(false);

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

  const [viewState, setViewState] = useState<OrthographicViewState>(initialViewState);

  // Reset view state when initial view state calculation changes
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
    // Prevent view state change while drawing lasso
    if (!isDrawingLasso || !isDragging.current) {
        setViewState(vs as OrthographicViewState);
    }
  }, [isDrawingLasso]);

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
      // Start drawing only on primary button click
      if (event.leftButton) {
          isDragging.current = true;
          setLassoPoints([info.coordinate as [number, number]]); // Start with the first point
          setSelectedLayerNames([]); // Clear previous selection
      }
  }, [isDrawingLasso]);

  const handleDrag = useCallback((info: PickingInfo, event: any) => {
      if (!isDrawingLasso || !isDragging.current) return;
      // Add point if dragging with primary button down
      if (event.leftButton) {
          setLassoPoints(prev => [...prev, info.coordinate as [number, number]]);
      }
  }, [isDrawingLasso]);

  const handleDragEnd = useCallback((info: PickingInfo, event: any) => {
      if (!isDrawingLasso || !isDragging.current) return;
      
      isDragging.current = false;
      
      // Finalize the polygon (close it)
      const finalPoints = [...lassoPoints, lassoPoints[0]]; // Add first point to end
      setLassoPoints(finalPoints);

      // Perform selection check
      if (finalPoints.length >= 4) { // Need at least 3 points + closing point
          const lassoPolygon = turfPolygon([finalPoints]);
          const newlySelectedLayers: string[] = [];
          
          Object.entries(rawLayerBoundaries).forEach(([layerName, boundary]) => {
              if (boundary && visibleLayers.has(layerName)) { // Only check visible layers with boundaries
                  // Check if any vertex of the layer boundary is inside the lasso
                  // More robust check: Check intersection of polygons (requires turf/intersect or similar)
                  // Simple check for now: is any boundary vertex inside lasso?
                  const boundaryPolygonFeature: GeoJsonFeature<Polygon> = turfPolygon([boundary]); // Create Feature for checking
                  const isIntersecting = boundary.some(coord => {
                      try {
                        // Use correct imports: turfPoint for point, lassoPolygon Feature
                        return booleanPointInPolygon(turfPoint(coord), lassoPolygon);
                      } catch (e) {
                        console.error("Error in booleanPointInPolygon:", e, coord, lassoPolygon);
                        return false;
                      }
                  });
                  
                  // Check if center of polygon is inside? May need centroid calculation
                  // const center = centroid(turfPolygon([boundary])); // Requires @turf/centroid
                  // const isCenterInside = booleanPointInPolygon(center, lassoPolygon);

                  if (isIntersecting) { // Or use a more robust intersection check
                      newlySelectedLayers.push(layerName);
                  }
              }
          });
          
          setSelectedLayerNames(newlySelectedLayers);
          console.log("[SpatialOverview] Lasso selected layers:", newlySelectedLayers);
          if (onLassoSelect) {
              onLassoSelect(newlySelectedLayers);
          }
      }
      // Don't clear points immediately, keep polygon visible until cleared
      // setIsDrawingLasso(false); // Optionally disable drawing after one lasso
  }, [isDrawingLasso, lassoPoints, rawLayerBoundaries, visibleLayers, onLassoSelect]);

  // --- Lasso UI Controls --- 
  const toggleLasso = () => {
      if (isDrawingLasso) {
          // If currently drawing, disable and clear
          clearLasso();
      }
      setIsDrawingLasso(!isDrawingLasso);
  };

  const clearLasso = () => {
      setIsDrawingLasso(false);
      setLassoPoints([]);
      setSelectedLayerNames([]);
      if (onLassoSelect) {
        onLassoSelect([]);
      }
  };

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
                  getFillColor: selectedLayerNames.includes(layerName) ? 
                                [255, 255, 0, 150] : // Yellow highlight
                                (layerColors[layerName] || [128, 128, 128, 50]), 
                  getLineColor: [255, 255, 255, 150], // White outline
                  getLineWidth: selectedLayerNames.includes(layerName) ? 3 : 1, // Thicker line for selected
                  lineWidthMinPixels: 1,
                  pickable: !isDrawingLasso, // Disable picking layer polygons while drawing lasso
                  autoHighlight: !isDrawingLasso, // Disable auto-highlight while drawing
                  highlightColor: [255, 255, 255, 100],
                  updateTriggers: { // Trigger update when selection changes
                    getFillColor: [selectedLayerNames],
                    getLineWidth: [selectedLayerNames]
                  }
              });
          })
          .filter(layer => layer !== null); // Filter out any null layers from map
          
       // Add Lasso Layer if drawing
       const lassoLayer = isDrawingLasso && lassoPoints.length > 1 ?
            new GeoJsonLayer({
                id: 'lasso-layer',
                // Data for GeoJsonLayer should be a Feature or FeatureCollection
                data: turfPolygon([lassoPoints.length >=3 ? [...lassoPoints, lassoPoints[0]] : lassoPoints ], { layerType: 'lasso' }), // Create a Turf polygon Feature directly
                // Example adding properties if needed later: turfPolygon([lassoPoints.length >=3 ? [...lassoPoints, lassoPoints[0]] : lassoPoints ], { name: 'My Lasso', isTemporary: true }),
                stroked: true,
                filled: lassoPoints.length >= 3, // Only fill if it's a polygon
                getLineColor: [255, 0, 0, 200], // Red line
                getFillColor: [255, 0, 0, 50], // Transparent red fill
                getLineWidth: 2,
                lineWidthMinPixels: 2,
                pickable: false // Lasso itself shouldn't be pickable
            }) : null;

       const finalLayers = [
           ...boundaryLayers,
           lassoLayer
       ].filter(Boolean); // Filter out null lasso layer
          
       console.log(`[SpatialOverviewVisualization] Rendering ${boundaryLayers.length} boundary layers. Lasso active: ${isDrawingLasso}`);
       return finalLayers as (PolygonLayer<any> | GeoJsonLayer<any>)[];

  }, [rawLayerBoundaries, visibleLayers, layerColors, isDrawingLasso, lassoPoints, selectedLayerNames]);

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  const onDeckError = useCallback((error: Error) => console.error('[DeckGL Error - SpatialOverview]', error), []);

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
    <div className={`${styles.spatialOverviewContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
        {/* Removed point count indicator */}
        {/* <div className={styles.resolutionIndicator}>...</div> */}
      <div className={styles.deckGlWrapper}>
        <DeckGL
          views={new OrthographicView({id: 'ortho-overview-view'})} // Unique view id
          viewState={viewState}
          onViewStateChange={onViewStateChange}
          controller={true}
          onError={onDeckError}
          layers={layers}
          getTooltip={getTooltipText} 
          // Add drag handlers for lasso
          onDragStart={handleDragStart}
          onDrag={handleDrag}
          onDragEnd={handleDragEnd}
          // Disable pan/zoom while drawing lasso?
          // controller={isDrawingLasso ? null : true}
          getCursor={({isDragging: deckIsDragging}) => isDrawingLasso ? 'crosshair' : (deckIsDragging ? 'grabbing' : 'grab')}
          style={{ width: '100%', height: '100%', position: 'relative' }}
        />
      </div>
      {/* Controls Overlay */}
      <div className={styles.controlsOverlay}>
          {/* Layer Visibility Controls */}
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
          {/* Lasso Controls */}
          <div className={`${styles.controlSection} ${styles.lassoControls}`}> 
                <h4>Selection Tool</h4>
                <button 
                    onClick={toggleLasso} 
                    disabled={!jobStatus?.results} 
                    // Add activeLasso class when drawing
                    className={isDrawingLasso ? styles.activeLasso : ''} 
                > 
                  {isDrawingLasso ? 'Cancel Lasso' : 'Start Lasso'}
                </button>
                {lassoPoints.length > 0 && (
                    <button onClick={clearLasso}> 
                      Clear Selection
                    </button>
                )}
                {selectedLayerNames.length > 0 && (
                    <div className={styles.selectionInfo}>
                      Selected: {selectedLayerNames.join(', ')}
                    </div>
                )}
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