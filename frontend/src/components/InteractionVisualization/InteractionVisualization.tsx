import React, { useState, useMemo, useCallback, useEffect, useRef } from 'react';
import { ScatterplotLayer, PolygonLayer, PathLayer } from '@deck.gl/layers';
import { ScreenGridLayer, HeatmapLayer, HexagonLayer } from '@deck.gl/aggregation-layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters, Color, PickingInfo, Layer } from '@deck.gl/core';
import { scaleOrdinal } from 'd3-scale';
import { schemeCategory10 } from 'd3-scale-chromatic';
import { Button, CircularProgress, Box, IconButton, Icon, Tooltip } from '@mui/material';
import styles from './InteractionVisualization.module.css';
import ScaleBar from '../VisualizationHelpers/ScaleBar';

// Interface for individual points from backend (now includes gene)
interface PointWithGene {
  x: number;
  y: number;
  gene: string; // Added gene property
}

// Interface for ligand points (doesn't need gene after processing)
interface LigandPoint {
  x: number;
  y: number;
}

// Type for the data prop from SharedDataStore
interface VisualizationData {
  ligand: LigandPoint[]; // Processed ligand points
  receptor: PointWithGene[]; // Raw receptor component points
  isComplex: boolean;
  receptorName: string; // Original receptor name (e.g., "R1_R2")
  warnings?: string[];
}

type ScopeType = 'whole_tissue' | 'layers' | 'custom';

interface InteractionVisualizationProps {
  data: VisualizationData | null; // Updated data prop type, allow null initially
  ligandName: string; // Keep ligandName for context
  // receptorName is now derived from data prop
  currentScope: 'whole_tissue' | 'layers';
  isLoading: boolean;
  cancelFetch: () => void;
  layerBoundaries?: Record<string, [number, number][] | null>;
}

// --- Centroid Calculation Logic (Ported from Python) ---
// Helper function to calculate distance squared
function distanceSq(p1: { x: number; y: number }, p2: { x: number; y: number }): number {
  const dx = p1.x - p2.x;
  const dy = p1.y - p2.y;
  return dx * dx + dy * dy;
}

// Calculates centroids for complex receptors
function calculateComplexCentroids(
  receptorComponentPoints: PointWithGene[],
  receptorName: string,
  proximityThresholdPixels: number = 100 // Default 100px = 50um
): { x: number; y: number }[] {
  if (!receptorName.includes('_') || receptorComponentPoints.length === 0) {
    // Not complex or no points, return empty centroids
    return [];
  }

  const components = receptorName.split('_');
  const proximityThresholdSq = proximityThresholdPixels * proximityThresholdPixels;

  // Separate points by component gene name
  const componentData: Record<string, { x: number; y: number }[]> = {};
  components.forEach(comp => componentData[comp] = []);
  receptorComponentPoints.forEach(p => {
    if (componentData[p.gene]) {
      componentData[p.gene].push({ x: p.x, y: p.y });
    }
  });

  // Check if all components are present
  if (!components.every(comp => componentData[comp].length > 0)) {
    console.warn(`[CentroidCalc] Cannot form complex ${receptorName}: component(s) missing points.`);
    return [];
  }

  // Determine anchor component (fewest points)
  let anchorComponentName: string | null = null;
  let minCount = Infinity;
  for (const comp of components) {
    if (componentData[comp].length < minCount) {
      minCount = componentData[comp].length;
      anchorComponentName = comp;
    }
  }

  if (!anchorComponentName) return []; // Should not happen if all components present

  const anchorCoords = componentData[anchorComponentName];
  const otherComponentNames = components.filter(c => c !== anchorComponentName);
  const validCentroids: { x: number; y: number }[] = [];

  console.log(`[CentroidCalc] Anchoring complex search on ${anchorComponentName} (${anchorCoords.length} points) for ${receptorName}`);

  for (const anchorPoint of anchorCoords) {
    const potentialClusterPoints: Record<string, { x: number; y: number }> = { [anchorComponentName]: anchorPoint };
    let isPotentialClusterValid = true;

    for (const otherCompName of otherComponentNames) {
      const otherCoords = componentData[otherCompName];
      let closestPoint: { x: number; y: number } | null = null;
      let minDistanceSq = proximityThresholdSq;

      for (const otherPoint of otherCoords) {
        const dSq = distanceSq(anchorPoint, otherPoint);
        if (dSq <= minDistanceSq) {
          minDistanceSq = dSq;
          closestPoint = otherPoint;
        }
      }

      if (closestPoint) {
        potentialClusterPoints[otherCompName] = closestPoint;
      } else {
        // If no point of this component is within threshold, invalidate this anchor
        isPotentialClusterValid = false;
        break;
      }
    }

    if (isPotentialClusterValid) {
      // Calculate centroid
      const clusterPointsArray = Object.values(potentialClusterPoints);
      let sumX = 0, sumY = 0;
      clusterPointsArray.forEach(p => { sumX += p.x; sumY += p.y; });
      const centroid = { x: sumX / clusterPointsArray.length, y: sumY / clusterPointsArray.length };
      validCentroids.push(centroid);
    }
  }

  console.log(`[CentroidCalc] Found ${validCentroids.length} potential complex centroids for ${receptorName}`);
  return validCentroids;
}
// --- End Centroid Calculation Logic ---

const getInitialViewState = (ligandData: LigandPoint[], receptorData: PointWithGene[], isComplex: boolean, complexCentroids: {x: number, y: number}[]) => {
    // Use centroids for bounds if complex, otherwise use component points
    const receptorPointsForBounds = isComplex ? complexCentroids : receptorData;
    const allPoints = [...ligandData, ...receptorPointsForBounds];
    if (allPoints.length === 0) {
        return { target: [0, 0, 0] as [number, number, number], zoom: 1 };
    }

    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    allPoints.forEach(p => {
        minX = Math.min(minX, p.x);
        minY = Math.min(minY, p.y);
        maxX = Math.max(maxX, p.x);
        maxY = Math.max(maxY, p.y);
    });

    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;
    const width = Math.max(1, maxX - minX);
    const height = Math.max(1, maxY - minY);
    
    const zoomX = Math.log2(window.innerWidth / width);
    const zoomY = Math.log2(window.innerHeight / height);
    const zoom = Math.min(zoomX, zoomY) - 1;

    return {
        target: [centerX, centerY, 0] as [number, number, number],
        zoom: zoom, 
        minZoom: -10,
        // Set maxZoom limit. Rationale:
        // - Coordinate units are typically ~2 units/µm.
        // - Orthographic zoom means pixels = worldUnits * 2^zoom.
        // - Scatterplot points have a base radius (e.g., 5 units = 2.5 µm).
        // - At zoom=6, 1 world unit (0.5 µm) = 2^6 = 64 pixels. A 5-unit radius point
        //   would appear very large (~320px radius before radiusScale/clamping).
        // - Zooming further overly magnifies the discrete coordinate space without 
        //   revealing more biological detail relative to point size.
        maxZoom: 6 
    };
};

type DensityScoringType = 'Off' | 'Ligand' | 'Receptor';
type AggregationLayerType = 'ScreenGrid' | 'Hexagon' | 'Heatmap';
type PointsDisplayType = 'off' | 'ligands' | 'receptors' | 'both';

const InteractionVisualization: React.FC<InteractionVisualizationProps> = ({ data, ligandName, currentScope, isLoading, cancelFetch, layerBoundaries }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [densityScoringType, setDensityScoringType] = useState<DensityScoringType>('Off'); 
  const [aggregationLayerType, setAggregationLayerType] = useState<AggregationLayerType>('ScreenGrid');
  const [pointsDisplayType, setPointsDisplayType] = useState<PointsDisplayType>('both');

  const [visibleLayers, setVisibleLayers] = useState<Set<string>>(new Set());
  const [layerColors, setLayerColors] = useState<Record<string, Color>>({});
  const [uniqueLayerNames, setUniqueLayerNames] = useState<string[]>([]);
  const [allLayersVisible, setAllLayersVisible] = useState(true);

  // --- Ruler State ---
  const [isRulerActive, setIsRulerActive] = useState(false);
  const [rulerPoints, setRulerPoints] = useState<[number, number][]>([]);
  const [rulerDistanceMicrons, setRulerDistanceMicrons] = useState<number | null>(null);
  const [rulerHoverCoord, setRulerHoverCoord] = useState<[number, number] | null>(null);

  const deckContainerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (layerBoundaries) {
      const validLayers = Object.keys(layerBoundaries).filter(l => layerBoundaries[l] !== null).sort();
      
      setUniqueLayerNames(validLayers);

      const colorScale = scaleOrdinal<string, string>(schemeCategory10);
      const colors: Record<string, Color> = {};
      validLayers.forEach((layer) => {
        const hexColor = colorScale(layer);
        const r = parseInt(hexColor.slice(1, 3), 16);
        const g = parseInt(hexColor.slice(3, 5), 16);
        const b = parseInt(hexColor.slice(5, 7), 16);
        colors[layer] = [r, g, b, 180];
      });
      setLayerColors(colors);

      setVisibleLayers(new Set(validLayers));
      setAllLayersVisible(true);
      console.log("[InteractionVisualization] Layers derived from boundaries:", validLayers);
    } else {
      console.warn("[InteractionVisualization] layerBoundaries prop not provided. Layer controls will be empty.");
      setUniqueLayerNames([]);
      setLayerColors({});
      setVisibleLayers(new Set());
      setAllLayersVisible(true);
    }
  }, [layerBoundaries]);

  const initialViewState = useMemo(() => {
      console.log("[InteractionVisualization] Received data:", data);
      if (!data) return { target: [0, 0, 0] as [number, number, number], zoom: 1 }; // Handle null data
      const ligandPts = data.ligand || [];
      const receptorPts = data.receptor || []; // Raw component points
      const state = getInitialViewState(ligandPts, receptorPts, data.isComplex, []);
      console.log("[InteractionVisualization] Calculated initialViewState:", state);
      return state;
  }, [data]);

  // Initialize state including zoom bounds
  const [viewState, setViewState] = useState<OrthographicViewState>({
    ...initialViewState,
    minZoom: initialViewState.minZoom ?? -10,
    maxZoom: initialViewState.maxZoom ?? 6,
  } as OrthographicViewState);

  useEffect(() => {
    // Update state when initialViewState changes, preserving bounds
    setViewState({
        ...initialViewState,
        minZoom: initialViewState.minZoom ?? -10,
        maxZoom: initialViewState.maxZoom ?? 6,
    } as OrthographicViewState);
  }, [initialViewState]);

  const onViewStateChange = useCallback(({ viewState: vs }: ViewStateChangeParameters<OrthographicViewState>) => {
    // Clamp zoom within bounds when updating state
    let newZoom = vs.zoom;
    // Clamp zoom only if it's a number, respecting the defined min/max bounds.
    if (typeof newZoom === 'number') {
        newZoom = Math.max(initialViewState.minZoom ?? -10, Math.min(initialViewState.maxZoom ?? 6, newZoom));
    }
    setViewState({...vs, zoom: newZoom, minZoom: initialViewState.minZoom ?? -10, maxZoom: initialViewState.maxZoom ?? 6 } as OrthographicViewState);
  }, [initialViewState.minZoom, initialViewState.maxZoom]);

  const handleZoom = useCallback((zoomIncrement: number) => {
    setViewState((currentViewState) => {
        const state = currentViewState as OrthographicViewState;
        const currentZoom = state.zoom ?? 0;
        let nextZoom: number = 0;

        if (typeof currentZoom === 'number') {
            nextZoom = currentZoom + zoomIncrement;
        } else {
            console.warn("Unexpected zoom format in Orthographic viewState:", currentZoom);
            nextZoom = zoomIncrement;
        }
        
        return {
            ...state,
            zoom: nextZoom
        } as OrthographicViewState;
    });
  }, []);

  const layers = useMemo(() => {
    // Check for null data or loading state AT THE START
    if (!data || isLoading) {
      return []; // Return empty if data is null or still loading
    }

    // --- Type Guard: 'data' is now guaranteed non-null --- 

    // Get receptor name from data 
    const receptorName = data.receptorName;

    // --- Calculate Centroids INSIDE this hook --- 
    const complexCentroids = data.isComplex
        ? calculateComplexCentroids(data.receptor, data.receptorName)
        : []; // Calculate only if complex, otherwise empty

    // --- Layer Boundary Logic (unchanged) --- 
    const boundaryLayers = layerBoundaries 
      ? Object.entries(layerBoundaries)
          // Filter for layers that are visible AND have non-null boundaries
          .filter(([layerName, boundary]) => boundary && visibleLayers.has(layerName))
          .map(([layerName, boundary]) => {
              if (!boundary) return null; // Should be filtered, but type guard
              
              // Define type for polygon data item
              interface PolygonDataItem { polygon: [number, number][]; layerName: string; }
              
              return new PolygonLayer<PolygonDataItem>({
                  id: `boundary-${layerName}`,
                  data: [{ polygon: boundary, layerName: layerName }], // Data format
                  getPolygon: (d: PolygonDataItem) => d.polygon,
                  // Use layerColors state for fill, provide default
                  getFillColor: layerColors[layerName] || [128, 128, 128, 50],
                  getLineColor: [255, 255, 255, 100], // White outline, slightly transparent
                  getLineWidth: 1,
                  lineWidthMinPixels: 1,
                  pickable: false, // Usually don't need to pick boundaries here
                  autoHighlight: false,
              });
          })
          .filter(layer => layer !== null) // Filter out nulls
      : []; // Empty array if no layerBoundaries prop

    // --- Determine Receptor Points to Plot --- 
    // Use calculated centroids if complex, otherwise the original component points
    const receptorPointsToPlot = data.isComplex ? complexCentroids : data.receptor;

    console.log(`[InteractionVisualization] Plotting receptors for ${receptorName}. Complex: ${data.isComplex}. Plotting ${receptorPointsToPlot.length} points.`);

    // --- Scatterplot Layers --- 
    const baseLayers = [];

    // Ligand Layer (always uses data.ligand)
    if (pointsDisplayType === 'ligands' || pointsDisplayType === 'both') {
        baseLayers.push(
          new ScatterplotLayer<LigandPoint>({
            id: 'ligand-layer',
            data: data.ligand || [], // data is non-null here
            getPosition: (d) => [d.x, d.y],
            getRadius: 5,
            getFillColor: [255, 0, 0, 180], // Red
            pickable: !isRulerActive, // Disable picking if ruler active
            radiusScale: 5,
            radiusMinPixels: 1,
            radiusMaxPixels: 50,
            visible: true, // Visibility controlled by inclusion in array
          })
        );
    }

    // Receptor Layer (uses centroids or original points)
    if (pointsDisplayType === 'receptors' || pointsDisplayType === 'both') {
        baseLayers.push(
          new ScatterplotLayer<PointWithGene | { x: number; y: number }>({ // Type union for data
            id: 'receptor-layer',
            data: receptorPointsToPlot, // Use determined points
            getPosition: (d) => [d.x, d.y],
            getRadius: 5,
            getFillColor: [0, 0, 255, 180], // Blue
            pickable: !isRulerActive, // Disable picking if ruler active
            radiusScale: 5,
            radiusMinPixels: 1,
            radiusMaxPixels: 50,
            visible: true, // Visibility controlled by inclusion in array
          })
        );
    }

    console.log(`[InteractionVisualization Layers] pointsDisplay: ${pointsDisplayType}, Aggregation: ${densityScoringType}, Ligands: ${data.ligand?.length || 0}, Receptor points plotted: ${receptorPointsToPlot.length}`);

    let aggregationLayer: ScreenGridLayer<any> | HeatmapLayer<any> | HexagonLayer<any> | null = null;
    const gridCellSizePixels = 20; 
    const hexagonRadius = 15;     
    const hexagonCoverage = 0.9;  
    const heatmapRadiusPixels = 40;
    const heatmapIntensity = 1;   
    const heatmapThreshold = 0.05;

    let layerData: any[] | null = null;
    let layerIdPrefix = '';
    if (densityScoringType === 'Ligand') {
        layerData = data.ligand || []; // data is non-null here
        layerIdPrefix = 'ligand';
    } else if (densityScoringType === 'Receptor') {
        // Use component points for density, not centroids
        layerData = data.receptor || []; // data is non-null here
        layerIdPrefix = 'receptor';
    }

    if (layerData && layerData.length > 0) { 
      const commonProps = {
          data: layerData,
          getPosition: (d: any) => [d.x, d.y] as [number, number],
          getWeight: 1,
          pickable: false,
      };
      const colorRange = densityScoringType === 'Receptor' ? [
          [237, 248, 251, 0],
          [179, 205, 227, 64],
          [140, 150, 198, 128],
          [136, 86, 167, 192],
          [129, 15, 124, 255]
      ] : [
          [255, 255, 178, 0],
          [254, 204, 92, 64],
          [253, 141, 60, 128],
          [240, 59, 32, 192],
          [189, 0, 38, 255]
      ];
       const heatmapColorRange = densityScoringType === 'Receptor' ? [
                [179, 205, 227, 64],
                [140, 150, 198, 128],
                [136, 86, 167, 192],
                [129, 15, 124, 255]
            ] : [
                [254, 204, 92, 64],
                [253, 141, 60, 128],
                [240, 59, 32, 192],
                [189, 0, 38, 255]
            ];

      if (aggregationLayerType === 'ScreenGrid') {
          console.log(`[InteractionVisualization Layers] Creating ${layerIdPrefix} ScreenGridLayer`);
          aggregationLayer = new ScreenGridLayer<any>({ 
            ...commonProps,
            id: `${layerIdPrefix}-screengrid-layer`, 
            cellSizePixels: gridCellSizePixels,
            colorRange: colorRange as Color[],
            gpuAggregation: false,
            aggregation: 'SUM' 
          });
      } else if (aggregationLayerType === 'Hexagon') {
          console.log(`[InteractionVisualization Layers] Creating ${layerIdPrefix} HexagonLayer`);
          aggregationLayer = new HexagonLayer<any>({
              ...commonProps,
              id: `${layerIdPrefix}-hexagon-layer`,
              radius: hexagonRadius,
              coverage: hexagonCoverage,
              colorRange: colorRange as Color[],
          });
      } else {
          console.log(`[InteractionVisualization Layers] Creating ${layerIdPrefix} HeatmapLayer`);
          aggregationLayer = new HeatmapLayer<any>({
            ...commonProps,
            id: `${layerIdPrefix}-heatmap-layer`,
            radiusPixels: heatmapRadiusPixels,
            intensity: heatmapIntensity,
            threshold: heatmapThreshold,
            colorRange: heatmapColorRange as Color[],
            aggregation: 'SUM' 
          });
      }
    }

    // --- Add Ruler Layers (Path and Points) ---
    let rulerLayers: Layer[] = [];
    let rulerPathData: [number, number][] | null = null;

    if (rulerPoints.length === 1 && rulerHoverCoord) {
        rulerPathData = [rulerPoints[0], rulerHoverCoord];
    } else if (rulerPoints.length === 2) {
        rulerPathData = rulerPoints;
    }

    if (rulerPathData) {
        rulerLayers.push(new PathLayer({
            id: 'ruler-path',
            data: [{ path: rulerPathData }],
            getPath: d => d.path,
            getColor: [0, 200, 200, 200], // Cyan with transparency
            getWidth: 2,
            widthMinPixels: 2,
            billboard: true, 
            pickable: false,
        }));
    }
    
    if (rulerPoints.length > 0) {
        rulerLayers.push(new ScatterplotLayer({
            id: 'ruler-points',
            data: rulerPoints.map(p => ({ coordinates: p })), 
            getPosition: d => d.coordinates,
            getFillColor: [0, 200, 200, 255],
            getRadius: 5, 
            radiusScale: 1, 
            radiusMinPixels: 3, 
            radiusMaxPixels: 10, 
            pickable: false,
        }));
    }

    // --- Combine Layers --- 
    const finalLayers = [
        ...(boundaryLayers as Layer[]), 
        aggregationLayer, 
        ...baseLayers, 
        ...rulerLayers 
    ].filter(Boolean); 

    return finalLayers as Layer[];
  }, [data, densityScoringType, pointsDisplayType, aggregationLayerType, isLoading, layerBoundaries, visibleLayers, layerColors, rulerPoints, isRulerActive, rulerHoverCoord]);

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  const onDeckError = useCallback((error: Error) => console.error('[DeckGL Error]', error), []);

  const toggleLayerVisibility = useCallback((layerName: string) => {
    setVisibleLayers(prev => {
        const next = new Set(prev);
        let allVisibleAfterToggle = false;
        if (next.has(layerName)) {
            next.delete(layerName);
        } else {
            next.add(layerName);
            if (next.size === uniqueLayerNames.length) {
                 allVisibleAfterToggle = true;
            }
        }
        setAllLayersVisible(allVisibleAfterToggle || next.size === uniqueLayerNames.length); 
        return next;
    });
  }, [uniqueLayerNames]);

  const toggleAllLayers = useCallback(() => {
    const nextVisibility = !allLayersVisible;
    setAllLayersVisible(nextVisibility);
    if (nextVisibility) {
      setVisibleLayers(new Set(uniqueLayerNames));
    } else {
      setVisibleLayers(new Set());
    }
  }, [allLayersVisible, uniqueLayerNames]);

  // --- Ruler Controls ---
  const toggleRuler = useCallback(() => {
    setIsRulerActive(prev => {
        const nextIsActive = !prev;
        if (!nextIsActive) { // Clear on deactivation
          setRulerPoints([]);
          setRulerDistanceMicrons(null);
        }
        return nextIsActive;
    });
  }, []);

  const clearRuler = useCallback(() => {
    setRulerPoints([]);
    setRulerDistanceMicrons(null);
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
        }
        return nextPoints;
      });
    }
  }, [isRulerActive]);

  const handleHover = useCallback((info: PickingInfo, event: any) => {
    if (isRulerActive && rulerPoints.length === 1 && info.coordinate) {
      setRulerHoverCoord([info.coordinate[0], info.coordinate[1]]);
    } else {
      if (rulerHoverCoord !== null) {
         setRulerHoverCoord(null);
      }
    }
  }, [isRulerActive, rulerPoints.length, rulerHoverCoord]);

  return (
    <div ref={deckContainerRef} className={`${styles.interactionVisualizationContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
      <div className={styles.deckGlWrapper}>
        {isLoading && (
          <Box
            sx={{
              position: 'absolute',
              top: 0,
              left: 0,
              right: 0,
              bottom: 0,
              backgroundColor: 'rgba(0, 0, 0, 0.5)',
              display: 'flex',
              flexDirection: 'column',
              justifyContent: 'center',
              alignItems: 'center',
              zIndex: 10,
              color: 'white'
            }}
          >
            <CircularProgress color="inherit" />
            <Button 
              variant="contained" 
              color="secondary" 
              onClick={cancelFetch} 
              sx={{ mt: 2 }}
            >
              Cancel Loading
            </Button>
          </Box>
        )}
        <DeckGL
          views={new OrthographicView({id: 'interaction-ortho-view'})}
          initialViewState={initialViewState as OrthographicViewState}
          viewState={viewState}
          onViewStateChange={onViewStateChange}
          controller={true}
          onError={onDeckError}
          layers={layers as any[]}
          onClick={handleDeckClick}
          onHover={handleHover}
          getCursor={({isDragging: deckIsDragging}) => isRulerActive ? 'crosshair' : (deckIsDragging ? 'grabbing' : 'grab')}
          style={{ width: '100%', height: '100%', position: 'relative' }}
        />
        {viewState && typeof viewState.zoom === 'number' && (
          <ScaleBar 
            currentZoom={viewState.zoom} 
            unitsPerMicron={2}
          />
        )}
        {rulerDistanceMicrons !== null && (
            <div className={styles.measurementDisplay}>
                <Icon sx={{ fontSize: 16, mr: 0.5, color: '#00c8c8' }}>straighten</Icon>
                <span>{`Distance: ${rulerDistanceMicrons.toFixed(1)} µm`}</span>
                <Tooltip title="Clear Measurement">
                    <IconButton size="small" onClick={clearRuler} sx={{ ml: 0.5, p: '2px' }}>
                        <Icon fontSize="inherit">close</Icon>
                    </IconButton>
                 </Tooltip>
            </div>
        )}
        {!isLoading && (
          <div className={styles.deckOverlayControls}>
              <div className={styles.deckButtons}>
                  <div className={styles.zoomControls}>
                      <Tooltip title="Zoom In">
                          <button onClick={() => handleZoom(1)} >+</button>
                      </Tooltip>
                      <Tooltip title="Zoom Out">
                           <button onClick={() => handleZoom(-1)} >-</button>
                      </Tooltip>
                  </div>
                  <Tooltip title={isRulerActive ? "Deactivate Ruler" : "Activate Ruler"}>
                       <IconButton onClick={toggleRuler} color={isRulerActive ? "primary" : "default"} size="small">
                           <Icon>straighten</Icon>
                       </IconButton>
                  </Tooltip>
                  <Tooltip title={isFullscreen ? "Exit fullscreen" : "Enter fullscreen"}>
                      <button 
                          className={styles.fullscreenButton}
                          onClick={toggleFullscreen}
                      >
                          {isFullscreen ? "✕" : "⛶"}
                      </button>
                   </Tooltip>
              </div>
              <div className={styles.legendAndControlsContainer}>
                  <div className={styles.legendContainer}>
                      {(pointsDisplayType === 'ligands' || pointsDisplayType === 'both') && (
                          <div className={styles.legendItem}>
                              <div className={styles.legendColor} style={{ backgroundColor: 'red' }} />
                              <span>{ligandName}</span>
                          </div>
                      )}
                      {(pointsDisplayType === 'receptors' || pointsDisplayType === 'both') && (
                          <div className={styles.legendItem}>
                              <div className={styles.legendColor} style={{ backgroundColor: 'blue' }} />
                              <span>{data ? data.receptorName : 'Receptor'}</span>
                          </div>
                      )}
                  </div>
                  <div className={styles.otherControlsContainer}>
                      <div className={styles.controlGroup}>
                          <span className={styles.controlLabel}>Density Scoring:</span>
                          <label>
                              <input type="radio" name="densityScoring" value="Off" checked={densityScoringType === 'Off'} onChange={() => setDensityScoringType('Off')} />
                              Off
                          </label>
                          <label>
                              <input type="radio" name="densityScoring" value="Ligand" checked={densityScoringType === 'Ligand'} onChange={() => setDensityScoringType('Ligand')} disabled={!data || data.ligand.length === 0} />
                              Ligand Density
                          </label>
                          <label>
                              <input type="radio" name="densityScoring" value="Receptor" checked={densityScoringType === 'Receptor'} onChange={() => setDensityScoringType('Receptor')} disabled={!data || data.receptor.length === 0} />
                              Receptor Density
                          </label>
                      </div>
                      <div className={styles.controlGroup} >
                           <span className={styles.controlLabel}>Density Style:</span>
                          <label>
                              <input type="radio" name="aggType" value="ScreenGrid" checked={aggregationLayerType === 'ScreenGrid'} onChange={() => setAggregationLayerType('ScreenGrid')} disabled={densityScoringType === 'Off'} />
                              Grid
                          </label>
                           <label>
                              <input type="radio" name="aggType" value="Hexagon" checked={aggregationLayerType === 'Hexagon'} onChange={() => setAggregationLayerType('Hexagon')} disabled={densityScoringType === 'Off'} />
                              Hexagon
                          </label>
                          <label>
                              <input type="radio" name="aggType" value="Heatmap" checked={aggregationLayerType === 'Heatmap'} onChange={() => setAggregationLayerType('Heatmap')} disabled={densityScoringType === 'Off'} />
                              Heatmap
                          </label>
                      </div>
                      <div className={styles.controlGroup}>
                          <span className={styles.controlLabel}>Points:</span>
                          <label>
                              <input type="radio" name="pointsDisplay" value="off" checked={pointsDisplayType === 'off'} onChange={() => setPointsDisplayType('off')} />
                              Off
                          </label>
                               <label>
                                  <input type="radio" name="pointsDisplay" value="ligands" checked={pointsDisplayType === 'ligands'} onChange={() => setPointsDisplayType('ligands')} disabled={!data || data.ligand.length === 0} />
                                  Ligands Only
                              </label>
                               <label>
                                  <input type="radio" name="pointsDisplay" value="receptors" checked={pointsDisplayType === 'receptors'} onChange={() => setPointsDisplayType('receptors')} disabled={!data || data.receptor.length === 0} />
                                  Receptors Only
                              </label>
                               <label>
                                  <input type="radio" name="pointsDisplay" value="both" checked={pointsDisplayType === 'both'} onChange={() => setPointsDisplayType('both')} disabled={!data || data.ligand.length === 0 || data.receptor.length === 0}/>
                                  Show Both
                              </label>
                      </div>
                  </div>
                  <div className={styles.controlSection}>
                      <div className={styles.layerHeader}>
                          <h4>Layers</h4>
                          <button 
                              onClick={toggleAllLayers} 
                              title={allLayersVisible ? "Hide All Layers" : "Show All Layers"}
                              className={styles.masterLayerToggle}
                              disabled={uniqueLayerNames.length === 0}
                          >
                              {allLayersVisible ? "Hide All" : "Show All"}
                          </button>
                      </div>
                      <div className={styles.layerLegendContainer}>
                        {uniqueLayerNames.length > 0 ? (
                          uniqueLayerNames.map((layerName) => (
                            <div 
                              key={layerName} 
                              className={styles.layerControlItem} 
                              onClick={() => toggleLayerVisibility(layerName)} 
                              title={`Toggle ${layerName}`}
                            >
                              <input
                                type="checkbox"
                                readOnly
                                className={styles.layerToggleCheckbox}
                                checked={visibleLayers.has(layerName)}
                              />
                              <span 
                                className={styles.legendColorSwatch} 
                                style={{ backgroundColor: `rgba(${(layerColors[layerName] || [128, 128, 128, 180]).join(',')})` }}
                              ></span>
                              <span className={styles.layerNameText}>{layerName}</span>
                            </div>
                          ))
                        ) : (
                          <p>No layers identified.</p>
                        )}
                      </div>
                  </div>
              </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default InteractionVisualization; 