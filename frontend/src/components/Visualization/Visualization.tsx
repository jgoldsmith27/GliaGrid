import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { ScatterplotLayer } from '@deck.gl/layers';
import { ScreenGridLayer } from '@deck.gl/aggregation-layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters } from '@deck.gl/core';
import styles from './Visualization.module.css';

interface Point {
  x: number;
  y: number;
}

interface VisualizationProps {
  data: {
    ligand: Point[];
    receptor: Point[];
  };
  ligandName: string;
  receptorName: string;
  currentScope: string;
}

// Helper function to calculate initial view state based on data bounds for OrthographicView
const getInitialViewState = (ligandData: Point[], receptorData: Point[]) => {
    const allPoints = [...ligandData, ...receptorData];
    if (allPoints.length === 0) {
        // Default view if no data
        return { target: [0, 0, 0] as [number, number, number], zoom: 1 }; // Use target and zoom, cast type
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
    const width = Math.max(1, maxX - minX); // Ensure width/height are at least 1 to avoid division by zero
    const height = Math.max(1, maxY - minY);
    
    // Calculate zoom level based on fitting the data extent into the viewport
    // This calculation might need further tuning depending on desired padding/behavior
    const zoomX = Math.log2(window.innerWidth / width);
    const zoomY = Math.log2(window.innerHeight / height);
    const zoom = Math.min(zoomX, zoomY) - 1; // Subtract 1 for some padding

    return {
        target: [centerX, centerY, 0] as [number, number, number], // Center the view, cast type explicitly
        zoom: zoom, 
        minZoom: -10, // Allow zooming out much further
        maxZoom: 10 // Adjust as needed
    };
};

// --- Add Heatmap Type Definition ---
type HeatmapType = 'None' | 'Ligand' | 'Receptor';

const Visualization: React.FC<VisualizationProps> = ({ data, ligandName, receptorName, currentScope }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  // --- Add State for Heatmap Type ---
  const [heatmapType, setHeatmapType] = useState<HeatmapType>('None'); 
  // --- Add State for Scatter Plot Visibility ---
  const [showScatterPlot, setShowScatterPlot] = useState(true);

  // Calculate initial view state memoized based on data
  const initialViewState = useMemo(() => {
      console.log("[Visualization] Received data:", data); // Log received data
      const state = getInitialViewState(data.ligand, data.receptor);
      console.log("[Visualization] Calculated initialViewState:", state);
      return state;
  }, [data]); // Recalculate only when data changes

  // --- Add State and Handlers for Zoom --- 
  const [viewState, setViewState] = useState<OrthographicViewState>(initialViewState as OrthographicViewState);

  // Update viewState when initialViewState changes (e.g., data updates)
  // This ensures view resets if data prop changes
  useEffect(() => {
    setViewState(initialViewState as OrthographicViewState);
  }, [initialViewState]);

  const onViewStateChange = useCallback(({ viewState: vs }: ViewStateChangeParameters<OrthographicViewState>) => {
    setViewState(vs as OrthographicViewState);
  }, []);

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
  // --- End Zoom Logic --- 

  // --- Dynamically create layers based on heatmap selection ---
  const layers = useMemo(() => {
    const baseLayers = [
      new ScatterplotLayer<Point>({
        id: 'ligand-layer',
        data: data.ligand,
        getPosition: (d) => [d.x, d.y],
        getRadius: 5,
        getFillColor: [255, 0, 0, 180], // Red for ligands
        pickable: true,
        radiusScale: 5,
        radiusMinPixels: 1,
        radiusMaxPixels: 50,
      }),
      new ScatterplotLayer<Point>({
        id: 'receptor-layer',
        data: data.receptor,
        getPosition: (d) => [d.x, d.y],
        getRadius: 5,
        getFillColor: [0, 0, 255, 180], // Blue for receptors
        pickable: true,
        radiusScale: 5,
        radiusMinPixels: 1,
        radiusMaxPixels: 50,
      }),
    ];

    console.log(`[Visualization Layers] heatmapType: ${heatmapType}, Ligands: ${data.ligand?.length || 0}, Receptors: ${data.receptor?.length || 0}`);

    let heatmapLayer = null;
    // Revert variable name slightly for clarity, still using ScreenGridLayer
    let typedAggregationLayer: ScreenGridLayer<Point> | null = null; 
    // --- Temporarily remove specific heatmap config for debugging ---
    const gridCellSizePixels = 20; // Cell size for ScreenGridLayer

    if (heatmapType === 'Ligand' && data.ligand.length > 0) {
      console.log(`[Visualization Layers] Creating Ligand ScreenGridLayer with ${data.ligand.length} points`);
      typedAggregationLayer = new ScreenGridLayer<Point>({ 
        id: 'ligand-screengrid-layer',
        data: data.ligand,
        getPosition: (d: Point) => [d.x, d.y],
        getWeight: 1, // Use uniform weight for now
        cellSizePixels: gridCellSizePixels,
        colorRange: [ // Use a similar Yellow -> Red range for now
          [255, 255, 178, 0], // Alpha 0 for lowest value
          [254, 204, 92, 64],
          [253, 141, 60, 128],
          [240, 59, 32, 192],
          [189, 0, 38, 255]  // Alpha 255 for highest value
        ],
        gpuAggregation: false,
        aggregation: 'SUM' // Aggregate by sum of weights (or COUNT)
      });
    } else if (heatmapType === 'Receptor' && data.receptor.length > 0) {
      console.log(`[Visualization Layers] Creating Receptor ScreenGridLayer with ${data.receptor.length} points`);
      typedAggregationLayer = new ScreenGridLayer<Point>({ 
        id: 'receptor-screengrid-layer',
        data: data.receptor,
        getPosition: (d: Point) => [d.x, d.y],
        getWeight: 1,
        cellSizePixels: gridCellSizePixels,
        colorRange: [ // Use a similar Light Blue -> Dark Blue range
          [237, 248, 251, 0],
          [179, 205, 227, 64],
          [140, 150, 198, 128],
          [136, 86, 167, 192],
          [129, 15, 124, 255]
        ],
        gpuAggregation: false,
        aggregation: 'SUM'
      });
    }

    // Conditionally include base layers based on showScatterPlot state
    const layersToShow = [];
    if (typedAggregationLayer) {
      layersToShow.push(typedAggregationLayer);
    }
    if (showScatterPlot) {
      layersToShow.push(...baseLayers);
    } else if (!typedAggregationLayer) {
        // If heatmap is off AND showScatterPlot is false, still show base layers? 
        // Or maybe just show nothing? Let's show baseLayers for now if heatmap is off.
        // This case might need refinement based on desired behavior when heatmap is None and toggle is off.
        // For now, if heatmap is None, always show base layers regardless of toggle.
       if (heatmapType === 'None') {
           layersToShow.push(...baseLayers);
       } 
    }
    
    const finalLayers = layersToShow;
    console.log('[Visualization Layers] Final layers array:', finalLayers.map(l => l?.id)); // Log layer IDs

    // Add aggregation layer *before* scatterplots so it's underneath
    return finalLayers;

  }, [data.ligand, data.receptor, heatmapType, showScatterPlot]); // Dependencies for useMemo

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
    // Note: DeckGL might need a resize notification after container size changes
    // This often happens automatically, but keep in mind if issues persist.
  };

  // --- Add onError handler for DeckGL ---
  const onDeckError = useCallback((error: Error) => console.error('[DeckGL Error]', error), []);

  return (
    <div className={`${styles.visualizationContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
      <div className={styles.deckGlWrapper}>
        <DeckGL
          views={new OrthographicView({id: 'ortho-view'})}
          viewState={viewState}
          onViewStateChange={onViewStateChange}
          controller={true}
          onError={onDeckError} // Added onError prop
          layers={layers}
          getTooltip={({object}) => object && `Point: (${object.x.toFixed(2)}, ${object.y.toFixed(2)})`}
          style={{ width: '100%', height: '100%', position: 'relative' }}
        />
        <div className={styles.deckOverlayControls}>
            <div className={styles.deckButtons}>
                <div className={styles.zoomControls}>
                    <button onClick={() => handleZoom(1)} title="Zoom In">+</button>
                    <button onClick={() => handleZoom(-1)} title="Zoom Out">-</button>
                </div>
                <button 
                    className={styles.fullscreenButton}
                    onClick={toggleFullscreen}
                    title={isFullscreen ? "Exit fullscreen" : "Enter fullscreen"}
                >
                    {isFullscreen ? "✕" : "⛶"}
                </button>
            </div>
            {/* --- Updated Legend Area --- */}
            <div className={styles.legendContainer}>
                {/* Scatterplot Legend */}
                <div className={styles.legendItem}>
                    <div className={styles.legendColor} style={{ backgroundColor: 'red' }} />
                    <span>{ligandName}</span>
                </div>
                <div className={styles.legendItem}>
                    <div className={styles.legendColor} style={{ backgroundColor: 'blue' }} />
                    <span>{receptorName}</span>
                </div>
                {/* Heatmap Selection */}
                <div className={styles.heatmapSelector}>
                    <span className={styles.heatmapLabel}>Heatmap:</span>
                    <label>
                        <input type="radio" name="heatmapType" value="None" checked={heatmapType === 'None'} onChange={() => setHeatmapType('None')} />
                        Off
                    </label>
                    <label>
                        <input type="radio" name="heatmapType" value="Ligand" checked={heatmapType === 'Ligand'} onChange={() => setHeatmapType('Ligand')} disabled={data.ligand.length === 0} />
                        Ligand Density
                    </label>
                    <label>
                        <input type="radio" name="heatmapType" value="Receptor" checked={heatmapType === 'Receptor'} onChange={() => setHeatmapType('Receptor')} disabled={data.receptor.length === 0} />
                        Receptor Density
                    </label>
                </div>
                {/* --- Add Scatter Plot Toggle --- */}
                <div className={styles.scatterToggle}>
                    <label>
                        <input 
                            type="checkbox" 
                            checked={showScatterPlot} 
                            onChange={(e) => setShowScatterPlot(e.target.checked)} 
                        />
                        Show Points
                    </label>
                </div>
            </div>
        </div>
      </div>
    </div>
  );
};

export default Visualization; 