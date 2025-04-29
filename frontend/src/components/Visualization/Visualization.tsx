import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { ScatterplotLayer } from '@deck.gl/layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters } from '@deck.gl/core';
import styles from './Visualization.module.css';
import DensityMapControls from './DensityMapControls';

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

const Visualization: React.FC<VisualizationProps> = ({ data, ligandName, receptorName, currentScope }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);

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

  const layers = [
    new ScatterplotLayer<Point>({
      id: 'ligand-layer',
      data: data.ligand,
      getPosition: (d) => [d.x, d.y],
      getRadius: 5,
      getFillColor: [255, 0, 0, 180], // Red for ligands (added alpha)
      pickable: true,
      radiusScale: 5, // Adjust size scaling if needed
      radiusMinPixels: 1,
      radiusMaxPixels: 50,
    }),
    new ScatterplotLayer<Point>({
      id: 'receptor-layer',
      data: data.receptor,
      getPosition: (d) => [d.x, d.y],
      getRadius: 5,
      getFillColor: [0, 0, 255, 180], // Blue for receptors (added alpha)
      pickable: true,
      radiusScale: 5,
      radiusMinPixels: 1,
      radiusMaxPixels: 50,
    }),
  ];

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
    // Note: DeckGL might need a resize notification after container size changes
    // This often happens automatically, but keep in mind if issues persist.
  };

  return (
    <div className={`${styles.visualizationContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
      <div className={styles.deckGlWrapper}>
        <DeckGL
          views={new OrthographicView({id: 'ortho-view'})}
          viewState={viewState}
          onViewStateChange={onViewStateChange}
          controller={true}
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
            <div className={styles.legend}>
                <div className={styles.legendItem}>
                    <div className={styles.legendColor} style={{ backgroundColor: 'red' }} />
                    <span>{ligandName}</span>
                </div>
                <div className={styles.legendItem}>
                    <div className={styles.legendColor} style={{ backgroundColor: 'blue' }} />
                    <span>{receptorName}</span>
                </div>
            </div>
        </div>
      </div>

      {!isFullscreen && (
        <DensityMapControls 
          ligandName={ligandName}
          receptorName={receptorName}
          currentScope={currentScope}
        />
      )}
    </div>
  );
};

export default Visualization; 