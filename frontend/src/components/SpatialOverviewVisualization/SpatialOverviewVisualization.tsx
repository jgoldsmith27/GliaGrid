import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { ScatterplotLayer } from '@deck.gl/layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters, Color, PickingInfo } from '@deck.gl/core';
import { scaleOrdinal } from 'd3-scale';
import { schemeCategory10 } from 'd3-scale-chromatic'; // Example color scheme
import styles from './SpatialOverviewVisualization.module.css';

interface AllPointsData {
    x: number;
    y: number;
    layer: string;
}

interface SpatialOverviewVisualizationProps {
  allPointsData: AllPointsData[];
  // Add other props if needed later, e.g., for lasso tool integration
}

// Helper to get initial view state based on all points
const getInitialViewStateForAllPoints = (allPoints: AllPointsData[]) => {
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
    // Adjust zoom calculation if needed for overview
    const zoomX = Math.log2(window.innerWidth / width);
    const zoomY = Math.log2(window.innerHeight / height);
    const zoom = Math.min(zoomX, zoomY) - 1; // Start with some padding
    return {
        target: [centerX, centerY, 0] as [number, number, number],
        zoom: zoom,
        minZoom: -10,
        maxZoom: 10
    };
};

// Helper to convert hex/named colors from d3 to RGBA array for DeckGL
const colorToRgba = (cssColor: string): Color => {
    const ctx = document.createElement('canvas').getContext('2d');
    if (!ctx) return [128, 128, 128, 255];
    ctx.fillStyle = cssColor;
    const colorStr = ctx.fillStyle; // Use the computed style
    let r = 128, g = 128, b = 128;

    // Check for rgba format first
    let match = colorStr.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*\d*\.?\d+)?\)/); // Corrected Regex
    if (match) {
         [r, g, b] = match.slice(1, 4).map(Number); // Get only R, G, B
    } else if (colorStr.startsWith('#')) {
        // Basic hex handling
        const hex = colorStr.replace('#', '');
        if (hex.length === 3) {
            [r, g, b] = [hex[0], hex[1], hex[2]].map(char => parseInt(char + char, 16));
        } else if (hex.length === 6) {
             r = parseInt(hex.substring(0, 2), 16);
             g = parseInt(hex.substring(2, 4), 16);
             b = parseInt(hex.substring(4, 6), 16);
        }
    }
    return [r, g, b, 200]; // Apply consistent alpha
};

const SpatialOverviewVisualization: React.FC<SpatialOverviewVisualizationProps> = ({ allPointsData }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [visibleLayers, setVisibleLayers] = useState<Set<string>>(new Set());
  const [layerColors, setLayerColors] = useState<Record<string, Color>>({});

  // Process data to find unique layers and assign colors
  useEffect(() => {
    if (allPointsData) {
        const uniqueLayers = [...new Set(allPointsData.map(p => p.layer))];
        // Use d3-scale for categorical colors
        const colorScale = scaleOrdinal<string, string>(schemeCategory10);
        const colors: Record<string, Color> = {};
        uniqueLayers.forEach((layer) => {
            colors[layer] = colorToRgba(colorScale(layer));
        });
        setLayerColors(colors);
        setVisibleLayers(new Set(uniqueLayers)); // Start with all layers visible
        console.log("[SpatialOverviewVisualization] Layers found:", uniqueLayers);
        console.log("[SpatialOverviewVisualization] Assigned colors:", colors);
    } else {
         setLayerColors({});
         setVisibleLayers(new Set());
    }
  }, [allPointsData]);

  const initialViewState = useMemo(() => {
      console.log("[SpatialOverviewVisualization] Received data:", allPointsData?.length);
      const state = getInitialViewStateForAllPoints(allPointsData || []);
      console.log("[SpatialOverviewVisualization] Calculated initialViewState:", state);
      return state;
  }, [allPointsData]);

  const [viewState, setViewState] = useState<OrthographicViewState>(initialViewState as OrthographicViewState);

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
        let nextZoom = 0;
        if (typeof currentZoom === 'number') {
            nextZoom = currentZoom + zoomIncrement;
        } else {
             console.warn("Unexpected zoom format:", currentZoom);
             nextZoom = zoomIncrement;
        }
        return { ...state, zoom: nextZoom } as OrthographicViewState;
    });
  }, []);

  const handleLayerToggle = useCallback((layerName: string) => {
    setVisibleLayers(prevVisibleLayers => {
        const newVisibleLayers = new Set(prevVisibleLayers);
        if (newVisibleLayers.has(layerName)) {
            newVisibleLayers.delete(layerName);
        } else {
            newVisibleLayers.add(layerName);
        }
        console.log("[SpatialOverviewVisualization] Visible layers updated:", newVisibleLayers);
        return newVisibleLayers;
    });
  }, []);

  const layers = useMemo(() => {
      if (!allPointsData) return [];

      const filteredData = allPointsData.filter(p => visibleLayers.has(p.layer));
      console.log(`[SpatialOverviewVisualization] Rendering ${filteredData.length} points for layers:`, visibleLayers);


      return [
          new ScatterplotLayer<AllPointsData>({
              id: 'all-points-layer',
              data: filteredData,
              getPosition: d => [d.x, d.y],
              getFillColor: d => layerColors[d.layer] || [128, 128, 128, 100], // Default gray
              getRadius: 3, // Potentially smaller radius for overview
              pickable: true, // Keep pickable for lasso tool later
              radiusScale: 3,
              radiusMinPixels: 1,
              radiusMaxPixels: 30, // Adjust max radius
              updateTriggers: { // Add update trigger for color changes
                  getFillColor: [layerColors]
              }
          }),
          // Add EditableGeoJsonLayer here later for lasso tool
      ];
  }, [allPointsData, visibleLayers, layerColors]);

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  const onDeckError = useCallback((error: Error) => console.error('[DeckGL Error - SpatialOverview]', error), []);

  // Type guard for tooltip object
  const isAllPointsData = (obj: any): obj is AllPointsData => {
      return typeof obj === 'object' && obj !== null && 'layer' in obj && 'x' in obj && 'y' in obj;
  };

  // Tooltip function using the type guard
  const getTooltipText = (info: PickingInfo | null | undefined): string | null => {
      if (info && info.object && isAllPointsData(info.object)) {
          const point = info.object;
          return `Layer: ${point.layer}\nPoint: (${point.x.toFixed(2)}, ${point.y.toFixed(2)})`;
      }
      return null;
  };

  return (
    <div className={`${styles.spatialOverviewContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
      <div className={styles.deckGlWrapper}>
        <DeckGL
          views={new OrthographicView({id: 'ortho-overview-view'})} // Unique view id
          viewState={viewState}
          onViewStateChange={onViewStateChange}
          controller={true}
          onError={onDeckError}
          layers={layers}
          getTooltip={getTooltipText} // Use the robust tooltip function
          style={{ width: '100%', height: '100%', position: 'relative' }}
        />
        <div className={styles.deckOverlayControls}>
            {/* Basic Controls: Zoom, Fullscreen */}
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

            {/* Layer Legend and Toggles */}
            <div className={styles.layerLegendContainer}>
                 <h4>Layers</h4>
                 {Object.entries(layerColors).sort(([a], [b]) => a.localeCompare(b)).map(([layerName, color]) => ( // Sort layers alphabetically
                    <div key={layerName} className={styles.layerLegendItem}>
                         <input
                            type="checkbox"
                            checked={visibleLayers.has(layerName)}
                            onChange={() => handleLayerToggle(layerName)}
                            id={`layer-toggle-${layerName}`}
                            className={styles.layerToggleCheckbox}
                         />
                         <label htmlFor={`layer-toggle-${layerName}`} className={styles.layerToggleLabel}>
                            <span
                                className={styles.legendColorSwatch}
                                style={{ backgroundColor: `rgba(${color[0]}, ${color[1]}, ${color[2]}, ${color[3] / 255})` }}
                             />
                            <span className={styles.layerNameText}>{layerName}</span>
                         </label>
                    </div>
                 ))}
            </div>
            {/* Add Lasso Tool controls here later */}
        </div>
      </div>
    </div>
  );
};

export default SpatialOverviewVisualization; 