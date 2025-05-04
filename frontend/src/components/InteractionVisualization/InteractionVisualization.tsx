import React, { useState, useMemo, useCallback, useEffect } from 'react';
import { ScatterplotLayer } from '@deck.gl/layers';
import { ScreenGridLayer, HeatmapLayer, HexagonLayer } from '@deck.gl/aggregation-layers';
import DeckGL from '@deck.gl/react';
import { OrthographicView, OrthographicViewState, ViewStateChangeParameters, Color } from '@deck.gl/core';
import { Button, CircularProgress, Box } from '@mui/material';
import styles from './InteractionVisualization.module.css';

interface Point {
  x: number;
  y: number;
}

interface AllPointsData {
    x: number;
    y: number;
    layer: string;
}

type ScopeType = 'whole_tissue' | 'layers' | 'custom';

interface InteractionVisualizationProps {
  data: {
    ligand: Point[];
    receptor: Point[];
  };
  ligandName: string;
  receptorName: string;
  currentScope: 'whole_tissue' | 'layers';
  isLoading: boolean;
  cancelFetch: () => void;
}

const getInitialViewState = (ligandData: Point[], receptorData: Point[]) => {
    const allPoints = [...ligandData, ...receptorData];
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
        maxZoom: 10
    };
};

type DensityScoringType = 'Off' | 'Ligand' | 'Receptor';
type AggregationLayerType = 'ScreenGrid' | 'Hexagon' | 'Heatmap';
type PointsDisplayType = 'off' | 'ligands' | 'receptors' | 'both';

const InteractionVisualization: React.FC<InteractionVisualizationProps> = ({ data, ligandName, receptorName, currentScope, isLoading, cancelFetch }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [densityScoringType, setDensityScoringType] = useState<DensityScoringType>('Off'); 
  const [aggregationLayerType, setAggregationLayerType] = useState<AggregationLayerType>('ScreenGrid');
  const [pointsDisplayType, setPointsDisplayType] = useState<PointsDisplayType>('both');

  const initialViewState = useMemo(() => {
      console.log("[InteractionVisualization] Received data:", data);
      const state = getInitialViewState(data.ligand, data.receptor);
      console.log("[InteractionVisualization] Calculated initialViewState:", state);
      return state;
  }, [data]);

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
    if (isLoading) return [];

    const baseLayers = [
      new ScatterplotLayer<Point>({
        id: 'ligand-layer',
        data: data.ligand,
        getPosition: (d) => [d.x, d.y],
        getRadius: 5,
        getFillColor: [255, 0, 0, 180],
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
        getFillColor: [0, 0, 255, 180],
        pickable: true,
        radiusScale: 5,
        radiusMinPixels: 1,
        radiusMaxPixels: 50,
      }),
    ];

    console.log(`[InteractionVisualization Layers] heatmapType: ${densityScoringType}, Ligands: ${data.ligand?.length || 0}, Receptors: ${data.receptor?.length || 0}`);

    let aggregationLayer: ScreenGridLayer<Point> | HeatmapLayer<Point> | HexagonLayer<Point> | null = null;
    const gridCellSizePixels = 20; 
    const hexagonRadius = 15;     
    const hexagonCoverage = 0.9;  
    const heatmapRadiusPixels = 40;
    const heatmapIntensity = 1;   
    const heatmapThreshold = 0.05;

    let layerData: Point[] | null = null;
    let layerIdPrefix = '';
    if (densityScoringType === 'Ligand') {
        layerData = data.ligand;
        layerIdPrefix = 'ligand';
    } else if (densityScoringType === 'Receptor') {
        layerData = data.receptor;
        layerIdPrefix = 'receptor';
    }

    if (layerData && layerData.length > 0) { 
      const commonProps = {
          data: layerData,
          getPosition: (d: Point) => [d.x, d.y] as [number, number],
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
          aggregationLayer = new ScreenGridLayer<Point>({ 
            ...commonProps,
            id: `${layerIdPrefix}-screengrid-layer`, 
            cellSizePixels: gridCellSizePixels,
            colorRange: colorRange as Color[],
            gpuAggregation: false,
            aggregation: 'SUM' 
          });
      } else if (aggregationLayerType === 'Hexagon') {
          console.log(`[InteractionVisualization Layers] Creating ${layerIdPrefix} HexagonLayer`);
          aggregationLayer = new HexagonLayer<Point>({
              ...commonProps,
              id: `${layerIdPrefix}-hexagon-layer`,
              radius: hexagonRadius,
              coverage: hexagonCoverage,
              colorRange: colorRange as Color[],
          });
      } else {
          console.log(`[InteractionVisualization Layers] Creating ${layerIdPrefix} HeatmapLayer`);
          aggregationLayer = new HeatmapLayer<Point>({
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

    const scatterLayersToShow = [];
    if (pointsDisplayType === 'ligands' || pointsDisplayType === 'both') {
        if (data.ligand.length > 0) scatterLayersToShow.push(baseLayers[0]);
    }
    if (pointsDisplayType === 'receptors' || pointsDisplayType === 'both') {
        if (data.receptor.length > 0) scatterLayersToShow.push(baseLayers[1]);
    }

    const finalLayers = [];
    if (!isLoading) {
        if (aggregationLayer) {
          finalLayers.push(aggregationLayer);
        }
        finalLayers.push(...scatterLayersToShow);
    }

    console.log('[InteractionVisualization Layers] Final layers array:', finalLayers.map(l => l?.id));

    return finalLayers;

  }, [data.ligand, data.receptor, densityScoringType, pointsDisplayType, aggregationLayerType, isLoading]);

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  const onDeckError = useCallback((error: Error) => console.error('[DeckGL Error]', error), []);

  return (
    <div className={`${styles.interactionVisualizationContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
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
          views={new OrthographicView({id: 'ortho-view'})}
          viewState={viewState}
          onViewStateChange={onViewStateChange}
          controller={true}
          onError={onDeckError}
          layers={layers}
          getTooltip={({object}) => object && `Point: (${object.x.toFixed(2)}, ${object.y.toFixed(2)})`}
          style={{ width: '100%', height: '100%', position: 'relative' }}
        />
        {!isLoading && (
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
                          <span>{receptorName}</span>
                      </div>
                  )}
                  <div className={styles.controlGroup}>
                      <span className={styles.controlLabel}>Density Scoring:</span>
                      <label>
                          <input type="radio" name="densityScoring" value="Off" checked={densityScoringType === 'Off'} onChange={() => setDensityScoringType('Off')} />
                          Off
                      </label>
                      <label>
                          <input type="radio" name="densityScoring" value="Ligand" checked={densityScoringType === 'Ligand'} onChange={() => setDensityScoringType('Ligand')} disabled={data.ligand.length === 0} />
                          Ligand Density
                      </label>
                      <label>
                          <input type="radio" name="densityScoring" value="Receptor" checked={densityScoringType === 'Receptor'} onChange={() => setDensityScoringType('Receptor')} disabled={data.receptor.length === 0} />
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
                              <input type="radio" name="pointsDisplay" value="ligands" checked={pointsDisplayType === 'ligands'} onChange={() => setPointsDisplayType('ligands')} disabled={data.ligand.length === 0} />
                              Ligands Only
                          </label>
                           <label>
                              <input type="radio" name="pointsDisplay" value="receptors" checked={pointsDisplayType === 'receptors'} onChange={() => setPointsDisplayType('receptors')} disabled={data.receptor.length === 0} />
                              Receptors Only
                          </label>
                           <label>
                              <input type="radio" name="pointsDisplay" value="both" checked={pointsDisplayType === 'both'} onChange={() => setPointsDisplayType('both')} disabled={data.ligand.length === 0 || data.receptor.length === 0}/>
                              Show Both
                          </label>
                  </div>
              </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default InteractionVisualization; 