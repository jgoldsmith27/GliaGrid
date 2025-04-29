import React, { useState } from 'react';
import { ScatterplotLayer } from '@deck.gl/layers';
import DeckGL from '@deck.gl/react';
import styles from './Visualization.module.css';

interface VisualizationProps {
  data: {
    ligand: {
      x: number;
      y: number;
    }[];
    receptor: {
      x: number;
      y: number;
    }[];
  };
  ligandName: string;
  receptorName: string;
}

const Visualization: React.FC<VisualizationProps> = ({ data, ligandName, receptorName }) => {
  const [isFullscreen, setIsFullscreen] = useState(false);

  const layers = [
    new ScatterplotLayer({
      id: 'ligand-layer',
      data: data.ligand,
      getPosition: (d: any) => [d.x, d.y],
      getRadius: 5,
      getFillColor: [255, 0, 0], // Red for ligands
      pickable: true,
    }),
    new ScatterplotLayer({
      id: 'receptor-layer',
      data: data.receptor,
      getPosition: (d: any) => [d.x, d.y],
      getRadius: 5,
      getFillColor: [0, 0, 255], // Blue for receptors
      pickable: true,
    }),
  ];

  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  return (
    <div className={`${styles.visualizationContainer} ${isFullscreen ? styles.fullscreen : ''}`}>
      <div className={styles.controls}>
        <button 
          className={styles.fullscreenButton}
          onClick={toggleFullscreen}
          title={isFullscreen ? "Exit fullscreen" : "Enter fullscreen"}
        >
          {isFullscreen ? "✕" : "⛶"}
        </button>
      </div>
      <DeckGL
        initialViewState={{
          longitude: 0,
          latitude: 0,
          zoom: 1,
        }}
        controller={true}
        layers={layers}
      />
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
  );
};

export default Visualization; 