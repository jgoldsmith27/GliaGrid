import React from 'react';
import styles from './LayerSelector.module.css';

interface LayerSelectorProps {
  availableLayers: string[];
  selectedLayers: string[];
  onLayersChange: (layers: string[]) => void;
}

const LayerSelector: React.FC<LayerSelectorProps> = ({
  availableLayers,
  selectedLayers,
  onLayersChange,
}) => {
  
  const handleRadioChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const layerName = event.target.value;
    // For radio buttons, we only need to set the selected layer
    onLayersChange([layerName]);
  };

  if (availableLayers.length === 0) {
    return (
        <div className={styles.layerSelector}>
            <h3 className={styles.title}>Layers</h3>
            <p className={styles.noLayers}>No layers found in results.</p>
        </div>
    );
  }

  return (
    <div className={styles.layerSelector}>
      <h3 className={styles.title}>Select Layer</h3>
      <div className={styles.optionsList}>
        {availableLayers.map(layer => (
          <label key={layer} className={styles.optionLabel}>
            <input
              type="radio"
              name="layer"
              value={layer}
              checked={selectedLayers.includes(layer)}
              onChange={handleRadioChange}
              className={styles.radioInput}
            />
            {layer}
          </label>
        ))}
      </div>
    </div>
  );
};

export default LayerSelector; 