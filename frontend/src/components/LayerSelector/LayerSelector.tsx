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
  
  const handleCheckboxChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const layerName = event.target.value;
    const isChecked = event.target.checked;

    if (isChecked) {
      // Add layer to selection
      onLayersChange([...selectedLayers, layerName]);
    } else {
      // Remove layer from selection
      onLayersChange(selectedLayers.filter(layer => layer !== layerName));
    }
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
      <h3 className={styles.title}>Select Layer(s)</h3>
      <div className={styles.optionsList}>
        {availableLayers.map(layer => (
          <label key={layer} className={styles.optionLabel}>
            <input
              type="checkbox"
              value={layer}
              checked={selectedLayers.includes(layer)}
              onChange={handleCheckboxChange}
              className={styles.checkboxInput}
            />
            {layer}
          </label>
        ))}
      </div>
      {/* Optional: Add Select All / Deselect All buttons later */}
    </div>
  );
};

export default LayerSelector; 