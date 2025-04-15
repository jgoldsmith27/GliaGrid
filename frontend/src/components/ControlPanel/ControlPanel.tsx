import React from 'react';
import styles from './ControlPanel.module.css';
import ScopeSelector from '../ScopeSelector/ScopeSelector'; // Placeholder
import LayerSelector from '../LayerSelector/LayerSelector'; // Placeholder
import AnalysisTypeSelector from '../AnalysisTypeSelector/AnalysisTypeSelector'; // Placeholder

// Match prop types from ResultsPage
type ScopeType = 'whole_tissue' | 'layers';
type AnalysisType = 'summary' | 'pathway_dominance' | 'module_context';

interface ControlPanelProps {
  selectedScope: ScopeType;
  onScopeChange: (scope: ScopeType) => void;
  availableLayers: string[];
  selectedLayers: string[];
  onLayersChange: (layers: string[]) => void;
  selectedAnalysisType: AnalysisType;
  onAnalysisTypeChange: (type: AnalysisType) => void;
  // TODO: Add filter props
}

const ControlPanel: React.FC<ControlPanelProps> = ({
  selectedScope,
  onScopeChange,
  availableLayers,
  selectedLayers,
  onLayersChange,
  selectedAnalysisType,
  onAnalysisTypeChange
}) => {
  return (
    <div className={styles.controlPanel}>
      <h2 className={styles.title}>Controls</h2>
      
      <div className={styles.controlGroup}>
        <ScopeSelector 
            selectedScope={selectedScope} 
            onScopeChange={onScopeChange} 
        />
      </div>

      {selectedScope === 'layers' && (
        <div className={styles.controlGroup}>
          <LayerSelector 
              availableLayers={availableLayers} 
              selectedLayers={selectedLayers} 
              onLayersChange={onLayersChange} 
          />
        </div>
      )}

      <div className={styles.controlGroup}>
        <AnalysisTypeSelector 
            selectedAnalysisType={selectedAnalysisType} 
            onAnalysisTypeChange={onAnalysisTypeChange} 
        />
      </div>

      {/* Placeholder for Filters */}
      <div className={styles.controlGroup}>
          <h3 className={styles.subTitle}>Filters (Coming Soon)</h3>
          {/* TODO: Add filter components based on selectedAnalysisType */}
      </div>

      {/* Placeholder for Comparison Controls */}
      <div className={styles.controlGroup}>
          <h3 className={styles.subTitle}>Comparison (Coming Soon)</h3>
          {/* TODO: Add comparison toggle/controls */}
      </div>

    </div>
  );
};

export default ControlPanel; 