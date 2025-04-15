import React from 'react';
import styles from './AnalysisTypeSelector.module.css';

type AnalysisType = 'summary' | 'pathway_dominance' | 'module_context';

interface AnalysisTypeSelectorProps {
  selectedAnalysisType: AnalysisType;
  onAnalysisTypeChange: (type: AnalysisType) => void;
}

// Helper to format the analysis type names for display
const formatAnalysisTypeName = (type: AnalysisType): string => {
    return type.replace('_', ' ').replace(/\b\w/g, l => l.toUpperCase());
};

const AnalysisTypeSelector: React.FC<AnalysisTypeSelectorProps> = ({
  selectedAnalysisType,
  onAnalysisTypeChange,
}) => {
  const analysisTypes: AnalysisType[] = ['summary', 'pathway_dominance', 'module_context'];

  return (
    <div className={styles.analysisTypeSelector}>
      <h3 className={styles.title}>Analysis Type</h3>
      <div className={styles.options}>
        {analysisTypes.map(type => (
            <label key={type} className={styles.optionLabel}>
                <input 
                    type="radio"
                    name="analysisType"
                    value={type}
                    checked={selectedAnalysisType === type}
                    onChange={() => onAnalysisTypeChange(type)}
                    className={styles.radioInput}
                />
                {formatAnalysisTypeName(type)}
            </label>
        ))}
      </div>
    </div>
  );
};

export default AnalysisTypeSelector; 