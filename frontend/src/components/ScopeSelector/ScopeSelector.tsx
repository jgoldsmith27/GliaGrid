import React from 'react';
import styles from './ScopeSelector.module.css';

// Export the type
export type ScopeType = 'whole_tissue' | 'layers' | 'custom';

interface ScopeSelectorProps {
  selectedScope: ScopeType;
  onScopeChange: (scope: ScopeType) => void;
}

const ScopeSelector: React.FC<ScopeSelectorProps> = ({
  selectedScope,
  onScopeChange,
}) => {
  return (
    <div className={styles.scopeSelector}>
      <h3 className={styles.title}>Scope</h3>
      <div className={styles.options}>
        <label className={styles.optionLabel}>
          <input
            type="radio"
            name="scope"
            value="whole_tissue"
            checked={selectedScope === 'whole_tissue'}
            onChange={() => onScopeChange('whole_tissue')}
            className={styles.radioInput}
          />
          Whole Tissue
        </label>
        <label className={styles.optionLabel}>
          <input
            type="radio"
            name="scope"
            value="layers"
            checked={selectedScope === 'layers'}
            onChange={() => onScopeChange('layers')}
            className={styles.radioInput}
          />
          Layers
        </label>
        <label className={styles.optionLabel}>
          <input
            type="radio"
            name="scope"
            value="custom"
            checked={selectedScope === 'custom'}
            onChange={() => onScopeChange('custom')}
            className={styles.radioInput}
          />
          Custom Selection
        </label>
      </div>
    </div>
  );
};

export default ScopeSelector; 