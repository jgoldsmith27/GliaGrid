import React from 'react';
import ToggleButton from '@mui/material/ToggleButton';
import ToggleButtonGroup from '@mui/material/ToggleButtonGroup';
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
  const handleChange = (
    event: React.MouseEvent<HTMLElement>,
    newScope: ScopeType | null, // ToggleButtonGroup can return null if all toggles are deselected
  ) => {
    if (newScope !== null) {
      onScopeChange(newScope);
    }
  };

  return (
    <div className={styles.scopeSelectorContainer}> {/* Renamed for clarity if needed */}
      <h3 className={styles.title}>SCOPE</h3> {/* Kept title, can be removed if ToggleButtonGroup is self-explanatory */}
      <ToggleButtonGroup
        value={selectedScope}
        exclusive // Ensures only one button can be active at a time
        onChange={handleChange}
        aria-label="Analysis Scope"
        className={styles.toggleButtonGroup} // Added class for styling
      >
        <ToggleButton value="whole_tissue" aria-label="Whole Tissue" className={styles.toggleButton}>
          Whole Tissue
        </ToggleButton>
        <ToggleButton value="layers" aria-label="Layers" className={styles.toggleButton}>
          Layers
        </ToggleButton>
        <ToggleButton value="custom" aria-label="Custom Selection" className={styles.toggleButton}>
          Custom Selection
        </ToggleButton>
      </ToggleButtonGroup>
    </div>
  );
};

export default ScopeSelector; 