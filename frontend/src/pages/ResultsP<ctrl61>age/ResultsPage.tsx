import React, { useState } from 'react';
import { TextField } from '@mui/material';
import styles from './ResultsPage.module.css';

const ResultsPage: React.FC = () => {
  const [ligandSearchTerm, setLigandSearchTerm] = useState('');
  const [receptorSearchTerm, setReceptorSearchTerm] = useState('');

  return (
    <>
      <div className={styles.controlGroup}>
        <TextField
          label="Search Ligand"
          variant="outlined"
          fullWidth
          value={ligandSearchTerm}
          onChange={(e) => setLigandSearchTerm(e.target.value)}
          sx={{
            backgroundColor: 'rgba(255, 0, 0, 0.1)', // Light red background for testing
            '.MuiInputBase-root': { 
              height: '36px',
              minHeight: '36px',
            },
            '.MuiOutlinedInput-root': { 
                height: '36px',
                minHeight: '36px', 
            },
            '.MuiInputBase-input': { // Directly overriding the problematic styles
              padding: '6px 8px !important', 
              height: 'auto !important',      // Let padding dictate internal height, override em-based height
              fontSize: '0.875rem',          // Keep font size adjustment
            },
            '.MuiInputLabel-root': {
              fontSize: '0.875rem',
            },
            '.MuiInputLabel-outlined.MuiInputLabel-shrink': {
              transform: 'translate(14px, -5px) scale(0.75)',
            },
          }}
        />
      </div>
      <div className={styles.controlGroup}>
        <TextField
          label="Search Receptor"
          variant="outlined"
          fullWidth
          value={receptorSearchTerm}
          onChange={(e) => setReceptorSearchTerm(e.target.value)}
          sx={{
            backgroundColor: 'rgba(255, 0, 0, 0.1)', // Light red background for testing
            '.MuiInputBase-root': {
              height: '36px',
              minHeight: '36px',
            },
            '.MuiOutlinedInput-root': {
                height: '36px',
                minHeight: '36px',
            },
            '.MuiInputBase-input': { // Directly overriding the problematic styles
              padding: '6px 8px !important',
              height: 'auto !important',
              fontSize: '0.875rem',
            },
            '.MuiInputLabel-root': {
              fontSize: '0.875rem',
            },
            '.MuiInputLabel-outlined.MuiInputLabel-shrink': {
              transform: 'translate(14px, -5px) scale(0.75)',
            },
          }}
        />
      </div>
    </>
  );
};

export default ResultsPage; 