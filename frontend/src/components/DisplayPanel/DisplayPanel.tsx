import React from 'react';
import styles from './DisplayPanel.module.css';
import SummaryStats from '../SummaryStats/SummaryStats';
import InteractionTable from '../InteractionTable/InteractionTable';

// Match prop types from ResultsPage
type ScopeType = 'whole_tissue' | 'layers';
type AnalysisType = 'summary' | 'pathway_dominance' | 'module_context';

interface DisplayPanelProps {
  analysisType: AnalysisType;
  data: any; // Data will vary based on analysis type
  scope: ScopeType;
  layers: string[]; // Selected layers if scope is 'layers'
}

const DisplayPanel: React.FC<DisplayPanelProps> = ({
  analysisType,
  data,
  scope,
  layers
}) => {

  const renderContent = () => {
    if (!data) {
      return <p className={styles.noData}>No data available for the current selection.</p>;
    }

    switch (analysisType) {
      case 'summary':
        return (
          <div>
            <h3 className={styles.subTitle}>Summary Statistics</h3>
            <SummaryStats data={data} />
          </div>
        );
      case 'pathway_dominance':
        return (
          <div>
            <h3 className={styles.subTitle}>Pathway Dominance</h3>
            <InteractionTable data={data} />
          </div>
        );
      case 'module_context':
        return (
          <div>
            <h3 className={styles.subTitle}>Module Context</h3>
            <InteractionTable data={data} />
          </div>
        );
      default:
        return <p className={styles.error}>Invalid analysis type selected.</p>;
    }
  };

  return (
    <div className={styles.displayPanel}>
       <h2 className={styles.title}>Results Display</h2>
       <p className={styles.contextInfo}>
         Showing results for: 
         <strong>{scope === 'whole_tissue' ? 'Whole Tissue' : `Layer(s): ${layers.join(', ') || 'None Selected'}`}</strong> - 
         <strong>{analysisType.replace('_', ' ').replace(/\b\w/g, l => l.toUpperCase())}</strong>
       </p>
       <div className={styles.contentArea}>
            {renderContent()}
       </div>
    </div>
  );
};

export default DisplayPanel; 