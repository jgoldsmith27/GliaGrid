import React from 'react';
import styles from './DisplayPanel.module.css';
import SummaryStats from '../SummaryStats/SummaryStats';
import PathwayDominanceTable from '../PathwayDominanceTable/PathwayDominanceTable';
import ModuleContextTable from '../ModuleContextTable/ModuleContextTable';
import { PathwayDominanceResult, ModuleContextResult, SummaryStatsData } from '../../types/analysisResults';

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
    if (data === null || data === undefined) {
      return <p className={styles.noData}>No data available for the current selection.</p>;
    }

    switch (analysisType) {
      case 'summary':
        if (typeof data !== 'object' || data === null || !('unique_ligands' in data && 'unique_receptors' in data)) {
            return <p className={styles.error}>Invalid data format for Summary Stats.</p>;
        }
        return (
          <div>
            <h3 className={styles.subTitle}>Summary Statistics</h3>
            <SummaryStats data={data as SummaryStatsData} />
          </div>
        );
      case 'pathway_dominance':
        if (!Array.isArray(data)) {
             return <p className={styles.error}>Invalid data format for Pathway Dominance.</p>;
         }
        return (
          <div>
            <h3 className={styles.subTitle}>Pathway Dominance</h3>
            <PathwayDominanceTable data={data as PathwayDominanceResult[]} />
          </div>
        );
      case 'module_context':
        if (!Array.isArray(data)) {
             return <p className={styles.error}>Invalid data format for Module Context.</p>;
         }
        return (
          <div>
            <h3 className={styles.subTitle}>Module Context</h3>
            <ModuleContextTable data={data as ModuleContextResult[]} />
          </div>
        );
      default:
        const _exhaustiveCheck: never = analysisType;
        console.warn(`Unhandled analysis type: ${_exhaustiveCheck}`);
        return <p className={styles.error}>Invalid analysis type selected.</p>;
    }
  };

  const capitalize = (s: string) => s.charAt(0).toUpperCase() + s.slice(1);
  const formatAnalysisType = (type: AnalysisType) => type.split('_').map(capitalize).join(' ');

  return (
    <div className={styles.displayPanel}>
       <h2 className={styles.title}>Results Display</h2>
       <p className={styles.contextInfo}>
         Showing results for: 
         <strong>{scope === 'whole_tissue' ? 'Whole Tissue' : `Layer(s): ${layers.join(', ') || 'None Selected'}`}</strong> - 
         <strong>{formatAnalysisType(analysisType)}</strong>
       </p>
       <div className={styles.contentArea}>
            {renderContent()}
       </div>
    </div>
  );
};

export default DisplayPanel; 