import React, { useState, useEffect } from 'react';
import styles from './DisplayPanel.module.css';
import SummaryStats from '../SummaryStats/SummaryStats';
import PathwayDominanceTable from '../PathwayDominanceTable/PathwayDominanceTable';
import ModuleContextTable from '../ModuleContextTable/ModuleContextTable';
import { PathwayDominanceResult, ModuleContextResult, SummaryStatsData } from '../../types/analysisResults';
import Visualization from '../Visualization/Visualization';

// Match prop types from ResultsPage
type ScopeType = 'whole_tissue' | 'layers';
type AnalysisType = 'summary' | 'pathway_dominance' | 'module_context';

interface DisplayPanelProps {
  analysisType: AnalysisType;
  data: any; // Data will vary based on analysis type
  scope: ScopeType;
  layers: string[]; // Selected layers if scope is 'layers'
}

interface VisualizationData {
  ligand: {
    x: number;
    y: number;
  }[];
  receptor: {
    x: number;
    y: number;
  }[];
}

const DisplayPanel: React.FC<DisplayPanelProps> = ({
  analysisType,
  data,
  scope,
  layers
}) => {
  const [selectedInteraction, setSelectedInteraction] = useState<PathwayDominanceResult | null>(null);
  const [visualizationData, setVisualizationData] = useState<VisualizationData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const fetchVisualizationData = async (ligand: string, receptor: string) => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch(`http://localhost:8000/api/visualization/${ligand}/${receptor}`);
      if (!response.ok) {
        throw new Error('Failed to fetch visualization data');
      }
      const data = await response.json();
      setVisualizationData(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An unknown error occurred');
    } finally {
      setLoading(false);
    }
  };

  const handleRowClick = (row: PathwayDominanceResult) => {
    setSelectedInteraction(row);
    fetchVisualizationData(row.ligand, row.receptor);
  };

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
            <PathwayDominanceTable 
              data={data as PathwayDominanceResult[]} 
              onRowClick={handleRowClick}
            />
            {loading && <p className={styles.loading}>Loading visualization...</p>}
            {error && <p className={styles.error}>{error}</p>}
            {visualizationData && selectedInteraction && (
              <div className={styles.visualizationContainer}>
                <Visualization
                  data={visualizationData}
                  ligandName={selectedInteraction.ligand}
                  receptorName={selectedInteraction.receptor}
                />
              </div>
            )}
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