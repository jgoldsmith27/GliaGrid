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

// Type guard to check for common interaction properties
interface BasicInteraction {
  ligand: string;
  receptor: string;
}
function isBasicInteraction(obj: any): obj is BasicInteraction {
    return obj && typeof obj.ligand === 'string' && typeof obj.receptor === 'string';
}

interface DisplayPanelProps {
  analysisType: AnalysisType;
  data: any; // Data will vary based on analysis type
  scope: ScopeType;
  layers: string[]; // Selected layers if scope is 'layers'
  jobId: string;
  selectedLayer: string;
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
  layers,
  jobId,
  selectedLayer
}) => {
  // Lifted state for selected interaction details and index
  const [selectedInteraction, setSelectedInteraction] = useState<BasicInteraction | null>(null);
  const [selectedRowIndex, setSelectedRowIndex] = useState<number | null>(null);
  const [visualizationData, setVisualizationData] = useState<VisualizationData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const fetchVisualizationData = async (ligand: string, receptor: string) => {
    setLoading(true);
    setError(null);
    setVisualizationData(null); // Clear previous visualization
    if (!jobId || !selectedLayer) {
      setError('Missing job ID or selected layer for visualization.');
      setLoading(false);
      return;
    }
    try {
      const url = `http://localhost:8000/api/visualization/${jobId}?ligand=${encodeURIComponent(ligand)}&receptor=${encodeURIComponent(receptor)}&layer=${encodeURIComponent(selectedLayer)}`;
      console.log(`[DisplayPanel] Fetching: ${url}`); // Add logging
      const response = await fetch(url);
      if (!response.ok) {
        const errMsg = await response.text();
        console.error(`[DisplayPanel] Fetch error ${response.status}:`, errMsg); // Log error details
        throw new Error(errMsg || 'Failed to fetch visualization data');
      }
      const data = await response.json();
      console.log(`[DisplayPanel] Fetch success:`, data); // Log success data
      setVisualizationData(data);
    } catch (err) {
      const message = err instanceof Error ? err.message : 'An unknown error occurred';
      console.error(`[DisplayPanel] Catch block error:`, message);
      setError(message);
      setVisualizationData(null); // Clear vis on error
    } finally {
      setLoading(false);
    }
  };

  // Generic handler updated to accept index and set state
  const handleInteractionRowClick = (row: unknown, index: number) => {
      if (isBasicInteraction(row)) {
        setSelectedInteraction(row); 
        setSelectedRowIndex(index);
        fetchVisualizationData(row.ligand, row.receptor);
      } else {
          console.warn('[DisplayPanel] Clicked row does not contain ligand/receptor:', row);
          setSelectedInteraction(null);
          setVisualizationData(null);
          setSelectedRowIndex(null);
      }
  };

  const renderContent = () => {
    if (data === null || data === undefined) {
      return <p className={styles.noData}>No data available for the current selection.</p>;
    }

    // Determine if the current analysis type should show visualization
    const showVisualization = analysisType === 'pathway_dominance' || analysisType === 'module_context';

    let tableComponent = null;
    switch (analysisType) {
      case 'summary':
        if (typeof data !== 'object' || data === null || !('unique_ligands' in data && 'unique_receptors' in data)) {
            return <p className={styles.error}>Invalid data format for Summary Stats.</p>;
        }
        tableComponent = (
          <div>
            <h3 className={styles.subTitle}>Summary Statistics</h3>
            <SummaryStats data={data as SummaryStatsData} />
          </div>
        );
        break; // Add break
      case 'pathway_dominance':
        if (!Array.isArray(data)) {
             return <p className={styles.error}>Invalid data format for Pathway Dominance.</p>;
         }
        tableComponent = (
          <div>
            <h3 className={styles.subTitle}>Pathway Dominance</h3>
            <PathwayDominanceTable 
              data={data as PathwayDominanceResult[]} 
              onRowClick={handleInteractionRowClick}
              loading={loading} 
              selectedRowIndex={selectedRowIndex}
            />
          </div>
        );
        break; // Add break
      case 'module_context':
        if (!Array.isArray(data)) {
             return <p className={styles.error}>Invalid data format for Module Context.</p>;
         }
        tableComponent = (
          <div>
            <h3 className={styles.subTitle}>Module Context</h3>
            <ModuleContextTable 
              data={data as ModuleContextResult[]} 
              onRowClick={handleInteractionRowClick}
              loading={loading} 
              selectedRowIndex={selectedRowIndex}
            />
          </div>
        );
        break; // Add break
      default:
        const _exhaustiveCheck: never = analysisType;
        console.warn(`Unhandled analysis type: ${_exhaustiveCheck}`);
        return <p className={styles.error}>Invalid analysis type selected.</p>;
    }

    // Render layout with table and optional visualization side-by-side
    return (
        <div className={styles.tableAndVizLayout}> 
            <div className={styles.tableArea}>
                {tableComponent}
            </div>
            {showVisualization && (
                <div className={styles.visualizationArea}> 
                    {loading && <p className={styles.loading}>Loading visualization...</p>}
                    {error && <p className={styles.error}>{error}</p>}
                    {/* Render visualization only if data exists and an interaction was selected */}
                    {visualizationData && selectedInteraction && (
                        <Visualization
                            data={visualizationData}
                            ligandName={selectedInteraction.ligand}
                            receptorName={selectedInteraction.receptor}
                        />
                    )}
                     {/* Placeholder or message if no interaction selected */}
                    {!loading && !error && !visualizationData && selectedInteraction && (
                        <p className={styles.noData}>Visualization data loaded but is empty.</p>
                    )}
                     {!loading && !error && !selectedInteraction && (
                        <p className={styles.noData}>Click on a row in the table to visualize interaction.</p>
                    )}
                </div>
            )}
        </div>
    );

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