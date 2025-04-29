import React, { useState, useEffect, useRef } from 'react';
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
  warnings?: string[]; // Optional warnings array
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
  const [vizWarnings, setVizWarnings] = useState<string[]>([]); // State for warnings
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  // Ref to hold the AbortController for the current fetch request
  const abortControllerRef = useRef<AbortController | null>(null);

  const fetchVisualizationData = async (ligand: string, receptor: string) => {
    // Cancel any previous ongoing request before starting a new one
    if (abortControllerRef.current) {
        abortControllerRef.current.abort();
        console.log("[DisplayPanel] Aborted previous fetch request.");
    }
    
    // Create a new AbortController for this request
    const controller = new AbortController();
    abortControllerRef.current = controller;

    setLoading(true);
    setError(null);
    setVisualizationData(null); 
    setVizWarnings([]);

    if (!jobId || !selectedLayer) {
      setError('Missing job ID or selected layer for visualization.');
      setLoading(false);
      abortControllerRef.current = null; // Clear controller if setup fails
      return;
    }

    try {
      const url = `http://localhost:8000/api/visualization/${jobId}?ligand=${encodeURIComponent(ligand)}&receptor=${encodeURIComponent(receptor)}&layer=${encodeURIComponent(selectedLayer)}`;
      console.log(`[DisplayPanel] Fetching: ${url}`);
      const response = await fetch(url, { 
          signal: controller.signal // Pass the signal here
      }); 

      // Check if the request was aborted after starting the fetch but before response
      if (controller.signal.aborted) {
          console.log("[DisplayPanel] Fetch aborted after starting.");
          // No need to set state here, finally block will handle loading
          return;
      }

      if (!response.ok) {
        const errMsg = await response.text();
        console.error(`[DisplayPanel] Fetch error ${response.status}:`, errMsg);
        throw new Error(errMsg || 'Failed to fetch visualization data');
      }

      const data: VisualizationData = await response.json();
      console.log(`[DisplayPanel] Fetch success:`, data);
      setVisualizationData(data);
      
      if (data.warnings && data.warnings.length > 0) {
          setVizWarnings(data.warnings);
      }

    } catch (err) {
      if ((err as Error).name === 'AbortError') {
          console.log("[DisplayPanel] Fetch request cancelled by user.");
          // Don't set a generic error message for cancellations
          setError('Request cancelled.'); // Optional: show a specific cancelled message
      } else {
          const message = err instanceof Error ? err.message : 'An unknown error occurred';
          console.error(`[DisplayPanel] Catch block error:`, message);
          setError(message);
          setVizWarnings([]); // Clear warnings on fetch error too
      }
      setVisualizationData(null); // Clear potentially partial data on error/abort
    } finally {
      // Ensure the controller associated with *this* request is cleared
      // Check if the current controller in the ref is the one we just used
      if (abortControllerRef.current === controller) {
          abortControllerRef.current = null;
      }
      setLoading(false);
    }
  };

  const handleCancelRequest = () => {
      if (abortControllerRef.current) {
          abortControllerRef.current.abort();
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
                    {/* Display Warnings */} 
                    {vizWarnings.length > 0 && (
                        <div className={styles.warningsContainer}>
                            <h4>Note:</h4>
                            <ul>
                                {vizWarnings.map((warning, index) => (
                                    <li key={index}>{warning}</li>
                                ))}
                            </ul>
                        </div>
                    )}

                    {/* Loading indicator with Cancel button */} 
                    {loading && (
                        <div className={styles.loadingContainer}> 
                            <p className={styles.loading}>Loading visualization...</p>
                            <button onClick={handleCancelRequest} className={styles.cancelButton}>Cancel</button>
                        </div>
                    )}
                    {error && <p className={styles.error}>{error}</p>}
                    
                    {/* Render visualization only if actual point data exists */} 
                    {visualizationData && (visualizationData.ligand.length > 0 || visualizationData.receptor.length > 0) && selectedInteraction && (
                        <Visualization
                            data={visualizationData}
                            ligandName={selectedInteraction.ligand}
                            receptorName={selectedInteraction.receptor}
                        />
                    )}
                    
                    {/* Improved Placeholder/Empty State Messages */}
                    {!loading && !error && !visualizationData && selectedInteraction && (
                         <p className={styles.noData}>Waiting for visualization data...</p> // Initial state before fetch finishes
                    )}
                    {!loading && !error && visualizationData && visualizationData.ligand.length === 0 && visualizationData.receptor.length === 0 && selectedInteraction && vizWarnings.length === 0 && (
                        <p className={styles.noData}>Visualization data loaded but is empty (no warnings received).</p> // Explicit empty state
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