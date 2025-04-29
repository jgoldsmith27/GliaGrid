import React, { useState, useEffect, useRef, useCallback } from 'react';
import styles from './SummaryTabContent.module.css'; // Create CSS module later
import AnalysisTable, { ColumnDefinition } from '../AnalysisTable/AnalysisTable';
import Visualization from '../Visualization/Visualization';
import { CombinedInteractionData } from '../../pages/ResultsPage/ResultsPage'; // Import shared type

// Copied from old DisplayPanel for fetch response
interface VisualizationData {
  ligand: { x: number; y: number }[];
  receptor: { x: number; y: number }[];
  warnings?: string[];
}

interface SummaryTabContentProps {
  jobId: string;
  combinedData: CombinedInteractionData[];
  selectedPair: [string, string] | null;
  onSelectPair: (pair: [string, string]) => void;
  selectedLayerName: string; // Name of layer/scope for API call
}

// Define columns for the combined table
const combinedColumns: ColumnDefinition<CombinedInteractionData>[] = [
    { key: 'ligand', header: 'Ligand' },
    { key: 'receptor', header: 'Receptor' },
    { key: 'score', header: 'Score (Pathway)', format: (v) => typeof v === 'number' ? v.toExponential(3) : 'N/A' },
    { key: 'ligand_norm_expr', header: 'Ligand Norm Expr', format: (v) => typeof v === 'number' ? v.toFixed(4) : 'N/A' },
    { key: 'receptor_avg_norm_expr', header: 'Receptor Norm Expr', format: (v) => typeof v === 'number' ? v.toFixed(4) : 'N/A' },
    { key: 'interaction_type', header: 'Interaction Type' },
    { key: 'ligand_module', header: 'Ligand Module' },
    { key: 'receptor_modules', header: 'Receptor Modules' }, // Array is handled by default formatter
    { key: 'is_same_module', header: 'Same Module?' }, // Boolean is handled by default formatter
];

const SummaryTabContent: React.FC<SummaryTabContentProps> = ({
  jobId,
  combinedData,
  selectedPair,
  onSelectPair,
  selectedLayerName,
}) => {
  const [visualizationData, setVisualizationData] = useState<VisualizationData | null>(null);
  const [vizWarnings, setVizWarnings] = useState<string[]>([]);
  const [loadingViz, setLoadingViz] = useState(false);
  const [vizError, setVizError] = useState<string | null>(null);
  const [selectedRowIndex, setSelectedRowIndex] = useState<number | null>(null); // For highlighting table row
  
  const abortControllerRef = useRef<AbortController | null>(null);

  const fetchVisualizationData = useCallback(async (ligand: string, receptor: string) => {
    if (abortControllerRef.current) {
        abortControllerRef.current.abort();
    }
    const controller = new AbortController();
    abortControllerRef.current = controller;

    setLoadingViz(true);
    setVizError(null);
    setVisualizationData(null); 
    setVizWarnings([]);

    if (!jobId || !selectedLayerName) {
      setVizError('Missing job ID or selected layer for visualization.');
      setLoadingViz(false);
      abortControllerRef.current = null;
      return;
    }

    try {
      const url = `http://localhost:8000/api/visualization/${jobId}?ligand=${encodeURIComponent(ligand)}&receptor=${encodeURIComponent(receptor)}&layer=${encodeURIComponent(selectedLayerName)}`;
      const response = await fetch(url, { signal: controller.signal }); 

      if (controller.signal.aborted) { return; }

      if (!response.ok) {
        const errMsg = await response.text();
        throw new Error(errMsg || 'Failed to fetch visualization data');
      }

      const data: VisualizationData = await response.json();
      setVisualizationData(data);
      if (data.warnings) setVizWarnings(data.warnings);

    } catch (err) {
      if ((err as Error).name === 'AbortError') {
          setVizError('Request cancelled.'); 
      } else {
          const message = err instanceof Error ? err.message : 'An unknown error occurred';
          setVizError(message);
          setVizWarnings([]);
      }
      setVisualizationData(null);
    } finally {
      if (abortControllerRef.current === controller) {
          abortControllerRef.current = null;
      }
      setLoadingViz(false);
    }
  }, [jobId, selectedLayerName]); // Dependencies for fetch

  // Fetch data when selectedPair changes
  useEffect(() => {
      if (selectedPair) {
          fetchVisualizationData(selectedPair[0], selectedPair[1]);
      } else {
          // Clear visualization if no pair is selected
          setVisualizationData(null);
          setVizWarnings([]);
          setVizError(null);
          setLoadingViz(false); 
          if (abortControllerRef.current) { // Cancel ongoing fetch if selection is cleared
              abortControllerRef.current.abort();
          }
      }
      // Cleanup function for aborting on unmount or before next fetch
      return () => {
         if (abortControllerRef.current) {
              abortControllerRef.current.abort();
          }
      };
  }, [selectedPair, fetchVisualizationData]);

  const handleTableRowClick = useCallback((row: CombinedInteractionData, index: number) => {
      setSelectedRowIndex(index);
      onSelectPair([row.ligand, row.receptor]); // Call parent handler to update selectedPair
      // Fetch logic is now handled by the useEffect hook above
  }, [onSelectPair]);

  const handleCancelVizRequest = useCallback(() => {
      if (abortControllerRef.current) {
          abortControllerRef.current.abort();
      }
  }, []);

  return (
    <div className={styles.summaryLayout}> {/* Use Flexbox/Grid */} 
        <div className={styles.tableArea}>
            <h3>Interaction Details</h3>
            <AnalysisTable 
                data={combinedData} 
                columns={combinedColumns} // Pass the defined columns
                onRowClick={handleTableRowClick} 
                selectedRowIndex={selectedRowIndex}
                loading={loadingViz} // Pass loading state to overlay table
            />
        </div>
        <div className={styles.visualizationArea}>
            <h3>Spatial Visualization</h3>
             {vizWarnings.length > 0 && (
                <div className={styles.warningsContainer}>
                    <h4>Note:</h4>
                    <ul>{vizWarnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
                </div>
            )}
            {loadingViz && (
                <div className={styles.loadingContainer}>
                    <p>Loading visualization...</p>
                    <button onClick={handleCancelVizRequest} className={styles.cancelButton}>Cancel</button>
                </div>
            )}
            {vizError && <p className={styles.error}>{vizError}</p>}
            
            {visualizationData && selectedPair ? (
                 <Visualization
                    data={visualizationData}
                    ligandName={selectedPair[0]}
                    receptorName={selectedPair[1]}
                    currentScope={selectedLayerName}
                />
            ) : (
                !loadingViz && <p>Select an interaction from the table to visualize.</p>
            )}
        </div>
    </div>
  );
};

export default SummaryTabContent; 