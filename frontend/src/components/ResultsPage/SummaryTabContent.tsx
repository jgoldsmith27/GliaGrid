import React, { useState, useEffect, useRef, useCallback } from 'react';
import styles from './SummaryTabContent.module.css'; // Create CSS module later
import AnalysisTable, { ColumnDefinition } from '../AnalysisTable/AnalysisTable';
import InteractionVisualization from '../InteractionVisualization/InteractionVisualization.tsx';
import SpatialOverviewVisualization from '../SpatialOverviewVisualization/SpatialOverviewVisualization';
import { CombinedInteractionData } from '../../pages/ResultsPage/ResultsPage'; // Import shared type
import { ScopeType } from '../ScopeSelector/ScopeSelector'; // Import ScopeType if defined there
import useInteractionData from '../../hooks/useInteractionData'; // Import the hook

// Type for the data returned by the new /api/points/{jobId}/all endpoint
interface AllPointsData {
    x: number;
    y: number;
    layer: string;
}

interface SummaryTabContentProps {
  jobId: string;
  combinedData: CombinedInteractionData[];
  selectedPair: [string, string] | null;
  onSelectPair: (pair: [string, string] | null) => void; // Allow null to clear selection
  currentScope: ScopeType; 
  apiScopeName: string | null; // ADDED: The actual scope name for the API call
  onLassoSelect?: (coords: [number, number][] | null) => void; // ADDED prop definition
}

// Define columns for the combined table
const combinedColumns: ColumnDefinition<CombinedInteractionData>[] = [
    { key: 'ligand', header: 'Ligand' },
    { key: 'receptor', header: 'Receptor' },
    { key: 'score', header: 'Score (Pathway)', format: (v) => typeof v === 'number' ? v.toExponential(3) : 'N/A' },
    { key: 'ligand_norm_expr', header: 'Ligand Norm Expr', format: (v) => typeof v === 'number' ? v.toFixed(4) : 'N/A' },
    { key: 'receptor_avg_norm_expr', header: 'Receptor Norm Expr', format: (v) => typeof v === 'number' ? v.toFixed(4) : 'N/A' },
    { 
        key: 'ligand_module', 
        header: 'Ligand Module',
        format: (v) => {
            if (v === null || v === undefined) return 'N/A';
            // Convert to string, removing trailing zeros for numbers
            return typeof v === 'number' ? v.toString() : String(v);
        }
    },
    { key: 'receptor_modules', header: 'Receptor Modules' }, // Array is handled by default formatter
];

const SummaryTabContent: React.FC<SummaryTabContentProps> = ({
  jobId,
  combinedData,
  selectedPair,
  onSelectPair,
  currentScope, 
  apiScopeName, // ADDED
  onLassoSelect, // ADDED prop destructuring
}) => {
  // State for SpatialOverviewVisualization
  const [allPointsData, setAllPointsData] = useState<AllPointsData[] | null>(null);
  const [loadingAllPoints, setLoadingAllPoints] = useState(false);
  const [allPointsError, setAllPointsError] = useState<string | null>(null);

  const [selectedRowIndex, setSelectedRowIndex] = useState<number | null>(null); // For highlighting table row
  
  const abortControllerRef = useRef<AbortController | null>(null);

  // Use the custom hook for interaction data
  const { 
    interactionVizData, 
    isLoading: isLoadingInteractionViz, // Rename hook output for consistency with existing JSX
    error: interactionVizError, 
    warnings: interactionVizWarnings,
    // fetchInteractionData, // Not directly needed unless manual trigger is implemented
    cancelFetch: cancelInteractionVizFetch // Rename for clarity
  } = useInteractionData(jobId, selectedPair, apiScopeName);

  // Fetch function for All Points data
  const fetchAllPoints = useCallback(async (currentJobId: string) => {
     if (abortControllerRef.current) {
        abortControllerRef.current.abort(); // Abort previous fetch
    }
    const controller = new AbortController();
    abortControllerRef.current = controller;

    setLoadingAllPoints(true);
    setAllPointsError(null);
    setAllPointsData(null);

    try {
        const url = `http://localhost:8000/api/points/${currentJobId}/all`;
        const response = await fetch(url, { signal: controller.signal });

        if (controller.signal.aborted) {
             console.log('All points fetch aborted');
             return;
        }

        if (!response.ok) { 
            throw new Error(await response.text() || 'Failed to fetch all points data'); 
        }
        const data: AllPointsData[] = await response.json();
        setAllPointsData(data);
    } catch (err) {
        if ((err as Error).name === 'AbortError') {
             setAllPointsError('Request cancelled.');
        } else {
            const message = err instanceof Error ? err.message : 'An unknown error occurred';
            setAllPointsError(message);
        }
        setAllPointsData(null);
    } finally {
         if (abortControllerRef.current === controller) {
             abortControllerRef.current = null;
        }
        setLoadingAllPoints(false);
    }
  }, []); // Dependency: jobId implicitly via usage

  // Effect to fetch data based on currentScope
  useEffect(() => {
    console.log(`[SummaryTabContent] useEffect triggered. Scope: ${currentScope}, Pair: ${selectedPair ? selectedPair.join('-') : 'null'}, apiScopeName: ${apiScopeName}`); 
    
    // Reset spatial overview state when dependencies change
    setAllPointsData(null);
    setAllPointsError(null);
    setLoadingAllPoints(false);
    // No need to reset interaction viz state, the hook manages it
    
    if (jobId) {
        if (currentScope === 'custom') {
             // Fetch all points data for custom scope
             fetchAllPoints(jobId);
             // Interaction viz hook will handle its state based on apiScopeName being 'custom'
        } 
        // No explicit fetch call needed for interaction viz here; the hook handles it based on its dependencies (jobId, selectedPair, apiScopeName)
    }

      // Cleanup function for aborting spatial overview fetch on unmount or before next effect run
      return () => {
         console.log("[SummaryTabContent] useEffect cleanup");
         if (abortControllerRef.current) { // Only cancel the spatial overview fetch
              console.log("[SummaryTabContent] Aborting spatial overview fetch.");
              abortControllerRef.current.abort();
              abortControllerRef.current = null;
          }
         // No need to cancel interaction viz fetch here; the hook handles its own cleanup
      };
  // Dependencies: Fetch when jobId, scope, or selected pair changes. apiScopeName is derived from scope, so it implicitly changes.
  // Removed fetchInteractionVisualizationData from dependencies
  // Keep fetchAllPoints
  }, [jobId, currentScope, apiScopeName, selectedPair, fetchAllPoints, onSelectPair]); 

  // Table row click handler - only relevant for non-custom scopes
  const handleTableRowClick = useCallback((row: CombinedInteractionData, index: number) => {
    console.log(`[SummaryTabContent] Row clicked: ${row.ligand}-${row.receptor}, Index: ${index}, Current scope: ${currentScope}`); // Log 1a
      if (currentScope !== 'custom') { // Only update pair if not in custom mode
         setSelectedRowIndex(index);
         const newPair: [string, string] = [row.ligand, row.receptor];
         console.log("[SummaryTabContent] Calling onSelectPair with:", newPair); // Log 1b
         onSelectPair(newPair); 
      }
  }, [onSelectPair, currentScope]);

  // Generic cancel handler - now only cancels the spatial overview fetch
  const handleCancelVizRequest = useCallback(() => {
      if (abortControllerRef.current) {
          console.log("[SummaryTabContent] Cancelling spatial overview fetch request.");
          abortControllerRef.current.abort();
          abortControllerRef.current = null; // Clear ref after abort
          setLoadingAllPoints(false); // Ensure loading stops
          setAllPointsError("Request cancelled by user."); // Set error state
      } else {
          console.log("[SummaryTabContent] No active spatial overview fetch request to cancel.");
      }
      // Note: To cancel the interaction viz, we would call cancelInteractionVizFetch()
      // but it needs a separate button or logic tied to the interaction viz state.
  }, []); // Removed dependency on abortControllerRef state, check directly

  // Helper to render the correct visualization component
  const renderVisualization = () => {
    if (currentScope === 'custom') {
        if (loadingAllPoints) {
            return (
                <div className={styles.loadingContainer}>
                  <p>Loading spatial overview...</p>
                  {/* This button now ONLY cancels the spatial overview fetch */}
                  <button onClick={handleCancelVizRequest} className={styles.cancelButton}>Cancel</button> 
                </div>
            );
        }
        if (allPointsError) {
            return <p className={styles.errorText}>Error loading spatial overview: {allPointsError}</p>;
        }
        if (allPointsData) {
            return <SpatialOverviewVisualization jobId={jobId} onLassoSelect={onLassoSelect} />;
        }
        return <p>No spatial overview data available.</p>; // Or initial state message

    } else { // 'whole_tissue' or 'layers'
         // Use state from the hook
        if (isLoadingInteractionViz) {
            return (
                <div className={styles.loadingContainer}>
                  <p>Loading interaction visualization...</p>
                   {/* FIX: Uncomment and enable the cancel button for interaction viz */}
                   <button onClick={cancelInteractionVizFetch} className={styles.cancelButton}>Cancel</button> 
                </div>
            );
        }
        if (interactionVizError) {
            return <p className={styles.errorText}>Error loading interaction visualization: {interactionVizError}</p>;
        }
        if (!selectedPair) {
             return <p>Select a Ligand-Receptor pair from the table to visualize interactions.</p>;
        }
        if (interactionVizData && selectedPair && apiScopeName && apiScopeName !== 'custom') {
            return (
                <> 
                   {/* Display warnings above the visualization if they exist */}
                   {interactionVizWarnings && interactionVizWarnings.length > 0 && (
                       <div className={styles.warningsContainer}>
                           <h4>Warnings:</h4>
                           <ul>{interactionVizWarnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
                       </div>
                   )}
                   <InteractionVisualization 
                       data={interactionVizData}
                       ligandName={selectedPair[0]}
                       receptorName={selectedPair[1]}
                       currentScope={currentScope}
                   />
                </>
            );
        }
         // Handle the case where there's no data AFTER loading and without error (e.g., API returned empty)
         if (!isLoadingInteractionViz && !interactionVizError && selectedPair) {
             return <p>No interaction visualization data found for {selectedPair.join('-')} in scope '{apiScopeName}'.</p>;
         }

        return <p>Select a Ligand-Receptor pair to begin.</p>; // Default initial state
    }
  };

  return (
    // Main container with flex layout
    <div className={styles.summaryLayout}> 
      {/* FIX: Conditionally render the table section only if scope is NOT custom */}
      {currentScope !== 'custom' && ( 
        <div className={styles.tableArea}> 
          <h3>Interaction Scores ({apiScopeName || currentScope})</h3> 
          <AnalysisTable
            data={combinedData}
            columns={combinedColumns}
            onRowClick={handleTableRowClick} 
            selectedRowIndex={selectedRowIndex} 
            loading={isLoadingInteractionViz}
          />
        </div>
      )} 

      {/* Always render the visualization area, it will expand if table is hidden */}
      <div className={styles.visualizationArea}> 
        <h3>
          {currentScope === 'custom' 
              ? 'Spatial Overview (Select Region)' // Updated title for clarity
              : `Interaction Visualization (${selectedPair ? selectedPair.join('-') : 'No Pair Selected'})`}
        </h3>
        {renderVisualization()}
      </div>
    </div>
  );
};

export default SummaryTabContent; 