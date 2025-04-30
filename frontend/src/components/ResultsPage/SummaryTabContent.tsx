import React, { useState, useEffect, useRef, useCallback } from 'react';
import styles from './SummaryTabContent.module.css'; // Create CSS module later
import AnalysisTable, { ColumnDefinition } from '../AnalysisTable/AnalysisTable';
import InteractionVisualization from '../InteractionVisualization/InteractionVisualization.tsx';
import SpatialOverviewVisualization from '../SpatialOverviewVisualization/SpatialOverviewVisualization';
import { CombinedInteractionData } from '../../pages/ResultsPage/ResultsPage'; // Import shared type
import { ScopeType } from '../ScopeSelector/ScopeSelector'; // Import ScopeType if defined there
import useInteractionData from '../../hooks/useInteractionData'; // Import the hook
import { PathwayDominanceResult, ModuleContextResult } from '../../types/analysisResults'; // ADDED import
import LoadingSpinner from '../LoadingSpinner/LoadingSpinner'; // ADDED import for loading state

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
  onAnalyzeSelection?: () => void;
  customAnalysisResults?: any | null;
  isLoadingCustomAnalysis?: boolean;
  customAnalysisError?: string | null;
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
  onAnalyzeSelection,
  customAnalysisResults,
  isLoadingCustomAnalysis,
  customAnalysisError,
}) => {
  // ADDED: Log received prop
  console.log("[SummaryTabContent] Received onAnalyzeSelection:", typeof onAnalyzeSelection);

  // State for SpatialOverviewVisualization
  // const [allPointsData, setAllPointsData] = useState<AllPointsData[] | null>(null);
  // const [loadingAllPoints, setLoadingAllPoints] = useState(false);
  // const [allPointsError, setAllPointsError] = useState<string | null>(null);

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
  } = useInteractionData(jobId, selectedPair, currentScope !== 'custom' ? apiScopeName : null); // Only fetch if not custom scope initially

  // Fetch function for All Points data
  const fetchAllPoints = useCallback(async (currentJobId: string) => {
     if (abortControllerRef.current) {
        abortControllerRef.current.abort(); // Abort previous fetch
    }
    const controller = new AbortController();
    abortControllerRef.current = controller;

    // setLoadingAllPoints(true);
    // setAllPointsError(null);
    // setAllPointsData(null);

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
        // setAllPointsData(data);
    } catch (err) {
        if ((err as Error).name === 'AbortError') {
             // setAllPointsError('Request cancelled.');
        } else {
            const message = err instanceof Error ? err.message : 'An unknown error occurred';
            // setAllPointsError(message);
        }
        // setAllPointsData(null);
    } finally {
         if (abortControllerRef.current === controller) {
             abortControllerRef.current = null;
        }
        // setLoadingAllPoints(false);
    }
  }, []); // Dependency: jobId implicitly via usage

  // Effect to fetch data based on currentScope
  useEffect(() => {
    console.log(`[SummaryTabContent] useEffect triggered. Scope: ${currentScope}, Pair: ${selectedPair ? selectedPair.join('-') : 'null'}, apiScopeName: ${apiScopeName}`); 
    
    // Reset spatial overview state when dependencies change
    // setAllPointsData(null);
    // setAllPointsError(null);
    // setLoadingAllPoints(false);
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
    // This handler applies to BOTH original and custom analysis tables
    console.log(`[SummaryTabContent] Row clicked: ${row.ligand}-${row.receptor}, Index: ${index}, Scope: ${currentScope}`); 
    setSelectedRowIndex(index);
    const newPair: [string, string] = [row.ligand, row.receptor];
    onSelectPair(newPair); // Call parent handler to manage selectedPair state
  }, [onSelectPair, currentScope]);

  // Generic cancel handler - now only cancels the spatial overview fetch
  const handleCancelVizRequest = useCallback(() => {
      if (abortControllerRef.current) {
          console.log("[SummaryTabContent] Cancelling spatial overview fetch request.");
          abortControllerRef.current.abort();
          abortControllerRef.current = null; // Clear ref after abort
          // setLoadingAllPoints(false); // Ensure loading stops
          // setAllPointsError("Request cancelled by user."); // Set error state
      } else {
          console.log("[SummaryTabContent] No active spatial overview fetch request to cancel.");
      }
      // Note: To cancel the interaction viz, we would call cancelInteractionVizFetch()
      // but it needs a separate button or logic tied to the interaction viz state.
  }, []); // Removed dependency on abortControllerRef state, check directly

  // Helper to render the correct component based on scope and custom results
  const renderContent = () => {
    if (currentScope === 'custom') {
        // Show loading spinner first if custom analysis is running
        if (isLoadingCustomAnalysis) {
            return (
                <div className={styles.loadingContainer}> 
                    <LoadingSpinner message="Running custom analysis..." />
                </div>
            );
        }
        // If not loading, check for errors
        if (customAnalysisError) {
            return <p className={styles.errorText}>Error during custom analysis: {customAnalysisError}</p>;
        }
        // If not loading and no error, check for results
        if (customAnalysisResults) {
            // Display results similar to non-custom scope
            // We need to process customAnalysisResults into CombinedInteractionData format
            // Assuming customAnalysisResults has { pathway_dominance: [], module_context: [] }
            const pathwayData = customAnalysisResults.pathway_dominance || [];
            const moduleData = customAnalysisResults.module_context || [];
            const moduleContextMap = new Map<string, ModuleContextResult>();
            moduleData.forEach((item: ModuleContextResult) => { moduleContextMap.set(`${item.ligand}-${item.receptor}`, item); });
            const customCombinedData: CombinedInteractionData[] = pathwayData.map((pathwayItem: PathwayDominanceResult) => {
                const moduleItem = moduleContextMap.get(`${pathwayItem.ligand}-${pathwayItem.receptor}`);
                return { ...pathwayItem, ...moduleItem }; // Simple merge
            });

            // Use the interaction hook for custom results visualization
            // NOTE: This re-uses the hook; ensure dependencies are correct
            const { interactionVizData: customVizData, isLoading: isLoadingCustomViz, error: customVizError } = 
                useInteractionData(jobId, selectedPair, null); // MODIFIED: Pass null for apiScopeName for custom results

            return (
              <div className={styles.summaryLayout}> {/* Reuse existing layout */}
                <div className={styles.tableArea}>
                  <h3>Interaction Scores (Custom Selection)</h3> 
                  <AnalysisTable
                    data={customCombinedData}
                    columns={combinedColumns}
                    onRowClick={handleTableRowClick} 
                    selectedRowIndex={selectedRowIndex} 
                    loading={isLoadingCustomViz} // Use custom loading state?
                  />
                </div>
                <div className={styles.visualizationArea}>
                  <h3>Interaction Visualization ({selectedPair ? selectedPair.join('-') : 'Select Pair'})</h3>
                  {isLoadingCustomViz && <p>Loading visualization...</p>}
                  {customVizError && <p className={styles.errorText}>Viz Error: {customVizError}</p>}
                  {!selectedPair && <p>Select a pair from the table.</p>}
                  {customVizData && selectedPair && (
                      <InteractionVisualization 
                         data={customVizData}
                         ligandName={selectedPair[0]}
                         receptorName={selectedPair[1]}
                         currentScope={'whole_tissue'}
                     />
                  )}
                </div>
              </div>
            );
        } else {
            // No custom results yet, show the spatial overview for selection
            return (
                <div className={styles.spatialOverviewOnlyArea}> {/* Optional: different style? */}
                   <h3>Spatial Overview (Select Region)</h3>
                   <SpatialOverviewVisualization 
                      jobId={jobId} 
                      onLassoSelect={onLassoSelect} 
                      onAnalyzeSelection={onAnalyzeSelection} // Pass down the trigger
                   />
                </div>
            );
        }
    } else { // 'whole_tissue' or 'layers'
        // Original logic using combinedData from props and interactionVizData from hook
        return (
          <div className={styles.summaryLayout}>
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
            <div className={styles.visualizationArea}>
              <h3>Interaction Visualization ({selectedPair ? selectedPair.join('-') : 'Select Pair'})</h3>
              {isLoadingInteractionViz && <p>Loading visualization...</p>}
              {interactionVizError && <p className={styles.errorText}>Viz Error: {interactionVizError}</p>}
              {!selectedPair && <p>Select a pair from the table.</p>}
              {interactionVizWarnings && interactionVizWarnings.length > 0 && (
                   <div className={styles.warningsContainer}><h4>Warnings:</h4><ul>{interactionVizWarnings.map((w, i) => <li key={i}>{w}</li>)}</ul></div>
              )}
              {interactionVizData && selectedPair && (
                  <InteractionVisualization 
                     data={interactionVizData}
                     ligandName={selectedPair[0]}
                     receptorName={selectedPair[1]}
                     currentScope={currentScope}
                 />
              )}
            </div>
          </div>
        );
    }
  };

  return renderContent(); // Render based on the logic above
};

export default SummaryTabContent; 