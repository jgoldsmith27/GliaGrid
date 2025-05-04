import React, { useState, useEffect, useCallback, useRef } from 'react';
import { useParams } from 'react-router-dom'; // Import useParams
import styles from './ResultsPage.module.css';
import LoadingSpinner from '../../components/LoadingSpinner/LoadingSpinner';
import SummaryTabContent from '../../components/ResultsPage/SummaryTabContent'; // New name
import ScopeSelector, { ScopeType } from '../../components/ScopeSelector/ScopeSelector'; // Import ScopeSelector and its exported type
import LayerSelector from '../../components/LayerSelector/LayerSelector'; // Assuming this exists
import { 
    PathwayDominanceResult, ModuleContextResult, 
    CustomAnalysisResultsBundle 
} from '../../types/analysisResults'; // Import types
// Import the new event-driven hook
import useJobStatus from '../../hooks/useJobStatus';
// Import the spatial visualization component
import SpatialOverviewVisualization from '../../components/SpatialOverviewVisualization/SpatialOverviewVisualization';
// Import the hook for spatial stream data - REMOVED
// import { useSpatialStreamData } from '../../services/data/SharedDataStore';

// Define the structure for combined data used in tables/visualizations
export interface CombinedInteractionData extends Partial<PathwayDominanceResult>, Partial<ModuleContextResult> {
  // Ensure required keys are present if needed, or keep all optional with Partial
  ligand: string; // Assume ligand/receptor are always present
  receptor: string;
}

// interface ResultsPageProps { // No longer needed
//   jobId: string | null;
// }

// ADDED: State for custom aggregation level (lifted from SummaryTabContent)
// const [customAggregationLevel, setCustomAggregationLevel] = useState<CustomAggregationLevel>('whole_custom');

const ResultsPage: React.FC = () => { // Define as standard functional component
  // ADDED: State for custom aggregation level (lifted from SummaryTabContent)
  const [customAggregationLevel, setCustomAggregationLevel] = useState<CustomAggregationLevel>('whole_custom');
  
  const { jobId } = useParams<{ jobId: string }>(); // Get jobId from URL
  
  // Use the event-driven hook to manage job status
  const { jobStatus, isLoading: isJobStatusLoading, error: jobStatusError, refreshStatus } = useJobStatus(jobId);
  
  // --- Get Spatial Stream State from SharedDataStore --- REMOVED ---
  // const { 
  //   points: spatialPoints, 
  //   isLoading: isSpatialLoading, 
  //   error: spatialError, 
  //   pointsReceived: spatialPointsReceived 
  // } = useSpatialStreamData();
  
  // --- UI/Data State ---
  // const [currentTab, setCurrentTab] = useState<ResultsTab>('Summary'); // Removed
  const [selectedScope, setSelectedScope] = useState<ScopeType>('whole_tissue');
  const [availableLayers, setAvailableLayers] = useState<string[]>([]);
  const [selectedLayers, setSelectedLayers] = useState<string[]>([]); 
  const [selectedPair, setSelectedPair] = useState<[string, string] | null>(null);
  const [combinedAnalysisData, setCombinedAnalysisData] = useState<CombinedInteractionData[]>([]);
  const [lassoCoords, setLassoCoords] = useState<[number, number][] | null>(null); // ADDED state for lasso coords
  // ADDED: State for custom analysis results
  // MODIFIED: Use the new bundle type for state
  const [customAnalysisResults, setCustomAnalysisResults] = useState<CustomAnalysisResultsBundle | null>(null); 
  const [isLoadingCustomAnalysis, setIsLoadingCustomAnalysis] = useState(false);
  const [customAnalysisError, setCustomAnalysisError] = useState<string | null>(null);

  // --- State for Spatial Overview Data (managed here for persistence) ---
  // REMOVED - Now managed by SharedDataStore via useSpatialStreamData hook
  // const [spatialPoints, setSpatialPoints] = useState<any[]>([]); 
  // const [isSpatialLoading, setIsSpatialLoading] = useState<boolean>(false);
  // const [spatialError, setSpatialError] = useState<string | null>(null);
  // const [spatialPointsReceived, setSpatialPointsReceived] = useState<number>(0);

  // Effect to extract available layers from jobStatus (updated from hook)
  useEffect(() => {
      // Access results from the 'outputs' field
      const outputs = jobStatus?.results?.outputs;
      if (outputs) {
         setAvailableLayers(Object.keys(outputs).filter(k => k !== 'whole_tissue'));
      } else {
         setAvailableLayers([]); // Clear if no results
      }
  }, [jobStatus]); // Depend on jobStatus from the hook

  // Effect to ensure a layer is selected when in 'layers' scope
  useEffect(() => {
    // Use the derived availableLayers state
    if (selectedScope === 'layers' && 
        (selectedLayers.length === 0 || !availableLayers.includes(selectedLayers[0]))) {
      if (availableLayers.length > 0) {
        setSelectedLayers([availableLayers[0]]);
      }
    }
  }, [selectedScope, availableLayers, selectedLayers]); // Depend on derived layers

  // --- Effect to Process Results and Combine Data ---
  useEffect(() => {
      // Access results from the 'outputs' field
      const outputs = jobStatus?.results?.outputs;
      if (!jobStatus || !outputs) { 
          setCombinedAnalysisData([]);
          setSelectedPair(null);
          return;
      }

      let scopeData: any = null;
      if (selectedScope === 'whole_tissue') {
          scopeData = outputs.whole_tissue;
      } else if (selectedScope === 'layers' && selectedLayers.length > 0) {
          scopeData = outputs[selectedLayers[0]];
      } else if (selectedScope === 'custom') {
          // Custom scope doesn't rely on combined data from results in the same way
          setCombinedAnalysisData([]);
          setSelectedPair(null);
          return;
      }

      if (!scopeData || !scopeData.pathway_dominance || !scopeData.module_context) {
          console.warn("Missing pathway_dominance or module_context data for scope:", selectedScope, selectedLayers[0]);
          setCombinedAnalysisData([]);
          setSelectedPair(null);
          return;
      }

      const pathwayData: PathwayDominanceResult[] = scopeData.pathway_dominance;
      const moduleData: ModuleContextResult[] = scopeData.module_context;

      // Create a map for faster module context lookup
      const moduleContextMap = new Map<string, ModuleContextResult>();
      moduleData.forEach(item => {
          moduleContextMap.set(`${item.ligand}-${item.receptor}`, item);
      });

      const combined: CombinedInteractionData[] = pathwayData.map(pathwayItem => {
          const moduleItem = moduleContextMap.get(`${pathwayItem.ligand}-${pathwayItem.receptor}`);
          
          return {
              ligand: pathwayItem.ligand,
              receptor: pathwayItem.receptor,
              score: pathwayItem.score,
              ligand_norm_expr: pathwayItem.ligand_norm_expr,
              receptor_avg_norm_expr: pathwayItem.receptor_avg_norm_expr,
              // Merge module context data if found
              interaction_type: moduleItem?.interaction_type,
              ligand_module: moduleItem?.ligand_module,
              receptor_modules: moduleItem?.receptor_modules,
              is_same_module: moduleItem?.is_same_module,
          };
      });

      setCombinedAnalysisData(combined);
      
      // If the currently selected pair is no longer in the combined data (e.g., scope changed), reset it
      if (selectedPair && !combined.some(item => item.ligand === selectedPair[0] && item.receptor === selectedPair[1])) {
           setSelectedPair(null);
      }
      // Ensure selectedPair is null if combined data is empty
      if (combined.length === 0) {
          setSelectedPair(null);
      }
      
  }, [jobStatus, selectedScope, selectedLayers]); // Re-run when results or scope change

  // --- Effect to Stream Spatial Data for Overview ---
  // REMOVED - Streaming is initiated from DataInputPage and managed by SharedDataStore
  // useEffect(() => { ... }, [...]); 

  const handleSelectPair = useCallback((pair: [string, string] | null) => {
    console.log("[ResultsPage] Setting selected pair:", pair);
    setSelectedPair(pair);
  }, []);

  const handleLassoSelect = useCallback((coords: [number, number][] | null) => {
    console.log("[ResultsPage] Lasso selection coords:", coords);
    setLassoCoords(coords);
    // Clear previous custom results when lasso selection changes/clears
    setCustomAnalysisResults(null);
    setCustomAnalysisError(null);
    setIsLoadingCustomAnalysis(false);
    // TODO: Add logic here later to use the coords, e.g., trigger a custom analysis
    // if (coords) { setSelectedScope('custom'); }
  }, []);

  // ADDED: Handler to trigger custom analysis
  const handleAnalyzeLasso = useCallback(async () => {
      if (!jobId || !lassoCoords) {
          console.error("Cannot analyze lasso without Job ID and coordinates.");
          setCustomAnalysisError("Missing job ID or lasso selection.");
          return;
      }
      // REMOVED: Resetting aggregation level - it now just controls view
      console.log(`[ResultsPage] Triggering custom analysis for job ${jobId} with coords:`, lassoCoords);
      setIsLoadingCustomAnalysis(true);
      setCustomAnalysisError(null);
      setCustomAnalysisResults(null);

      try {
          const apiUrl = `/api/analysis/custom/${jobId}`; 
          const response = await fetch(apiUrl, {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              // MODIFIED: Remove aggregation level from body
              body: JSON.stringify({ 
                  polygon: lassoCoords 
                  // REMOVED: aggregation: customAggregationLevel 
              }),
          });

          if (!response.ok) {
              const errorText = await response.text();
              throw new Error(`Custom analysis failed: ${response.status} ${errorText || response.statusText}`);
          }

          // Store the raw results - expecting the bundle structure now
          const results: CustomAnalysisResultsBundle = await response.json(); 
          console.log("[ResultsPage] Received custom analysis results bundle:", results);
          setCustomAnalysisResults(results); // Store the bundle
          setSelectedPair(null); // Clear selected pair for the new custom results view
          // Set view to whole results by default after new analysis
          setCustomAggregationLevel('whole_custom'); 

      } catch (error) {
          console.error("[ResultsPage] Custom analysis error:", error);
          setCustomAnalysisError(error instanceof Error ? error.message : String(error));
          setCustomAnalysisResults(null);
      } finally {
          setIsLoadingCustomAnalysis(false);
      }
  // MODIFIED: Remove customAggregationLevel from dependencies - call is independent of view state
  }, [jobId, lassoCoords]); 
  
  // --- Scope Change Handlers ---
  const handleScopeChange = useCallback((scope: ScopeType) => {
    console.log("[ResultsPage] Scope changed to:", scope);
    setSelectedScope(scope);
    // Clear custom analysis state when changing scope
    setLassoCoords(null);
    setCustomAnalysisResults(null);
    setCustomAnalysisError(null);
    setIsLoadingCustomAnalysis(false);

    if (scope === 'whole_tissue') {
        setSelectedLayers([]);
    } else if (scope === 'custom') {
        setSelectedLayers([]);
        setSelectedPair(null); 
    }
  }, []);

  const handleLayersChange = useCallback((layers: string[]) => {
    setSelectedLayers(layers);
  }, []);

  // --- UI Rendering --- 

  // Display initial loading message
  if (isJobStatusLoading && !jobStatus) {
    return (
      <div className={styles.container}>
        <LoadingSpinner message="Loading Analysis Status..." />
      </div>
    );
  }
  
  if (jobStatus && jobStatus.status === 'running') {
    return (
      <div className={styles.container}>
        <LoadingSpinner message={jobStatus?.message || "Processing Analysis..."} />
        {jobStatus.progress !== null && jobStatus.progress !== undefined && (
          <progress value={jobStatus.progress} max="1"></progress>
        )}
      </div>
    );
  }

  // Show waiting message if we're still waiting for results
  if (!jobStatus || (jobStatus.status === 'pending' && !isJobStatusLoading)) {
    return (
      <div className={styles.container}>
        <LoadingSpinner message="Waiting for analysis to start..." />
      </div>
    );
  }

  if (jobStatusError) {
    return (
      <div className={styles.container}>
        <h1>Error Loading Job Status</h1>
        <p className={styles.errorMessage}>{jobStatusError}</p>
      </div>
    );
  }

  if (jobStatus.status === 'failed') {
    return (
      <div className={styles.container}>
        <h1>Analysis Failed</h1>
        <p className={styles.errorMessage}>Error: {jobStatus.message || 'Unknown error occurred during analysis.'}</p>
      </div>
    );
  }

  // We have successful results here
  // const results = jobStatus.results; // Keep results accessible if needed

  // Determine which layer name to pass to SummaryTabContent for API calls
  // For custom scope, we don't need a specific layer name for the interaction fetch (it won't run)
  // For layer scope, use the first selected layer.
  // For whole_tissue, use 'whole_tissue' (or whatever the API expects).
  const scopeForApi = selectedScope === 'layers' && selectedLayers.length > 0 
                      ? selectedLayers[0] 
                      : selectedScope; // Fallback to 'whole_tissue' or 'custom'

  // --- Render Main Layout --- 
  console.log("[ResultsPage] Passing onAnalyzeSelection:", typeof handleAnalyzeLasso); // LOGGING
  return (
    <div className={styles.resultsPageLayout}> {/* New overall layout class */} 
      {/* Control Area */}
      <div className={styles.controlArea}>
        <h3>Scope & Settings</h3>
        <div className={styles.controlGroup}>
           {/* Remove this label */}
           {/* <label className={styles.controlLabel}>Analysis Scope:</label> */}
            <ScopeSelector 
                selectedScope={selectedScope} 
                onScopeChange={handleScopeChange} 
            />
        </div>
        
        {/* Conditionally render the Layer Selector container */} 
        {selectedScope === 'layers' && (
            <div className={styles.layerSelectorContainer}> {/* Wrapper Div */} 
                <LayerSelector 
                    availableLayers={availableLayers} 
                    selectedLayers={selectedLayers} 
                    onLayersChange={handleLayersChange} 
                />
                {/* Keep warning inside the container if relevant */} 
                {selectedLayers.length > 1 && (
                    <p className={styles.warning}>Warning: Multiple layers selected, analysis uses only the first ({selectedLayers[0]}).</p>
                )}
            </div>
        )}
        
         {/* Add other filters/settings here later */}
      </div>

      {/* Content Area - Render based on selected scope */}
      <div className={styles.contentArea}>
          {selectedScope === 'custom' ? (
              // Render Spatial Overview PLUS Loading/Results/Error Area
              <div className={styles.customScopeLayout}> {/* Use a specific layout class */}
                  <div className={styles.spatialOverviewAreaCustom}> {/* Specific class for spatial in custom */} 
                      {jobId && (
                          <SpatialOverviewVisualization 
                              jobId={jobId}
                              onLassoSelect={handleLassoSelect}
                              onAnalyzeSelection={handleAnalyzeLasso}
                          />
                      )}
                  </div>
                  <div className={styles.customResultsArea}> {/* Area below spatial for results/loading */} 
                      {isLoadingCustomAnalysis ? (
                          <div className={styles.loadingContainer}> 
                              <LoadingSpinner message="Running custom analysis..." />
                          </div>
                      ) : customAnalysisError ? (
                          <p className={styles.errorText}>Error during custom analysis: {customAnalysisError}</p>
                      ) : customAnalysisResults ? (
                           // Pass custom results to SummaryTabContent for display
                           // Note: SummaryTabContent might need internal adaptation 
                           // to handle receiving custom results this way.
                           <SummaryTabContent 
                                jobId={jobId || ''} 
                                combinedData={[]} // Pass empty for original data
                                selectedPair={selectedPair}
                                onSelectPair={handleSelectPair}
                                apiScopeName={null} // No API scope for custom results display
                                currentScope={selectedScope} 
                                // Pass custom analysis state and handlers
                                customAnalysisResults={customAnalysisResults} 
                                isLoadingCustomAnalysis={false} // Already handled loading above
                                customAnalysisError={null} // Already handled error above
                                // ADDED: Pass aggregation state and setter
                                customAggregationLevel={customAggregationLevel}
                                setCustomAggregationLevel={setCustomAggregationLevel}
                            />
                      ) : (
                          // Initial state before analysis is run
                          <p className={styles.placeholderText}>Select a region in the spatial view and click 'Analyze Selection'.</p>
                      )}
                  </div>
              </div>
          ) : (
              // Render Summary Table/Interaction Viz Area for other scopes
              <div className={styles.summaryInteractionArea}> 
                  <SummaryTabContent 
                        jobId={jobId || ''} 
                        combinedData={combinedAnalysisData}
                        selectedPair={selectedPair}
                        onSelectPair={handleSelectPair}
                        apiScopeName={scopeForApi} 
                        currentScope={selectedScope} 
                        // Pass lasso/analyze handlers only for non-custom scopes implicitly via SummaryTab
                        // onLassoSelect={handleLassoSelect} // Removed, SummaryTab doesn't need it directly?
                        // onAnalyzeSelection={handleAnalyzeLasso} // Removed, SummaryTab doesn't trigger this
                        // Pass null for custom state when not in custom scope
                        customAnalysisResults={null} 
                        isLoadingCustomAnalysis={false} 
                        customAnalysisError={null} 
                        // ADDED: Pass aggregation state and setter (even if null/default)
                        customAggregationLevel={customAggregationLevel} 
                        setCustomAggregationLevel={setCustomAggregationLevel}
                    />
              </div>
          )}
      </div>
    </div>
  );
};

// Type for custom aggregation level (can be moved to a shared types file)
type CustomAggregationLevel = 'whole_custom' | 'custom_by_layer';

// REMOVED Simplified Type definitions at the end
// interface AnalysisResultItem { // Simplified for example
//     ligand?: string | null;
//     receptor?: string | null;
//     score?: number | null;
//     [key: string]: any; // Allow extra fields
// }
// 
// interface CustomAnalysisResponse {
//     pathway_dominance: AnalysisResultItem[];
//     module_context: AnalysisResultItem[];
// }
// 
// interface LayeredCustomAnalysisResponse {
//     results_by_layer: { [layerName: string]: CustomAnalysisResponse };
// }

export default ResultsPage; 