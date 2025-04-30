import React, { useState, useEffect, useCallback } from 'react';
import { useParams } from 'react-router-dom'; // Import useParams
import styles from './ResultsPage.module.css';
import LoadingSpinner from '../../components/LoadingSpinner/LoadingSpinner';
import SummaryTabContent from '../../components/ResultsPage/SummaryTabContent'; // New name
import ScopeSelector, { ScopeType } from '../../components/ScopeSelector/ScopeSelector'; // Import ScopeSelector and its exported type
import LayerSelector from '../../components/LayerSelector/LayerSelector'; // Assuming this exists
import { PathwayDominanceResult, ModuleContextResult } from '../../types/analysisResults'; // Import types
// Import the new event-driven hook
import useJobStatus from '../../hooks/useJobStatus';

// Define the structure for combined data
export interface CombinedInteractionData {
  ligand: string;
  receptor: string;
  score?: number; // Pathway score
  ligand_norm_expr?: number;
  receptor_avg_norm_expr?: number;
  interaction_type?: string;
  ligand_module?: string;
  receptor_modules?: string[];
  is_same_module?: boolean;
}

// interface ResultsPageProps { // No longer needed
//   jobId: string | null;
// }

const ResultsPage: React.FC = () => { // Define as standard functional component
  const { jobId } = useParams<{ jobId: string }>(); // Get jobId from URL
  
  // Use the event-driven hook to manage job status
  const { jobStatus, isLoading, error, refreshStatus } = useJobStatus(jobId);
  
  // --- UI/Data State ---
  // const [currentTab, setCurrentTab] = useState<ResultsTab>('Summary'); // Removed
  const [selectedScope, setSelectedScope] = useState<ScopeType>('whole_tissue');
  const [availableLayers, setAvailableLayers] = useState<string[]>([]);
  const [selectedLayers, setSelectedLayers] = useState<string[]>([]); 
  const [selectedPair, setSelectedPair] = useState<[string, string] | null>(null);
  const [combinedAnalysisData, setCombinedAnalysisData] = useState<CombinedInteractionData[]>([]);

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

  const handleSelectPair = useCallback((pair: [string, string] | null) => {
    console.log("[ResultsPage] Setting selected pair:", pair);
    setSelectedPair(pair);
  }, []);
  
  // --- Scope Change Handlers ---
  const handleScopeChange = useCallback((scope: ScopeType) => {
    console.log("[ResultsPage] Scope changed to:", scope);
    setSelectedScope(scope);
    if (scope === 'whole_tissue') {
        setSelectedLayers([]);
        // Optionally reset pair: setSelectedPair(null);
    } else if (scope === 'custom') {
        setSelectedLayers([]);
        setSelectedPair(null); // Reset pair for custom view
    }
  }, []);

  const handleLayersChange = useCallback((layers: string[]) => {
    setSelectedLayers(layers);
  }, []);

  // --- UI Rendering --- 

  // Display initial loading message
  if (isLoading && !jobStatus) {
    return (
      <div className={styles.container}>
        <LoadingSpinner message="Loading Analysis Results..." />
      </div>
    );
  }
  
  if (isLoading && jobStatus && jobStatus.status === 'running') {
    return (
      <div className={styles.container}>
        <LoadingSpinner message={jobStatus?.message || "Processing Analysis..."} />
        {jobStatus && jobStatus.progress !== null && jobStatus.progress !== undefined && (
          <progress value={jobStatus.progress} max="1"></progress>
        )}
      </div>
    );
  }

  // Show waiting message if we're still waiting for results
  if (!jobStatus || (jobStatus.status === 'pending' && !isLoading)) {
    return (
      <div className={styles.container}>
        <LoadingSpinner message="Waiting for analysis to start..." />
      </div>
    );
  }

  if (error) {
    return (
      <div className={styles.container}>
        <h1>Error Loading Results</h1>
        <p className={styles.errorMessage}>{error}</p>
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

      {/* Content Area */} 
      <div className={styles.contentArea}>
          <SummaryTabContent 
                jobId={jobId || ''} 
                combinedData={combinedAnalysisData}
                selectedPair={selectedPair}
                onSelectPair={handleSelectPair}
                // Pass the calculated scope name for API calls
                apiScopeName={scopeForApi} 
                // Keep passing the general scope type for conditional rendering if needed
                currentScope={selectedScope} 
            />
      </div>
    </div>
  );
};

export default ResultsPage; 