import React, { useState, useEffect, useCallback } from 'react';
import { useParams } from 'react-router-dom'; // Import useParams
import styles from './ResultsPage.module.css';
import LoadingSpinner from '../../components/LoadingSpinner/LoadingSpinner';
import SummaryTabContent from '../../components/ResultsPage/SummaryTabContent'; // New name
import ScopeSelector from '../../components/ScopeSelector/ScopeSelector'; // Assuming this exists
import LayerSelector from '../../components/LayerSelector/LayerSelector'; // Assuming this exists
import { PathwayDominanceResult, ModuleContextResult } from '../../types/analysisResults'; // Import types

// Define types for the selections
type ScopeType = 'whole_tissue' | 'layers';
// type AnalysisType = 'summary' | 'pathway_dominance' | 'module_context'; // No longer needed
// type ResultsTab = 'Summary' | 'Custom'; // Removed

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
  
  const [jobStatus, setJobStatus] = useState<any>(null); // Store full status/results
  const [isConnected, setIsConnected] = useState<boolean>(false); // Track WebSocket connection
  const [isLoading, setIsLoading] = useState<boolean>(true); // Still useful for initial load/connection
  const [error, setError] = useState<string | null>(null);
  
  // --- UI/Data State ---
  // const [currentTab, setCurrentTab] = useState<ResultsTab>('Summary'); // Removed
  const [selectedScope, setSelectedScope] = useState<ScopeType>('whole_tissue');
  const [availableLayers, setAvailableLayers] = useState<string[]>([]);
  const [selectedLayers, setSelectedLayers] = useState<string[]>([]); 
  const [selectedPair, setSelectedPair] = useState<[string, string] | null>(null);
  const [combinedAnalysisData, setCombinedAnalysisData] = useState<CombinedInteractionData[]>([]);

  useEffect(() => {
    if (!jobId) {
      setError('No analysis job ID provided.');
      setIsLoading(false);
      return;
    }

    // --- WebSocket Logic ---
    const wsUrl = `ws://localhost:8000/ws/analysis/status/${jobId}`;
    const ws = new WebSocket(wsUrl);

    ws.onopen = () => {
      console.log(`[WebSocket] Connected to ${wsUrl}`);
      setIsConnected(true);
      setIsLoading(true); // Reset loading state on connect/reconnect
      setError(null);
      // No need to fetch initial status, endpoint sends it on connect
    };

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data);
        console.log('[WebSocket] Message received:', data);
        setJobStatus(data);
        if (data.results) {
            setAvailableLayers(Object.keys(data.results).filter(k => k !== 'whole_tissue'));
        }

        // Update loading state based on received status
        if (data.status === 'success' || data.status === 'failed' || data.status === 'error') {
          setIsLoading(false);
        } else {
          setIsLoading(true); // Still processing
        }
        
        // Handle specific error messages from WebSocket
        if (data.status === 'error') {
            setError(data.message || 'An error occurred via WebSocket.');
        }

      } catch (parseError) {
        console.error('[WebSocket] Error parsing message:', parseError, 'Data:', event.data);
        setError('Received invalid data from server.');
        setIsLoading(false);
      }
    };

    ws.onerror = (event) => {
      console.error('[WebSocket] Error:', event);
      setError('WebSocket connection error. Please try refreshing the page.');
      setIsConnected(false);
      setIsLoading(false);
    };

    ws.onclose = (event) => {
      console.log(`[WebSocket] Disconnected from ${wsUrl}. Code: ${event.code}, Reason: ${event.reason}`);
      setIsConnected(false);
      // Optionally set loading to false only if the final status wasn't success/failed
      if (jobStatus?.status !== 'success' && jobStatus?.status !== 'failed') {
         setIsLoading(false);
         // Keep existing status unless it was clearly an error during connection
         if (!error && event.code !== 1000) { // 1000 = Normal closure
            setError('Disconnected. Analysis might be incomplete. Refresh?');
         }
      }
    };

    // Cleanup function to close WebSocket when component unmounts
    return () => {
      if (ws.readyState === WebSocket.OPEN || ws.readyState === WebSocket.CONNECTING) {
          console.log(`[WebSocket] Closing connection to ${wsUrl}`);
          ws.close(1000, "Component unmounting"); // Normal closure
      }
    };
    // --- End WebSocket Logic ---

    /* Remove Polling Logic
    const fetchStatus = async () => {
      // ... removed ...
    };
    fetchStatus();
    */

  }, [jobId]); // Re-run effect if jobId changes

  // Effect to ensure a layer is selected when in 'layers' scope
  useEffect(() => {
    // Extract available layers from results
    const availableLayers = jobStatus?.results ? 
      Object.keys(jobStatus.results).filter(k => k !== 'whole_tissue') : [];
    
    // If we're in layer scope and no layer is selected, select the first available layer
    if (selectedScope === 'layers' && 
        (selectedLayers.length === 0 || !availableLayers.includes(selectedLayers[0]))) {
      if (availableLayers.length > 0) {
        setSelectedLayers([availableLayers[0]]);
      }
    }
  }, [selectedScope, jobStatus, selectedLayers]);

  // --- Effect to Process Results and Combine Data ---
  useEffect(() => {
      if (!jobStatus || !jobStatus.results) {
          setCombinedAnalysisData([]);
          setSelectedPair(null);
          return;
      }

      let scopeData = selectedScope === 'whole_tissue' 
          ? jobStatus.results.whole_tissue 
          : selectedLayers.length > 0 
              ? jobStatus.results[selectedLayers[0]]
              : null;

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
           setSelectedPair(null); // Reset to null instead of first item
      }
      // Ensure selectedPair is null if combined data is empty
      if (combined.length === 0) {
          setSelectedPair(null);
      }
      
  }, [jobStatus, selectedScope, selectedLayers]); // Re-run when results or scope change

  const handleSelectPair = useCallback((pair: [string, string]) => {
    setSelectedPair(pair);
  }, []);
  
  // --- Scope Change Handlers ---
  const handleScopeChange = useCallback((scope: ScopeType) => {
    setSelectedScope(scope);
    // Reset selected layers when switching to whole_tissue
    if (scope === 'whole_tissue') {
        setSelectedLayers([]);
    }
  }, []);

  const handleLayersChange = useCallback((layers: string[]) => {
    setSelectedLayers(layers);
  }, []);

  // --- UI Rendering --- 

  // Display connection status or initial loading message
  if (!isConnected && isLoading && !error) {
     return (
      <div className={styles.container}>
        <LoadingSpinner message="Connecting to Analysis Stream..." />
      </div>
    );
  }
  
  if (isLoading && isConnected) {
    return (
      <div className={styles.container}>
        <LoadingSpinner message={jobStatus?.message || "Loading Analysis Results..."} />
        {jobStatus && jobStatus.progress !== null && jobStatus.progress !== undefined && (
          <progress value={jobStatus.progress} max="1"></progress>
        )}
      </div>
    );
  }

  // Always show loading until we have valid results
  if (!jobStatus || !jobStatus.results) {
    return (
      <div className={styles.container}>
        <LoadingSpinner message="Waiting for analysis results..." />
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

  // Extract available layers from results for the LayerSelector - REMOVE THIS
  // const availableLayers = results ? Object.keys(results).filter(k => k !== 'whole_tissue') : [];

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
                // Pass the specific layer string needed for the API call
                selectedLayerName={selectedScope === 'layers' && selectedLayers.length > 0 ? selectedLayers[0] : 'whole_tissue'} 
            />
      </div>
    </div>
  );
};

export default ResultsPage; 