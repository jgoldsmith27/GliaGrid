import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom'; // Import useParams
import styles from './ResultsPage.module.css';
import ControlPanel from '../../components/ControlPanel/ControlPanel';
import DisplayPanel from '../../components/DisplayPanel/DisplayPanel';

// Define types for the selections
type ScopeType = 'whole_tissue' | 'layers';
type AnalysisType = 'summary' | 'pathway_dominance' | 'module_context';

// interface ResultsPageProps { // No longer needed
//   jobId: string | null;
// }

const ResultsPage: React.FC = () => { // Define as standard functional component
  const { jobId } = useParams<{ jobId: string }>(); // Get jobId from URL
  
  const [jobStatus, setJobStatus] = useState<any>(null); // Store full status/results
  const [isConnected, setIsConnected] = useState<boolean>(false); // Track WebSocket connection
  const [isLoading, setIsLoading] = useState<boolean>(true); // Still useful for initial load/connection
  const [error, setError] = useState<string | null>(null);
  
  // State for UI controls
  const [selectedScope, setSelectedScope] = useState<ScopeType>('whole_tissue');
  const [selectedLayers, setSelectedLayers] = useState<string[]>([]); // e.g., ['layer_1', 'layer_2']
  const [selectedAnalysisType, setSelectedAnalysisType] = useState<AnalysisType>('summary');
  // TODO: Add state for filters

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

  // --- UI Rendering --- 

  // Display connection status or initial loading message
  if (!isConnected && isLoading && !error) {
     return (
      <div className={styles.container}>
        <h1>Connecting to Analysis Stream...</h1>
      </div>
    );
  }
  
  if (isLoading && isConnected) {
    return (
      <div className={styles.container}>
        <h1>Loading Analysis Results...</h1>
        {jobStatus && jobStatus.message && <p>{jobStatus.message}</p>}
        {/* Optional: Add a spinner or refined progress indicator here */}
        {jobStatus && jobStatus.progress !== null && jobStatus.progress !== undefined && (
          <progress value={jobStatus.progress} max="1"></progress>
        )}
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

  if (!jobStatus) {
     return (
      <div className={styles.container}>
        <h1>No Status Information</h1>
        <p>Could not retrieve status for job ID: {jobId}</p>
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

  if (jobStatus.status !== 'success' || !jobStatus.results) {
      return (
      <div className={styles.container}>
        <h1>Analysis In Progress or Incomplete</h1>
        <p>Status: {jobStatus.status}</p>
        <p>{jobStatus.message}</p>
        {/* Re-render loading or provide a manual refresh? */}
      </div>
    );
  }

  // We have successful results here
  const results = jobStatus.results;

  // TODO: Extract available layers from results for the LayerSelector
  const availableLayers = results ? Object.keys(results).filter(k => k !== 'whole_tissue') : [];

  // TODO: Filter/select data based on state (selectedScope, selectedLayers, selectedAnalysisType)
  const getDataForDisplay = () => {
      if (!results) return null;
      
      let scopeData = selectedScope === 'whole_tissue' 
          ? results.whole_tissue 
          : selectedLayers.length > 0 
              ? results[selectedLayers[0]] // Simplistic: show first selected layer for now
              : null;

      if (!scopeData) return null;

      switch(selectedAnalysisType) {
          case 'summary': return scopeData.ligand_receptor_counts;
          case 'pathway_dominance': return scopeData.pathway_dominance;
          case 'module_context': return scopeData.module_context;
          default: return null;
      }
  };

  const displayData = getDataForDisplay();

  return (
    <div className={styles.container}>
      <h1 className={styles.pageTitle}>Analysis Results</h1>
      
      <div className={styles.resultsLayout}>
          <div className={styles.controlPanelContainer}>
              <ControlPanel 
                  selectedScope={selectedScope}
                  onScopeChange={setSelectedScope}
                  availableLayers={availableLayers}
                  selectedLayers={selectedLayers}
                  onLayersChange={setSelectedLayers}
                  selectedAnalysisType={selectedAnalysisType}
                  onAnalysisTypeChange={setSelectedAnalysisType}
              />
          </div>
          <div className={styles.displayPanelContainer}>
              <DisplayPanel 
                  analysisType={selectedAnalysisType}
                  data={displayData}
                  scope={selectedScope}
                  layers={selectedLayers}
              />
          </div>
          {/* Optional Comparison Panel - can add later */}
          {/* <div className={styles.comparisonPanelContainer}>
              <h2>Comparison (Placeholder)</h2>
          </div> */}
      </div>

    </div>
  );
};

export default ResultsPage; 