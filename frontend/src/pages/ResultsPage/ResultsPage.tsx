import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { useParams, useNavigate } from 'react-router-dom'; // Import useParams and useNavigate
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
// ADDED: Import types/hooks needed for decoupled viz fetch
import { SharedDataStore, useSharedData, DataRequestOptions } from '../../services/data/SharedDataStore';
import { InteractionVisualizationData } from '../../hooks/useInteractionData'; // Assuming type is exported
import { Tab, Tabs, Box, Typography, CircularProgress, Alert, Button } from '@mui/material';
// ADD TextField for search inputs
import TextField from '@mui/material/TextField';

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
  
  const { jobId } = useParams<{ jobId: string }>();
  const navigate = useNavigate(); // Get navigate function
  
  // Use the event-driven hook to manage job status
  const { jobStatus, isLoading: isJobStatusLoading, error: jobStatusError, refreshStatus } = useJobStatus(jobId);
  
  // ADDED: Extract boundaries directly here
  const layerBoundaries = jobStatus?.results?.outputs?.layer_boundaries;
  
  // --- Get Spatial Stream State from SharedDataStore --- REMOVED ---
  // const { 
  //   points: spatialPoints, 
  //   isLoading: isSpatialLoading, 
  //   error: spatialError, 
  //   pointsReceived: spatialPointsReceived 
  // } = useSpatialStreamData();
  
  // --- UI/Data State ---
  const [currentTab, setCurrentTab] = useState(0);
  const [selectedScope, setSelectedScope] = useState<ScopeType>('whole_tissue');
  const [availableLayers, setAvailableLayers] = useState<string[]>([]);
  const [selectedLayer, setSelectedLayer] = useState<string | null>(null);
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

  // ADDED: State for the currently displayed visualization
  const [displayedVizPair, setDisplayedVizPair] = useState<[string, string] | null>(null);
  const [displayedVizData, setDisplayedVizData] = useState<InteractionVisualizationData | null>(null);
  const [isLoadingDisplayedViz, setIsLoadingDisplayedViz] = useState(false);
  const [displayedVizError, setDisplayedVizError] = useState<string | null>(null);

  // ADDED: State for search terms
  const [ligandSearchTerm, setLigandSearchTerm] = useState<string>('');
  const [receptorSearchTerm, setReceptorSearchTerm] = useState<string>('');

  // ADDED: Get data store instance
  const dataStore = useSharedData();

  // Effect to extract available layers from jobStatus (updated from hook)
  useEffect(() => {
      // Access results from the 'outputs' field
      const outputs = jobStatus?.results?.outputs;
      if (outputs) {
         // Filter out both 'whole_tissue' and 'layer_boundaries'
         setAvailableLayers(Object.keys(outputs).filter(k => k !== 'whole_tissue' && k !== 'layer_boundaries'));
      } else {
         setAvailableLayers([]); // Clear if no results
      }
  }, [jobStatus]); // Depend on jobStatus from the hook

  // Effect to ensure a layer is selected when in 'layers' scope
  useEffect(() => {
    // MODIFIED: Logic for single selectedLayer
    if (selectedScope === 'layers') {
      // If no layer is selected OR the selected layer isn't available anymore
      if (!selectedLayer || !availableLayers.includes(selectedLayer)) {
        // Select the first available layer, or null if none exist
        setSelectedLayer(availableLayers.length > 0 ? availableLayers[0] : null);
      }
    } else {
       // Clear selected layer if not in layers scope
       if (selectedLayer !== null) {
           setSelectedLayer(null);
       }
    }
  }, [selectedScope, availableLayers, selectedLayer]); // MODIFIED dependencies

  // --- Effect to Process Results and Combine Data ---
  useEffect(() => {
      const outputs = jobStatus?.results?.outputs;
      // MODIFIED: Log selectedLayer
      console.log(`[ResultsPage] Processing results. Scope: ${selectedScope}, Selected Layer: ${selectedLayer}, Outputs available: ${!!outputs}`);
      
      if (!jobStatus || !outputs) { 
          setCombinedAnalysisData([]);
          setSelectedPair(null);
          return;
      }

      let scopeData: any = null;
      if (selectedScope === 'whole_tissue') {
          scopeData = outputs.whole_tissue;
      } else if (selectedScope === 'layers' && selectedLayer) { // MODIFIED: Check selectedLayer
          const layerKey = selectedLayer;
          console.log(`[ResultsPage] Accessing layer data for key: ${layerKey}`);
          scopeData = outputs[layerKey];
          console.log(`[ResultsPage] Retrieved scopeData for layer ${layerKey}:`, scopeData ? 'Data found' : 'Data NOT found', scopeData);
      } else if (selectedScope === 'custom') {
          // Custom scope doesn't rely on combined data from results in the same way
          setCombinedAnalysisData([]);
          setSelectedPair(null);
          return;
      }

      if (!scopeData || !scopeData.pathway_dominance || !scopeData.module_context) {
          // MODIFIED: Log selectedLayer
          console.warn(`[ResultsPage] Data structure invalid or missing for scope: ${selectedScope}, Layer Key: ${selectedLayer}. Setting empty table data.`, scopeData);
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

      // ADDED log here
      console.log(`[ResultsPage] Setting combinedAnalysisData with ${combined.length} rows for scope ${selectedScope}${selectedLayer ? ` (${selectedLayer})` : ''}.`);
      setCombinedAnalysisData(combined);
      
      // If the currently selected pair is no longer in the combined data (e.g., scope changed), reset it
      if (selectedPair && !combined.some(item => item.ligand === selectedPair[0] && item.receptor === selectedPair[1])) {
           setSelectedPair(null);
      }
      // Ensure selectedPair is null if combined data is empty
      if (combined.length === 0) {
          setSelectedPair(null);
      }
      
  }, [jobStatus, selectedScope, selectedLayer]); // MODIFIED: Depend on selectedLayer

  // --- Effect to Stream Spatial Data for Overview ---
  // REMOVED - Streaming is initiated from DataInputPage and managed by SharedDataStore
  // useEffect(() => { ... }, [...]); 

  // --- Filtering logic for search terms ---
  const filteredAnalysisData = useMemo(() => {
    if (!combinedAnalysisData) return [];

    let filtered = combinedAnalysisData;

    const ligandQuery = ligandSearchTerm.toLowerCase().trim();
    const receptorQuery = receptorSearchTerm.toLowerCase().trim();

    if (ligandQuery) {
      filtered = filtered.filter(item =>
        item.ligand?.toLowerCase().includes(ligandQuery)
      );
    }

    if (receptorQuery) {
      filtered = filtered.filter(item =>
        item.receptor?.toLowerCase().includes(receptorQuery)
      );
    }

    return filtered;
  }, [combinedAnalysisData, ligandSearchTerm, receptorSearchTerm]);

  // Determine scopeForApi (needed for handleSelectPair)
  const scopeForApi = 
      selectedScope === 'layers' 
          ? selectedLayer // Pass the selected layer string or null
      : selectedScope === 'custom' 
          ? null // Pass null for custom scope, viz fetch handled differently? Check custom logic
          : 'whole_tissue'; // Default to whole_tissue for that scope

  // MODIFIED: handleSelectPair triggers viz fetch for non-custom scopes
  const handleSelectPair = useCallback(async (pair: [string, string] | null) => {
    // Update the pair state regardless of scope (so UI reflects selection)
    setDisplayedVizPair(pair);

    if (!pair) {
        console.log("[ResultsPage] Pair deselected.");
        // Decide if clearing the viz is desired on deselect
        // setDisplayedVizData(null);
        // setIsLoadingDisplayedViz(false);
        // setDisplayedVizError(null);
        return;
    }

    // Only fetch main viz data if NOT in custom scope
    if (selectedScope === 'custom') {
        console.log("[ResultsPage] In custom scope, skipping main viz fetch.");
        // Clear main viz state to avoid showing stale data from other scopes
        setDisplayedVizData(null);
        setIsLoadingDisplayedViz(false);
        setDisplayedVizError(null);
        return;
    }

    // Determine the scope for the API call (whole_tissue or layer name)
    const currentScopeForViz = scopeForApi; // Use the calculated scopeForApi

    // Check if we should fetch (valid scope and jobId)
    if (!jobId || currentScopeForViz === null) {
         // This case handles 'layers' scope before a layer is selected
         console.log(`[ResultsPage] Skipping visualization fetch. Scope is null or no JobId.`);
         setDisplayedVizData(null);
         setIsLoadingDisplayedViz(false);
         setDisplayedVizError('Select a specific Layer to visualize interactions.');
         return;
    }

    // Proceed with fetching for whole_tissue or selected layer
    console.log(`[ResultsPage] Triggering main viz fetch for pair: ${pair.join('-')}, scope: ${currentScopeForViz}`);
    setIsLoadingDisplayedViz(true);
    setDisplayedVizData(null);
    setDisplayedVizError(null);

    try {
        const options: DataRequestOptions = {
            ligand: pair[0],
            receptor: pair[1],
            layer: currentScopeForViz ?? undefined,
            // polygon: lassoCoords || undefined,
        };
        
        const data = await dataStore.requestData(jobId, 'interactionPoints', options);

        console.log("[ResultsPage] Main viz fetch success:", data);
        setDisplayedVizData(data as InteractionVisualizationData);
        setDisplayedVizError(null);

    } catch (error) {
        console.error("[ResultsPage] Main viz fetch error:", error);
        setDisplayedVizError(error instanceof Error ? error.message : String(error));
        setDisplayedVizData(null);
    } finally {
        setIsLoadingDisplayedViz(false);
    }
  }, [jobId, selectedScope, selectedLayer, scopeForApi, dataStore]);

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
    // Clear custom analysis state
    setLassoCoords(null);
    setCustomAnalysisResults(null);
    setCustomAnalysisError(null);
    setIsLoadingCustomAnalysis(false);

    // Clear layer selection if needed
    if (scope !== 'layers') {
        setSelectedLayer(null);
    }
    // Clear custom scope pair selection
    if (scope === 'custom') {
       // setSelectedPair(null); // Removed, table selection not directly tied here
    }
    // *** ADDED: Clear visualization state on scope change ***
    setDisplayedVizPair(null);
    setDisplayedVizData(null);
    setIsLoadingDisplayedViz(false);
    setDisplayedVizError(null);

  }, []);

  // MODIFIED: Simplified handler for single layer change
  const handleLayerChange = useCallback((layer: string | null) => {
    setSelectedLayer(layer);
    // ADDED: Clear visualization when layer changes
    setDisplayedVizPair(null);
    setDisplayedVizData(null);
    setIsLoadingDisplayedViz(false);
    setDisplayedVizError(null);
  }, []);

  const handleTabChange = (event: React.SyntheticEvent, newValue: number) => {
    setCurrentTab(newValue);
  };

  const handleBackClick = () => {
    navigate('/'); // Navigate to the Data Input Page (assuming it's the root)
  };

  // --- UI Rendering --- 

  // ADDED: Define table loading state
  // Considered loading if job is loading AND we don't have results yet,
  // OR if we are in layers scope but haven't determined the specific layer data yet.
  const isTableDataLoading = (isJobStatusLoading && !jobStatus?.results?.outputs) || 
                             (selectedScope === 'layers' && !selectedLayer && availableLayers.length > 0); // Add other conditions if needed

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

  // --- Render Main Layout --- 
  console.log("[ResultsPage] Passing onAnalyzeSelection:", typeof handleAnalyzeLasso); // LOGGING

  const searchBarInputs = (
    <div style={{ display: 'flex', gap: '16px', alignItems: 'baseline' }}> {/* Changed to baseline alignment */}
      <div className={styles.controlGroup} style={{ marginBottom: '0px' }}> {/* Adjust margin if needed */}
        <TextField
          label="Search Ligand"
          variant="outlined"
          value={ligandSearchTerm}
          onChange={(e) => setLigandSearchTerm(e.target.value)}
          sx={{
            width: '240px', // Adjusted width
            '.MuiOutlinedInput-root': { 
                height: '36px', 
                display: 'flex', 
                alignItems: 'center',
                // Removed overflow: 'hidden' to allow natural input scroll if text is too long
            },
            '.MuiInputBase-input': { 
              padding: '6px 10px', 
              fontSize: '0.875rem',
              // Ensure default overflow behavior for text (usually scroll)
            },
            '.MuiInputLabel-root': {
              fontSize: '0.875rem',
              lineHeight: '1', 
            },
            '.MuiInputLabel-outlined.MuiInputLabel-shrink': {
              transform: 'translate(14px, -5px) scale(0.75)', 
            },
          }}
        />
      </div>
      <div className={styles.controlGroup} style={{ marginBottom: '0px'}}> {/* Adjust margin if needed */}
        <TextField
          label="Search Receptor"
          variant="outlined"
          value={receptorSearchTerm}
          onChange={(e) => setReceptorSearchTerm(e.target.value)}
          sx={{
            width: '240px', // Adjusted width
            '.MuiOutlinedInput-root': { 
                height: '36px',
                display: 'flex',
                alignItems: 'center',
                // Removed overflow: 'hidden'
            },
            '.MuiInputBase-input': { 
              padding: '6px 10px',
              fontSize: '0.875rem',
            },
            '.MuiInputLabel-root': {
              fontSize: '0.875rem',
              lineHeight: '1',
            },
            '.MuiInputLabel-outlined.MuiInputLabel-shrink': {
              transform: 'translate(14px, -5px) scale(0.75)',
            },
          }}
        />
      </div>
    </div>
  );

  // --- ADDED: Handler to initiate comparison ---
  const handleInitiateComparison = () => {
    if (!jobStatus?.results?.inputs) {
        console.error("[ResultsPage] Job status or input data not available for comparison.");
        // TODO: Optionally show a user-facing error (e.g., using a snackbar)
        alert("Cannot initiate comparison: Job data is missing.");
        return;
    }

    const inputs = jobStatus.results.inputs;
    // Assuming 'spatialFileId' is the key for the primary data file used in the initial analysis
    const spatialFileId = inputs.files?.spatialFileId; 
    const spatialMapping = inputs.mappings?.spatialMapping;

    if (!spatialFileId || !spatialMapping) { // Keeping this check for the primary spatial file/mapping
        console.error("[ResultsPage] Spatial File ID or Spatial Column Mappings not found in job status for comparison.", { inputs });
        // TODO: Optionally show a user-facing error
        alert("Cannot initiate comparison: Essential spatial file ID or column mappings are missing.");
        return;
    }
    
    // Get all file IDs and all mappings
    const allFileIds = inputs.files; 
    const allMappings = inputs.mappings;

    if (!allFileIds || !allMappings) {
        console.error("[ResultsPage] Complete File IDs or Column Mappings sets not found in job status.", { inputs });
        alert("Cannot initiate comparison: Complete file or mapping information is missing from job status.");
        return;
    }

    let definition: Record<string, any> = {}; // Ensure definition is Record<string, any> or a more specific type
    if (selectedScope === 'layers') {
        if (!selectedLayer) {
            console.error("[ResultsPage] Layer not selected for comparison.");
            alert("Please select a layer to compare.");
            return; 
        }
        definition = { layer_name: selectedLayer };
    } else if (selectedScope === 'custom') {
        if (!lassoCoords || lassoCoords.length === 0) {
            console.error("[ResultsPage] Lasso selection not defined for comparison.");
            alert("Please make a custom selection (lasso) to compare.");
            return; 
        }
        definition = { polygon_coords: lassoCoords };
    }
    // For 'whole_tissue', definition remains {} as per plan

    const selection1Data = {
        source_job_id: jobId, // From useParams
        // file_id: spatialFileId, // Replaced by files object
        files: allFileIds, // Store all file IDs
        type: selectedScope === 'layers' 
              ? 'layer' 
              : selectedScope === 'custom' 
              ? 'lasso' 
              : 'whole_tissue', // Correctly map to backend expected types
        definition: definition,
        // column_mappings: spatialMapping // Replaced by mappings object
        mappings: allMappings // Store all mapping objects
    };

    console.log("[ResultsPage] Initiating comparison with Selection 1:", selection1Data);
    // TODO: Navigate to Comparison Tool UI, passing selection1Data
    // Example: navigate('/comparison-tool', { state: { selection1: selection1Data } });
    navigate('/comparison', { state: { selection1: selection1Data } });
  };
  // --- End of ADDED Handler ---

  return (
    <div className={styles.resultsPageLayout}> {/* New overall layout class */} 
      {/* Back Button Container - Added */}
      <div className={styles.backButtonContainer}>
        <button onClick={handleBackClick} className={styles.backButton}>
          &larr; Back to Data Input
        </button>
      </div>
      
      {/* Control Area */}
      <div className={styles.controlArea}>
        <h3>Scope & Settings</h3>
        <div className={styles.leftPanel}>
          <ScopeSelector selectedScope={selectedScope} onScopeChange={handleScopeChange} />
          
          {/* --- ADDED: Comparison Buttons --- */}
          <Box sx={{ pt: 1, pb: 1, display: 'flex', flexDirection: 'column', gap: 1, alignItems: 'flex-start' }}>
            {selectedScope === 'whole_tissue' && jobStatus?.results?.inputs?.files?.spatialFileId && jobStatus?.results?.inputs?.mappings?.spatialMapping && (
                <Button onClick={handleInitiateComparison}>
                    Compare Whole Tissue
                </Button>
            )}
            {selectedScope === 'layers' && selectedLayer && jobStatus?.results?.inputs?.files?.spatialFileId && jobStatus?.results?.inputs?.mappings?.spatialMapping && (
                <Button onClick={handleInitiateComparison}>
                    Compare Layer: {selectedLayer}
                </Button>
            )}
            {selectedScope === 'custom' && lassoCoords && lassoCoords.length > 0 && jobStatus?.results?.inputs?.files?.spatialFileId && jobStatus?.results?.inputs?.mappings?.spatialMapping && (
                 <Button onClick={handleInitiateComparison}>
                    Compare Custom Selection
                </Button>
            )}
          </Box>
          {/* --- End of ADDED Comparison Buttons --- */}
          
          {selectedScope === 'layers' && (
            <div className={styles.layerSelectorContainer}> {/* Wrapper Div */} 
                <LayerSelector 
                    availableLayers={availableLayers} 
                    // MODIFIED: Adapt single layer state to expected props
                    selectedLayers={selectedLayer ? [selectedLayer] : []} 
                    onLayersChange={(newLayers: string[]) => {
                        // Extract first layer or null and call our single-layer handler
                        handleLayerChange(newLayers.length > 0 ? newLayers[0] : null);
                    }} 
                />
            </div>
          )}
          
          {/* Conditionally render search bars in Control Area for non-custom scopes */}
          {selectedScope !== 'custom' && searchBarInputs}
          
           {/* Add other filters/settings here later */}
        </div>
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
                  {/* Render search bars here for custom scope */}
                  {searchBarInputs}
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
                                combinedData={filteredAnalysisData} // MODIFIED: Pass filtered data
                                onSelectPair={handleSelectPair}
                                apiScopeName={null} // No API scope for custom results display
                                currentScope={selectedScope} 
                                // Pass custom analysis state and handlers
                                customAnalysisResults={customAnalysisResults} 
                                isLoadingCustomAnalysis={isLoadingCustomAnalysis}
                                customAnalysisError={customAnalysisError}
                                // ADDED: Pass aggregation state and setter
                                customAggregationLevel={customAggregationLevel}
                                setCustomAggregationLevel={setCustomAggregationLevel}
                                lassoCoords={lassoCoords} // ADDED: Pass lassoCoords
                                layerBoundaries={layerBoundaries}
                                isTableDataLoading={false} // Table data isn't loading from jobStatus here
                                displayedVizPair={displayedVizPair}
                                displayedVizData={displayedVizData}
                                isLoadingDisplayedViz={isLoadingDisplayedViz}
                                displayedVizError={displayedVizError}
                                // ADDED: Pass search terms for custom results filtering
                                ligandSearchTerm={ligandSearchTerm}
                                receptorSearchTerm={receptorSearchTerm}
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
                        combinedData={filteredAnalysisData}
                        onSelectPair={handleSelectPair}
                        apiScopeName={scopeForApi}
                        currentScope={selectedScope} 
                        // Pass null for custom state when not in custom scope
                        customAnalysisResults={null} 
                        // Pass calculated table loading state
                        isTableDataLoading={isTableDataLoading}
                        // REMOVED: isLoadingCustomAnalysis/customAnalysisError (use table loading)
                        customAnalysisError={null}
                        // ADDED: Pass aggregation state and setter (even if null/default)
                        customAggregationLevel={customAggregationLevel} 
                        setCustomAggregationLevel={setCustomAggregationLevel}
                        lassoCoords={null} // ADDED: Pass null lassoCoords for non-custom scopes
                        layerBoundaries={layerBoundaries}
                        // ADDED: Pass decoupled viz state
                        displayedVizPair={displayedVizPair}
                        displayedVizData={displayedVizData}
                        isLoadingDisplayedViz={isLoadingDisplayedViz}
                        displayedVizError={displayedVizError}
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

// Added named export for potential use in tests or other components, and default export
export { ResultsPage };
export default ResultsPage; 