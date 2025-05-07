import React from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import styles from './ComparisonToolPage.module.css';
import { Box, Typography, Paper, Button, CircularProgress } from '@mui/material';

// Import hooks and components needed for Selection 2 (same project)
import useJobStatus from '../../hooks/useJobStatus'; // To get details of selection1's project
import ScopeSelector, { ScopeType } from '../../components/ScopeSelector/ScopeSelector';
import LayerSelector from '../../components/LayerSelector/LayerSelector';
import SpatialOverviewVisualization from '../../components/SpatialOverviewVisualization/SpatialOverviewVisualization';

// Define a type for Selection Data (based on what's sent from ResultsPage)
// This should align with the structure in ComparisonFeature.md for selection1/selection2
interface SelectionData {
  source_job_id?: string; // Optional if not always present for every selection type
  file_id: string;
  type: 'whole_tissue' | 'layer' | 'lasso'; // Or your ScopeType
  definition: {
    layer_name?: string;
    polygon_coords?: [number, number][];
  };
  column_mappings: { // Assuming a generic object, refine if specific keys are known
    geneCol?: string;
    xCol?: string;
    yCol?: string;
    layerCol?: string;
    [key: string]: any; // Allow other mappings
  };
}

interface ComparisonLocationState {
  selection1?: SelectionData;
  newSelectionFor?: 'selection2';
  selectionData?: SelectionData; // Data coming back from DataInputPage for selection2
}

const ComparisonToolPage: React.FC = () => {
  const location = useLocation();
  const navigate = useNavigate();
  const state = location.state as ComparisonLocationState | null;

  const [selection1, setSelection1] = React.useState<SelectionData | null>(null);
  const [selection2, setSelection2] = React.useState<SelectionData | null>(null);
  const [selection2ProjectOrigin, setSelection2ProjectOrigin] = React.useState<'same' | 'different' | null>(null);
  
  // --- State for defining Selection 2 (if 'same' project) ---
  const [jobStatusForSelection1Project, setJobStatusForSelection1Project] = React.useState<any | null>(null); // Using any for now, refine with JobStatusData type from useJobStatus
  const [selectedScopeForSelection2, setSelectedScopeForSelection2] = React.useState<ScopeType>('whole_tissue');
  const [availableLayersForSelection2, setAvailableLayersForSelection2] = React.useState<string[]>([]);
  const [selectedLayerForSelection2, setSelectedLayerForSelection2] = React.useState<string | null>(null);
  const [lassoCoordsForSelection2, setLassoCoordsForSelection2] = React.useState<[number, number][] | null>(null); // Placeholder
  // TODO: Add states for loading, errors, comparison results etc.

  // --- Effect to set initial selection1 and handle data coming back for selection2 --- 
  React.useEffect(() => {
    if (state?.selection1 && !selection1) { // Only set initial selection1 once
      setSelection1(state.selection1);
    }
    // Logic to handle data coming back for selection2 from DataInputPage
    if (state?.newSelectionFor === 'selection2' && state?.selectionData) {
      setSelection2(state.selectionData);
      setSelection2ProjectOrigin(null); // Reset origin choice, selection 2 is now defined
      // Optionally clear the state from navigation to prevent re-processing on refresh if not desired
      // navigate(location.pathname, { replace: true, state: { ...state, newSelectionFor: undefined, selectionData: undefined } });
    }
  }, [state, selection1, navigate, location.pathname]);

  // --- End of Effect ---

  // --- Effect to fetch job details for selection1's project if 'same project' is chosen for selection2 ---
  // We use selection1.source_job_id to fetch the job details for the current project context
  const { jobStatus: currentProjectJobStatus, isLoading: isLoadingCurrentProjectJobStatus, error: currentProjectJobStatusError } = useJobStatus(selection1?.source_job_id && selection2ProjectOrigin === 'same' ? selection1.source_job_id : null);

  React.useEffect(() => {
    if (selection2ProjectOrigin === 'same' && currentProjectJobStatus) {
      setJobStatusForSelection1Project(currentProjectJobStatus);
      const outputs = currentProjectJobStatus?.results?.outputs;
      if (outputs) {
        setAvailableLayersForSelection2(Object.keys(outputs).filter(k => k !== 'whole_tissue' && k !== 'layer_boundaries'));
      } else {
        setAvailableLayersForSelection2([]);
      }
    } else {
      // Reset if not in 'same project' mode or if selection1 is not yet defined
      setJobStatusForSelection1Project(null);
      setAvailableLayersForSelection2([]);
      setSelectedLayerForSelection2(null);
    }
  }, [selection2ProjectOrigin, currentProjectJobStatus]);

  // --- Effect to ensure a layer is selected for Selection 2 when in 'layers' scope and available layers change ---
  React.useEffect(() => {
    if (selection2ProjectOrigin === 'same' && selectedScopeForSelection2 === 'layers') {
      if (!selectedLayerForSelection2 || !availableLayersForSelection2.includes(selectedLayerForSelection2)) {
        setSelectedLayerForSelection2(availableLayersForSelection2.length > 0 ? availableLayersForSelection2[0] : null);
      }
    } else {
      if (selectedLayerForSelection2 !== null) {
        setSelectedLayerForSelection2(null); // Clear if not in 'layers' scope for sel2
      }
    }
  }, [selectedScopeForSelection2, availableLayersForSelection2, selectedLayerForSelection2, selection2ProjectOrigin]);
  // --- End of Effects ---

  const handleBackToResults = () => {
    if (selection1?.source_job_id) {
      navigate(`/analysis/${selection1.source_job_id}`);
    } else {
      navigate('/'); // Fallback to data input or another sensible page
    }
  };

  const handleChooseSelection2Origin = (origin: 'same' | 'different') => {
    setSelection2ProjectOrigin(origin);
    if (origin === 'different') {
      // Navigate to DataInputPage to select/load a different project for selection 2
      // Pass current selection1 so DataInputPage knows where to return and potentially other context.
      // const returnState = {
      //   returnTo: '/comparison',
      //   selectionPurpose: 'selection2',
      //   currentSelection1: selection1 // Pass selection1 along
      // };
      // navigate('/data-input', { state: returnState });
      // For now, mark as coming soon
      console.log("Select Different Project for Selection 2 - Marked as Coming Soon");
    }
  };

  // --- Callback Handlers for Selection 2 Scope/Layer changes ---
  const handleScopeChangeForSelection2 = (scope: ScopeType) => {
    setSelectedScopeForSelection2(scope);
    // Potentially clear layer/lasso for selection 2 if scope changes
    if (scope !== 'layers') {
      setSelectedLayerForSelection2(null);
    }
    if (scope !== 'custom') {
      setLassoCoordsForSelection2(null);
    }
  };

  const handleLayerChangeForSelection2 = (newLayers: string[]) => {
    // Assuming LayerSelector might return multiple, but we only use one for now
    setSelectedLayerForSelection2(newLayers.length > 0 ? newLayers[0] : null);
  };

  const handleLassoSelectForSelection2 = React.useCallback((coords: [number, number][] | null) => {
    console.log("[ComparisonToolPage] Lasso for Selection 2 coords:", coords);
    setLassoCoordsForSelection2(coords);
  }, []);

  const handleConfirmSelection2SameProject = () => {
    if (!selection1) return; // Should not happen if this UI is visible

    let definitionForSel2: SelectionData['definition'] = {};
    if (selectedScopeForSelection2 === 'layers') {
      if (!selectedLayerForSelection2) {
        alert("Please select a layer for Selection 2.");
        return;
      }
      definitionForSel2 = { layer_name: selectedLayerForSelection2 };
    } else if (selectedScopeForSelection2 === 'custom') {
      if (!lassoCoordsForSelection2 || lassoCoordsForSelection2.length === 0) {
        alert("Please define a lasso selection for Selection 2.");
        // For now, we'll allow proceeding without lasso for placeholder UI
        // return;
        console.warn("Proceeding with empty custom selection for Selection 2 (dev placeholder)")
      }
      definitionForSel2 = { polygon_coords: lassoCoordsForSelection2 || [] };
    }
    // For 'whole_tissue', definition remains {}

    const sel2Data: SelectionData = {
      source_job_id: selection1.source_job_id, // Same project
      file_id: selection1.file_id,             // Same file
      type: selectedScopeForSelection2 === 'layers' 
            ? 'layer' 
            : selectedScopeForSelection2 === 'custom' 
            ? 'lasso' 
            : 'whole_tissue',
      definition: definitionForSel2,
      column_mappings: selection1.column_mappings // Same mappings
    };
    setSelection2(sel2Data);
    // After confirming, we might want to hide the selection2ProjectOrigin choice section
    // or change the UI flow. For now, setting selection2 will make its display box appear.
  };
  // --- End of Callbacks ---

  if (!selection1) {
    // This might happen if navigated directly without state or state is lost.
    // Could show a message and a button to go back or select a primary dataset.
    return (
      <Paper className={styles.container} elevation={3}>
        <Typography variant="h5" gutterBottom>Comparison Tool</Typography>
        <Typography>
          Selection 1 data is missing. Please initiate a comparison from a results page.
        </Typography>
        <Button variant="outlined" onClick={() => navigate('/')} sx={{ mt: 2 }}>
          Go to Data Input
        </Button>
      </Paper>
    );
  }

  return (
    <Paper className={styles.container} elevation={3}>
      <Typography variant="h4" gutterBottom className={styles.pageTitle}>
        Comparison Configuration
      </Typography>

      {/* Display Selection 1 */}
      <Box className={styles.selectionBox} sx={{ mb: 4 }}>
        <Typography variant="h6" gutterBottom>Selection 1 (From Job: {selection1.source_job_id || 'N/A'})</Typography>
        <Typography><strong>File ID:</strong> {selection1.file_id}</Typography>
        <Typography><strong>Type:</strong> {selection1.type}</Typography>
        <Typography><strong>Definition:</strong></Typography>
        <pre className={styles.preformattedText}>{JSON.stringify(selection1.definition, null, 2)}</pre>
        <Typography><strong>Column Mappings:</strong></Typography>
        <pre className={styles.preformattedText}>{JSON.stringify(selection1.column_mappings, null, 2)}</pre>
      </Box>

      {/* Choose Origin for Selection 2 */}
      {!selection2 && !selection2ProjectOrigin && (
        <Box className={styles.selectionOriginBox} sx={{ mb: 4 }}>
          <Typography variant="h6" gutterBottom>Define Selection 2</Typography>
          <Button 
            variant="contained" 
            onClick={() => handleChooseSelection2Origin('same')} 
            sx={{ mr: 2 }}
            className={styles.actionButton}
          >
            Use Same Project for Selection 2
          </Button>
          <Button 
            variant="contained" 
            onClick={() => handleChooseSelection2Origin('different')}
            className={styles.actionButton}
          >
            Select Different Project for Selection 2
          </Button>
        </Box>
      )}

      {/* UI for "Different Project" if chosen - Display Coming Soon */}
      {selection2ProjectOrigin === 'different' && !selection2 && (
        <Box className={styles.selectionBox} sx={{ mb: 4, p:3, textAlign: 'center' }}>
           <Typography variant="h6" gutterBottom>Select Different Project</Typography>
           <Typography sx={{ fontStyle: 'italic', color: 'gray' }}>
            This feature (selecting a different project for the second comparison selection) is coming soon!
          </Typography>
           <Button variant="outlined" onClick={() => setSelection2ProjectOrigin(null)} sx={{mt: 2}}>
            Back to Selection 2 Choice
          </Button>
        </Box>
      )}

      {/* UI for defining Selection 2 if 'same' project is chosen */}
      {selection2ProjectOrigin === 'same' && !selection2 && (
        <Box className={styles.selectionBox} sx={{ mb: 4 }}>
          <Typography variant="h6" gutterBottom>Define Selection 2 (From Project: {selection1.source_job_id})</Typography>
          
          {isLoadingCurrentProjectJobStatus && <CircularProgress size={24} sx={{mr: 1}}/>}
          {currentProjectJobStatusError && <Typography color="error">Error loading project details: {currentProjectJobStatusError}</Typography>}
          
          {jobStatusForSelection1Project && (
            <>
              <ScopeSelector 
                selectedScope={selectedScopeForSelection2} 
                onScopeChange={handleScopeChangeForSelection2} 
              />

              {selectedScopeForSelection2 === 'layers' && (
                <Box sx={{mt: 2}} className={styles.layerSelectorContainer}>
                  <LayerSelector 
                    availableLayers={availableLayersForSelection2}
                    selectedLayers={selectedLayerForSelection2 ? [selectedLayerForSelection2] : []}
                    onLayersChange={handleLayerChangeForSelection2}
                    // Removed multiSelect if LayerSelector supports a direct single select mode or if we adapt
                  />
                  {availableLayersForSelection2.length === 0 && <Typography sx={{mt:1, fontSize: '0.9em'}}>No layers found in this project.</Typography>}
                </Box>
              )}

              {selectedScopeForSelection2 === 'custom' && (
                <Box sx={{mt: 2}}>
                  <Typography sx={{ fontStyle: 'italic', color: 'gray', mb:1 }}>
                    Draw a lasso on the spatial overview below to define your custom selection for Selection 2.
                  </Typography>
                  {selection1?.source_job_id && (
                    <Box className={styles.spatialVizContainerForSel2}> {/* Optional: Add specific styling for container */}
                        <SpatialOverviewVisualization 
                            jobId={selection1.source_job_id} // Use current project's job ID
                            onLassoSelect={handleLassoSelectForSelection2}
                            // onAnalyzeSelection can be a no-op or adapted if needed later for live feedback
                            onAnalyzeSelection={() => console.log("Analyze Selection clicked in Comparison (Sel2) - No-op for now")}
                            showAnalyzeButton={false}
                            // We might need to pass layerBoundaries if relevant for this viz context
                            // layerBoundaries={jobStatusForSelection1Project?.results?.outputs?.layer_boundaries}
                        />
                    </Box>
                  )}
                  {!selection1?.source_job_id && <Typography color="error">Source Job ID for spatial view is missing.</Typography>}
                </Box>
              )}

              <Button variant="outlined" sx={{ mt: 2 }} onClick={handleConfirmSelection2SameProject}>
                Confirm Selection 2 (Same Project)
              </Button>
            </>
          )}
        </Box>
      )}
      
      {/* Display Selection 2 if defined */}
      {selection2 && (
         <Box className={styles.selectionBox} sx={{ mb: 4, backgroundColor: '#e3f2fd' /* Light blue to differentiate */ }}>
          <Typography variant="h6" gutterBottom>Selection 2 (From Job: {selection2.source_job_id || 'N/A'})</Typography>
          <Typography><strong>File ID:</strong> {selection2.file_id}</Typography>
          <Typography><strong>Type:</strong> {selection2.type}</Typography>
          <Typography><strong>Definition:</strong></Typography>
          <pre className={styles.preformattedText}>{JSON.stringify(selection2.definition, null, 2)}</pre>
          <Typography><strong>Column Mappings:</strong></Typography>
          <pre className={styles.preformattedText}>{JSON.stringify(selection2.column_mappings, null, 2)}</pre>
        </Box>
      )}

      {/* Placeholder for Analysis Configuration and Trigger */}
      {selection1 && selection2 && (
        <Box sx={{ mt: 4, pt: 2, borderTop: '1px solid #ccc' }}>
          <Typography variant="h6" gutterBottom>Configure Analysis</Typography>
          {/* Dropdowns/inputs for analysis types, parameters (FDR, etc.) */}
          <Typography sx={{ fontStyle: 'italic', color: 'gray', mb:2 }}>
            (UI for selecting analysis types like 'differential_expression' and their parameters will go here.)
          </Typography>
          <Button 
            variant="contained" 
            color="primary" 
            size="large"
            className={styles.actionButton}
            onClick={() => alert("Run Comparison - TBD")}
          >
            Run Comparison
          </Button>
        </Box>
      )}

      <Box sx={{ mt: 4, pt: 2, borderTop: '1px solid #ccc' }}>
        <Button variant="outlined" onClick={handleBackToResults} className={styles.utilityButton}>
          Back to Results Page (of Selection 1)
        </Button>
      </Box>

    </Paper>
  );
};

export default ComparisonToolPage; 