import React, { useState, useCallback, useMemo, useEffect } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import styles from './ComparisonToolPage.module.css';
import { Box, Typography, Paper, Button, CircularProgress, TextField, Alert } from '@mui/material';
import InfoIcon from '@mui/icons-material/Info';

// Import hooks and components needed for Selection 2 (same project)
import useJobStatus from '../../hooks/useJobStatus'; // To get details of selection1's project
import ScopeSelector, { ScopeType } from '../../components/ScopeSelector/ScopeSelector';
import LayerSelector from '../../components/LayerSelector/LayerSelector';
import SpatialOverviewVisualization from '../../components/SpatialOverviewVisualization/SpatialOverviewVisualization';
import InteractionVisualization from '../../components/InteractionVisualization/InteractionVisualization';
import useInteractionData from '../../hooks/useInteractionData'; // InteractionVisualizationData might be local to this file or a shared type

import ComparisonResultsDisplay from '../../components/ComparisonResultsDisplay/ComparisonResultsDisplay'; // NEW IMPORT

// --- Revised Types for Selection Data (Keep here if used by parts of ComparisonToolPage not moved out) ---
interface FileSet {
  spatialFileId?: string;
  interactionsFileId?: string;
  modulesFileId?: string;
}

interface MappingSet {
  spatialMapping?: { geneCol?: string; xCol?: string; yCol?: string; layerCol?: string; [key: string]: any };
  interactionsMapping?: { ligandCol?: string; receptorCol?: string; [key: string]: any };
  modulesMapping?: { geneCol?: string; moduleCol?: string; [key: string]: any };
}

interface SelectionData {
  source_job_id?: string; 
  files: FileSet;
  type: 'whole_tissue' | 'layer' | 'lasso'; 
  definition: {
    layer_name?: string;
    polygon_coords?: [number, number][];
  };
  mappings: MappingSet;
}
// --- End of Revised Types ---

// --- Frontend-specific types matching backend response/request shapes (Keep here if used by parts not moved out) ---
interface DifferentialExpressionResultFE {
    molecule_id: string;
    type: string;
    log2_fold_change?: number | null;
    p_value?: number | null;
    q_value?: number | null;
    mean_selection1?: number | null;
    mean_selection2?: number | null;
    ligand_id?: string;
    receptor_id?: string;
}

// Interface for the request payload (matching backend ComparisonRequest)
interface ComparisonRequestPayload {
    comparison_name?: string | null;
    selection1: SelectionData;
    selection2: SelectionData;
    fdr_threshold: number;
}

// Interface for the job response from /api/analysis/compare
interface ComparisonJobResponseFE {
    status: string; 
    message: string;
    job_id?: string | null;
}

// --- End of Frontend Types ---

interface ComparisonLocationState {
  selection1?: SelectionData;
  newSelectionFor?: 'selection2';
  selectionData?: SelectionData;
}

const ComparisonToolPage: React.FC = () => {
  const location = useLocation();
  const navigate = useNavigate();
  const state = location.state as ComparisonLocationState | null;

  const [selection1, setSelection1] = React.useState<SelectionData | null>(null);
  const [selection2, setSelection2] = React.useState<SelectionData | null>(null);
  const [selection2ProjectOrigin, setSelection2ProjectOrigin] = React.useState<'same' | 'different' | null>(null);
  
  const [jobStatusForSelection1Project, setJobStatusForSelection1Project] = React.useState<any | null>(null);
  const [selectedScopeForSelection2, setSelectedScopeForSelection2] = React.useState<ScopeType>('whole_tissue');
  const [availableLayersForSelection2, setAvailableLayersForSelection2] = React.useState<string[]>([]);
  const [selectedLayerForSelection2, setSelectedLayerForSelection2] = React.useState<string | null>(null);
  const [lassoCoordsForSelection2, setLassoCoordsForSelection2] = React.useState<[number, number][] | null>(null);
  
  const [fdrThreshold, setFdrThreshold] = React.useState<number>(0.05);
  
  const [comparisonError, setComparisonError] = React.useState<string | null>(null); // For job submission error
  
  const [comparisonJobId, setComparisonJobId] = React.useState<string | null>(null);

  const [isStartingComparison, setIsStartingComparison] = React.useState<boolean>(false);
  const [startComparisonError, setStartComparisonError] = React.useState<string | null>(null);

  const [filterMinLog2FC, setFilterMinLog2FC] = React.useState<number>(0);
  const [filterMaxQValue, setFilterMaxQValue] = React.useState<number>(0.05);
  const [filterMinMeanValue, setFilterMinMeanValue] = React.useState<number>(0);
  const [showLigands, setShowLigands] = React.useState<boolean>(true);
  const [showReceptors, setShowReceptors] = React.useState<boolean>(true);
  const [showLRPairs, setShowLRPairs] = React.useState<boolean>(true);

  const [selectedMolecule, setSelectedMolecule] = React.useState<{
    id: string;
    type: string;
    ligandId?: string;
    receptorId?: string;
  } | null>(null);

  const [isSelection1Fullscreen, setIsSelection1Fullscreen] = useState<boolean>(false);
  const [isSelection2Fullscreen, setIsSelection2Fullscreen] = useState<boolean>(false);

  const {
    jobStatus: comparisonJobDetails,
    isLoading: isLoadingComparisonJobStatus, 
    error: comparisonJobErrorFromHook, // Renamed to avoid conflict with submission error
  } = useJobStatus(comparisonJobId);
  
  // Combine submission error and hook error if needed, or decide which one to prioritize
  const displayableComparisonJobError = comparisonJobErrorFromHook || comparisonError;

  const selectedPairMemo = React.useMemo(() => {
    if (!selectedMolecule) return null;
    return [
      selectedMolecule.type === 'ligand_receptor_pair' ? selectedMolecule.ligandId || '' : selectedMolecule.id, 
      selectedMolecule.type === 'ligand_receptor_pair' ? selectedMolecule.receptorId || '' : selectedMolecule.id
    ] as [string, string];
  }, [selectedMolecule]);
  
  const selection1Scope = React.useMemo(() => selection1?.type === 'layer' ? 'layers' : 'whole_tissue', [selection1]);
  const selection1Polygon = React.useMemo(() => selection1?.type === 'lasso' ? selection1.definition.polygon_coords || null : null, [selection1]);
  const selection2Scope = React.useMemo(() => selection2?.type === 'layer' ? 'layers' : 'whole_tissue', [selection2]);
  const selection2Polygon = React.useMemo(() => selection2?.type === 'lasso' ? selection2.definition.polygon_coords || null : null, [selection2]);
  
  const { 
    interactionVizData: selection1VizData, 
    isLoading: isLoadingSelection1Viz, 
    error: selection1VizError,
    cancelFetch: cancelSelection1VizFetch
  } = useInteractionData(selection1?.source_job_id || null, selectedPairMemo, selection1Scope, selection1Polygon);
  
  const { 
    interactionVizData: selection2VizData, 
    isLoading: isLoadingSelection2Viz, 
    error: selection2VizError,
    cancelFetch: cancelSelection2VizFetch
  } = useInteractionData(selection2?.source_job_id || null, selectedPairMemo, selection2Scope, selection2Polygon);
  
  const selection1VizDataMemo = React.useMemo(() => {
    if (!selection1VizData || !selectedMolecule) return null;
    return {
      ligand: selection1VizData.ligand || [],
      receptor: selection1VizData.receptor.map((p: { layer?: string; x: number; y: number }) => ({ ...p, gene: p.layer || '' })) || [],
      isComplex: true,
      receptorName: selectedMolecule.type === 'ligand_receptor_pair' ? selectedMolecule.receptorId || '' : selectedMolecule.id
    };
  }, [selection1VizData, selectedMolecule]);
  
  const selection2VizDataMemo = React.useMemo(() => {
    if (!selection2VizData || !selectedMolecule) return null;
    return {
      ligand: selection2VizData.ligand || [],
      receptor: selection2VizData.receptor.map((p: { layer?: string; x: number; y: number }) => ({ ...p, gene: p.layer || '' })) || [],
      isComplex: true,
      receptorName: selectedMolecule.type === 'ligand_receptor_pair' ? selectedMolecule.receptorId || '' : selectedMolecule.id
    };
  }, [selection2VizData, selectedMolecule]);

  const renderSelection1Visualization = React.useMemo(() => {
    if (!selectedMolecule || !selection1VizDataMemo) return null;
    return (
      <InteractionVisualization
        data={selection1VizDataMemo}
        ligandName={selectedMolecule.type === 'ligand_receptor_pair' ? selectedMolecule.ligandId || '' : selectedMolecule.id}
        currentScope={selection1Scope}
        isLoading={isLoadingSelection1Viz}
        cancelFetch={cancelSelection1VizFetch}
        layerBoundaries={undefined}
      />
    );
  }, [selection1VizDataMemo, selectedMolecule, selection1Scope, isLoadingSelection1Viz, cancelSelection1VizFetch]);

  const renderSelection2Visualization = React.useMemo(() => {
    if (!selectedMolecule || !selection2VizDataMemo) return null;
    return (
      <InteractionVisualization
        data={selection2VizDataMemo}
        ligandName={selectedMolecule.type === 'ligand_receptor_pair' ? selectedMolecule.ligandId || '' : selectedMolecule.id}
        currentScope={selection2Scope}
        isLoading={isLoadingSelection2Viz}
        cancelFetch={cancelSelection2VizFetch}
        layerBoundaries={undefined}
      />
    );
  }, [selection2VizDataMemo, selectedMolecule, selection2Scope, isLoadingSelection2Viz, cancelSelection2VizFetch]);

  const renderSelection1Loading = React.useMemo(() => {
    if (!isLoadingSelection1Viz) return null;
    return (<Box sx={{ display: 'flex', justifyContent: 'center', p: 4 }}><CircularProgress /></Box>);
  }, [isLoadingSelection1Viz]);

  const renderSelection2Loading = React.useMemo(() => {
    if (!isLoadingSelection2Viz) return null;
    return (<Box sx={{ display: 'flex', justifyContent: 'center', p: 4 }}><CircularProgress /></Box>);
  }, [isLoadingSelection2Viz]);
  
  const toggleSelection1Fullscreen = useCallback(() => setIsSelection1Fullscreen(!isSelection1Fullscreen), [isSelection1Fullscreen]);
  const toggleSelection2Fullscreen = useCallback(() => setIsSelection2Fullscreen(!isSelection2Fullscreen), [isSelection2Fullscreen]);

  // Helper function to format selection definition for display (REMAINS HERE)
  const formatSelectionDefinition = useCallback((selectionType: string, definition: any): string => {
    if (selectionType === 'whole_tissue' || (!definition || Object.keys(definition).length === 0)) return 'Whole Tissue';
    if (selectionType === 'layer' && definition.layer_name) return `Layer: ${definition.layer_name}`;
    if (selectionType === 'lasso' && definition.polygon_coords) return `Custom Selection (${definition.polygon_coords.length} points)`;
    return selectionType;
  }, []);

  React.useEffect(() => {
    if (state?.selection1 && !selection1) setSelection1(state.selection1);
    if (state?.newSelectionFor === 'selection2' && state?.selectionData) {
      setSelection2(state.selectionData);
      setSelection2ProjectOrigin(null);
    }
  }, [state, selection1, navigate, location.pathname]);

  const { jobStatus: currentProjectJobStatus, isLoading: isLoadingCurrentProjectJobStatus, error: currentProjectJobStatusError } = useJobStatus(selection1?.source_job_id && selection2ProjectOrigin === 'same' ? selection1.source_job_id : null);

  React.useEffect(() => {
    if (selection2ProjectOrigin === 'same' && currentProjectJobStatus) {
      setJobStatusForSelection1Project(currentProjectJobStatus);
      const outputs = currentProjectJobStatus?.results?.outputs;
      setAvailableLayersForSelection2(outputs ? Object.keys(outputs).filter(k => k !== 'whole_tissue' && k !== 'layer_boundaries') : []);
    } else {
      setJobStatusForSelection1Project(null);
      setAvailableLayersForSelection2([]);
      setSelectedLayerForSelection2(null);
    }
  }, [selection2ProjectOrigin, currentProjectJobStatus]);

  React.useEffect(() => {
    if (selection2ProjectOrigin === 'same' && selectedScopeForSelection2 === 'layers') {
      if (!selectedLayerForSelection2 || !availableLayersForSelection2.includes(selectedLayerForSelection2)) {
        setSelectedLayerForSelection2(availableLayersForSelection2.length > 0 ? availableLayersForSelection2[0] : null);
      }
    } else {
      if (selectedLayerForSelection2 !== null) setSelectedLayerForSelection2(null);
    }
  }, [selectedScopeForSelection2, availableLayersForSelection2, selectedLayerForSelection2, selection2ProjectOrigin]);

  const handleBackToResults = () => {
    if (selection1?.source_job_id) navigate(`/analysis/${selection1.source_job_id}`);
    else navigate('/');
  };

  const handleChooseSelection2Origin = (origin: 'same' | 'different') => {
    setSelection2ProjectOrigin(origin);
    if (origin === 'different') console.log("Select Different Project for Selection 2 - Marked as Coming Soon");
  };

  const handleScopeChangeForSelection2 = (scope: ScopeType) => {
    setSelectedScopeForSelection2(scope);
    if (scope !== 'layers') setSelectedLayerForSelection2(null);
    if (scope !== 'custom') setLassoCoordsForSelection2(null);
  };

  const handleLayerChangeForSelection2 = (newLayers: string[]) => {
    setSelectedLayerForSelection2(newLayers.length > 0 ? newLayers[0] : null);
  };

  const handleLassoSelectForSelection2 = React.useCallback((coords: [number, number][] | null) => {
    console.log("[ComparisonToolPage] Lasso for Selection 2 coords:", coords);
    setLassoCoordsForSelection2(coords);
  }, []);

  const handleConfirmSelection2SameProject = () => {
    if (!selection1) return;
    let definitionForSel2: SelectionData['definition'] = {};
    if (selectedScopeForSelection2 === 'layers') {
      if (!selectedLayerForSelection2) { alert("Please select a layer for Selection 2."); return; }
      definitionForSel2 = { layer_name: selectedLayerForSelection2 };
    } else if (selectedScopeForSelection2 === 'custom') {
      if (!lassoCoordsForSelection2 || lassoCoordsForSelection2.length === 0) console.warn("Proceeding with empty custom selection for Selection 2 (dev placeholder)")
      definitionForSel2 = { polygon_coords: lassoCoordsForSelection2 || [] };
    }
    const sel2Data: SelectionData = {
      source_job_id: selection1.source_job_id,
      files: selection1.files,
      type: selectedScopeForSelection2 === 'layers' ? 'layer' : selectedScopeForSelection2 === 'custom' ? 'lasso' : 'whole_tissue',
      definition: definitionForSel2,
      mappings: selection1.mappings
    };
    setSelection2(sel2Data);
  };

  const handleRunComparison = async () => {
    if (!selection1 || !selection2) { alert("Both selections must be defined."); return; }
    setComparisonJobId(null);
    setStartComparisonError(null);
    setComparisonError(null); // Clear previous hook error too
    setIsStartingComparison(true);
    const requestPayload: ComparisonRequestPayload = { selection1, selection2, fdr_threshold: fdrThreshold };
    try {
      const response = await fetch('/api/analysis/compare', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(requestPayload),
      });
      const responseData: ComparisonJobResponseFE = await response.json(); 
      if (!response.ok) {
        if (response.status === 422) console.error("FastAPI Validation Error Details:", responseData);
        throw new Error(responseData.message || `Failed to start job (status ${response.status})`);
      }
      setComparisonJobId(responseData.job_id ?? null);
    } catch (error) {
      console.error("Error starting comparison job:", error);
      const errorMsg = error instanceof Error ? error.message : String(error);
      setStartComparisonError(errorMsg); // Use for immediate feedback on submission
      setComparisonError(errorMsg); // Also set general comparison error
      setComparisonJobId(null);
    } finally {
      setIsStartingComparison(false);
    }
  };

  const handleCloseVisualization = () => setSelectedMolecule(null);

  // All hooks are now declared above this point.
  // Now, we can have the early return.
  // if (!selection1 || !selection1.files || !selection1.mappings) {
  //   return (<Paper><Typography>Loading or data missing...</Typography></Paper>);
  // }

  const handleRowClick = (item: DifferentialExpressionResultFE) => {
    console.log("ComparisonToolPage - Row clicked:", item);
    setSelectedMolecule(item.type === 'ligand_receptor_pair' ? 
      { id: `${item.ligand_id}-${item.receptor_id}`, type: 'ligand_receptor_pair', ligandId: item.ligand_id, receptorId: item.receptor_id } :
      { id: item.molecule_id, type: item.type }
    );
  };
  
  useEffect(() => {
    console.log("selectedMolecule changed:", selectedMolecule);
  }, [selectedMolecule]);
  
  return (
    <Paper elevation={2} className={styles.mainContainer}>
      {/* Conditionally render content based on selection1 being ready */}
      {selection1 && selection1.files && selection1.mappings ? (
        <>
          <Typography variant="h4" gutterBottom className={styles.pageTitle}>
            Comparison Configuration
          </Typography>

          {/* Display Selection 1 Box */}
          <Box className={styles.selectionBox} sx={{ mb: 4 }}>
            <Typography variant="h6" gutterBottom>Selection 1</Typography>
            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
              <Typography><strong>Source:</strong> {selection1.source_job_id || 'N/A'}</Typography>
              <Typography><strong>Type:</strong> {formatSelectionDefinition(selection1.type, selection1.definition)}</Typography>
            </Box>
          </Box>

          {!selection2 && !selection2ProjectOrigin && (
            <Box className={styles.selectionOriginBox} sx={{ mb: 4 }}>
              <Typography variant="h6" gutterBottom>Define Selection 2</Typography>
              <Button variant="contained" onClick={() => handleChooseSelection2Origin('same')} sx={{ mr: 2 }} className={styles.actionButton}>
                Use Same Project for Selection 2
              </Button>
              <Button variant="contained" onClick={() => handleChooseSelection2Origin('different')} className={styles.actionButton}>
                Select Different Project for Selection 2
              </Button>
            </Box>
          )}

          {selection2ProjectOrigin === 'different' && !selection2 && (
            <Box className={styles.selectionBox} sx={{ mb: 4, p:3, textAlign: 'center' }}>
              <Typography variant="h6" gutterBottom>Select Different Project</Typography>
              <Typography sx={{ fontStyle: 'italic', color: 'gray' }}>Feature coming soon!</Typography>
              <Button variant="outlined" onClick={() => setSelection2ProjectOrigin(null)} sx={{mt: 2}}>
                Back to Selection 2 Choice
              </Button>
            </Box>
          )}

          {selection2ProjectOrigin === 'same' && !selection2 && (
            <Box className={styles.selectionBox} sx={{ mb: 4 }}>
              <Typography variant="h6" gutterBottom>Define Selection 2 (From Project: {selection1.source_job_id})</Typography>
              {isLoadingCurrentProjectJobStatus && <CircularProgress size={24} sx={{mr: 1}}/>}
              {currentProjectJobStatusError && <Typography color="error">Error loading project details: {currentProjectJobStatusError}</Typography>}
              {jobStatusForSelection1Project && (
                <>
                  <ScopeSelector selectedScope={selectedScopeForSelection2} onScopeChange={handleScopeChangeForSelection2} />
                  {selectedScopeForSelection2 === 'layers' && (
                    <Box sx={{mt: 2}} className={styles.layerSelectorContainer}>
                      <LayerSelector availableLayers={availableLayersForSelection2} selectedLayers={selectedLayerForSelection2 ? [selectedLayerForSelection2] : []} onLayersChange={handleLayerChangeForSelection2} />
                      {availableLayersForSelection2.length === 0 && <Typography sx={{mt:1, fontSize: '0.9em'}}>No layers found.</Typography>}
                    </Box>
                  )}
                  {selectedScopeForSelection2 === 'custom' && (
                    <Box sx={{mt: 2}}>
                      <Typography sx={{ fontStyle: 'italic', color: 'gray', mb:1 }}>Draw on spatial overview for Selection 2.</Typography>
                      {selection1?.source_job_id && (
                        <Box className={styles.spatialVizContainerForSel2}>
                            <SpatialOverviewVisualization jobId={selection1.source_job_id} onLassoSelect={handleLassoSelectForSelection2} showAnalyzeButton={false} />
                        </Box>
                      )}
                      {!selection1?.source_job_id && <Typography color="error">Source Job ID missing.</Typography>}
                    </Box>
                  )}
                  <Button variant="outlined" sx={{ mt: 2 }} onClick={handleConfirmSelection2SameProject}>Confirm Selection 2</Button>
                </>
              )}
            </Box>
          )}

          {selection2 && (
            <Box className={styles.selectionBox} sx={{ mb: 4, backgroundColor: '#e3f2fd' }}>
              <Typography variant="h6" gutterBottom>Selection 2</Typography>
              <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
                <Typography><strong>Source:</strong> {selection2.source_job_id || 'N/A'}</Typography>
                <Typography><strong>Type:</strong> {formatSelectionDefinition(selection2.type, selection2.definition)}</Typography>
              </Box>
            </Box>
          )}

          {selection1 && selection2 && (
            <Box sx={{ mt: 4, pt: 2, borderTop: '1px solid #ccc' }} className={styles.analysisConfigSection}>
              <Typography variant="h6" gutterBottom>Comparison Parameters</Typography>
              <TextField 
                label="Significance Threshold (FDR)" type="number" value={fdrThreshold}
                onChange={(e) => { const val = parseFloat(e.target.value); if (!isNaN(val) && val >= 0 && val <= 1) setFdrThreshold(val); }}
                inputProps={{ step: 0.01, min: 0, max: 1 }} sx={{ maxWidth: 300, mb: 3 }} variant="outlined" size="small"
                helperText="False Discovery Rate threshold (e.g., 0.05)."
              />
              <Button variant="contained" color="primary" size="large" className={styles.actionButton} onClick={handleRunComparison} 
                disabled={isStartingComparison || isLoadingComparisonJobStatus || (comparisonJobDetails?.status === 'running' || comparisonJobDetails?.status === 'pending')}>
                {(isStartingComparison || isLoadingComparisonJobStatus || comparisonJobDetails?.status === 'running' || comparisonJobDetails?.status === 'pending') 
                  ? <CircularProgress size={24} color="inherit" /> : 'Run Comparison'} 
              </Button>
              {startComparisonError && (<Alert severity="error" sx={{ mt: 2 }}>Failed to start: {startComparisonError}</Alert>)}
            </Box>
          )}

          {/* Results Display Component */}
          {comparisonJobId && selection1 && selection2 && (
            <ComparisonResultsDisplay
              comparisonJobId={comparisonJobId}
              comparisonJobDetails={comparisonJobDetails}
              isLoadingComparisonJobStatus={isLoadingComparisonJobStatus}
              comparisonJobError={displayableComparisonJobError}
              fdrThreshold={fdrThreshold}
              selection1={selection1}
              selection2={selection2}
              filterMinLog2FC={filterMinLog2FC}
              setFilterMinLog2FC={setFilterMinLog2FC}
              filterMaxQValue={filterMaxQValue}
              setFilterMaxQValue={setFilterMaxQValue}
              filterMinMeanValue={filterMinMeanValue}
              setFilterMinMeanValue={setFilterMinMeanValue}
              showLigands={showLigands}
              setShowLigands={setShowLigands}
              showReceptors={showReceptors}
              setShowReceptors={setShowReceptors}
              showLRPairs={showLRPairs}
              setShowLRPairs={setShowLRPairs}
              selectedMolecule={selectedMolecule}
              handleRowClick={handleRowClick}
              formatSelectionDefinition={formatSelectionDefinition}
              styles={styles}
            />
          )}

          {/* Interaction Visualization Section */}
          {selectedMolecule && comparisonJobDetails?.status === 'success' && selection1 && selection2 && (
            <Box sx={{ 
              mt: 4, p: 2, border: '1px solid #e0e0e0', borderRadius: 1, position: 'relative',
              ...(isSelection1Fullscreen || isSelection2Fullscreen ? {
                position: 'fixed', top: 0, left: 0, right: 0, bottom: 0, zIndex: 1300,
                backgroundColor: 'white', padding: 4, overflowY: 'auto'
              } : {})
            }}>
              <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
                <Typography variant="h6">
                  Interaction Visualization: {selectedMolecule.type === 'ligand_receptor_pair' ? 
                    `${selectedMolecule.ligandId} → ${selectedMolecule.receptorId}` : selectedMolecule.id}
                </Typography>
                <Box>
                  {(isSelection1Fullscreen || isSelection2Fullscreen) && (
                    <Button variant="outlined" size="small" onClick={() => { setIsSelection1Fullscreen(false); setIsSelection2Fullscreen(false); }} sx={{ mr: 1 }}>
                      Exit Fullscreen
                    </Button>
                  )}
                  <Button variant="outlined" size="small" onClick={handleCloseVisualization}>Close</Button>
                </Box>
              </Box>
              <Box sx={{ display: 'flex', flexDirection: { xs: 'column', md: isSelection1Fullscreen || isSelection2Fullscreen ? 'column' : 'row' }, gap: 2 }}>
                {!isSelection2Fullscreen && (
                  <Box sx={{ flex: 1, border: '1px solid #e0e0e0', borderRadius: 1, p: 1, position: 'relative', height: isSelection1Fullscreen ? 'calc(100vh - 180px)' : 'auto' }}>
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', borderBottom: '1px solid #e0e0e0', pb: 1, mb: 1 }}>
                      <Typography variant="subtitle1">Selection 1: {formatSelectionDefinition(selection1.type, selection1.definition)}</Typography>
                      <Button variant="text" size="small" onClick={toggleSelection1Fullscreen}>{isSelection1Fullscreen ? 'Exit Fullscreen' : 'Fullscreen'}</Button>
                    </Box>
                    {selection1VizError && <Typography color="error">{selection1VizError}</Typography>}
                    {renderSelection1Visualization}
                    {renderSelection1Loading}
                  </Box>
                )}
                {!isSelection1Fullscreen && (
                  <Box sx={{ flex: 1, border: '1px solid #e0e0e0', borderRadius: 1, p: 1, position: 'relative', height: isSelection2Fullscreen ? 'calc(100vh - 180px)' : 'auto' }}>
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', borderBottom: '1px solid #e0e0e0', pb: 1, mb: 1 }}>
                      <Typography variant="subtitle1">Selection 2: {formatSelectionDefinition(selection2.type, selection2.definition)}</Typography>
                      <Button variant="text" size="small" onClick={toggleSelection2Fullscreen}>{isSelection2Fullscreen ? 'Exit Fullscreen' : 'Fullscreen'}</Button>
                    </Box>
                    {selection2VizError && <Typography color="error">{selection2VizError}</Typography>}
                    {renderSelection2Visualization}
                    {renderSelection2Loading}
                  </Box>
                )}
              </Box>
            </Box>
          )}
        
          {/* Back Button Section */}
          <Box sx={{ mt: 4, pt: 2, borderTop: '1px solid #ccc' }}>
            <Button variant="outlined" onClick={handleBackToResults} className={styles.utilityButton}>
              Back to Results Page (of Selection 1)
            </Button>
          </Box>
        </>
      ) : (
        <Typography>Loading or data missing...</Typography> // Loading state if selection1 not ready
      )}
    </Paper>
  );
}; 

export default ComparisonToolPage;
