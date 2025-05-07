import React from 'react';
import { Box, Typography, Paper, Button, CircularProgress, TextField, Alert, Table, TableBody, TableCell, TableContainer, TableHead, TableRow, Tooltip, tableCellClasses, FormGroup, FormControlLabel, Checkbox, Accordion, AccordionSummary, AccordionDetails } from '@mui/material';
import { styled } from '@mui/material/styles';
import InfoIcon from '@mui/icons-material/Info';
import ArrowUpward from '@mui/icons-material/ArrowUpward';
import ArrowDownward from '@mui/icons-material/ArrowDownward';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';

// --- Copied Type Definitions (can be moved to a types.ts file later) ---
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
// --- End of Copied Type Definitions ---

// --- Helper functions moved from ComparisonToolPage ---
const getSignificanceStars = (qValue: number | null | undefined): string => {
  if (qValue === null || qValue === undefined) return '';
  if (qValue < 0.001) return '***';
  if (qValue < 0.01) return '**';
  if (qValue < 0.05) return '*';
  return '';
};

const normalizeExpressionValue = (value: number | null | undefined): { display: string, raw: number } => {
  if (value === null || value === undefined) {
    return { display: '-', raw: 0 };
  }
  const raw = value;
  let display = '';
  const relativePercentage = (value / 1e-3) * 100;
  if (relativePercentage < 0.01) {
    display = '<0.01%';
  } else if (relativePercentage < 0.1) {
    display = relativePercentage.toFixed(2) + '%';
  } else if (relativePercentage < 10) {
    display = relativePercentage.toFixed(1) + '%';
  } else {
    display = Math.round(relativePercentage) + '%';
  }
  return { display, raw };
};

const formatQValue = (qValue: number | null | undefined): { display: string, significance: string } => {
  if (qValue === null || qValue === undefined) {
    return { display: '-', significance: '' };
  }
  let significance = '';
  if (qValue < 0.001) significance = '***';
  else if (qValue < 0.01) significance = '**';
  else if (qValue < 0.05) significance = '*';
  let display = '';
  if (qValue < 0.001) {
    display = '<0.1%';
  } else if (qValue < 0.01) {
    display = '<1%';
  } else if (qValue < 0.05) {
    display = '<5%';
  } else {
    display = (qValue * 100).toFixed(0) + '%';
  }
  return { display, significance };
};

const safeFormatLog2FC = (log2fc: number | null | undefined): { display: string, foldChange: string, isPositive: boolean } => {
  if (log2fc === null || log2fc === undefined) {
    return { display: '-', foldChange: '-', isPositive: false };
  }
  const absLog2FC = Math.abs(log2fc);
  const foldChange = Math.pow(2, absLog2FC).toFixed(1);
  const isPositive = log2fc > 0;
  return {
    display: absLog2FC.toFixed(2),
    foldChange: `${foldChange}×`,
    isPositive
  };
};

const getDisplayNameAndRole = (item: DifferentialExpressionResultFE): { displayName: string, roleName: string } => {
  if (item.type === 'ligand_receptor_pair') {
    return {
      displayName: `${item.ligand_id || '?'} → ${item.receptor_id || '?'}`,
      roleName: "L-R Pair"
    };
  } else if (item.type === 'single_ligand') {
    return {
      displayName: item.molecule_id,
      roleName: "Ligand"
    };
  } else if (item.type === 'single_receptor') {
    return {
      displayName: item.molecule_id,
      roleName: "Receptor"
    };
  } else {
    return {
      displayName: item.molecule_id,
      roleName: item.type || "Unknown"
    };
  }
};

const StyledTableCell = styled(TableCell)(({ theme }) => ({
  [`&.${tableCellClasses.head}`]: {
    // Optional: Adjust head cell styling if needed
  },
  [`&.${tableCellClasses.body}`]: {
    fontSize: 14,
    padding: '6px 10px',
  },
}));

// const Stars = ({ significance }: { significance: number }) => { // significance prop seems unused here in original code, qValue is used for stars
//   // Re-evaluate if this Stars component is still needed or if getSignificanceStars is used directly
//   return null; 
// };


// --- Props Interface ---
interface ComparisonResultsDisplayProps {
  comparisonJobId: string | null;
  comparisonJobDetails: any; // Replace 'any' with a more specific type for job details
  isLoadingComparisonJobStatus: boolean;
  comparisonJobError: string | null;
  fdrThreshold: number;
  selection1: SelectionData | null;
  selection2: SelectionData | null;
  
  filterMinLog2FC: number;
  setFilterMinLog2FC: (value: number) => void;
  filterMaxQValue: number;
  setFilterMaxQValue: (value: number) => void;
  filterMinMeanValue: number;
  setFilterMinMeanValue: (value: number) => void;
  showLigands: boolean;
  setShowLigands: (value: boolean) => void;
  showReceptors: boolean;
  setShowReceptors: (value: boolean) => void;
  showLRPairs: boolean;
  setShowLRPairs: (value: boolean) => void;
  
  selectedMolecule: { id: string; type: string; ligandId?: string; receptorId?: string } | null;
  handleRowClick: (item: DifferentialExpressionResultFE) => void;
  
  formatSelectionDefinition: (selectionType: string, definition: any) => string;
  styles: { [key: string]: string }; // For ComparisonToolPage.module.css styles
}

const ComparisonResultsDisplay: React.FC<ComparisonResultsDisplayProps> = ({
  comparisonJobId,
  comparisonJobDetails,
  isLoadingComparisonJobStatus,
  comparisonJobError,
  fdrThreshold,
  selection1,
  selection2,
  filterMinLog2FC,
  setFilterMinLog2FC,
  filterMaxQValue,
  setFilterMaxQValue,
  filterMinMeanValue,
  setFilterMinMeanValue,
  showLigands,
  setShowLigands,
  showReceptors,
  setShowReceptors,
  showLRPairs,
  setShowLRPairs,
  selectedMolecule,
  handleRowClick,
  formatSelectionDefinition,
  styles,
}) => {

  // If selection1 or selection2 is null, we might not be able to render parts of this.
  // Consider adding checks or returning a fallback UI.
  if (!selection1 || !selection2) {
      // This case should ideally be handled by the parent, 
      // or this component should show a specific loading/error state.
      // For now, if this happens, it implies an issue in prop passing.
      return <Typography>Selection data missing for results display.</Typography>;
  }
    
  return (
    <Box sx={{ mt: 4, pt: 2, borderTop: '1px solid #ccc' }} className={styles.resultsDisplaySection}>
      <Typography variant="h6" gutterBottom>Comparison Results (Job: {comparisonJobId})</Typography>
      
      {isLoadingComparisonJobStatus && (
        <Box sx={{display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100px'}}>
            <CircularProgress />
            <Typography sx={{ml: 2}}>Loading comparison status...</Typography>
        </Box>
      )}

      {!isLoadingComparisonJobStatus && (comparisonJobDetails?.status === 'pending' || comparisonJobDetails?.status === 'running') && (
        <Box sx={{display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100px'}}>
            <CircularProgress />
            <Typography sx={{ml: 2}}>{comparisonJobDetails.message || `Job is ${comparisonJobDetails.status}...`} (Progress: {comparisonJobDetails?.progress !== null && comparisonJobDetails?.progress !== undefined ? `${(comparisonJobDetails.progress * 100).toFixed(0)}%` : 'N/A'})</Typography>
        </Box>
      )}

      {comparisonJobError && (
        <Alert severity="error" sx={{mt: 2}}>
           Error tracking comparison job: {comparisonJobError}
        </Alert>
      )}
      
      {!isLoadingComparisonJobStatus && comparisonJobDetails?.status === 'failed' && (
        <Alert severity="error" sx={{mt: 2}}>
           Comparison Failed: {comparisonJobDetails.message || 'Unknown error'}
           {comparisonJobDetails.results?.errors && comparisonJobDetails.results.errors.length > 0 && (
             <ul style={{ marginTop: '8px', marginBottom: 0 }}>
               {comparisonJobDetails.results.errors.map((err: string, i: number) => <li key={i}>{err}</li>)}
             </ul>
           )}
        </Alert>
      )}

      {!isLoadingComparisonJobStatus && comparisonJobDetails?.status === 'success' && (
        <Box sx={{mt: 2}}>
          {comparisonJobDetails.results?.results?.differential_expression && comparisonJobDetails.results.results.differential_expression.length > 0 ? (
            <>
              <Accordion sx={{ mb: 3 }}>
                <AccordionSummary
                  expandIcon={<ExpandMoreIcon />}
                  aria-controls="explanation-content"
                  id="explanation-header"
                >
                  <Typography variant="h6">What do these results mean? (Click to expand)</Typography>
                </AccordionSummary>
                <AccordionDetails>
                  <Typography variant="subtitle1" gutterBottom>Understanding Molecular Comparison Results:</Typography>
                  <Box sx={{ ml: 2, mb: 2 }}>
                    <Typography variant="body1" paragraph>
                      <strong>Normalized Mean Values</strong> - Average expression level for each molecule in the selected regions, 
                      normalized to account for differences in total counts. Very small values may appear as 0.000 but have tooltips showing the actual value.
                    </Typography>
                    <Typography variant="body1" paragraph>
                      <strong>Log₂FC (Log2 Fold Change)</strong> - Shows how much expression differs between selections:
                      <Box component="ul" sx={{ mt: 1, mb: 1 }}>
                        <li>Positive values (🔼) indicate <strong>higher expression in Selection 2</strong></li>
                        <li>Negative values (🔽) indicate <strong>higher expression in Selection 1</strong></li>
                        <li>The value is shown alongside the actual fold change (e.g., "2×" means "twice as much")</li>
                      </Box>
                    </Typography>
                    <Typography variant="body1" paragraph>
                      <strong>q-value</strong> - The statistical significance (false discovery rate corrected p-value):
                      <Box component="ul" sx={{ mt: 1, mb: 1 }}>
                        <li>* (q &lt; 0.05) - Significant</li>
                        <li>** (q &lt; 0.01) - Highly significant</li>
                        <li>*** (q &lt; 0.001) - Extremely significant</li>
                      </Box>
                    </Typography>
                    <Typography variant="body1" paragraph>
                      <strong>Types</strong> - Molecules are classified as:
                      <Box component="ul" sx={{ mt: 1, mb: 1 }}>
                        <li><strong>Ligands</strong> - Signaling molecules that are released by cells</li>
                        <li><strong>Receptors</strong> - Cell surface proteins that bind to ligands</li>
                        <li><strong>L-R Pairs</strong> - Ligand-receptor interactions that may indicate cell-cell communication</li>
                      </Box>
                    </Typography>
                    <Typography variant="body1">
                      <strong>Biological Interpretation:</strong> Differential expression may suggest region-specific functions, 
                      cell type differences, or responses to environmental cues.
                    </Typography>
                  </Box>
                </AccordionDetails>
              </Accordion>

              <Box sx={{ mb: 4 }}>
                <Typography variant="h6" gutterBottom sx={{ 
                  borderBottom: '2px solid #3f51b5', 
                  pb: 1, 
                  display: 'flex', 
                  alignItems: 'center' 
                }}>
                  <span style={{ marginRight: '8px' }}>🔍</span> Top Changed Interactions
                </Typography>
                
                <Box sx={{ 
                  display: 'flex', 
                  justifyContent: 'space-between', 
                  mb: 2, 
                  p: 2, 
                  backgroundColor: '#f5f5f5', 
                  borderRadius: 1 
                }}>
                  <Box sx={{ 
                    flex: 1, 
                    p: 1, 
                    backgroundColor: 'rgba(0, 0, 255, 0.1)', 
                    borderRadius: 1, 
                    mr: 1 
                  }}>
                    <Typography variant="subtitle1" sx={{ fontWeight: 'bold' }}>
                      Selection 1: {formatSelectionDefinition(selection1.type, selection1.definition)}
                    </Typography>
                    <Typography variant="body2">
                      {selection1.source_job_id}
                    </Typography>
                  </Box>
                  <Box sx={{ 
                    flex: 1, 
                    p: 1, 
                    backgroundColor: 'rgba(255, 0, 0, 0.1)', 
                    borderRadius: 1, 
                    ml: 1 
                  }}>
                    <Typography variant="subtitle1" sx={{ fontWeight: 'bold' }}>
                       Selection 2: {formatSelectionDefinition(selection2.type, selection2.definition)}
                    </Typography>
                    <Typography variant="body2">
                      {selection2.source_job_id}
                    </Typography>
                  </Box>
                </Box>

                <Box>
                  <Typography variant="subtitle1" sx={{ 
                    mt: 2, 
                    mb: 1, 
                    fontWeight: 'bold', 
                    color: '#d32f2f',
                    display: 'flex',
                    alignItems: 'center'
                  }}>
                    <span style={{ marginRight: '8px' }}>🔺</span> More Active in Selection 2
                  </Typography>
                  <Box sx={{ 
                    mb: 3, 
                    maxHeight: '250px', 
                    overflowY: 'auto', 
                    border: '1px solid #e0e0e0', 
                    borderRadius: 1 
                  }}>
                    {comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => {
                        const absLog2FC = Math.abs(item.log2_fold_change || 0);
                        const maxMean = Math.max(item.mean_selection1 || 0, item.mean_selection2 || 0);
                        const typeFilter = 
                          (item.type === 'single_ligand' && showLigands) || 
                          (item.type === 'single_receptor' && showReceptors) || 
                          (item.type === 'ligand_receptor_pair' && showLRPairs);
                          
                        return absLog2FC >= filterMinLog2FC &&
                               (item.q_value || 1) <= filterMaxQValue &&
                               maxMean >= filterMinMeanValue &&
                               typeFilter && (item.log2_fold_change || 0) > 0; // Added condition for more active in Sel 2
                      }).sort((a: DifferentialExpressionResultFE, b: DifferentialExpressionResultFE) => (b.log2_fold_change || 0) - (a.log2_fold_change || 0))
                      .slice(0, 10)
                      .map((item: DifferentialExpressionResultFE, index: number) => {
                        const nameInfo = getDisplayNameAndRole(item);
                        const fcInfo = safeFormatLog2FC(item.log2_fold_change);
                        const stars = getSignificanceStars(item.q_value);
                        return (
                          <Box key={`${item.molecule_id}-${index}-sel2up`} sx={{ p: 1.5, borderBottom: '1px solid #e0e0e0', backgroundColor: index % 2 === 0 ? 'white' : '#fafafa', '&:hover': { backgroundColor: '#f5f5f5' } }}>
                            <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                              <Box sx={{ display: 'flex', alignItems: 'center' }}>
                                <Typography sx={{ mr: 1, color: '#666' }}>🔼</Typography>
                                <Typography sx={{ fontWeight: 'bold' }}>{nameInfo.displayName} {stars}</Typography>
                                <Typography variant="caption" sx={{ ml: 1, color: '#666' }}>({nameInfo.roleName})</Typography>
                              </Box>
                              <Typography variant="body2" sx={{ color: '#d32f2f', fontWeight: 'bold', backgroundColor: 'rgba(211, 47, 47, 0.1)', px: 1, py: 0.5, borderRadius: 1 }}>
                                {fcInfo.display} ({fcInfo.foldChange})
                              </Typography>
                            </Box>
                            <Box sx={{ display: 'flex', mt: 1, justifyContent: 'space-between' }}>
                              <Typography variant="caption" sx={{ color: '#666' }}>Mean expr: {item.mean_selection1?.toFixed(3)} → {item.mean_selection2?.toFixed(3)}</Typography>
                              <Typography variant="caption" sx={{ color: '#666' }}>q-value: {item.q_value?.toExponential(1)}</Typography>
                            </Box>
                          </Box>
                        );
                      })}
                      {comparisonJobDetails.results.results.differential_expression.filter((item: DifferentialExpressionResultFE) => (item.log2_fold_change || 0) > 0 && Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC && (item.q_value || 1) <= filterMaxQValue && Math.max(item.mean_selection1 || 0, item.mean_selection2 || 0) >= filterMinMeanValue && ((item.type === 'single_ligand' && showLigands) || (item.type === 'single_receptor' && showReceptors) || (item.type === 'ligand_receptor_pair' && showLRPairs))).length === 0 && (
                        <Box sx={{ p: 2, textAlign: 'center', color: '#666' }}><Typography>No significant upregulated molecules found in Selection 2 with current filters.</Typography></Box>
                      )}
                  </Box>
                  
                  <Typography variant="subtitle1" sx={{ mt: 3, mb: 1, fontWeight: 'bold', color: '#1976d2', display: 'flex', alignItems: 'center' }}>
                    <span style={{ marginRight: '8px' }}>🔻</span> More Active in Selection 1
                  </Typography>
                  <Box sx={{ mb: 3, maxHeight: '250px', overflowY: 'auto', border: '1px solid #e0e0e0', borderRadius: 1 }}>
                    {comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => {
                        const absLog2FC = Math.abs(item.log2_fold_change || 0);
                        const maxMean = Math.max(item.mean_selection1 || 0, item.mean_selection2 || 0);
                        const typeFilter = 
                          (item.type === 'single_ligand' && showLigands) || 
                          (item.type === 'single_receptor' && showReceptors) || 
                          (item.type === 'ligand_receptor_pair' && showLRPairs);
                        return absLog2FC >= filterMinLog2FC &&
                               (item.q_value || 1) <= filterMaxQValue &&
                               maxMean >= filterMinMeanValue &&
                               typeFilter && (item.log2_fold_change || 0) < 0; // Added condition for more active in Sel 1
                      }).sort((a: DifferentialExpressionResultFE, b: DifferentialExpressionResultFE) => (a.log2_fold_change || 0) - (b.log2_fold_change || 0))
                      .slice(0, 10)
                      .map((item: DifferentialExpressionResultFE, index: number) => {
                        const nameInfo = getDisplayNameAndRole(item);
                        const fcInfo = safeFormatLog2FC(item.log2_fold_change);
                        const stars = getSignificanceStars(item.q_value);
                        return (
                          <Box key={`${item.molecule_id}-${index}-sel1up`} sx={{ p: 1.5, borderBottom: '1px solid #e0e0e0', backgroundColor: index % 2 === 0 ? 'white' : '#fafafa', '&:hover': { backgroundColor: '#f5f5f5' } }}>
                            <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                              <Box sx={{ display: 'flex', alignItems: 'center' }}>
                                <Typography sx={{ mr: 1, color: '#666' }}>🔽</Typography>
                                <Typography sx={{ fontWeight: 'bold' }}>{nameInfo.displayName} {stars}</Typography>
                                <Typography variant="caption" sx={{ ml: 1, color: '#666' }}>({nameInfo.roleName})</Typography>
                              </Box>
                              <Typography variant="body2" sx={{ color: '#1976d2', fontWeight: 'bold', backgroundColor: 'rgba(25, 118, 210, 0.1)', px: 1, py: 0.5, borderRadius: 1 }}>
                                {fcInfo.display} ({fcInfo.foldChange})
                              </Typography>
                            </Box>
                            <Box sx={{ display: 'flex', mt: 1, justifyContent: 'space-between' }}>
                              <Typography variant="caption" sx={{ color: '#666' }}>Mean expr: {item.mean_selection1?.toFixed(3)} → {item.mean_selection2?.toFixed(3)}</Typography>
                              <Typography variant="caption" sx={{ color: '#666' }}>q-value: {item.q_value?.toExponential(1)}</Typography>
                            </Box>
                          </Box>
                        );
                      })}
                      {comparisonJobDetails.results.results.differential_expression.filter((item: DifferentialExpressionResultFE) => (item.log2_fold_change || 0) < 0 && Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC && (item.q_value || 1) <= filterMaxQValue && Math.max(item.mean_selection1 || 0, item.mean_selection2 || 0) >= filterMinMeanValue && ((item.type === 'single_ligand' && showLigands) || (item.type === 'single_receptor' && showReceptors) || (item.type === 'ligand_receptor_pair' && showLRPairs))).length === 0 && (
                        <Box sx={{ p: 2, textAlign: 'center', color: '#666' }}><Typography>No significant upregulated molecules found in Selection 1 with current filters.</Typography></Box>
                      )}
                  </Box>
                </Box>
              </Box>
              
              <Box sx={{ mb: 4 }}>
                <Typography variant="h6" gutterBottom sx={{ borderBottom: '2px solid #3f51b5', pb: 1, display: 'flex', alignItems: 'center' }}>
                  <span style={{ marginRight: '8px' }}>🧠</span> What It Might Mean
                </Typography>
                <Box sx={{ p: 2, backgroundColor: '#f5f5f5', borderRadius: 1 }}>
                  {(() => {
                    const lrPairs = comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => 
                        item.type === 'ligand_receptor_pair' && 
                        (item.q_value || 1) <= filterMaxQValue &&
                        Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC
                      );
                    const ligandsUpInSel2 = comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => 
                        item.type === 'single_ligand' && (item.log2_fold_change || 0) > 0 &&
                        (item.q_value || 1) <= filterMaxQValue && Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC
                      );
                    const ligandsUpInSel1 = comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => 
                        item.type === 'single_ligand' && (item.log2_fold_change || 0) < 0 &&
                        (item.q_value || 1) <= filterMaxQValue && Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC
                      );
                    const receptorsUpInSel2 = comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => 
                        item.type === 'single_receptor' && (item.log2_fold_change || 0) > 0 &&
                        (item.q_value || 1) <= filterMaxQValue && Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC
                      );
                    const receptorsUpInSel1 = comparisonJobDetails.results.results.differential_expression
                      .filter((item: DifferentialExpressionResultFE) => 
                        item.type === 'single_receptor' && (item.log2_fold_change || 0) < 0 &&
                        (item.q_value || 1) <= filterMaxQValue && Math.abs(item.log2_fold_change || 0) >= filterMinLog2FC
                      );
                    const topPair = [...lrPairs].sort((a, b) => Math.abs(b.log2_fold_change || 0) - Math.abs(a.log2_fold_change || 0))[0];
                    const sel1Name = formatSelectionDefinition(selection1.type, selection1.definition);
                    const sel2Name = formatSelectionDefinition(selection2.type, selection2.definition);
                    
                    return (
                      <>
                        <Typography variant="body1" sx={{ mb: 2 }}>
                          Comparison between <strong>{sel1Name}</strong> and <strong>{sel2Name}</strong> reveals:
                        </Typography>
                        <Box sx={{ mb: 2 }}>
                          <Typography variant="body1" sx={{ mb: 1, fontWeight: 'medium' }}>Cell Communication Patterns:</Typography>
                          <ul style={{ marginTop: 0, paddingLeft: '1.5rem' }}>
                            {lrPairs.length > 0 ? (
                              <li><Typography variant="body2">
                                {lrPairs.length} significantly different ligand-receptor interactions.
                                {topPair && ` Strongest: ${topPair.ligand_id} → ${topPair.receptor_id} (more active in ${(topPair.log2_fold_change || 0) > 0 ? sel2Name : sel1Name}).`}
                              </Typography></li>
                            ) : ( <li><Typography variant="body2">No significantly different L-R pairs at current thresholds.</Typography></li>)}
                            {ligandsUpInSel2.length > 0 && (<li><Typography variant="body2">{ligandsUpInSel2.length} ligands more active in {sel2Name}.</Typography></li>)}
                            {ligandsUpInSel1.length > 0 && (<li><Typography variant="body2">{ligandsUpInSel1.length} ligands more active in {sel1Name}.</Typography></li>)}
                            {receptorsUpInSel2.length > 0 && (<li><Typography variant="body2">{receptorsUpInSel2.length} receptors more active in {sel2Name}.</Typography></li>)}
                            {receptorsUpInSel1.length > 0 && (<li><Typography variant="body2">{receptorsUpInSel1.length} receptors more active in {sel1Name}.</Typography></li>)}
                          </ul>
                        </Box>
                        <Typography variant="body2" color="text.secondary">Note: Automated interpretation. Biological validation needed.</Typography>
                      </>
                    );
                  })()}
                </Box>
              </Box>
              
              <Box sx={{ mb: 4 }}>
                <Typography variant="h6" gutterBottom sx={{ borderBottom: '2px solid #3f51b5', pb: 1, display: 'flex', alignItems: 'center' }}>
                  <span style={{ marginRight: '8px' }}>📊</span> Complete Results Table
                </Typography>
                <TableContainer component={Paper} sx={{ mt: 2, maxHeight: '600px' }}> 
                  <Table stickyHeader className={styles.resultsTable} aria-label="comparison results table">
                    <TableHead>
                      <TableRow>
                        <TableCell>Molecule <Tooltip title="Unique molecule identifier or ligand-receptor pair"><InfoIcon sx={{ml: 0.5, fontSize: '0.9rem', color: 'text.secondary'}} /></Tooltip></TableCell>
                        <TableCell>Rel. Expression (Sel 1) <Tooltip title="Relative expression level in the first selection"><InfoIcon sx={{ml: 0.5, fontSize: '0.9rem', color: 'text.secondary'}} /></Tooltip></TableCell>
                        <TableCell>Rel. Expression (Sel 2) <Tooltip title="Relative expression level in the second selection"><InfoIcon sx={{ml: 0.5, fontSize: '0.9rem', color: 'text.secondary'}} /></Tooltip></TableCell>
                        <TableCell>Log₂FC <Tooltip title="Log2 Fold Change. Positive: higher in Sel 2."><InfoIcon sx={{ml: 0.5, fontSize: '0.9rem', color: 'text.secondary'}} /></Tooltip></TableCell>
                        <TableCell>FDR <Tooltip title="False discovery rate as a percentage"><InfoIcon sx={{ml: 0.5, fontSize: '0.9rem', color: 'text.secondary'}} /></Tooltip></TableCell>
                        <TableCell>Type <Tooltip title="Molecule classification"><InfoIcon sx={{ml: 0.5, fontSize: '0.9rem', color: 'text.secondary'}} /></Tooltip></TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {comparisonJobDetails.results.results.differential_expression
                        .filter((item: DifferentialExpressionResultFE) => {
                          const absLog2FC = Math.abs(item.log2_fold_change || 0);
                          const maxMean = Math.max(item.mean_selection1 || 0, item.mean_selection2 || 0);
                          const typeFilter = 
                            (item.type === 'single_ligand' && showLigands) || 
                            (item.type === 'single_receptor' && showReceptors) || 
                            (item.type === 'ligand_receptor_pair' && showLRPairs);
                          return absLog2FC >= filterMinLog2FC && (item.q_value || 1) <= filterMaxQValue && maxMean >= filterMinMeanValue && typeFilter;
                        })
                        .map((item: DifferentialExpressionResultFE, index: number) => {
                          const nameInfo = getDisplayNameAndRole(item);
                          const fcInfo = safeFormatLog2FC(item.log2_fold_change);
                          const qValueInfo = formatQValue(item.q_value);
                          const sel1Norm = normalizeExpressionValue(item.mean_selection1);
                          const sel2Norm = normalizeExpressionValue(item.mean_selection2);
                          
                          return (
                            <TableRow 
                              key={`${item.molecule_id}-${index}-full`}
                              onClick={() => handleRowClick(item)}
                              sx={{ cursor: 'pointer', '&:hover': { backgroundColor: 'rgba(0, 0, 0, 0.04)' }, ...(selectedMolecule?.id === (item.type === 'ligand_receptor_pair' ? `${item.ligand_id}-${item.receptor_id}` : item.molecule_id) ? { backgroundColor: 'rgba(25, 118, 210, 0.08)' } : {})}}
                            >
                              <StyledTableCell>{nameInfo.displayName}</StyledTableCell>
                              <StyledTableCell><Tooltip title={`Exact: ${item.mean_selection1?.toExponential(4) || 'N/A'}`}><span>{sel1Norm.display}</span></Tooltip></StyledTableCell>
                              <StyledTableCell><Tooltip title={`Exact: ${item.mean_selection2?.toExponential(4) || 'N/A'}`}><span>{sel2Norm.display}</span></Tooltip></StyledTableCell>
                              <StyledTableCell>
                                <Box sx={{display: 'flex', alignItems: 'center'}}>
                                  {fcInfo.isPositive ? <ArrowUpward sx={{color: 'success.main', mr: 0.5, fontSize: '1rem'}} /> : <ArrowDownward sx={{color: 'error.main', mr: 0.5, fontSize: '1rem'}} />}
                                  {fcInfo.display}
                                  <Tooltip title={`${fcInfo.isPositive ? 'Increase' : 'Decrease'} of ${fcInfo.foldChange}`}><Box component="span" sx={{ml: 0.5, color: 'text.secondary', fontSize: '0.8rem'}}>({fcInfo.foldChange})</Box></Tooltip>
                                </Box>
                              </StyledTableCell>
                              <StyledTableCell><Tooltip title={`Exact q: ${item.q_value?.toExponential(4) || 'N/A'}`}><span>{qValueInfo.display} <span style={{ color: '#1976d2', marginLeft: 4 }}>{qValueInfo.significance}</span></span></Tooltip></StyledTableCell>
                              <StyledTableCell>{nameInfo.roleName}</StyledTableCell>
                            </TableRow>
                          );
                        })}
                    </TableBody>
                  </Table>
                </TableContainer>
              </Box>
            </>
          ) : (
            <Typography sx={{mt: 2, fontStyle: 'italic'}}>Comparison successful, but no significant differential expression found at FDR &lt; {fdrThreshold}.</Typography>
          )}
          {comparisonJobDetails.results?.errors && comparisonJobDetails.results.errors.length > 0 && (
            <Alert severity="warning" sx={{ mt: 2 }}>
               Comparison completed with issues:
               <ul style={{ marginTop: '8px', marginBottom: 0 }}>
                 {comparisonJobDetails.results.errors.map((err: string, i: number) => <li key={i}>{err}</li>)}
               </ul>
            </Alert>
           )}
      </Box>
    )}
    
    {!isLoadingComparisonJobStatus && !comparisonJobDetails && !comparisonJobError && (
        <Typography sx={{mt: 2, fontStyle: 'italic'}}>Comparison job initiated. Waiting for status updates...</Typography>
    )}
    
    {comparisonJobDetails?.status === 'success' && (
      <Box sx={{ mt: 3, mb: 2, p: 2, border: '1px solid #e0e0e0', borderRadius: 1 }}>
        <Typography variant="h6" gutterBottom>Filter Results</Typography>
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
          <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 2 }}>
            <Box sx={{ flex: '1 1 250px', minWidth: '200px' }}>
              <Tooltip title="Minimum absolute Log2 Fold Change value to display">
                <TextField label="Min |Log2FC|" type="number" fullWidth size="small" InputProps={{ inputProps: { min: 0, step: 0.1 } }} value={filterMinLog2FC} onChange={(e) => setFilterMinLog2FC(Math.max(0, parseFloat(e.target.value) || 0))} />
              </Tooltip>
            </Box>
            <Box sx={{ flex: '1 1 250px', minWidth: '200px' }}>
              <Tooltip title="Maximum adjusted p-value (q-value) to display">
                <TextField label="Max q-value" type="number" fullWidth size="small" InputProps={{ inputProps: { min: 0, max: 0.1, step: 0.001 } }} value={filterMaxQValue} onChange={(e) => setFilterMaxQValue(Math.max(0, Math.min(0.1, parseFloat(e.target.value) || 0.05)))} />
              </Tooltip>
            </Box>
            <Box sx={{ flex: '1 1 250px', minWidth: '200px' }}>
              <Tooltip title="Minimum mean value in either selection to display">
                <TextField label="Min Mean Value" type="number" fullWidth size="small" InputProps={{ inputProps: { min: 0, step: 0.0001 } }} value={filterMinMeanValue} onChange={(e) => setFilterMinMeanValue(Math.max(0, parseFloat(e.target.value) || 0))} />
              </Tooltip>
            </Box>
          </Box>
          <Box>
            <FormGroup row>
              <FormControlLabel control={<Checkbox checked={showLigands} onChange={(e) => setShowLigands(e.target.checked)}/>} label="Show Ligands" />
              <FormControlLabel control={<Checkbox checked={showReceptors} onChange={(e) => setShowReceptors(e.target.checked)}/>} label="Show Receptors" />
              <FormControlLabel control={<Checkbox checked={showLRPairs} onChange={(e) => setShowLRPairs(e.target.checked)}/>} label="Show L-R Pairs" />
            </FormGroup>
          </Box>
          <Typography variant="body2" color="text.secondary">
            {comparisonJobDetails.results?.results?.differential_expression
              .filter((item: DifferentialExpressionResultFE) => {
                const absLog2FC = Math.abs(item.log2_fold_change || 0);
                const maxMean = Math.max(item.mean_selection1 || 0, item.mean_selection2 || 0);
                const typeFilter = (item.type === 'single_ligand' && showLigands) || (item.type === 'single_receptor' && showReceptors) || (item.type === 'ligand_receptor_pair' && showLRPairs);
                return absLog2FC >= filterMinLog2FC && (item.q_value || 1) <= filterMaxQValue && maxMean >= filterMinMeanValue && typeFilter;
              }).length
            } of {comparisonJobDetails.results?.results?.differential_expression.length} significant results displayed
          </Typography>
        </Box>
      </Box>
    )}
  </Box>
  );
};

export default ComparisonResultsDisplay; 