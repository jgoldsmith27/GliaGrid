import React, { useState, useEffect, useRef, useCallback } from 'react';
import styles from './SummaryTabContent.module.css'; // Create CSS module later
import AnalysisTable, { ColumnDefinition } from '../AnalysisTable/AnalysisTable';
import InteractionVisualization from '../InteractionVisualization/InteractionVisualization.tsx';
import SpatialOverviewVisualization from '../SpatialOverviewVisualization/SpatialOverviewVisualization';
import { CombinedInteractionData } from '../../pages/ResultsPage/ResultsPage'; // Import shared type
import { ScopeType } from '../ScopeSelector/ScopeSelector'; // Import ScopeType if defined there
import useInteractionData from '../../hooks/useInteractionData'; // Import the hook
import { 
    PathwayDominanceResult, ModuleContextResult, 
    CustomAnalysisScopeResult, CustomAnalysisResultsBundle 
} from '../../types/analysisResults'; // ADDED import
import LoadingSpinner from '../LoadingSpinner/LoadingSpinner'; // ADDED import for loading state

// Type for the data returned by the new /api/points/{jobId}/all endpoint
interface AllPointsData {
    x: number;
    y: number;
    layer: string;
}

// MODIFIED: Add new props
interface SummaryTabContentProps {
  jobId: string;
  combinedData: CombinedInteractionData[];
  selectedPair: [string, string] | null;
  onSelectPair: (pair: [string, string] | null) => void; // Allow null to clear selection
  currentScope: ScopeType; 
  apiScopeName: string | null; // ADDED: The actual scope name for the API call
  onLassoSelect?: (coords: [number, number][] | null) => void; // ADDED prop definition
  onAnalyzeSelection?: () => void;
  // MODIFIED: Use the new bundle type
  customAnalysisResults?: CustomAnalysisResultsBundle | null; 
  isLoadingCustomAnalysis?: boolean;
  customAnalysisError?: string | null;
  // ADDED: Props for aggregation level state (lifted up)
  customAggregationLevel: CustomAggregationLevel;
  setCustomAggregationLevel: React.Dispatch<React.SetStateAction<CustomAggregationLevel>>;
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
  // ADDED: Destructure new props
  customAggregationLevel,
  setCustomAggregationLevel,
}) => {
  // REMOVED: Internal state for custom aggregation level
  // const [customAggregationLevel, setCustomAggregationLevel] = useState<CustomAggregationLevel>('whole_custom');

  // ADDED: Log received prop
  console.log("[SummaryTabContent] Received onAnalyzeSelection:", typeof onAnalyzeSelection);

  // State for SpatialOverviewVisualization
  // const [allPointsData, setAllPointsData] = useState<AllPointsData[] | null>(null);
  // const [loadingAllPoints, setLoadingAllPoints] = useState(false);
  // const [allPointsError, setAllPointsError] = useState<string | null>(null);

  const [selectedRowIndex, setSelectedRowIndex] = useState<number | null>(null); // For highlighting table row
  
  // Use the custom hook for interaction data (non-custom scopes)
  const { interactionVizData, isLoading: isLoadingInteractionViz, error: interactionVizError, warnings: interactionVizWarnings, cancelFetch: cancelInteractionVizFetch } = 
    useInteractionData(jobId, selectedPair, currentScope !== 'custom' ? apiScopeName : null); 

  // Effect to fetch data based on currentScope (REMOVED fetchAllPoints call)
  useEffect(() => {
    console.log(`[SummaryTabContent] useEffect triggered. Scope: ${currentScope}, Pair: ${selectedPair ? selectedPair.join('-') : 'null'}, apiScopeName: ${apiScopeName}`); 
    
    if (jobId) {
        if (currentScope === 'custom') {
             // REMOVED: fetchAllPoints(jobId);
             // Interaction viz hook will handle its state based on apiScopeName being 'custom'
        } 
        // No explicit fetch call needed for interaction viz here; the hook handles it based on its dependencies
    }

      // REMOVED: Cleanup logic related to fetchAllPoints abortController
      return () => {
         console.log("[SummaryTabContent] useEffect cleanup");
         // No cleanup needed here now for fetchAllPoints
      };
  // MODIFIED: Removed fetchAllPoints from dependency array
  }, [jobId, currentScope, apiScopeName, selectedPair, onSelectPair]); 

  // Table row click handler
  const handleTableRowClick = useCallback((row: CombinedInteractionData, index: number) => {
    // This handler applies to BOTH original and custom analysis tables
    console.log(`[SummaryTabContent] Row clicked: ${row.ligand}-${row.receptor}, Index: ${index}, Scope: ${currentScope}`); 
    setSelectedRowIndex(index);
    const newPair: [string, string] = [row.ligand, row.receptor];
    onSelectPair(newPair); // Call parent handler to manage selectedPair state
  }, [onSelectPair, currentScope]);

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

            // Use the interaction hook for custom results visualization
            // NOTE: This re-uses the hook; ensure dependencies are correct
            const { interactionVizData: customVizData, isLoading: isLoadingCustomViz, error: customVizError } =
                useInteractionData(jobId, selectedPair, null);

            // --- Extract data based on aggregation level from the bundle --- 
            const layeredData = customAnalysisResults.layered_results || {};
            const isLayeredAvailable = Object.keys(layeredData).length > 0;
            const availableLayers = Object.keys(layeredData);
            
            // State for selecting a layer if results are layered
            // Initialize based on available layers, only if layered view is selected AND layers exist
            const [selectedCustomLayer, setSelectedCustomLayer] = useState<string | null>(() => {
                return customAggregationLevel === 'custom_by_layer' && isLayeredAvailable ? availableLayers[0] : null;
            });

            // Update selected layer if aggregation level changes or results update
            useEffect(() => {
                if (customAggregationLevel === 'custom_by_layer' && isLayeredAvailable && !availableLayers.includes(selectedCustomLayer || '')) {
                    setSelectedCustomLayer(availableLayers[0]);
                } else if (customAggregationLevel === 'whole_custom') {
                    setSelectedCustomLayer(null); // Clear layer selection when switching to whole view
                }
            }, [customAggregationLevel, isLayeredAvailable, availableLayers, selectedCustomLayer]);
            
            // Get the specific data scope to display based on aggregation level and selected layer
            let currentCustomData: CustomAnalysisScopeResult | null = null;
            if (customAggregationLevel === 'whole_custom') {
                currentCustomData = customAnalysisResults.whole_results;
            } else if (customAggregationLevel === 'custom_by_layer' && isLayeredAvailable && selectedCustomLayer) {
                 currentCustomData = layeredData[selectedCustomLayer];
            }

            // Process the *current* data for the table
            const pathwayData = currentCustomData?.pathway_dominance || [];
            const moduleData = currentCustomData?.module_context || [];
            const moduleContextMap = new Map<string, ModuleContextResult>();
            moduleData.forEach((item: ModuleContextResult) => { 
                if (item.ligand && item.receptor) { 
                    moduleContextMap.set(`${item.ligand}-${item.receptor}`, item); 
                }
            });
            const customCombinedData: CombinedInteractionData[] = pathwayData.map((pathwayItem: PathwayDominanceResult) => {
                const moduleItem = pathwayItem.ligand && pathwayItem.receptor ? moduleContextMap.get(`${pathwayItem.ligand}-${pathwayItem.receptor}`) : undefined;
                return { ...pathwayItem, ...(moduleItem || {}) }; 
            });

            return (
              <>
                <div className={styles.customAggregationSelector}>
                  <label>
                    <input
                      type="radio"
                      name="customAggregation"
                      value="whole_custom"
                      checked={customAggregationLevel === 'whole_custom'}
                      onChange={() => setCustomAggregationLevel('whole_custom')}
                    />
                    Whole Custom Selection
                  </label>
                  <label style={{ marginLeft: '10px' }}>
                    <input
                      type="radio"
                      name="customAggregation"
                      value="custom_by_layer"
                      checked={customAggregationLevel === 'custom_by_layer'}
                      onChange={() => setCustomAggregationLevel('custom_by_layer')}
                      disabled={!isLayeredAvailable} 
                    />
                    Custom Points by Layer
                  </label>
                  {/* Layer selector dropdown (only shows if layered view active AND layers available) */}
                  {customAggregationLevel === 'custom_by_layer' && isLayeredAvailable && (
                      <select 
                          value={selectedCustomLayer || ''} 
                          onChange={(e) => setSelectedCustomLayer(e.target.value)} 
                          className={styles.layerSelectorDropdown} 
                          disabled={availableLayers.length === 0}
                      >
                          {availableLayers.length === 0 && <option value="">No Layers Found</option>}
                          {availableLayers.map(layer => (
                              <option key={layer} value={layer}>{layer}</option>
                          ))}
                      </select>
                  )}
                </div>
                
                <div className={styles.summaryLayout}> 
                  {/* Conditional Rendering based on aggregation level */}
                  {customAggregationLevel === 'whole_custom' ? (
                    <>
                      <div className={styles.tableArea}>
                        <h3>Interaction Scores (Whole Custom Selection)</h3>
                        <AnalysisTable
                          data={customCombinedData} 
                          columns={combinedColumns}
                          onRowClick={handleTableRowClick}
                          selectedRowIndex={selectedRowIndex}
                          loading={isLoadingCustomViz} 
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
                              currentScope={'whole_tissue'} // Representing the whole custom scope
                          />
                        )}
                      </div>
                    </>
                  ) : (
                    // MODIFIED: Display layered results or placeholder
                    isLayeredAvailable && selectedCustomLayer && currentCustomData ? (
                      <>
                        <div className={styles.tableArea}>
                          <h3>Interaction Scores (Layer: {selectedCustomLayer})</h3>
                          <AnalysisTable
                            data={customCombinedData} 
                            columns={combinedColumns}
                            onRowClick={handleTableRowClick}
                            selectedRowIndex={selectedRowIndex}
                            loading={isLoadingCustomViz}
                          />
                        </div>
                        <div className={styles.visualizationArea}>
                          <h3>Interaction Visualization ({selectedPair ? `${selectedPair.join('-')} [${selectedCustomLayer}]` : `Select Pair [${selectedCustomLayer}]`})</h3>
                          {isLoadingCustomViz && <p>Loading visualization...</p>}
                          {customVizError && <p className={styles.errorText}>Viz Error: {customVizError}</p>}
                          {!selectedPair && <p>Select a pair from the table.</p>}
                          {customVizData && selectedPair && (
                              <InteractionVisualization
                                data={customVizData} 
                                ligandName={selectedPair[0]}
                                receptorName={selectedPair[1]}
                                currentScope={'layers'} // Representing layer scope
                            />
                          )}
                        </div>
                      </>
                    ) : (
                      // Placeholder if layer view selected but results/layer not available/selected
                      <div className={styles.placeholderArea}>
                        <h3>Interaction Scores (Custom Points by Layer)</h3>
                        {!isLayeredAvailable && <p>Layer-specific results are not available for this selection.</p>}
                        {isLayeredAvailable && !selectedCustomLayer && <p>Select a layer to view results.</p>}
                      </div>
                    )
                  )}
                </div>
              </>
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
                     currentScope={currentScope === 'layers' ? 'layers' : 'whole_tissue'}
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

// ADDED Types (copy from SummaryTabContent for now, consider moving to shared types file)
type CustomAggregationLevel = 'whole_custom' | 'custom_by_layer';

// Define response types (consider moving to shared types file)
interface AnalysisResultItem { // Simplified for example
    ligand?: string | null;
    receptor?: string | null;
    score?: number | null;
    [key: string]: any; // Allow extra fields
}

interface CustomAnalysisResponse {
    pathway_dominance: AnalysisResultItem[];
    module_context: AnalysisResultItem[];
}

interface LayeredCustomAnalysisResponse {
    results_by_layer: { [layerName: string]: CustomAnalysisResponse };
} 