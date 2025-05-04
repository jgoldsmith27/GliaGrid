// Defines the structure of data returned by specific analysis types

/**
 * Represents a single row in the Pathway Dominance results.
 */
export interface PathwayDominanceResult {
  ligand: string;
  receptor: string;
  score: number;
  pathway?: string; // Add pathway name (optional if not always present)
  ligand_norm_expr: number; // Add normalized expression
  receptor_avg_norm_expr: number; // Add normalized expression
  // Add other potential fields if known, e.g., pathway name
  // pathway?: string; 
}

/**
 * Represents a single row in the Module Context results.
 */
export interface ModuleContextResult {
  ligand: string;
  receptor: string;
  interaction_type: string;
  ligand_module: string;
  receptor_modules: string[];
  is_same_module: boolean;
}

/**
 * Represents the structure for Summary Statistics (Ligand/Receptor Counts).
 */
export interface SummaryStatsData {
    unique_ligands: number;
    unique_receptors: number;
    // Add any other summary stats if they exist
}

// Could potentially define a union type if useful elsewhere
// export type InteractionResult = PathwayDominanceResult | ModuleContextResult; 

// Combined type (potentially useful, but not strictly necessary if using individual results)
// export interface CombinedResultItem extends PathwayDominanceResult, ModuleContextResult {}


// --- Custom Analysis Response Structures (Refactored) ---

// Structure for a single scope's result (either whole custom OR one layer)
export interface CustomAnalysisScopeResult {
    pathway_dominance: PathwayDominanceResult[];
    module_context: ModuleContextResult[];
}

// Response bundle containing both whole and layered results
export interface CustomAnalysisResultsBundle {
    whole_results: CustomAnalysisScopeResult; // Always present
    layered_results: { [layerName: string]: CustomAnalysisScopeResult }; // Optional/Can be empty
}

// DEPRECATED - Replaced by CustomAnalysisResultsBundle
// export interface CustomAnalysisResponse extends CustomAnalysisScopeResult {}
// export interface LayeredCustomAnalysisResponse {
//     results_by_layer: { [layerName: string]: CustomAnalysisScopeResult };
// } 