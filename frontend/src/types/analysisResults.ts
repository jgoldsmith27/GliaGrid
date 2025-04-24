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
  interaction_type: 'intra-module' | 'inter-module' | string; // Use 'interaction_type' for consistency, allow flexibility
  ligand_module?: string; // Add module info (optional if not always present)
  receptor_modules?: string[]; // Add receptor modules info (optional array)
  is_same_module?: boolean; // Add same module flag (optional boolean)
  // num_receptor_components_in_modules: number; // Likely not needed for display
  // Add other potential fields if known, e.g., module names
  // ligand_module?: string;
  // receptor_module?: string;
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