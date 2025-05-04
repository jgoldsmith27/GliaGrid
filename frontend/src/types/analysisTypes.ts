// Defines shared types related to analysis scope and results

export type AnalysisScope = 'layers' | 'pairs' | 'custom';

// Combined data structure used in ResultsPage and SummaryTabContent
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

// Add other relevant types here later, e.g., CombinedInteractionData if needed elsewhere 