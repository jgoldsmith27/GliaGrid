"""Core analysis functions for calculating counts, pathway dominance, and module context."""

import pandas as pd
import numpy as np
from typing import Dict, Any, List, Set
from sklearn.preprocessing import MinMaxScaler # For pathway dominance normalization
import time # For placeholder simulation - REMOVE if not needed


# --- Helper Functions --- 

def _get_present_genes(spatial_df_subset: pd.DataFrame) -> Set[str]:
    """Get the set of unique gene names present in the spatial subset."""
    if 'gene' not in spatial_df_subset.columns or spatial_df_subset.empty:
        return set()
    return set(spatial_df_subset['gene'].unique())

def _split_receptor_complex(receptor_name: str) -> List[str]:
    """Split a receptor name into components if it's a complex (contains '_')."""
    return receptor_name.split('_')

# --- Stage 1: Ligand-Receptor Counts --- 

def _calculate_lr_counts_for_scope(spatial_df_subset: pd.DataFrame, interactions_df: pd.DataFrame) -> Dict[str, int]:
    """Calculates unique ligand/receptor counts for a given spatial data subset."""
    print(f"    Calculating counts for scope with {len(spatial_df_subset)} spots...")
    
    present_genes = _get_present_genes(spatial_df_subset)
    if not present_genes:
        return {"unique_ligands": 0, "unique_receptors": 0}

    unique_ligands_found = set()
    unique_receptors_found = set()

    for _, interaction in interactions_df.iterrows():
        ligand = interaction['ligand']
        receptor_complex = str(interaction['receptor']) # Ensure string
        
        # Check ligand presence
        if ligand in present_genes:
            unique_ligands_found.add(ligand)
            
        # Check receptor components presence
        receptor_components = _split_receptor_complex(receptor_complex)
        for component in receptor_components:
            if component in present_genes:
                unique_receptors_found.add(component)
                # Optimization: if one component is found, we count the gene, 
                # no need to check other components of the *same* complex for *this* gene count.
                # However, we add the specific component name found to the set.

    return {"unique_ligands": len(unique_ligands_found), "unique_receptors": len(unique_receptors_found)}

def run_stage1_counts_pipeline(spatial_df: pd.DataFrame, interactions_df: pd.DataFrame) -> Dict[str, Dict[str, int]]:
    """Runs the LR count calculation for whole tissue and each layer."""
    print("Running Stage 1 Counts Pipeline...")
    results = {}
    # Ensure 'layer' column exists and handle potential non-string layer names gracefully
    if 'layer' not in spatial_df.columns:
        print("ERROR: 'layer' column not found in spatial data.")
        return {'whole_tissue': _calculate_lr_counts_for_scope(spatial_df, interactions_df)} # Calculate for whole tissue only
        
    try:
        # Attempt to convert layer names to string to ensure consistency
        layers = list(spatial_df['layer'].astype(str).unique())
    except Exception as e:
        print(f"ERROR converting layer names to string: {e}. Proceeding with original layer types.")
        layers = list(spatial_df['layer'].unique())
        
    scopes = ['whole_tissue'] + layers

    for scope in scopes:
        print(f"  Processing scope: {scope}")
        if scope == 'whole_tissue':
            subset_df = spatial_df
        else:
            # Filter using the potentially converted scope name
            subset_df = spatial_df[spatial_df['layer'] == scope] 
        
        if subset_df.empty:
            print(f"    Skipping empty scope: {scope}")
            results[str(scope)] = {"unique_ligands": 0, "unique_receptors": 0}
            continue
            
        results[str(scope)] = _calculate_lr_counts_for_scope(subset_df, interactions_df)
        
    print("Stage 1 Counts Pipeline Complete.")
    return results

# --- Stage 2a: Pathway Dominance --- 

def _calculate_pathway_dominance_for_scope(spatial_df_subset: pd.DataFrame, interactions_df: pd.DataFrame) -> List[Dict[str, Any]]:
    """Calculates pathway dominance scores for a given spatial data subset.
    1. Finds all ligands/receptors present in the scope.
    2. Filters interactions to only those possible within the scope.
    3. Calculates relative frequency for present L/R as expression proxy.
    4. Normalizes expression proxy using MinMaxScaler.
    5. Scores possible interactions based on normalized expression.
    """
    print(f"    Calculating pathway dominance for scope with {len(spatial_df_subset)} spots...")
    
    if spatial_df_subset.empty or 'gene' not in spatial_df_subset.columns:
        return []

    # 1. Find all L/R present in the scope
    present_genes = _get_present_genes(spatial_df_subset)
    if not present_genes:
        return []
        
    # Also identify which specific genes are ligands or receptors based on the full interactions list
    all_ligands_in_db = set(interactions_df['ligand'].unique())
    all_receptor_components_in_db = set(comp for receptor in interactions_df['receptor'].astype(str) for comp in _split_receptor_complex(receptor))
    
    present_ligands = present_genes.intersection(all_ligands_in_db)
    present_receptor_components = present_genes.intersection(all_receptor_components_in_db)
    
    if not present_ligands or not present_receptor_components:
        print("    Skipping pathway dominance: No ligands or receptor components present in scope.")
        return []

    # 2. Filter interactions_df to possible interactions within the scope
    possible_interactions = interactions_df[
        interactions_df['ligand'].isin(present_ligands) & 
        interactions_df['receptor'].astype(str).apply(lambda r: any(comp in present_receptor_components for comp in _split_receptor_complex(r)))
    ].copy() # Create a copy to avoid SettingWithCopyWarning

    if possible_interactions.empty:
        print("    Skipping pathway dominance: No possible interactions found for present L/R.")
        return []
    print(f"    Found {len(possible_interactions)} possible interactions for the scope.")

    # 3. Calculate relative frequency for *all* present genes (ligands and receptors)
    gene_counts = spatial_df_subset['gene'].value_counts()
    total_spots = len(spatial_df_subset)
    # Focus on genes relevant to the possible interactions for normalization
    relevant_genes_for_norm = present_ligands.union(present_receptor_components)
    gene_expression = {
        gene: count / total_spots 
        for gene, count in gene_counts.items() 
        if gene in relevant_genes_for_norm
    }

    # 4. Normalize expression proxy (0-1 scale) for relevant genes
    if not gene_expression:
         print("    Skipping pathway dominance: No relevant genes found for expression calculation.")
         return [] 
         
    expression_genes = list(gene_expression.keys())
    expression_values = np.array([gene_expression[g] for g in expression_genes]).reshape(-1, 1)
    
    # Handle case with only one relevant gene (MinMaxScaler needs >1 value range)
    if len(expression_values) == 1:
        gene_normalized_expression = {expression_genes[0]: 1.0} # Assign max normalized value
    elif len(expression_values) > 1:
        scaler = MinMaxScaler()
        normalized_values = scaler.fit_transform(expression_values).flatten()
        gene_normalized_expression = dict(zip(expression_genes, normalized_values))
    else: # Should not happen if gene_expression check passed, but defensively...
        gene_normalized_expression = {}

    # 5. Analyze and score the *possible* interaction pairs
    scored_pairs = []
    for _, interaction in possible_interactions.iterrows():
        ligand = interaction['ligand']
        receptor_complex = str(interaction['receptor']) 
        pathway = interaction.get('pathway_name', 'Unknown')
        
        # Ligand must be in normalized expression dict (was present and relevant)
        if ligand not in gene_normalized_expression:
            continue 
            
        ligand_normalized = gene_normalized_expression[ligand]
        
        # Calculate average normalized expression for receptor components *present AND relevant*
        receptor_components = _split_receptor_complex(receptor_complex)
        receptor_norm_expr_values = []
        for component in receptor_components:
            if component in gene_normalized_expression: # Check if component was present and included in normalization
                receptor_norm_expr_values.append(gene_normalized_expression[component])
        
        # Only score if at least one receptor component has a normalized expression value
        if not receptor_norm_expr_values:
            continue
            
        avg_receptor_normalized = sum(receptor_norm_expr_values) / len(receptor_norm_expr_values)
        
        normalized_score = ligand_normalized * avg_receptor_normalized
        
        scored_pairs.append({
            'ligand': ligand,
            'receptor': receptor_complex,
            'pathway': pathway,
            'ligand_norm_expr': ligand_normalized,
            'receptor_avg_norm_expr': avg_receptor_normalized,
            'score': normalized_score
        })

    # 6. Sort by score and return
    if not scored_pairs:
        return []
        
    pairs_df = pd.DataFrame(scored_pairs)
    pairs_df = pairs_df.sort_values('score', ascending=False)
    
    print(f"    Scored {len(pairs_df)} interactions for the scope.")
    return pairs_df.to_dict('records')

def run_pathway_dominance_pipeline(spatial_df: pd.DataFrame, interactions_df: pd.DataFrame) -> Dict[str, List[Dict[str, Any]]]:
    """Runs the pathway dominance calculation for whole tissue and each layer."""
    print("Running Pathway Dominance Pipeline...")
    results = {}
    if 'layer' not in spatial_df.columns:
        print("ERROR: 'layer' column not found in spatial data.")
        return {'whole_tissue': _calculate_pathway_dominance_for_scope(spatial_df, interactions_df)}
        
    try:
        layers = list(spatial_df['layer'].astype(str).unique())
    except Exception as e:
        print(f"ERROR converting layer names to string: {e}. Proceeding with original layer types.")
        layers = list(spatial_df['layer'].unique())
        
    scopes = ['whole_tissue'] + layers

    for scope in scopes:
        print(f"  Processing scope: {scope}")
        if scope == 'whole_tissue':
            subset_df = spatial_df
        else:
            subset_df = spatial_df[spatial_df['layer'] == scope]
        
        if subset_df.empty:
            print(f"    Skipping empty scope: {scope}")
            results[str(scope)] = []
            continue
            
        results[str(scope)] = _calculate_pathway_dominance_for_scope(subset_df, interactions_df)
        
    print("Pathway Dominance Pipeline Complete.")
    return results

# --- Stage 2b: Module Context --- 

def _calculate_module_context_for_scope(modules_df: pd.DataFrame, pathway_results_subset: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Calculates module context for significant pairs in a given spatial data subset.
    Assumes pathway_results_subset is a list of dicts with 'ligand' and 'receptor' keys.
    """
    print(f"    Calculating module context for scope...") # No subset needed here directly

    if not pathway_results_subset or modules_df.empty or 'gene' not in modules_df.columns or 'module' not in modules_df.columns:
        print("    Skipping module context: Missing pathway results or invalid modules data.")
        return []
        
    # Create module lookup for faster access (gene -> module)
    # Ensure gene column is string type for reliable lookup
    try:
         modules_df['gene'] = modules_df['gene'].astype(str)
         module_lookup = modules_df.set_index('gene')['module'].to_dict()
    except Exception as e:
        print(f"ERROR creating module lookup: {e}. Check modules file columns ('gene', 'module').")
        return []

    module_context_results = []
    for pair in pathway_results_subset: # Iterate through the provided pathway pairs
        ligand = pair.get('ligand')
        receptor_complex = str(pair.get('receptor', '')) # Ensure string
        
        if not ligand or not receptor_complex:
            continue # Skip if pair is missing essential info

        receptors = _split_receptor_complex(receptor_complex)
        ligand_module = module_lookup.get(ligand)
        
        receptor_modules = []
        num_receptor_components_in_modules = 0
        for r in receptors:
            r_module = module_lookup.get(r)
            if r_module is not None:
                receptor_modules.append(r_module)
                num_receptor_components_in_modules += 1

        # Determine interaction type based on module presence
        interaction_type = 'Unknown'
        is_same_module = None 
        if ligand_module is not None and receptor_modules: # Both ligand and receptor(s) have modules
            # All found receptor components must match the ligand module
            is_same_module = all(m == ligand_module for m in receptor_modules)
            interaction_type = 'Intra-module' if is_same_module else 'Inter-module'
        elif ligand_module is not None: # Only ligand module known
             interaction_type = 'Ligand Module Known'
        elif receptor_modules: # Only receptor module(s) known
             interaction_type = 'Receptor Module Known'
             # Can check if all receptor components are in the *same* module, even if ligand unknown
             if len(set(receptor_modules)) == 1:
                  is_same_module = True # All receptor components are in the same module
             else:
                  is_same_module = False # Receptor components span different modules
        # else: Both unknown, remains 'Unknown'

        context = pair.copy() # Start with original pair data (like score)
        context.update({
            'ligand_module': ligand_module, 
            'receptor_modules': receptor_modules if receptor_modules else None, # Use None if empty
            'num_receptor_components_in_modules': num_receptor_components_in_modules,
            'is_same_module': is_same_module, # Can be True, False, or None
            'interaction_type': interaction_type
        })
        module_context_results.append(context)

    return module_context_results

def run_module_context_pipeline(spatial_df: pd.DataFrame, interactions_df: pd.DataFrame, modules_df: pd.DataFrame, full_pathway_results: Dict[str, List[Dict[str, Any]]]) -> Dict[str, List[Dict[str, Any]]]:
    """Runs the module context calculation for whole tissue and each layer."""
    print("Running Module Context Pipeline...")
    results = {}
    if 'layer' not in spatial_df.columns:
        print("ERROR: 'layer' column not found in spatial data. Cannot run module context per layer.")
        # Calculate only for whole tissue if pathway results exist
        pathway_subset = full_pathway_results.get('whole_tissue', [])
        return {'whole_tissue': _calculate_module_context_for_scope(modules_df, pathway_subset)} 
        
    try:
        layers = list(spatial_df['layer'].astype(str).unique())
    except Exception as e:
        print(f"ERROR converting layer names to string: {e}. Proceeding with original layer types.")
        layers = list(spatial_df['layer'].unique())

    scopes = ['whole_tissue'] + layers

    for scope in scopes:
        print(f"  Processing scope: {scope}")
        # Get the relevant pathway results for this scope
        pathway_subset = full_pathway_results.get(str(scope), []) # Use str(scope) for lookup
        
        # Module context calculation doesn't directly need the spatial subset,
        # only the pathway results for the scope and the global modules_df.
        if not pathway_subset:
             print(f"    Skipping empty/missing pathway results for scope: {scope}")
             results[str(scope)] = []
             continue
             
        results[str(scope)] = _calculate_module_context_for_scope(modules_df, pathway_subset)
        
    print("Module Context Pipeline Complete.")
    return results 