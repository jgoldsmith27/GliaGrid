"""Core analysis functions for calculating counts, pathway dominance, and module context."""

import pandas as pd
import numpy as np
from typing import Dict, Any, List, Set
from sklearn.preprocessing import MinMaxScaler # For pathway dominance normalization
import time # For placeholder simulation - REMOVE if not needed
import logging # ADD logging import

# --- Setup Logger ---
logger = logging.getLogger(__name__) # ADD logger instance

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
    logger.debug(f"    Calculating counts for scope with {len(spatial_df_subset)} spots...")
    
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
    logger.info("Running Stage 1 Counts Pipeline...")
    results = {}
    # Ensure 'layer' column exists and handle potential non-string layer names gracefully
    if 'layer' not in spatial_df.columns:
        logger.error("'layer' column not found in spatial data.")
        return {'whole_tissue': _calculate_lr_counts_for_scope(spatial_df, interactions_df)} # Calculate for whole tissue only
        
    try:
        # Attempt to convert layer names to string to ensure consistency
        layers = list(spatial_df['layer'].astype(str).unique())
    except Exception as e:
        logger.error(f"Error converting layer names to string: {e}. Proceeding with original layer types.")
        layers = list(spatial_df['layer'].unique())
        
    scopes = ['whole_tissue'] + layers

    for scope in scopes:
        logger.info(f"  Processing Counts scope: {scope}")
        try:
            if scope == 'whole_tissue':
                subset_df = spatial_df
            else:
                # Filter using the potentially converted scope name
                subset_df = spatial_df[spatial_df['layer'] == scope] 
            
            if subset_df.empty:
                logger.warning(f"    Skipping empty scope: {scope}")
                results[str(scope)] = {"unique_ligands": 0, "unique_receptors": 0}
                continue
            
            results[str(scope)] = _calculate_lr_counts_for_scope(subset_df, interactions_df)
        except Exception as e:
            logger.exception(f"    ERROR processing counts for scope {scope}: {e}")
            results[str(scope)] = {"unique_ligands": -1, "unique_receptors": -1} # Indicate error
        
    logger.info("Stage 1 Counts Pipeline Complete.")
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
    logger.debug(f"    Calculating pathway dominance for scope with {len(spatial_df_subset)} spots...")
    
    if spatial_df_subset.empty or 'gene' not in spatial_df_subset.columns:
        logger.warning("    Skipping pathway dominance: Missing required columns (gene) or empty subset.")
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
        logger.warning("    Skipping pathway dominance: No ligands or receptor components present in scope.")
        return []

    # 2. Filter interactions_df to possible interactions within the scope
    possible_interactions = interactions_df[
        interactions_df['ligand'].isin(present_ligands) & 
        interactions_df['receptor'].astype(str).apply(lambda r: any(comp in present_receptor_components for comp in _split_receptor_complex(r)))
    ].copy() # Create a copy to avoid SettingWithCopyWarning

    if possible_interactions.empty:
        logger.warning("    Skipping pathway dominance: No possible interactions found for present L/R.")
        return []
    logger.debug(f"    Found {len(possible_interactions)} possible interactions for the scope.")

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
         logger.warning("    Skipping pathway dominance: No relevant genes found for expression calculation.")
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
        
        if ligand not in gene_normalized_expression:
            continue 
            
        ligand_normalized = gene_normalized_expression[ligand]
        
        receptor_components = _split_receptor_complex(receptor_complex)
        receptor_norm_expr_values = []
        for component in receptor_components:
            if component in gene_normalized_expression:
                receptor_norm_expr_values.append(gene_normalized_expression[component])
        
        if not receptor_norm_expr_values:
            continue
            
        avg_receptor_normalized = sum(receptor_norm_expr_values) / len(receptor_norm_expr_values)
        
        normalized_score = ligand_normalized * avg_receptor_normalized
        
        scored_pairs.append({
            'ligand': ligand,
            'receptor': receptor_complex,
            'ligand_norm_expr': ligand_normalized,
            'receptor_avg_norm_expr': avg_receptor_normalized,
            'score': normalized_score
        })

    # Sort by score and return
    if not scored_pairs:
        return []
        
    pairs_df = pd.DataFrame(scored_pairs)
    pairs_df = pairs_df.sort_values('score', ascending=False)
    
    logger.debug(f"    Scored {len(pairs_df)} interactions for the scope.")
    return pairs_df.to_dict('records') # Return no longer includes points

def run_pathway_dominance_pipeline(spatial_df: pd.DataFrame, interactions_df: pd.DataFrame) -> Dict[str, List[Dict[str, Any]]]:
    """Runs the pathway dominance calculation for whole tissue and each layer."""
    logger.info("Running Pathway Dominance Pipeline...")
    results = {}
    if 'layer' not in spatial_df.columns:
        logger.error("'layer' column not found in spatial data.")
        return {'whole_tissue': _calculate_pathway_dominance_for_scope(spatial_df, interactions_df)}
        
    try:
        layers = list(spatial_df['layer'].astype(str).unique())
    except Exception as e:
        logger.error(f"Error converting layer names to string: {e}. Proceeding with original layer types.")
        layers = list(spatial_df['layer'].unique())
        
    scopes = ['whole_tissue'] + layers

    for scope in scopes:
        logger.info(f"  Processing Pathway Dominance scope: {scope}")
        try:
            if scope == 'whole_tissue':
                subset_df = spatial_df
            else:
                subset_df = spatial_df[spatial_df['layer'] == scope]
            
            if subset_df.empty:
                logger.warning(f"    Skipping empty scope: {scope}")
                results[str(scope)] = []
                continue
            
            results[str(scope)] = _calculate_pathway_dominance_for_scope(subset_df, interactions_df)
        except Exception as e:
            logger.exception(f"    ERROR processing pathway dominance for scope {scope}: {e}")
            results[str(scope)] = [] # Return empty list on error for this scope
        
    logger.info("Pathway Dominance Pipeline Complete.")
    return results

# --- Stage 2b: Module Context --- 

def _calculate_module_context_for_scope(modules_df: pd.DataFrame, pathway_results_subset: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Calculates module context for significant pairs in a given spatial data subset.
    Assumes pathway_results_subset is a list of dicts with 'ligand' and 'receptor' keys.
    """
    logger.debug("    Calculating module context for scope...")

    if not pathway_results_subset or modules_df.empty or 'gene' not in modules_df.columns or 'module' not in modules_df.columns:
        logger.warning("    Skipping module context: Missing pathway results or invalid modules data.")
        return []
        
    # Create module lookup for faster access (gene -> module)
    # Ensure gene column is string type for reliable lookup
    try:
         modules_df['gene'] = modules_df['gene'].astype(str)
         module_lookup = modules_df.set_index('gene')['module'].to_dict()
    except Exception as e:
        logger.error(f"ERROR creating module lookup: {e}. Check modules file columns ('gene', 'module').")
        return []

    module_context_results = []
    for pair in pathway_results_subset: # Iterate through the provided pathway pairs
        ligand = pair.get('ligand')
        receptor_complex = str(pair.get('receptor', '')) # Ensure string
        
        if not ligand or not receptor_complex:
            continue # Skip if pair is missing essential info

        receptors = _split_receptor_complex(receptor_complex)
        ligand_module = module_lookup.get(ligand)
        
        # Collect receptor modules more simply, without counting components
        receptor_modules = []
        for r in receptors:
            r_module = module_lookup.get(r)
            if r_module is not None:
                receptor_modules.append(r_module)

        context = pair.copy() # Start with original pair data (like score)
        context.update({
            'ligand_module': ligand_module, 
            'receptor_modules': receptor_modules if receptor_modules else None, # Use None if empty
            # Removed 'num_receptor_components_in_modules' field
            # Removed 'is_same_module' field
            # Removed 'interaction_type' field
        })
        module_context_results.append(context)

    return module_context_results

def run_module_context_pipeline(spatial_df: pd.DataFrame, interactions_df: pd.DataFrame, modules_df: pd.DataFrame, full_pathway_results: Dict[str, Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """Runs the module context calculation for whole tissue and each layer."""
    logger.info("Running Module Context Pipeline...")
    results = {}
    if 'layer' not in spatial_df.columns:
        logger.error("'layer' column not found in spatial data. Cannot run module context per layer.")
        # Calculate only for whole tissue if pathway results exist
        pathway_subset = full_pathway_results.get('whole_tissue', {}).get('pathway_dominance', [])
        return {'whole_tissue': _calculate_module_context_for_scope(modules_df, pathway_subset)} 
        
    try:
        layers = list(spatial_df['layer'].astype(str).unique())
    except Exception as e:
        logger.error(f"Error converting layer names to string: {e}. Proceeding with original layer types.")
        layers = list(spatial_df['layer'].unique())

    scopes = ['whole_tissue'] + layers

    for scope in scopes:
        logger.info(f"  Processing Module Context scope: {scope}")
        try:
            # Get the relevant pathway results list for this scope
            pathway_subset = full_pathway_results.get(str(scope), {}).get('pathway_dominance', []) 
            
            # Module context calculation doesn't directly need the spatial subset,
            # only the pathway results for the scope and the global modules_df.
            if not pathway_subset:
                 logger.warning(f"    Skipping empty/missing pathway results for scope: {scope}")
                 results[str(scope)] = []
                 continue
                 
            results[str(scope)] = _calculate_module_context_for_scope(modules_df, pathway_subset)
        except Exception as e:
            logger.exception(f"    ERROR processing module context for scope {scope}: {e}")
            results[str(scope)] = [] # Return empty list on error for this scope
        
    logger.info("Module Context Pipeline Complete.")
    return results 


# --- NEW: Adapted Pipeline for Custom DataFrames ---

async def run_analysis_pipeline_from_dataframes(
    all_spatial_df: pd.DataFrame, 
    interactions_df: pd.DataFrame, 
    modules_df: pd.DataFrame
) -> Dict[str, Dict[str, Any]]:
    """Runs the analysis pipeline stages directly on provided DataFrames for a custom scope."""
    logger.info("Running Analysis Pipeline from DataFrames (Custom Selection)...")
    scope_name = "custom_selection" # Define the scope name for results
    final_results = {scope_name: {}}

    # Simulate async nature if needed, though core funcs might be CPU-bound
    # await asyncio.sleep(0.1)

    try:
        # 1. Calculate Counts
        logger.info(f"  Calculating counts for {scope_name}...")
        count_results = _calculate_lr_counts_for_scope(all_spatial_df, interactions_df)
        final_results[scope_name]['ligand_receptor_counts'] = count_results
        logger.info(f"  Counts complete for {scope_name}.")

        # 2. Calculate Pathway Dominance
        logger.info(f"  Calculating pathway dominance for {scope_name}...")
        pathway_results_list = _calculate_pathway_dominance_for_scope(all_spatial_df, interactions_df)
        final_results[scope_name]['pathway_dominance'] = pathway_results_list
        logger.info(f"  Pathway dominance complete for {scope_name}.")

        # 3. Calculate Module Context
        logger.info(f"  Calculating module context for {scope_name}...")
        # Pass the pathway results calculated for this specific scope
        module_context_list = _calculate_module_context_for_scope(modules_df, pathway_results_list)
        final_results[scope_name]['module_context'] = module_context_list
        logger.info(f"  Module context complete for {scope_name}.")

        logger.info(f"Analysis Pipeline from DataFrames Complete for {scope_name}.")
        return final_results

    except Exception as e:
        # Log the specific error from the pipeline stage
        logger.exception(f"ERROR during analysis pipeline from DataFrames for {scope_name}: {e}")
        # Re-raise the exception so the service layer can catch it and update job status
        raise e 

# --- Make sure analysis_logic/__init__.py exports the new function --- 