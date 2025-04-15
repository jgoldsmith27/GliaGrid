import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import argparse
import gc

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent

def load_module_data(module_file_path: Path) -> pd.DataFrame:
    """Load the cell module data."""
    if not module_file_path.exists():
        raise FileNotFoundError(f"Module data not found at {module_file_path}")
    
    logger.info(f"Loading module data from {module_file_path}...")
    # Assuming header=None and specific column names based on original script
    module_data = pd.read_csv(module_file_path, 
                              header=None, 
                              names=['gene', 'ensembl', 'score1', 'score2', 'module'])
    
    # Validate required columns
    required_columns = ['gene', 'module']
    missing_columns = [col for col in required_columns if col not in module_data.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in module data: {missing_columns}")
    
    # Ensure gene column is suitable for merging (string type)
    module_data['gene'] = module_data['gene'].astype(str)
    
    logger.info(f"Loaded {len(module_data)} module entries")
    return module_data

def load_layer_pathway_data(layer_pathway_dir: Path) -> pd.DataFrame:
    """Load ligand-receptor pair data for a specific layer's pathway analysis."""
    pairs_path = layer_pathway_dir / "ligand_receptor_pairs.csv"

    if not pairs_path.exists():
        logger.warning(f"Ligand-receptor pair data not found for layer {layer_pathway_dir.parent.name} at {pairs_path}")
        return None

    logger.info(f"Loading pair data for layer {layer_pathway_dir.parent.name}...")
    pairs_data = pd.read_csv(pairs_path)

    # Basic validation
    required_columns = ['ligand', 'receptor', 'normalized_score']
    missing_columns = [col for col in required_columns if col not in pairs_data.columns]
    if missing_columns:
        logger.warning(f"Pairs data for layer {layer_pathway_dir.parent.name} is missing columns: {missing_columns}")
        return None
    if pairs_data.empty:
        logger.warning(f"Pairs data is empty for layer {layer_pathway_dir.parent.name}")
        return None
        
    return pairs_data

def analyze_module_context_for_layer(pairs_data: pd.DataFrame, module_data: pd.DataFrame, layer_name: str, num_top_pairs: int):
    """Analyze module context for the top ligand-receptor pairs of a layer."""
    logger.info(f"Analyzing module context for top {num_top_pairs} pairs in {layer_name}...")
    
    # Get top N pairs based on normalized_score
    top_pairs = pairs_data.sort_values('normalized_score', ascending=False).head(num_top_pairs)
    
    if top_pairs.empty:
        logger.warning(f"No pairs data available to analyze module context for {layer_name}.")
        return None, None # Return None for both context_df and pivot_df

    # Create module lookup for faster access
    module_lookup = module_data.set_index('gene')['module'].to_dict()

    module_context = []
    interactions = []
    for _, pair in top_pairs.iterrows():
        ligand = pair['ligand']
        # Handle potential non-string receptors before splitting
        receptor_complex = str(pair['receptor']) 
        receptors = receptor_complex.split('_')
        score = pair['normalized_score']
        
        ligand_module = module_lookup.get(ligand)
        
        receptor_modules = []
        valid_receptor_components_in_module_data = []
        for r in receptors:
            r_module = module_lookup.get(r)
            if r_module is not None:
                receptor_modules.append(r_module)
                valid_receptor_components_in_module_data.append(r)
                # Record interaction for heatmap
                if ligand_module is not None:
                    interactions.append((ligand_module, r_module, score))

        # Determine interaction type
        is_same_module = None
        interaction_type = 'Unknown'
        if ligand_module is not None and receptor_modules:
            is_same_module = all(m == ligand_module for m in receptor_modules)
            interaction_type = 'Intra-module' if is_same_module else 'Inter-module'
        elif ligand_module is None and receptor_modules: 
            interaction_type = 'Receptor Module Known' # Ligand module unknown
        elif ligand_module is not None and not receptor_modules:
            interaction_type = 'Ligand Module Known' # Receptor module(s) unknown
            
        module_context.append({
            'ligand': ligand,
            'receptor': receptor_complex,
            'normalized_score': score,
            'ligand_module': ligand_module,
            'receptor_modules': receptor_modules if receptor_modules else None, # Store None if empty list
            'num_receptor_components_in_modules': len(receptor_modules),
            'is_same_module': is_same_module,
            'interaction_type': interaction_type
        })

    # Convert context results to DataFrame
    context_df = pd.DataFrame(module_context)

    # Create interaction pivot table for heatmap
    pivot_df = None
    if interactions:
        interaction_df = pd.DataFrame(interactions, columns=['ligand_module', 'receptor_module', 'score'])
        # Aggregate scores (sum) for pairs involving the same modules
        pivot_df = interaction_df.pivot_table(values='score', 
                                            index='ligand_module', 
                                            columns='receptor_module', 
                                            aggfunc='sum', 
                                            fill_value=0) # Fill NaN with 0
    else:
        logger.warning(f"No module interactions found to create heatmap for {layer_name}.")

    return context_df, pivot_df

def create_module_visualizations(context_df: pd.DataFrame, pivot_df: pd.DataFrame, layer_output_dir: Path, layer_name: str):
    """Create and save visualizations for module context analysis."""
    logger.info(f"Creating module context visualizations for {layer_name}...")
    layer_output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Module interaction heatmap
    plt.figure(figsize=(12, 10))
    if pivot_df is not None and not pivot_df.empty:
        sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='.2f', linewidths=.5)
        plt.title(f'Module Interaction Strength (Top Pairs - {layer_name})')
        plt.xlabel('Receptor Module')
        plt.ylabel('Ligand Module')
    else:
        plt.text(0.5, 0.5, 'No module interactions found for heatmap', ha='center', va='center')
    plt.tight_layout()
    plt.savefig(layer_output_dir / "module_interactions_heatmap.png", dpi=300)
    plt.close()

    # 2. Interaction type distribution
    if not context_df.empty:
        plt.figure(figsize=(8, 6))
        interaction_counts = context_df['interaction_type'].value_counts()
        if not interaction_counts.empty:
            sns.barplot(x=interaction_counts.index, y=interaction_counts.values, palette="Set2")
            plt.title(f'Distribution of Interaction Types (Top Pairs - {layer_name})')
            plt.xticks(rotation=45, ha='right')
            plt.ylabel('Count')
        else:
             plt.text(0.5, 0.5, 'No interaction type data to plot', ha='center', va='center')
        plt.tight_layout()
        plt.savefig(layer_output_dir / "interaction_types_distribution.png", dpi=300)
        plt.close()
    else:
         logger.warning(f"No context data available for interaction type plot for {layer_name}")

    logger.info(f"Module context visualizations saved to {layer_output_dir}")

def main(args):
    """Main function to run module context analysis for each layer."""
    try:
        pathway_results_base_dir = PROJECT_ROOT / args.pathway_results_base_dir
        module_file_path = PROJECT_ROOT / args.module_file
        results_base_dir = PROJECT_ROOT / args.results_base_dir
        
        module_data = load_module_data(module_file_path)

        # Find all layer directories in the pathway results path
        # Assumes structure like: .../results/[LayerName]/pathway_analysis/
        layer_pathway_dirs = [d for d in pathway_results_base_dir.glob('*/pathway_analysis') if d.is_dir()]
        
        if not layer_pathway_dirs:
            logger.error(f"No layer pathway analysis directories found in {pathway_results_base_dir}/*/pathway_analysis/. Did you run analyze_pathway_dominance_by_layer.py first?")
            return

        logger.info(f"Found {len(layer_pathway_dirs)} layer directories to process: {[d.parent.name for d in layer_pathway_dirs]}")

        # --- Loop through each layer directory ---
        for layer_pathway_dir in layer_pathway_dirs:
            layer_name = layer_pathway_dir.parent.name # Get layer name from parent dir
            logger.info(f"\n{'='*20} Processing Layer: {layer_name} {'='*20}")

            # Load pathway pairs data for this layer
            pairs_data = load_layer_pathway_data(layer_pathway_dir)
            if pairs_data is None:
                continue # Skip layer if pathway data loading failed

            # Analyze module context for this layer
            context_df, pivot_df = analyze_module_context_for_layer(pairs_data, module_data, layer_name, args.num_top_pairs)
            
            if context_df is None:
                logger.warning(f"Module context analysis failed or yielded no results for layer {layer_name}. Skipping save and visualization.")
                continue # Skip if analysis failed

            # Define layer-specific output directory for module analysis
            layer_module_results_dir = results_base_dir / layer_name / "module_analysis"
            layer_module_results_dir.mkdir(parents=True, exist_ok=True)

            # Save results
            context_output_path = layer_module_results_dir / "module_context_analysis.csv"
            pivot_output_path = layer_module_results_dir / "module_interaction_matrix.csv"
            summary_stats_path = layer_module_results_dir / "module_analysis_summary.csv"
            
            context_df.to_csv(context_output_path, index=False)
            logger.info(f"Module context results saved for {layer_name}: {context_output_path}")
            
            if pivot_df is not None:
                pivot_df.to_csv(pivot_output_path)
                logger.info(f"Module interaction matrix saved for {layer_name}: {pivot_output_path}")
            else:
                # Create an empty file to indicate no matrix was generated
                pivot_output_path.touch()
                logger.info(f"No module interaction matrix generated for {layer_name}, empty file created at: {pivot_output_path}")

            # Calculate and save summary statistics
            if not context_df.empty:
                 summary_stats = {
                     'Metric': list(context_df['interaction_type'].value_counts().index),
                     'Count': list(context_df['interaction_type'].value_counts().values)
                 }
                 summary_df = pd.DataFrame(summary_stats)
                 summary_df.to_csv(summary_stats_path, index=False)
                 logger.info(f"Summary statistics saved for {layer_name}: {summary_stats_path}")
            else:
                 summary_stats_path.touch()
                 logger.info(f"No context data for summary stats for {layer_name}, empty file created at: {summary_stats_path}")
                 
            # Create visualizations
            create_module_visualizations(context_df, pivot_df, layer_module_results_dir, layer_name)
            
            # Memory management
            del pairs_data, context_df, pivot_df
            gc.collect()

        logger.info("\n=== Layer-specific module context analysis complete ===")

    except FileNotFoundError as e:
        logger.error(f"Input file not found: {e}")
    except ValueError as e:
        logger.error(f"Data validation error: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred during layer-specific module context analysis: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze module context for top LR pairs in each spatial layer.')
    parser.add_argument('--pathway-results-base-dir', type=str, 
                        default='results',
                        help='Base directory relative to project root containing layer-specific pathway analysis results.')
    parser.add_argument('--module-file', type=str, 
                        default='data/cell_modules/human_cerebral_cortex_condition_control_integrated_zpower_matrix.csv', # Adjusted path assumption
                        help='Path relative to project root for the cell module data file.')
    parser.add_argument('--results-base-dir', type=str, 
                        default='results',
                        help='Base directory relative to project root where layer-specific module analysis results will be saved.')
    parser.add_argument('--num-top-pairs', type=int, default=15,
                        help='Number of top ligand-receptor pairs (by normalized score) to analyze for module context.')
    
    args = parser.parse_args()
    main(args)
