import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
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

def load_layer_data(layer_processed_dir: Path):
    """Load ligand and receptor data for a specific layer."""
    ligand_path = layer_processed_dir / "ligand_spots.csv"
    receptor_path = layer_processed_dir / "receptor_spots.csv"

    if not ligand_path.exists():
        logger.warning(f"Ligand data not found for layer {layer_processed_dir.name} at {ligand_path}")
        return None, None
    if not receptor_path.exists():
        logger.warning(f"Receptor data not found for layer {layer_processed_dir.name} at {receptor_path}")
        return None, None

    logger.info(f"Loading data for layer {layer_processed_dir.name}...")
    ligand_data = pd.read_csv(ligand_path)
    receptor_data = pd.read_csv(receptor_path)

    # Basic validation
    if ligand_data.empty or 'gene' not in ligand_data.columns or 'percent_spots' not in ligand_data.columns:
        logger.warning(f"Ligand data invalid or empty for layer {layer_processed_dir.name}")
        return None, None
    if receptor_data.empty or 'gene' not in receptor_data.columns or 'percent_spots' not in receptor_data.columns:
        logger.warning(f"Receptor data invalid or empty for layer {layer_processed_dir.name}")
        return None, None
        
    return ligand_data, receptor_data

def load_interaction_data() -> pd.DataFrame:
    """Load the interaction data from CellChatDB."""
    interaction_path = PROJECT_ROOT / "data/cellchat/processed/interaction.csv"
    if not interaction_path.exists():
        raise FileNotFoundError(f"Interaction data not found at {interaction_path}")
    
    logger.info("Loading global interaction data...")
    interactions = pd.read_csv(interaction_path)
    # Add a check for is_complex_receptor if it exists, otherwise create it
    if 'is_complex_receptor' not in interactions.columns:
         logger.info("'is_complex_receptor' column not found in interactions. Adding it based on '_' in receptor name.")
         interactions['is_complex_receptor'] = interactions['receptor'].str.contains('_')
    return interactions

def analyze_pathway_for_layer(ligand_data, receptor_data, interactions, layer_name, min_expression_threshold=0.001):
    """Perform pathway dominance analysis for a single layer."""
    logger.info(f"Filtering and normalizing expression data for {layer_name}...")
    
    # 1. Filter expression data
    ligand_data_filtered = ligand_data[ligand_data['percent_spots'] > min_expression_threshold].copy()
    receptor_data_filtered = receptor_data[receptor_data['percent_spots'] > min_expression_threshold].copy()

    logger.info(f"Filtered ligands ({layer_name}): {len(ligand_data_filtered)}/{len(ligand_data)} ({len(ligand_data_filtered)/len(ligand_data):.1%})" if len(ligand_data) > 0 else f"Filtered ligands ({layer_name}): 0/0")
    logger.info(f"Filtered receptors ({layer_name}): {len(receptor_data_filtered)}/{len(receptor_data)} ({len(receptor_data_filtered)/len(receptor_data):.1%})" if len(receptor_data) > 0 else f"Filtered receptors ({layer_name}): 0/0")

    if ligand_data_filtered.empty or receptor_data_filtered.empty:
        logger.warning(f"Insufficient filtered data for {layer_name}. Skipping pathway analysis.")
        return None, None

    # 2. Normalize expression values (0-1 scale)
    l_scaler = MinMaxScaler()
    r_scaler = MinMaxScaler()
    ligand_data_filtered['normalized_expression'] = l_scaler.fit_transform(ligand_data_filtered[['percent_spots']])
    receptor_data_filtered['normalized_expression'] = r_scaler.fit_transform(receptor_data_filtered[['percent_spots']])

    # Create lookups
    ligand_expr = {row['gene']: {'percent': row['percent_spots'], 'normalized': row['normalized_expression']} 
                   for _, row in ligand_data_filtered.iterrows()}
    receptor_expr = {row['gene']: {'percent': row['percent_spots'], 'normalized': row['normalized_expression']} 
                     for _, row in receptor_data_filtered.iterrows()}

    logger.info(f"Analyzing interaction pairs for {layer_name}...")
    # 3. Analyze interaction pairs
    valid_pairs = []
    for _, row in interactions.iterrows():
        ligand = row['ligand']
        receptor_complex = row['receptor']
        pathway = row['pathway_name']
        is_complex = row.get('is_complex_receptor', '_' in receptor_complex) # Use existing or derive
        
        if ligand not in ligand_expr:
            continue
            
        receptors = receptor_complex.split('_')
        receptor_scores = {'percent': [], 'normalized': []}
        valid_receptors_in_pair = []
        
        for receptor in receptors:
            if receptor in receptor_expr:
                receptor_scores['percent'].append(receptor_expr[receptor]['percent'])
                receptor_scores['normalized'].append(receptor_expr[receptor]['normalized'])
                valid_receptors_in_pair.append(receptor)
        
        if not valid_receptors_in_pair:
            continue
        
        avg_receptor_percent = sum(receptor_scores['percent']) / len(receptor_scores['percent'])
        avg_receptor_normalized = sum(receptor_scores['normalized']) / len(receptor_scores['normalized'])
        
        raw_score = ligand_expr[ligand]['percent'] * avg_receptor_percent
        normalized_score = ligand_expr[ligand]['normalized'] * avg_receptor_normalized
        
        valid_pairs.append({
            'ligand': ligand,
            'receptor': receptor_complex,
            'valid_receptors': '_'.join(valid_receptors_in_pair),
            'pathway': pathway,
            'is_complex_receptor': is_complex, # Add complex flag
            'ligand_percent': ligand_expr[ligand]['percent'],
            'receptor_percent': avg_receptor_percent,
            'ligand_normalized': ligand_expr[ligand]['normalized'],
            'receptor_normalized': avg_receptor_normalized,
            'raw_score': raw_score,
            'normalized_score': normalized_score
        })

    if not valid_pairs:
        logger.warning(f"No valid ligand-receptor pairs found for {layer_name} after filtering.")
        return None, None
        
    logger.info(f"Found {len(valid_pairs)} valid ligand-receptor pairs for {layer_name}")

    # 4. Convert to DataFrame and sort
    pairs_df = pd.DataFrame(valid_pairs)
    pairs_df = pairs_df.sort_values('normalized_score', ascending=False)

    logger.info(f"Analyzing pathways for {layer_name}...")
    # 5. Group by pathway
    pathway_summary = pairs_df.groupby('pathway').agg(
        avg_score=('normalized_score', 'mean'),
        total_score=('normalized_score', 'sum'),
        pair_count=('normalized_score', 'count'),
        unique_ligands=('ligand', 'nunique'),
        # Use a lambda to count unique genes within the 'valid_receptors' strings
        unique_receptors=('valid_receptors', lambda x: len(set('_'.join(x).split('_'))))
    ).reset_index()
    
    # Calculate weighted score
    pathway_summary['weighted_score'] = pathway_summary['total_score'] * np.sqrt(pathway_summary['pair_count'])
    pathway_summary = pathway_summary.sort_values('weighted_score', ascending=False)

    return pairs_df, pathway_summary

def create_visualizations(pairs_df, pathway_summary, layer_output_dir, layer_name):
    """Create and save visualizations for a layer's pathway analysis."""
    logger.info(f"Creating visualizations for {layer_name}...")
    layer_output_dir.mkdir(parents=True, exist_ok=True)

    # Top pairs bar chart
    plt.figure(figsize=(12, 8))
    top_pairs = pairs_df.head(15)
    # Assign y to hue and set legend=False as per FutureWarning
    sns.barplot(x='normalized_score', y=top_pairs['ligand'] + ' - ' + top_pairs['receptor'], 
                hue=top_pairs['ligand'] + ' - ' + top_pairs['receptor'], 
                data=top_pairs, palette="viridis", legend=False)
    plt.title(f'Top 15 Ligand-Receptor Pairs by Normalized Score ({layer_name})')
    plt.tight_layout()
    plt.savefig(layer_output_dir / "top_pairs.png", dpi=300)
    plt.close()

    # Top pathways bar chart
    plt.figure(figsize=(12, 8))
    top_pathways = pathway_summary.head(10)
    # Assign y to hue and set legend=False as per FutureWarning
    sns.barplot(x='weighted_score', y='pathway', hue='pathway', 
                data=top_pathways, palette="magma", legend=False)
    plt.title(f'Top 10 Pathways by Weighted Score ({layer_name})')
    plt.tight_layout()
    plt.savefig(layer_output_dir / "top_pathways.png", dpi=300)
    plt.close()

    # --- REMOVED Pathway Metrics Heatmap ---
    # logger.info(f"Skipping Pathway Metrics Heatmap for {layer_name}.")

    # --- REMOVED Top Pathway Composition Chart ---
    # logger.info(f"Skipping Top Pathway Composition Chart for {layer_name}.")
        
    logger.info(f"Visualizations saved to {layer_output_dir}")

def main(args):
    """Main function to run pathway dominance analysis for each layer."""
    try:
        processed_data_dir = PROJECT_ROOT / args.processed_data_base_dir
        results_base_dir = PROJECT_ROOT / args.results_base_dir
        interactions = load_interaction_data()

        # Find all layer directories in the processed data path
        layer_dirs = [d for d in processed_data_dir.iterdir() if d.is_dir()]
        if not layer_dirs:
            logger.error(f"No layer directories found in {processed_data_dir}. Did you run analyze_spatial_expression_by_layer.py first?")
            return

        logger.info(f"Found {len(layer_dirs)} layer directories to process: {[d.name for d in layer_dirs]}")

        # --- Loop through each layer directory ---
        for layer_dir in layer_dirs:
            layer_name = layer_dir.name
            logger.info(f"\n{'='*20} Processing Layer: {layer_name} {'='*20}")

            # Load data for this layer
            ligand_data, receptor_data = load_layer_data(layer_dir)
            if ligand_data is None or receptor_data is None:
                continue # Skip layer if data loading failed

            # Analyze pathways for this layer
            pairs_df, pathway_summary = analyze_pathway_for_layer(ligand_data, receptor_data, interactions, layer_name, args.min_expression)
            
            if pairs_df is None or pathway_summary is None:
                logger.warning(f"Pathway analysis failed or yielded no results for layer {layer_name}. Skipping save and visualization.")
                continue # Skip if analysis failed

            # Define layer-specific output directory
            layer_results_dir = results_base_dir / layer_name / "pathway_analysis"
            layer_results_dir.mkdir(parents=True, exist_ok=True)

            # Save results
            pairs_output_path = layer_results_dir / "ligand_receptor_pairs.csv"
            pathway_output_path = layer_results_dir / "pathway_summary.csv"
            
            pairs_df.to_csv(pairs_output_path, index=False)
            pathway_summary.to_csv(pathway_output_path, index=False)
            logger.info(f"Results for {layer_name} saved:")
            logger.info(f"- Pairs: {pairs_output_path}")
            logger.info(f"- Pathway Summary: {pathway_output_path}")

            # Create visualizations
            create_visualizations(pairs_df, pathway_summary, layer_results_dir, layer_name)
            
            # Memory management
            del ligand_data, receptor_data, pairs_df, pathway_summary
            gc.collect()

        logger.info("\n=== Layer-specific pathway analysis complete ===")

    except FileNotFoundError as e:
        logger.error(f"Input file not found: {e}")
    except ValueError as e:
        logger.error(f"Data validation error: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred during layer-specific pathway analysis: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze pathway dominance for each spatial layer.')
    parser.add_argument('--processed-data-base-dir', type=str, 
                        default='data/spatial/processed',
                        help='Base directory relative to project root containing layer-specific processed data (ligand/receptor spots).')
    parser.add_argument('--results-base-dir', type=str, 
                        default='results',
                        help='Base directory relative to project root where layer-specific results will be saved.')
    parser.add_argument('--min-expression', type=float, default=0.001,
                        help='Minimum percentage of spots a gene must be present in to be considered expressed (0.001 = 0.1%%).')
    
    args = parser.parse_args()
    main(args)
