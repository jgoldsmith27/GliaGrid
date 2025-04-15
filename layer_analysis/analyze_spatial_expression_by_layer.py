import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Set
import scipy.sparse as sp
import gc
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent

def load_interaction_data() -> pd.DataFrame:
    """Load the interaction data from CellChatDB."""
    interaction_path = PROJECT_ROOT / "data/cellchat/processed/interaction.csv"
    if not interaction_path.exists():
        raise FileNotFoundError(f"Interaction data not found at {interaction_path}")
    
    logger.info("Loading interaction data...")
    interactions = pd.read_csv(interaction_path)
    
    # Validate required columns
    required_columns = ['ligand', 'receptor']
    missing_columns = [col for col in required_columns if col not in interactions.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in interaction data: {missing_columns}")
    
    # Validate no empty values in key columns
    if interactions['ligand'].isna().any() or interactions['receptor'].isna().any():
        raise ValueError("Found empty values in ligand or receptor columns")
    
    logger.info(f"Loaded {len(interactions)} interactions")
    return interactions

def load_spatial_data(spatial_csv_path: Path) -> pd.DataFrame:
    """Load the spatial data from the specified CSV file."""
    if not spatial_csv_path.exists():
        raise FileNotFoundError(f"Spatial data not found at {spatial_csv_path}")
    
    logger.info(f"Loading spatial data from {spatial_csv_path}...")
    # Read the CSV file
    df = pd.read_csv(spatial_csv_path)
    
    # Validate required columns
    required_columns = ['geneID', 'bin1_ID', 'layer'] # Added 'layer'
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in spatial data: {missing_columns}")
    
    # Validate no empty values in key columns
    if df['geneID'].isna().any() or df['bin1_ID'].isna().any() or df['layer'].isna().any():
        raise ValueError("Found empty values in geneID, bin1_ID, or layer columns")
        
    # Check for and handle potential NaN string values in 'layer'
    if df['layer'].astype(str).str.contains('nan', case=False).any():
        logger.warning("Found potential 'nan' strings in layer column. Attempting to handle.")
        # Convert to string, replace 'nan' (case-insensitive) with actual NaN, then drop rows with NaN layer
        df['layer'] = df['layer'].astype(str)
        df['layer'] = df['layer'].replace('nan', np.nan, regex=True)
        original_rows = len(df)
        df.dropna(subset=['layer'], inplace=True)
        rows_dropped = original_rows - len(df)
        if rows_dropped > 0:
            logger.warning(f"Dropped {rows_dropped} rows with NaN or 'nan' string in layer column.")
            
    # Print available columns
    logger.info("Available columns in spatial data:")
    for col in df.columns:
        logger.info(f"- {col}")
    
    logger.info(f"Loaded {len(df)} rows of spatial data")
    return df

def extract_unique_genes(interactions: pd.DataFrame) -> Dict[str, Set[str]]:
    """Extract unique ligands and receptors from interaction data."""
    logger.info("\nExtracting unique ligands and receptors from CellChatDB interactions...")
    
    # Get unique ligands
    ligands = set(interactions['ligand'].str.upper())
    logger.info(f"Total unique ligands in CellChatDB: {len(ligands)}")
    logger.info("Sample of ligands (first 10):")
    for ligand in sorted(list(ligands))[:10]:
        logger.info(f"- {ligand}")
    
    # Get unique receptors (split complexes)
    all_receptors = interactions['receptor'].str.upper()
    receptors = set()
    complex_receptors = set()
    
    for receptor in all_receptors:
        if '_' in receptor:
            complex_receptors.add(receptor)
            # Split complex into individual genes
            receptors.update(receptor.split('_'))
        else:
            receptors.add(receptor)
    
    logger.info(f"\nTotal unique receptor genes in CellChatDB: {len(receptors)}")
    logger.info("Sample of receptor genes (first 10):")
    for receptor in sorted(list(receptors))[:10]:
        logger.info(f"- {receptor}")
    
    logger.info(f"\nTotal receptor complexes: {len(complex_receptors)}")
    logger.info("Sample of receptor complexes (first 10):")
    for complex_r in sorted(list(complex_receptors))[:10]:
        logger.info(f"- {complex_r}")
    
    # Validate no empty sets
    if not ligands or not receptors:
        raise ValueError("No ligands or receptors found in interaction data")
    
    return {
        'ligands': ligands,
        'receptors': receptors
    }

def analyze_gene_expression(
    df_layer: pd.DataFrame,
    genes: Set[str],
    gene_type: str,
    current_layer: str,
    chunk_size: int = 1000
) -> pd.DataFrame:
    """Analyze presence of genes in spatial spots for a specific layer."""
    logger.info(f"\nAnalyzing {gene_type} presence in spatial spots for {current_layer}...")
    logger.info(f"Total {gene_type} to analyze: {len(genes)}")
    
    # Validate input
    if not genes:
        raise ValueError(f"No {gene_type} provided for analysis")
    if df_layer.empty:
        logger.warning(f"Empty spatial data provided for layer {current_layer}. Skipping analysis for {gene_type}.")
        return pd.DataFrame() # Return empty DataFrame if no data for this layer
    
    results = []
    
    # Convert gene names to uppercase for case-insensitive matching
    genes_upper = {g.upper() for g in genes}
    # Ensure geneID is string before applying .str
    df_layer['geneID'] = df_layer['geneID'].astype(str) 
    df_layer['geneID_upper'] = df_layer['geneID'].str.upper()
    
    # Find valid genes present in this layer's data
    valid_genes = set(df_layer['geneID_upper'].unique()) & genes_upper
    missing_genes = genes_upper - valid_genes
    
    logger.info(f"\nGene matching summary for {gene_type} in {current_layer}:")
    logger.info(f"- Total {gene_type} in CellChatDB: {len(genes)}")
    logger.info(f"- Found in {current_layer} spatial data: {len(valid_genes)}")
    logger.info(f"- Missing in {current_layer} spatial data: {len(missing_genes)}")
    
    if len(valid_genes) == 0:
        logger.warning(f"No {gene_type} found in the dataset for {current_layer}!")
        return pd.DataFrame()
    
    # Get total number of unique spots in this layer
    total_spots = len(df_layer['bin1_ID'].unique())
    logger.info(f"\nTotal number of spatial spots in {current_layer}: {total_spots}")
    
    # Validate total spots
    if total_spots == 0:
        logger.warning(f"No spatial spots found in {current_layer}. Skipping gene analysis.")
        return pd.DataFrame()
    
    # Process genes in chunks
    valid_genes = list(valid_genes)
    n_genes = len(valid_genes)
    n_chunks = (n_genes + chunk_size - 1) // chunk_size
    
    logger.info(f"\nProcessing {n_genes} {gene_type} in {n_chunks} chunks for {current_layer}...")
    
    for chunk_idx in range(n_chunks):
        chunk_start = chunk_idx * chunk_size
        chunk_end = min(chunk_start + chunk_size, n_genes)
        chunk_genes = valid_genes[chunk_start:chunk_end]
        
        logger.info(f"\nProcessing chunk {chunk_idx + 1}/{n_chunks} ({chunk_start+1}-{chunk_end} of {n_genes} {gene_type}) for {current_layer}")
        
        for i, gene in enumerate(chunk_genes, 1):
            # Get spots where this gene is present within the current layer
            gene_spots = df_layer[df_layer['geneID_upper'] == gene]['bin1_ID'].unique()
            n_spots = len(gene_spots)
            percent_spots = (n_spots / total_spots) * 100 if total_spots > 0 else 0
            
            # Validate spot counts
            if n_spots > total_spots:
                raise ValueError(f"Gene {gene} ({current_layer}) has more spots ({n_spots}) than total spots ({total_spots})")
            
            # Find original case of gene name (handle potential missing gene after filtering)
            orig_gene_series = df_layer[df_layer['geneID_upper'] == gene]['geneID']
            orig_gene = orig_gene_series.iloc[0] if not orig_gene_series.empty else gene # Fallback to upper if original case missing
            
            # Log progress for each gene
            # logger.debug(f"Processing {gene_type} {i}/{len(chunk_genes)} ({current_layer}): {orig_gene} - {n_spots} spots ({percent_spots:.2f}%)")
            
            results.append({
                'gene': orig_gene,  # Use original case
                'spots_present': int(n_spots),
                'percent_spots': float(percent_spots)
            })
        
        # Clean up memory
        gc.collect()
    
    # Convert to DataFrame and sort
    if results:
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('spots_present', ascending=False)
        
        # Validate results
        if len(results_df) != len(valid_genes):
            # This might happen if a gene existed in the layer but had 0 spots? Check logic.
            # For now, just log a warning.
            logger.warning(f"Result count ({len(results_df)}) mismatch for {gene_type} in {current_layer}. Expected {len(valid_genes)}.")
        
        logger.info(f"\nTop 10 most prevalent {gene_type} in {current_layer} (by number of spots):")
        for _, row in results_df.head(10).iterrows():
            logger.info(f"- {row['gene']}: {row['spots_present']} spots ({row['percent_spots']:.2f}%)")
    else:
        results_df = pd.DataFrame()
    
    # Drop temporary column
    df_layer.drop(columns=['geneID_upper'], inplace=True, errors='ignore')

    return results_df

def main(args):
    """Main function to run the layer-specific spatial expression analysis."""
    try:
        # Define paths using PROJECT_ROOT
        spatial_data_path = PROJECT_ROOT / args.spatial_data_file
        base_output_dir = PROJECT_ROOT / "data/spatial/processed"
        
        # Load data
        interactions = load_interaction_data()
        df_all = load_spatial_data(spatial_data_path)
        
        # Extract unique genes from CellChatDB
        genes = extract_unique_genes(interactions)
        
        # Get unique layers from the data
        unique_layers = sorted(df_all['layer'].unique())
        logger.info(f"\nFound unique layers: {unique_layers}")
        
        # --- Loop through each layer ---
        for layer in unique_layers:
            logger.info(f"\n{'='*20} Processing Layer: {layer} {'='*20}")
            
            # Filter data for the current layer
            df_layer = df_all[df_all['layer'] == layer].copy() # Use copy to avoid SettingWithCopyWarning
            
            if df_layer.empty:
                logger.warning(f"No data found for layer {layer}. Skipping.")
                continue
                
            # Create layer-specific output directory
            layer_output_dir = base_output_dir / str(layer) # Ensure layer name is string for path
            layer_output_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Output directory for {layer}: {layer_output_dir}")

            # Analyze ligand presence for this layer
            logger.info(f"\n--- Starting Ligand Analysis for {layer} ---")
            ligand_results = analyze_gene_expression(df_layer, genes['ligands'], "ligands", layer)
            
            # Analyze receptor presence for this layer
            logger.info(f"\n--- Starting Receptor Analysis for {layer} ---")
            receptor_results = analyze_gene_expression(df_layer, genes['receptors'], "receptors", layer)
            
            # Save results for this layer
            ligand_path = layer_output_dir / "ligand_spots.csv"
            receptor_path = layer_output_dir / "receptor_spots.csv"
            
            # Save only if results are not empty
            if not ligand_results.empty:
                ligand_results.to_csv(ligand_path, index=False)
                logger.info(f"Ligand results saved for {layer}: {ligand_path}")
            else:
                logger.warning(f"No ligand results to save for {layer}.")

            if not receptor_results.empty:
                receptor_results.to_csv(receptor_path, index=False)
                logger.info(f"Receptor results saved for {layer}: {receptor_path}")
            else:
                 logger.warning(f"No receptor results to save for {layer}.")

            # Print summary for the layer
            logger.info(f"\nSummary for {layer}:")
            logger.info(f"- Total ligands analyzed: {len(ligand_results)}")
            logger.info(f"- Total receptors analyzed: {len(receptor_results)}")
            
            # Optional: Clear memory if processing many large layers
            del df_layer, ligand_results, receptor_results
            gc.collect()

        logger.info("\n=== Layer-specific analysis complete ===")
        
    except FileNotFoundError as e:
        logger.error(f"Input file not found: {e}")
    except ValueError as e:
        logger.error(f"Data validation error: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred during layer-specific spatial expression analysis: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze spatial gene expression layer by layer.')
    parser.add_argument('--spatial-data-file', type=str, 
                        default='data/spatial/raw/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop_bin100_clustered.csv',
                        help="Path relative to project root for the clustered spatial data CSV file (must contain a 'layer' column).")
    
    args = parser.parse_args()
    main(args)
