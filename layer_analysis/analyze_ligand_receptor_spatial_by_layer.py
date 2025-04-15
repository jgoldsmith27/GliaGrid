#!/usr/bin/env python3
"""
Analyze ligand-receptor pairs in spatial data layer by layer, 
calculating distance/density scores and visualizing interactions within each layer.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.stats import gaussian_kde
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
from pathlib import Path
import argparse
from tqdm import tqdm
import logging
import networkx as nx
from scipy import stats
import gc

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent

def load_layer_lr_pairs(layer_pathway_dir: Path, num_pairs: int):
    """Load the top N ligand-receptor pairs for a specific layer."""
    pairs_path = layer_pathway_dir / "ligand_receptor_pairs.csv"
    if not pairs_path.exists():
        logger.warning(f"Ligand-receptor pair data not found for layer {layer_pathway_dir.parent.name} at {pairs_path}")
        return None
    
    logger.info(f"Loading top {num_pairs} pairs for layer {layer_pathway_dir.parent.name}...")
    lr_pairs = pd.read_csv(pairs_path)
    
    # Validate necessary columns
    required_cols = ['ligand', 'receptor', 'normalized_score', 'is_complex_receptor', 'valid_receptors']
    missing_cols = [col for col in required_cols if col not in lr_pairs.columns]
    if missing_cols:
         logger.warning(f"Pairs file {pairs_path} missing required columns: {missing_cols}. Cannot perform spatial analysis.")
         return None
         
    # Sort by score and take top N
    lr_pairs = lr_pairs.sort_values('normalized_score', ascending=False).head(num_pairs)
    logger.info(f"Loaded {len(lr_pairs)} pairs for layer {layer_pathway_dir.parent.name}.")
    return lr_pairs

def load_spatial_data(file_path: Path):
    """Load the main clustered spatial data file."""
    logger.info(f"Loading main spatial data from {file_path}")
    if not file_path.exists():
        raise FileNotFoundError(f"Main spatial data file not found at {file_path}")
    df = pd.read_csv(file_path)
    # Validate required columns
    required_columns = ['geneID', 'x', 'y', 'layer'] 
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in main spatial data: {missing_columns}")
    # Handle potential NaN layers
    if df['layer'].isna().any():
        original_rows = len(df)
        df.dropna(subset=['layer'], inplace=True)
        logger.warning(f"Dropped {original_rows - len(df)} rows with NaN layer in main spatial data.")
    return df

def get_gene_spatial_data_for_layer(spatial_data_layer: pd.DataFrame, gene_name: str):
    """Extract spatial coordinates for a specific gene within a given layer's DataFrame."""
    gene_data = spatial_data_layer[spatial_data_layer['geneID'] == gene_name]
    if len(gene_data) == 0:
        # This is expected if a gene from the top pair isn't in this specific layer
        # logger.debug(f"No spatial data found for gene {gene_name} within the current layer.")
        return None
    
    # Return x, y coordinates
    return gene_data[['x', 'y']].values

def calculate_distance_density_score(ligand_coords, receptor_coords, ligand_name=None, receptor_name=None):
    """
    Calculate an interaction density score based on the product of 
    Ligand KDE and Receptor KDE evaluated at each receptor's location.
    Higher score = higher co-localization density of ligands and receptors.
    """
    if ligand_coords is None or receptor_coords is None:
        logger.error(f"Missing ligand or receptor coordinates for {ligand_name}-{receptor_name}.")
        return None
    
    # KDE requires at least 2 points for meaningful calculation
    if len(ligand_coords) < 2:
        logger.warning(f"Cannot calculate KDE score for {ligand_name}-{receptor_name}: requires >= 2 ligand points, found {len(ligand_coords)}.")
        return {"average_score": 0.0, "receptor_scores": {}} 
    if len(receptor_coords) < 2:
        logger.warning(f"Cannot calculate KDE score for {ligand_name}-{receptor_name}: requires >= 2 receptor points/centroids, found {len(receptor_coords)}.")
        return {"average_score": 0.0, "receptor_scores": {}} 
        
    if len(receptor_coords) == 0:
        logger.warning(f"Cannot calculate KDE score for {ligand_name}-{receptor_name}: no receptor points found.")
        return {"average_score": 0.0, "receptor_scores": {}}

    try:
        ligand_kde = gaussian_kde(ligand_coords.T)
        receptor_kde = gaussian_kde(receptor_coords.T)
        
        ligand_density_at_receptors = ligand_kde(receptor_coords.T)
        receptor_density_at_receptors = receptor_kde(receptor_coords.T)
        
        ligand_density_at_receptors[ligand_density_at_receptors < 0] = 0
        receptor_density_at_receptors[receptor_density_at_receptors < 0] = 0
        
        interaction_scores = ligand_density_at_receptors * receptor_density_at_receptors
        average_score = np.nanmean(interaction_scores) 
        
        receptor_scores_detail = { 
            f"receptor_{i}": {
                "interaction_density_score": interaction_scores[i],
                "ligand_density": ligand_density_at_receptors[i],
                "receptor_density": receptor_density_at_receptors[i]
                } 
            for i in range(len(receptor_coords))
        }
        
        logger.info(f"Calculated Interaction Density score for {ligand_name}-{receptor_name}. Avg score: {average_score:.4g}")
        
        return {
            "average_score": average_score,
            "receptor_scores": receptor_scores_detail
        }

    except Exception as e:
        logger.error(f"Failed to calculate Interaction Density score for {ligand_name}-{receptor_name}: {e}")
        return {"average_score": np.nan, "receptor_scores": {}}

def visualize_interaction(ligand_coords, receptor_clusters, score_data, ligand_name, receptor_name, layer_name, output_dir):
    """
    Visualize the spatial interaction density for a specific layer.
    """
    receptor_coords = np.array([cluster['centroid'] for cluster in receptor_clusters]) if receptor_clusters else np.array([])
    
    if ligand_coords is None or not isinstance(ligand_coords, np.ndarray) or ligand_coords.ndim != 2 or ligand_coords.shape[1] != 2:
        logger.warning(f"Invalid/missing ligand coords for {ligand_name}-{receptor_name} in {layer_name}. Skipping visualization.")
        return
    if receptor_coords is None or not isinstance(receptor_coords, np.ndarray) or receptor_coords.ndim != 2 or receptor_coords.shape[1] != 2:
         logger.warning(f"Invalid/missing receptor coords for {ligand_name}-{receptor_name} in {layer_name}. Skipping visualization.")
         return   
    if len(ligand_coords) == 0 or len(receptor_coords) == 0 or score_data is None or score_data.get('average_score') is None or np.isnan(score_data['average_score']):
        logger.warning(f"Cannot visualize density for {ligand_name}-{receptor_name} in {layer_name}: missing data/coords/valid score.")
        return
    
    avg_score_value = score_data['average_score']
    normalized_score = score_data.get('normalized_score_within_layer', None) # Use layer-specific normalized score
    score_display = f'Avg. Score: {avg_score_value:.4g}'
    if normalized_score is not None:
        score_display += f' (Layer Norm: {normalized_score:.3f})'

    plt.figure(figsize=(12, 10))
    ax_density = plt.gca()
    
    # Define grid boundaries based on layer data
    all_points = np.vstack([ligand_coords, receptor_coords])
    x_min, x_max = np.min(all_points[:, 0]), np.max(all_points[:, 0])
    y_min, y_max = np.min(all_points[:, 1]), np.max(all_points[:, 1])
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_padding = x_range * 0.05 if x_range > 0 else 10
    y_padding = y_range * 0.05 if y_range > 0 else 10
    x_min -= x_padding; x_max += x_padding
    y_min -= y_padding; y_max += y_padding
    
    grid_res = 150j
    x_grid, y_grid = np.mgrid[x_min:x_max:grid_res, y_min:y_max:grid_res]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    
    ligand_density = np.zeros(x_grid.shape)
    receptor_density = np.zeros(x_grid.shape)
    interaction_density = np.zeros(x_grid.shape)
    plotted_interaction_density = False

    # Calculate Densities (KDE)
    can_calc_ligand_density = False
    if len(ligand_coords) >= 2:
        try:
            ligand_kernel = gaussian_kde(ligand_coords.T)
            ligand_density = np.reshape(ligand_kernel(positions), x_grid.shape)
            can_calc_ligand_density = True
        except Exception as e:
            logger.warning(f"Could not compute ligand density for {ligand_name} ({layer_name}): {e}")

    can_calc_receptor_density = False
    if len(receptor_coords) >= 2:
        try:
            receptor_kernel = gaussian_kde(receptor_coords.T)
            receptor_density = np.reshape(receptor_kernel(positions), x_grid.shape)
            can_calc_receptor_density = True
        except Exception as e:
            logger.warning(f"Could not compute receptor density for {receptor_name} ({layer_name}) centroids: {e}")

    # Plot Interaction Density
    if can_calc_ligand_density and can_calc_receptor_density:
        try:
            interaction_density = ligand_density * receptor_density
            max_interaction = np.max(interaction_density)
            
            if max_interaction > 1e-15: # Lowered threshold
                interaction_density_normalized = interaction_density / max_interaction 
                im_interaction = ax_density.imshow(np.rot90(interaction_density_normalized), cmap='Greens', 
                                              extent=[x_min, x_max, y_min, y_max], 
                                              alpha=0.7, aspect='auto', vmin=0.01, vmax=1.0, zorder=1)
                plotted_interaction_density = True
                
                density_for_contours = interaction_density_normalized[interaction_density_normalized > 0.01]
                if len(density_for_contours) > 0:
                    density_thresholds = np.percentile(density_for_contours, [60, 80, 95]) 
                    contour_levels = sorted(list(set(d for d in density_thresholds if d > 1e-3)))
                    
                    if len(contour_levels) > 0:
                        contour = ax_density.contour(x_grid, y_grid, interaction_density_normalized, 
                                                levels=contour_levels, colors='purple', alpha=0.8, linewidths=0.8, zorder=2)
                        ax_density.clabel(contour, inline=True, fontsize=6, fmt='%.2f')
                else:
                     logger.info(f"Norm. interaction density low for contours ({ligand_name}-{receptor_name}, {layer_name})")
            else:
                 logger.info(f"Interaction density below plotting threshold (1e-15) ({ligand_name}-{receptor_name}, {layer_name})")
        except Exception as e:
            logger.warning(f"Could not compute/plot interaction density ({ligand_name}-{receptor_name}, {layer_name}): {e}")
            
    # Plot Raw Points
    point_size = 6
    ligand_alpha = 0.7
    receptor_alpha = 0.8
    ax_density.scatter(ligand_coords[:, 0], ligand_coords[:, 1], s=point_size, color='blue', alpha=ligand_alpha, label=f'{ligand_name} (Ligand)', marker='.', zorder=3)

    max_receptors_to_plot = 5000 
    indices_to_plot = np.random.choice(len(receptor_clusters), size=min(len(receptor_clusters), max_receptors_to_plot), replace=False) if len(receptor_clusters) > max_receptors_to_plot else range(len(receptor_clusters))
    plotted_receptor_label = False
    for i in indices_to_plot:
        centroid = receptor_clusters[i]['centroid']
        receptor_label = f'{receptor_name} (Receptor/Centroid)' if not plotted_receptor_label else ""
        ax_density.scatter(centroid[0], centroid[1], s=point_size+2, color='red', marker='.', label=receptor_label, zorder=4)
        plotted_receptor_label = True

    # Final plot settings
    plt.title(f'Interaction Density: {ligand_name} ({len(ligand_coords)}) - {receptor_name} ({len(receptor_coords)}) in {layer_name}\n'
              f'{score_display} (Higher = Denser Co-localization)')
    ax_density.set_xlabel('X Coordinate')
    ax_density.set_ylabel('Y Coordinate')
    ax_density.set_xlim(x_min, x_max)
    ax_density.set_ylim(y_min, y_max)
    ax_density.set_aspect('equal', adjustable='box')
    ax_density.grid(True, linestyle=':', alpha=0.3)
    
    handles, labels = ax_density.get_legend_handles_labels()
    if plotted_interaction_density:
        handles.append(plt.Rectangle((0,0),1,1,fc="lightgreen", alpha=0.7)); labels.append('Interaction Density')
        handles.append(plt.Line2D([0], [0], color='purple', lw=1)); labels.append('Interaction Contours')
    legend_fontsize = 'small' if len(handles) < 8 else 'x-small'
    ax_density.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=legend_fontsize)

    # Add scale bar (assuming 1 pixel = 0.5 micrometers)
    scale_bar_length_pixels = 100 # 50 um
    scale_x = x_min + (x_max - x_min) * 0.05
    scale_y = y_min + (y_max - y_min) * 0.05
    ax_density.plot([scale_x, scale_x + scale_bar_length_pixels], [scale_y, scale_y], 'k-', linewidth=3)
    ax_density.text(scale_x + scale_bar_length_pixels/2, scale_y - (y_max - y_min) * 0.03, '50 μm', ha='center', fontsize=10)
    
    density_output_path = output_dir / f"{ligand_name}_{receptor_name}_density_map.png"
    plt.savefig(density_output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Density map for {layer_name} saved to {density_output_path}")

def handle_complex_receptors_for_layer(receptor_name, is_complex, valid_receptors, spatial_data_layer):
    """Handle complex receptors within a specific layer's spatial data."""
    if not is_complex:
        coords = get_gene_spatial_data_for_layer(spatial_data_layer, receptor_name)
        if coords is None:
            return None, 0, {}
        clusters = [{'centroid': c, 'component_points': [c], 'component_names': [receptor_name]} for c in coords]
        return clusters, len(coords), {receptor_name: len(coords)}

    # Complex Receptor Logic
    components = valid_receptors.split('_') if isinstance(valid_receptors, str) else []
    if not components:
         logger.warning(f"Complex receptor {receptor_name} has no valid components listed.")
         return None, 0, {}
         
    # logger.info(f"Processing complex receptor {receptor_name} ({components}) for layer...")

    component_data = {}
    component_counts = {}
    all_components_found_in_layer = True
    for component in components:
        coords = get_gene_spatial_data_for_layer(spatial_data_layer, component)
        if coords is None or len(coords) == 0:
            # logger.debug(f"Component {component} not found in current layer for complex {receptor_name}.")
            all_components_found_in_layer = False
            component_counts[component] = 0
            component_data[component] = {'coords': np.array([])} 
        else:
            count = len(coords)
            component_counts[component] = count
            component_data[component] = {'coords': coords}
            
    if not all_components_found_in_layer:
         # Decide if *any* presence is enough or if *all* must be present *in the layer*.
         # Current logic: If any component is missing *in this layer*, cannot form complex *in this layer*.
         # logger.debug(f"Cannot form clusters for {receptor_name} in layer: component(s) missing.")
         return None, 0, component_counts
         
    # Determine anchor component (fewest points *in this layer*)
    min_count_in_layer = float('inf')
    anchor_component_name = None
    for comp, count in component_counts.items():
        if count > 0 and count < min_count_in_layer:
             min_count_in_layer = count
             anchor_component_name = comp
             
    if anchor_component_name is None: # Should not happen if all_components_found_in_layer is True
         logger.warning(f"Could not determine anchor component for {receptor_name} in layer.")
         return None, 0, component_counts
         
    other_component_names = [c for c in components if c != anchor_component_name]
    # logger.debug(f"Anchoring complex search on {anchor_component_name} ({component_counts[anchor_component_name]} points) for layer...")

    proximity_threshold = 100 # 50 micrometers = 100 pixels
    valid_clusters = []
    anchor_coords = component_data[anchor_component_name]['coords']

    for i, anchor_point in enumerate(anchor_coords):
        potential_cluster_points = {anchor_component_name: anchor_point}
        is_valid_cluster = True

        for other_comp_name in other_component_names:
            other_coords = component_data[other_comp_name]['coords']
            if len(other_coords) == 0: # Should be caught earlier, but double-check
                 is_valid_cluster = False; break
            distances = np.sqrt(np.sum((other_coords - anchor_point)**2, axis=1))
            within_threshold_indices = np.where(distances <= proximity_threshold)[0]

            if len(within_threshold_indices) == 0:
                is_valid_cluster = False; break

            closest_index = within_threshold_indices[np.argmin(distances[within_threshold_indices])]
            potential_cluster_points[other_comp_name] = other_coords[closest_index]

        if is_valid_cluster:
            cluster_points_array = np.array(list(potential_cluster_points.values()))
            centroid = np.mean(cluster_points_array, axis=0)
            valid_clusters.append({
                'centroid': centroid,
                'component_points': cluster_points_array,
                'component_names': list(potential_cluster_points.keys())
            })

    if valid_clusters:
        # logger.info(f"Found {len(valid_clusters)} complex locations for {receptor_name} in layer.")
        return valid_clusters, len(valid_clusters), component_counts
    else:
        # logger.debug(f"No valid component clusters found within threshold for {receptor_name} in layer.")
        return None, 0, component_counts

def normalize_interaction_scores_within_layer(scores):
    """Normalize interaction scores within a single layer using min-max scaling."""
    valid_scores = [s for s in scores if s is not None and not np.isnan(s)]
    if not valid_scores or len(valid_scores) < 2: # Need at least 2 scores to normalize meaningfully
        # Return 0.5 for all if normalization not possible (or keep as None?)
        return [0.5 if s is not None and not np.isnan(s) else None for s in scores] 
    
    score_array = np.array(valid_scores)
    min_val = np.min(score_array)
    max_val = np.max(score_array)
    
    normalized_scores = [None] * len(scores)
    valid_idx = 0
    for i, score in enumerate(scores):
         if score is not None and not np.isnan(score):
             if max_val > min_val:
                 normalized_scores[i] = float((score - min_val) / (max_val - min_val))
             else:
                 normalized_scores[i] = 1.0 # All scores were identical, set to max
             valid_idx += 1
             
    return normalized_scores

def main(args):
    """Main function to run spatial interaction analysis for each layer."""
    try:
        pathway_results_base_dir = PROJECT_ROOT / args.pathway_results_base_dir
        spatial_data_path = PROJECT_ROOT / args.spatial_data_file
        results_base_dir = PROJECT_ROOT / args.results_base_dir
        
        # Load the full spatial data once
        spatial_data_all = load_spatial_data(spatial_data_path)
        unique_layers = sorted(spatial_data_all['layer'].unique())
        logger.info(f"Found layers in main spatial data: {unique_layers}")

        # Find layer pathway analysis directories
        layer_pathway_dirs = [d for d in pathway_results_base_dir.glob('*/pathway_analysis') if d.is_dir()]
        if not layer_pathway_dirs:
            logger.error(f"No layer pathway analysis directories found in {pathway_results_base_dir}/*/pathway_analysis/.")
            return

        logger.info(f"Found {len(layer_pathway_dirs)} layer pathway directories to process.")

        # --- Loop through each layer --- 
        for layer_pathway_dir in layer_pathway_dirs:
            layer_name = layer_pathway_dir.parent.name
            logger.info(f"\n{'='*20} Processing Layer: {layer_name} {'='*20}")

            # Filter main spatial data for the current layer
            spatial_data_layer = spatial_data_all[spatial_data_all['layer'] == layer_name].copy()
            if spatial_data_layer.empty:
                logger.warning(f"No spatial data rows found for layer {layer_name} in main file. Skipping.")
                continue
            logger.info(f"Filtered spatial data for {layer_name}: {len(spatial_data_layer)} rows.")
            
            # Load the top LR pairs identified for this layer
            lr_pairs_layer = load_layer_lr_pairs(layer_pathway_dir, args.num_pairs)
            if lr_pairs_layer is None or lr_pairs_layer.empty:
                logger.warning(f"No valid LR pairs found for layer {layer_name}. Skipping spatial analysis.")
                continue

            # --- Stage 1: Calculate raw scores for this layer's pairs --- 
            raw_results_layer = []
            # Cache coordinates *within this layer* to avoid repeated filtering
            spatial_coords_cache_layer = {}

            logger.info(f"Calculating raw spatial scores for {len(lr_pairs_layer)} pairs in {layer_name}...")
            for index, row in tqdm(lr_pairs_layer.iterrows(), total=len(lr_pairs_layer), desc=f"Spatial Scores {layer_name}"):
                ligand = row['ligand']
                receptor = row['receptor']
                is_complex = row['is_complex_receptor']
                valid_receptors = row['valid_receptors']
                
                # Get ligand coordinates *from the layer-specific dataframe*
                if ligand not in spatial_coords_cache_layer:
                    spatial_coords_cache_layer[ligand] = get_gene_spatial_data_for_layer(spatial_data_layer, ligand)
                ligand_coords = spatial_coords_cache_layer[ligand]
                num_ligands_found = len(ligand_coords) if ligand_coords is not None else 0

                # Get receptor coordinates/clusters *from the layer-specific dataframe*
                receptor_clusters, num_centroids, component_counts = handle_complex_receptors_for_layer(
                    receptor, is_complex, valid_receptors, spatial_data_layer
                )

                receptor_coords = None
                if receptor_clusters is not None and len(receptor_clusters) > 0:
                    receptor_coords = np.array([cluster['centroid'] for cluster in receptor_clusters])

                avg_score_value = np.nan
                score_details = {}
                if ligand_coords is not None and len(ligand_coords) > 0 and receptor_coords is not None and len(receptor_coords) > 0:
                    score_data = calculate_distance_density_score(
                        ligand_coords, receptor_coords, ligand_name=ligand, receptor_name=receptor
                    )
                    if score_data:
                        avg_score_value = score_data.get('average_score', np.nan)
                        score_details = score_data.get('receptor_scores', {})
                else:
                    logger.debug(f"Skipping score calculation for {ligand}-{receptor} in {layer_name}: insufficient data.")
                    
                component_counts_ordered = []
                if is_complex and valid_receptors:
                    component_names = valid_receptors.split('_')
                    component_counts_ordered = [component_counts.get(comp, 0) for comp in component_names]
                elif not is_complex:
                    component_counts_ordered = [component_counts.get(receptor, 0)]
                else:
                    component_counts_ordered = [] # Should not happen if valid_receptors is handled
                
                raw_results_layer.append({
                    'ligand': ligand,
                    'receptor': receptor,
                    'is_complex': is_complex,
                    'valid_receptors': valid_receptors,
                    'avg_score': avg_score_value,
                    'num_ligands_found_layer': num_ligands_found,
                    'num_receptors_found_layer': num_centroids,
                    'receptor_scores_detail': score_details,
                    'component_counts_ordered_layer': component_counts_ordered,
                    'ligand_coords_ref': ligand_coords, # Store refs needed for visualization
                    'receptor_clusters_ref': receptor_clusters
                })

            # --- Stage 2: Normalize scores *within this layer* --- 
            logger.info(f"Normalizing scores within {layer_name}...")
            if not raw_results_layer:
                 logger.warning(f"No raw results generated for {layer_name}, skipping normalization and saving.")
                 continue
                 
            results_layer_df = pd.DataFrame(raw_results_layer)
            raw_scores_list = results_layer_df['avg_score'].tolist()
            normalized_scores_layer = normalize_interaction_scores_within_layer(raw_scores_list)
            results_layer_df['score_normalized_within_layer'] = normalized_scores_layer

            # --- Stage 3: Save results and visualize for this layer --- 
            logger.info(f"Saving results and generating visualizations for {layer_name}...")
            layer_spatial_results_dir = results_base_dir / layer_name / "spatial_interaction_analysis"
            summary_dir = layer_spatial_results_dir / "summary"
            density_maps_dir = layer_spatial_results_dir / "density_maps"
            summary_dir.mkdir(parents=True, exist_ok=True)
            density_maps_dir.mkdir(parents=True, exist_ok=True)

            if not results_layer_df.empty:
                cols_to_save = [
                    'ligand', 'receptor', 'is_complex', 'valid_receptors', 
                    'avg_score', 'score_normalized_within_layer', 
                    'num_ligands_found_layer', 'num_receptors_found_layer', 
                    'component_counts_ordered_layer'
                ]
                csv_df = results_layer_df[[col for col in cols_to_save if col in results_layer_df.columns]]
                csv_path = summary_dir / 'lr_spatial_analysis_summary.csv'
                csv_df.to_csv(csv_path, index=False)
                logger.info(f"Spatial analysis summary saved for {layer_name}: {csv_path}")

                # Generate density map visualizations
                for index, row in tqdm(results_layer_df.iterrows(), total=len(results_layer_df), desc=f"Visualizing {layer_name}"):
                    if row['avg_score'] is not None and not np.isnan(row['avg_score']):
                        viz_score_data = {
                            'average_score': row['avg_score'],
                            'normalized_score_within_layer': row['score_normalized_within_layer'],
                            'receptor_scores': row['receptor_scores_detail'] 
                        }
                        visualize_interaction(
                            row['ligand_coords_ref'], row['receptor_clusters_ref'], 
                            viz_score_data, 
                            row['ligand'], row['receptor'], layer_name,
                            density_maps_dir
                        )
                    else:
                         logger.debug(f"Skipping visualization for {row['ligand']}-{row['receptor']} ({layer_name}) due to NaN score.")
            else:
                logger.warning(f"No results generated for {layer_name} to save or visualize.")
                
            # --- Stage 4: Create summary plots *for this layer* ---
            logger.info(f"Generating summary score plots for {layer_name}...")
            if not results_layer_df.empty:
                plot_df = results_layer_df.dropna(subset=['avg_score']).copy()
                plot_df = plot_df.sort_values('score_normalized_within_layer', ascending=False)
                
                if not plot_df.empty:
                    plt.figure(figsize=(12, max(6, len(plot_df) * 0.4))) # Adjusted size
                    sns.barplot(data=plot_df, x='score_normalized_within_layer', y='ligand', hue='receptor', dodge=True)
                    plt.title(f'Top {len(plot_df)} Pairs by Normalized Interaction Score ({layer_name})')
                    plt.xlabel('Normalized Score (Within Layer, 0-1)')
                    plt.ylabel('Ligand')
                    plt.legend(title='Receptor', bbox_to_anchor=(1.02, 1), loc='upper left')
                    plt.tight_layout(rect=[0, 0, 0.85, 1])
                    plt.savefig(summary_dir / 'lr_scores_summary_normalized.png', dpi=300)
                    plt.close()
                    logger.info(f"Layer summary plot saved: {summary_dir / 'lr_scores_summary_normalized.png'}")
                else: 
                    logger.warning(f"No valid data to plot summary scores for {layer_name}.")
            else:
                logger.warning(f"No results generated, skipping summary plots for {layer_name}.")
                
            # Memory management
            del spatial_data_layer, lr_pairs_layer, raw_results_layer, results_layer_df, spatial_coords_cache_layer
            gc.collect()

        logger.info("\n=== Layer-specific spatial interaction analysis complete ===")

    except FileNotFoundError as e:
        logger.error(f"Input file not found: {e}")
    except ValueError as e:
        logger.error(f"Data validation error: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred during layer-specific spatial interaction analysis: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze spatial LR interactions layer by layer.')
    parser.add_argument('--pathway-results-base-dir', type=str, 
                        default='results',
                        help="Base directory containing layer-specific pathway analysis results (e.g., results/Layer1/pathway_analysis/).")
    parser.add_argument('--spatial-data-file', type=str, 
                        default='data/spatial/raw/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop_bin100_clustered.csv',
                        help="Path relative to project root for the main clustered spatial data CSV file (must contain 'layer', 'x', 'y' columns).")
    parser.add_argument('--results-base-dir', type=str, 
                        default='results',
                        help="Base directory where layer-specific spatial analysis results will be saved.")
    parser.add_argument('--num-pairs', type=int, default=15,
                        help="Number of top ligand-receptor pairs from pathway analysis to analyze spatially per layer.")
    # parser.add_argument('--k-nearest', type=int, default=50, help='k for k-NN density (Not used in current KDE method).')
    
    args = parser.parse_args()
    main(args)
