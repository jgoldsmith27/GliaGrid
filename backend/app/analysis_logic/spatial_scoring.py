#!/usr/bin/env python3
"""
Spatial analysis functions including KDE-based scoring and visualization.
Migrated and adapted from layer_analysis/analyze_ligand_receptor_spatial_by_layer.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.stats import gaussian_kde
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
from pathlib import Path
import logging
import networkx as nx
from scipy import stats
import gc
from typing import Dict, Any, List

# Configure logger (consider integrating with FastAPI logging if applicable)
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


# --- Helper Functions (potentially move or refactor if shared) ---

def get_gene_spatial_data_for_layer(spatial_data_layer: pd.DataFrame, gene_name: str):
    """Extract spatial coordinates for a specific gene within a given layer DataFrame."""
    if spatial_data_layer is None or not isinstance(spatial_data_layer, pd.DataFrame) or 'gene' not in spatial_data_layer.columns:
        logger.warning(f"Invalid spatial data provided for layer when searching for {gene_name}.")
        return None
    
    gene_data = spatial_data_layer[spatial_data_layer['gene'] == gene_name]
    if gene_data.empty:
        # logger.debug(f"Gene {gene_name} not found in the provided spatial data layer.")
        return None
    
    # Assuming columns 'x' and 'y' exist for coordinates
    if 'x' not in gene_data.columns or 'y' not in gene_data.columns:
        logger.warning(f"Missing 'x' or 'y' coordinate columns in spatial data for gene {gene_name}.")
        return None
        
    return gene_data[['x', 'y']].to_numpy()

# --- KDE Score Calculation ---

def calculate_kde_scores(ligand_coords: np.ndarray, receptor_coords: np.ndarray, ligand_name: str = None, receptor_name: str = None) -> Dict[str, Any]:
    """
    Calculate spatial interaction scores based on Kernel Density Estimation (KDE).

    Computes three types of scores evaluated at receptor locations:
    1.  **Co-occurrence Score:** Product of Ligand KDE and Receptor KDE.
        Indicates co-localization density.
    2.  **Ligand Density Score:** Ligand KDE evaluated at receptor locations.
        Indicates ligand density near receptors.
    3.  **Receptor Density Score:** Receptor KDE evaluated at receptor locations.
        Indicates receptor density near other receptors (reflects clustering).

    Args:
        ligand_coords: Numpy array of shape (n_ligands, 2) for ligand coordinates.
        receptor_coords: Numpy array of shape (n_receptors, 2) for receptor coordinates.
        ligand_name: Optional name of the ligand for logging.
        receptor_name: Optional name of the receptor for logging.

    Returns:
        A dictionary containing average scores and detailed scores for each receptor,
        or a dictionary with NaN/empty values if calculation fails.
        Example:
        {
            "avg_co_occurrence_score": float,
            "avg_ligand_density_score": float,
            "avg_receptor_density_score": float,
            "receptor_scores": {
                "receptor_0": {
                    "co_occurrence_score": float,
                    "ligand_density_score": float,
                    "receptor_density_score": float
                }, ...
            }
        }
    """
    pair_id = f"{ligand_name}-{receptor_name}" if ligand_name and receptor_name else "pair"

    # Basic validation
    if ligand_coords is None or not isinstance(ligand_coords, np.ndarray) or ligand_coords.ndim != 2 or ligand_coords.shape[1] != 2:
        logger.error(f"Invalid ligand coordinates provided for {pair_id}.")
        return {"avg_co_occurrence_score": np.nan, "avg_ligand_density_score": np.nan, "avg_receptor_density_score": np.nan, "receptor_scores": {}}
    if receptor_coords is None or not isinstance(receptor_coords, np.ndarray) or receptor_coords.ndim != 2 or receptor_coords.shape[1] != 2:
        logger.error(f"Invalid receptor coordinates provided for {pair_id}.")
        return {"avg_co_occurrence_score": np.nan, "avg_ligand_density_score": np.nan, "avg_receptor_density_score": np.nan, "receptor_scores": {}}

    n_ligands = len(ligand_coords)
    n_receptors = len(receptor_coords)

    # Fallback for insufficient points for KDE
    fallback_return = {
        "avg_co_occurrence_score": 0.0,
        "avg_ligand_density_score": 0.0,
        "avg_receptor_density_score": 0.0,
        "receptor_scores": {}
    }
    if n_ligands < 2:
        logger.warning(f"Cannot calculate KDE scores for {pair_id}: requires >= 2 ligand points, found {n_ligands}.")
        return fallback_return
    if n_receptors < 2:
        logger.warning(f"Cannot calculate KDE scores for {pair_id}: requires >= 2 receptor points, found {n_receptors}.")
        return fallback_return
    if n_receptors == 0: # Should be caught by n_receptors < 2, but defensive check
        logger.warning(f"Cannot calculate KDE scores for {pair_id}: no receptor points found.")
        return fallback_return

    try:
        # Calculate KDEs
        ligand_kde = gaussian_kde(ligand_coords.T)
        receptor_kde = gaussian_kde(receptor_coords.T)

        # Evaluate densities at receptor locations
        ligand_density_at_receptors = ligand_kde(receptor_coords.T)
        receptor_density_at_receptors = receptor_kde(receptor_coords.T)

        # Ensure non-negative densities (KDE can sometimes yield small negatives)
        ligand_density_at_receptors[ligand_density_at_receptors < 0] = 0
        receptor_density_at_receptors[receptor_density_at_receptors < 0] = 0

        # Calculate scores per receptor
        co_occurrence_scores = ligand_density_at_receptors * receptor_density_at_receptors

        # Calculate average scores, handling potential NaNs if all scores were zero etc.
        avg_co_occurrence = np.nanmean(co_occurrence_scores) if len(co_occurrence_scores) > 0 else 0.0
        avg_ligand_density = np.nanmean(ligand_density_at_receptors) if len(ligand_density_at_receptors) > 0 else 0.0
        avg_receptor_density = np.nanmean(receptor_density_at_receptors) if len(receptor_density_at_receptors) > 0 else 0.0

        # Structure detailed scores per receptor
        receptor_scores_detail = {
            f"receptor_{i}": {
                "co_occurrence_score": float(co_occurrence_scores[i]), # Ensure standard float type
                "ligand_density_score": float(ligand_density_at_receptors[i]),
                "receptor_density_score": float(receptor_density_at_receptors[i])
            }
            for i in range(n_receptors)
        }

        logger.info(f"Calculated KDE scores for {pair_id}. Avg Co-occurrence: {avg_co_occurrence:.4g}, Avg Ligand Density: {avg_ligand_density:.4g}, Avg Receptor Density: {avg_receptor_density:.4g}")

        return {
            "avg_co_occurrence_score": float(avg_co_occurrence),
            "avg_ligand_density_score": float(avg_ligand_density),
            "avg_receptor_density_score": float(avg_receptor_density),
            "receptor_scores": receptor_scores_detail
        }

    except Exception as e:
        logger.error(f"Failed to calculate KDE scores for {pair_id}: {e}")
        # Use traceback logger here if more detail is needed: logger.exception(...)
        return {
            "avg_co_occurrence_score": np.nan,
            "avg_ligand_density_score": np.nan,
            "avg_receptor_density_score": np.nan,
            "receptor_scores": {}
        }

# --- Visualization ---

def visualize_spatial_density(
    ligand_coords: np.ndarray, 
    receptor_coords: np.ndarray, 
    score_data: Dict[str, Any], 
    ligand_name: str, 
    receptor_name: str, 
    scope_name: str, 
    output_dir: Path, 
    score_type: str = 'co_occurrence',
    custom_plot_description: str = None, # Added for UI customization
    custom_legend_label: str = None     # Added for UI customization
) -> Path:
    """
    Visualize the spatial density for a ligand-receptor pair based on the specified score type.

    Generates and saves a plot showing the selected density (co-occurrence, ligand, or receptor)
    as a heatmap, with ligand/receptor points overlaid.

    Args:
        ligand_coords: Numpy array of ligand coordinates (n_ligands, 2).
        receptor_coords: Numpy array of receptor coordinates (n_receptors, 2).
        score_data: Dictionary containing calculated KDE scores from calculate_kde_scores.
        ligand_name: Name of the ligand.
        receptor_name: Name of the receptor.
        scope_name: Name of the scope (e.g., layer name) for context in title/logging.
        output_dir: Path to the directory where the plot image will be saved.
        score_type: The type of score to visualize ('co_occurrence', 'ligand_density', 
                    'receptor_density'). Defaults to 'co_occurrence'.
        custom_plot_description: Optional custom description for the plot title.
        custom_legend_label: Optional custom label for the plot legend.

    Returns:
        The Path object pointing to the saved image file, or None if visualization failed.
    """
    pair_id = f"{ligand_name}-{receptor_name}"
    logger.info(f"Visualizing {score_type} for {pair_id} in {scope_name}...")

    # --- Input Validation ---
    if ligand_coords is None or not isinstance(ligand_coords, np.ndarray) or ligand_coords.ndim != 2 or ligand_coords.shape[1] != 2:
        logger.warning(f"Invalid/missing ligand coords for {pair_id} in {scope_name}. Skipping visualization.")
        return None
    if receptor_coords is None or not isinstance(receptor_coords, np.ndarray) or receptor_coords.ndim != 2 or receptor_coords.shape[1] != 2:
         logger.warning(f"Invalid/missing receptor coords for {pair_id} in {scope_name}. Skipping visualization.")
         return None   

    n_ligands = len(ligand_coords)
    n_receptors = len(receptor_coords)
    
    # Check if score data is valid for the requested type
    avg_score_key = f"avg_{score_type}_score"
    if avg_score_key not in score_data:
        logger.warning(f"Requested score_type '{score_type}' not found in score_data for {pair_id}. Skipping visualization.")
        return None
    avg_score_value = score_data[avg_score_key]
    if avg_score_value is None or np.isnan(avg_score_value):
        logger.warning(f"Cannot visualize {score_type} density for {pair_id} in {scope_name}: missing or invalid average score.")
        return None
        
    # Check if enough points exist for the required KDE calculations
    min_points_for_kde = 2
    req_ligand = score_type in ['co_occurrence', 'ligand_density']
    req_receptor = score_type in ['co_occurrence', 'receptor_density']
    
    if req_ligand and n_ligands < min_points_for_kde:
        logger.warning(f"Cannot visualize {score_type} density for {pair_id} in {scope_name}: requires >= {min_points_for_kde} ligand points, found {n_ligands}.")
        return None
    if req_receptor and n_receptors < min_points_for_kde:
         logger.warning(f"Cannot visualize {score_type} density for {pair_id} in {scope_name}: requires >= {min_points_for_kde} receptor points, found {n_receptors}.")
         return None
    if n_ligands == 0 or n_receptors == 0: # Basic check
         logger.warning(f"Cannot visualize density for {pair_id} in {scope_name}: zero ligand or receptor points.")
         return None

    # --- Plotting Setup ---
    plt.figure(figsize=(12, 10))
    ax_density = plt.gca()
    output_dir.mkdir(parents=True, exist_ok=True) # Ensure output dir exists

    # Define grid boundaries based on combined data
    all_points = np.vstack([ligand_coords, receptor_coords])
    x_min, x_max = np.min(all_points[:, 0]), np.max(all_points[:, 0])
    y_min, y_max = np.min(all_points[:, 1]), np.max(all_points[:, 1])
    x_range = x_max - x_min
    y_range = y_max - y_min
    padding_factor = 0.05
    x_padding = x_range * padding_factor if x_range > 0 else 10
    y_padding = y_range * padding_factor if y_range > 0 else 10
    x_min -= x_padding; x_max += x_padding
    y_min -= y_padding; y_max += y_padding
    
    # Create grid for KDE evaluation
    grid_res_complex = 150j # Use complex number for mgrid resolution
    x_grid, y_grid = np.mgrid[x_min:x_max:grid_res_complex, y_min:y_max:grid_res_complex]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])

    # --- Calculate Densities on Grid ---
    ligand_density_grid = np.zeros(x_grid.shape)
    receptor_density_grid = np.zeros(x_grid.shape)
    can_calc_ligand = False
    can_calc_receptor = False

    if n_ligands >= min_points_for_kde:
        try:
            ligand_kernel = gaussian_kde(ligand_coords.T)
            ligand_density_grid = np.reshape(ligand_kernel(positions), x_grid.shape)
            ligand_density_grid[ligand_density_grid < 0] = 0 # Ensure non-negative
            can_calc_ligand = True
        except Exception as e:
            logger.warning(f"Could not compute ligand density grid for {ligand_name} ({scope_name}): {e}")

    if n_receptors >= min_points_for_kde:
        try:
            receptor_kernel = gaussian_kde(receptor_coords.T)
            receptor_density_grid = np.reshape(receptor_kernel(positions), x_grid.shape)
            receptor_density_grid[receptor_density_grid < 0] = 0 # Ensure non-negative
            can_calc_receptor = True
        except Exception as e:
            logger.warning(f"Could not compute receptor density grid for {receptor_name} ({scope_name}): {e}")

    # --- Select Density and Plot Parameters based on score_type ---
    density_to_plot = None
    cmap = 'viridis' # Default colormap
    plot_title_main = ""
    plot_title_suffix = ""
    legend_label = ""
    density_description = ""
    
    # Determine default descriptions based on score_type
    if score_type == 'co_occurrence':
        default_density_description = "Co-occurrence"
        default_cmap = 'Greens'
        default_legend_label = 'Co-occurrence Density'
    elif score_type == 'ligand_density':
        default_density_description = "Ligand"
        default_cmap = 'Blues'
        default_legend_label = 'Ligand Density'
    elif score_type == 'receptor_density':
        default_density_description = "Receptor"
        default_cmap = 'Reds'
        default_legend_label = 'Receptor Density'
    else:
        logger.error(f"Invalid score_type '{score_type}' provided for visualization. Cannot proceed.")
        plt.close() # Close the figure
        return None
        
    # Use custom labels if provided, otherwise use defaults
    density_description = custom_plot_description if custom_plot_description else default_density_description
    cmap = default_cmap # Colormap is tied to the calculation type, not custom label
    legend_label = custom_legend_label if custom_legend_label else default_legend_label

    # Select the density grid based on the actual score_type calculation needed
    if score_type == 'co_occurrence' and can_calc_ligand and can_calc_receptor:
        density_to_plot = ligand_density_grid * receptor_density_grid
    elif score_type == 'ligand_density' and can_calc_ligand:
        density_to_plot = ligand_density_grid
    elif score_type == 'receptor_density' and can_calc_receptor:
        density_to_plot = receptor_density_grid
        
    # Log warnings if the required density for the *calculation* couldn't be made
    if density_to_plot is None:
         logger.warning(f"Cannot calculate required density for score_type '{score_type}' ({pair_id}, {scope_name}). Skipping heatmap.")

    # Construct plot titles using the potentially customized density_description
    plot_title_main = f"{density_description} Density: {ligand_name} ({n_ligands}) - {receptor_name} ({n_receptors})"
    plot_title_suffix = f"in {scope_name}"

    # --- Plot Density Heatmap & Contours ---
    plotted_density_heatmap = False
    if density_to_plot is not None:
        try:
            max_density = np.max(density_to_plot)
            density_threshold = 1e-9 # Threshold to avoid plotting near-zero density
            
            if max_density > density_threshold: 
                # Normalize density for visualization (0 to 1)
                density_normalized = density_to_plot / max_density 
                
                # Plot heatmap
                im_density = ax_density.imshow(np.rot90(density_normalized), cmap=cmap, 
                                              extent=[x_min, x_max, y_min, y_max], 
                                              alpha=0.75, aspect='auto', 
                                              vmin=0.01, vmax=1.0, # Start color range slightly above zero
                                              zorder=1)
                plotted_density_heatmap = True
                
                # Plot contours (optional, based on normalized density)
                contour_threshold = 0.05 # Minimum normalized density for contour calc
                density_for_contours = density_normalized[density_normalized > contour_threshold]
                if len(density_for_contours) > 0:
                    # Define contour levels (e.g., percentiles of the density above threshold)
                    contour_levels_perc = [60, 80, 95] 
                    density_thresholds = np.percentile(density_for_contours, contour_levels_perc) 
                    contour_levels = sorted(list(set(d for d in density_thresholds if d > contour_threshold))) # Ensure levels are unique and above threshold
                    
                    if len(contour_levels) > 0:
                        contour = ax_density.contour(x_grid, y_grid, density_normalized, 
                                                levels=contour_levels, 
                                                colors='purple', # Consistent contour color
                                                alpha=0.8, linewidths=0.8, zorder=2)
                        # Add labels to contours
                        ax_density.clabel(contour, inline=True, fontsize=6, fmt='%.2f') 
                else:
                     logger.debug(f"Normalized {score_type} density too low for contours ({pair_id}, {scope_name})")
            else:
                 logger.info(f"{density_description} density below plotting threshold ({density_threshold:.1g}) for {pair_id} in {scope_name}.")
        except Exception as e:
            logger.warning(f"Could not plot {score_type} density heatmap/contours for {pair_id} in {scope_name}: {e}")
            plotted_density_heatmap = False # Ensure flag is false if plotting fails
            
    # --- Plot Raw Points ---
    point_size = 6
    ligand_alpha = 0.7
    receptor_alpha = 0.8
    
    # Plot Ligands
    if n_ligands > 0:
        ax_density.scatter(ligand_coords[:, 0], ligand_coords[:, 1], 
                           s=point_size, color='mediumblue', alpha=ligand_alpha, 
                           label=f'{ligand_name} (Ligand)', marker='.', zorder=3)
    
    # Plot Receptors (potentially subsample if very dense)
    max_receptors_to_plot = 5000 
    if n_receptors > 0:
        if n_receptors > max_receptors_to_plot:
            indices_to_plot = np.random.choice(n_receptors, size=max_receptors_to_plot, replace=False)
            receptor_coords_to_plot = receptor_coords[indices_to_plot, :]
            receptor_label = f'{receptor_name} (Receptor, subsampled)'
        else:
            receptor_coords_to_plot = receptor_coords
            receptor_label = f'{receptor_name} (Receptor)'
            
        ax_density.scatter(receptor_coords_to_plot[:, 0], receptor_coords_to_plot[:, 1], 
                           s=point_size + 2, color='firebrick', alpha=receptor_alpha, 
                           label=receptor_label, marker='.', zorder=4)

    # --- Final Plot Styling ---
    
    # Title
    score_display = f'Avg. Score ({default_density_description}): {avg_score_value:.4g}'
    plt.title(f'{plot_title_main} {plot_title_suffix}\n{score_display}')
    
    # Labels and Limits
    ax_density.set_xlabel('X Coordinate')
    ax_density.set_ylabel('Y Coordinate')
    ax_density.set_xlim(x_min, x_max)
    ax_density.set_ylim(y_min, y_max)
    
    # Aspect Ratio and Grid
    ax_density.set_aspect('equal', adjustable='box')
    ax_density.grid(True, linestyle=':', alpha=0.3, zorder=0)
    
    # Legend
    handles, labels = ax_density.get_legend_handles_labels()
    if plotted_density_heatmap:
        # Add patch for heatmap color
        legend_color = plt.get_cmap(cmap)(0.6) # Sample color from the heatmap's cmap
        # Use the potentially customized legend_label here
        handles.append(plt.Rectangle((0,0), 1, 1, fc=legend_color, alpha=0.7))
        labels.append(legend_label) 
        # Add line for contours if they were potentially plotted
        if 'contour' in locals() and contour is not None:
            handles.append(plt.Line2D([0], [0], color='purple', lw=1))
            # Contour label is always normalized density
            labels.append('Density Contours (Norm.)') 
            
    legend_fontsize = 'small' if len(handles) < 8 else 'x-small'
    ax_density.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.02, 1), 
                      fontsize=legend_fontsize, frameon=False)

    # Scale Bar (Example: Assuming 1 unit = 0.5 micrometers)
    # Make scale bar relative to the plot range
    scale_bar_fraction = 0.1 # Fraction of x-range for scale bar
    scale_bar_length_units = (x_max - x_min) * scale_bar_fraction
    pixels_per_unit = 2 # Example: if 1 unit = 0.5 um, then 2 units = 1 um
    scale_bar_microns = scale_bar_length_units / pixels_per_unit 
    scale_bar_label = f'{scale_bar_microns:.0f} μm' if scale_bar_microns >= 1 else f'{scale_bar_microns*1000:.0f} nm'
    
    scale_y_pos_fraction = 0.05 # Position from bottom
    scale_x_pos_fraction = 0.05 # Position from left
    scale_y = y_min + (y_max - y_min) * scale_y_pos_fraction
    scale_x = x_min + (x_max - x_min) * scale_x_pos_fraction
    
    ax_density.plot([scale_x, scale_x + scale_bar_length_units], [scale_y, scale_y], 
                    'k-', linewidth=3)
    ax_density.text(scale_x + scale_bar_length_units / 2, scale_y - (y_max - y_min) * 0.03, 
                    scale_bar_label, ha='center', va='top', fontsize=9)
    
    # --- Save Figure ---
    # Sanitize names for filename
    safe_ligand = "".join(c if c.isalnum() else "_" for c in ligand_name)
    safe_receptor = "".join(c if c.isalnum() else "_" for c in receptor_name)
    safe_scope = "".join(c if c.isalnum() else "_" for c in scope_name)
    
    # Use the potentially customized description in the filename for clarity
    safe_description = "".join(c if c.isalnum() else "_" for c in density_description.lower().replace(' ','_'))
    output_filename = f"{safe_ligand}_{safe_receptor}_{safe_scope}_{safe_description}_density.png"
    output_path = output_dir / output_filename
    
    try:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        # Use the potentially customized description in the log message
        logger.info(f"{density_description} density map saved to {output_path}") 
    except Exception as e:
        logger.error(f"Failed to save density map to {output_path}: {e}")
        output_path = None # Indicate failure
    finally:
        plt.close() # Ensure figure is closed to free memory
        
    return output_path


# --- Complex Receptor Handling ---
# (Skipped for now - assuming simple receptors)
# Placeholder function for structure
def get_receptor_coordinates(receptor_name: str, spatial_data_layer: pd.DataFrame) -> np.ndarray:
    """Get coordinates for a simple receptor (complex handling omitted for now)."""
    logger.debug(f"Getting coordinates for simple receptor: {receptor_name}")
    # Uses the existing helper, assuming gene column is named 'gene' in this context
    # TODO: Ensure column name consistency ('gene' vs 'geneID')
    # For now, let's assume the helper function uses the correct column name internally
    # or the DataFrame passed has been adapted.
    return get_gene_spatial_data_for_layer(spatial_data_layer, receptor_name)


# --- Normalization ---

# (Removing normalize_scores function as requested)

# (End of file or add API integration functions below later) 