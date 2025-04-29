#!/usr/bin/env python3
"""API Endpoint for Spatial Analysis Visualization."""

from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from pathlib import Path
import pandas as pd # Placeholder - data loading needs actual implementation
import logging
import numpy as np

# Import analysis functions (adjust path if needed based on final structure)
from ..analysis_logic import (
    calculate_kde_scores, 
    visualize_spatial_density, 
    get_receptor_coordinates, 
    get_gene_spatial_data_for_layer
)
from ..core.config import settings # Assuming settings define data paths

logger = logging.getLogger(__name__)
router = APIRouter()

# --- Models --- 

class SpatialVisualizationRequest(BaseModel):
    ligand_name: str = Field(..., description="Name of the ligand gene.")
    receptor_name: str = Field(..., description="Name of the receptor gene/complex.")
    scope_name: str = Field(..., description="Scope of analysis (e.g., 'Layer 1', 'whole_tissue').")
    visualization_type: str = Field("Co-occurrence Density", description="Desired visualization type (e.g., 'Co-occurrence Density', 'Ligand Density', 'Receptor Density', 'Ligand Density at Receptors', 'Receptor Density at Ligands').")

class SpatialVisualizationResponse(BaseModel):
    image_path: Optional[str] = Field(None, description="Relative path to the generated visualization image within static files.")
    average_scores: Dict[str, Optional[float]] = Field({}, description="Dictionary of calculated average scores (co-occurrence, ligand, receptor).")
    error: Optional[str] = Field(None, description="Error message if visualization failed.")

# --- Helper Functions (Placeholder - Replace with actual data loading) ---

def get_spatial_data() -> pd.DataFrame:
    # !!! Placeholder: Implement actual loading from configured path !!!
    # Example: 
    # try:
    #     data_path = Path(settings.SPATIAL_DATA_PATH) 
    #     df = pd.read_csv(data_path)
    #     # Validate columns: 'gene', 'x', 'y', 'layer' 
    #     required_cols = ['gene', 'x', 'y', 'layer'] # Adjust based on actual file
    #     if not all(col in df.columns for col in required_cols):
    #         raise ValueError(f"Spatial data missing required columns: {required_cols}")
    #     return df
    # except Exception as e:
    #     logger.error(f"Failed to load spatial data: {e}")
    #     raise HTTPException(status_code=500, detail="Failed to load spatial data.")
    logger.warning("Using placeholder spatial data.")
    # Create minimal placeholder data
    data = {
        'gene': ['L1','L1','L1', 'R1', 'R1', 'R1', 'L1', 'R1', 'L2', 'R2'], 
        'x': np.random.rand(10) * 100, 
        'y': np.random.rand(10) * 100, 
        'layer': ['Layer1']*5 + ['Layer2']*5
    }
    return pd.DataFrame(data)

# --- Endpoint --- 

@router.post("/spatial/visualize/", response_model=SpatialVisualizationResponse)
def run_spatial_visualization(
    request: SpatialVisualizationRequest,
    spatial_data: pd.DataFrame = Depends(get_spatial_data)
):
    """
    Generates a spatial density visualization for a ligand-receptor pair within a specific scope.
    Allows selection from five visualization types influencing the heatmap and/or interpretation.
    """
    logger.info(f"Received spatial visualization request: {request.dict()}")

    # 1. Filter data for the requested scope
    if request.scope_name == 'whole_tissue':
        scope_data = spatial_data
    else:
        scope_data = spatial_data[spatial_data['layer'] == request.scope_name].copy()

    if scope_data.empty:
        logger.warning(f"No spatial data found for scope: {request.scope_name}")
        return SpatialVisualizationResponse(error=f"No spatial data found for scope: {request.scope_name}")

    # 2. Get Coordinates (Using placeholder simple receptor logic)
    ligand_coords = get_gene_spatial_data_for_layer(scope_data, request.ligand_name)
    # Using placeholder simple receptor function - update when complex logic is added
    receptor_coords = get_receptor_coordinates(request.receptor_name, scope_data) 

    if ligand_coords is None or receptor_coords is None or len(ligand_coords) == 0 or len(receptor_coords) == 0:
        err_msg = f"Could not find sufficient coordinate data for pair {request.ligand_name}-{request.receptor_name} in scope {request.scope_name}."
        logger.warning(err_msg)
        return SpatialVisualizationResponse(error=err_msg)
        
    logger.info(f"Found {len(ligand_coords)} ligand and {len(receptor_coords)} receptor points for {request.ligand_name}-{request.receptor_name} in {request.scope_name}.")

    # 3. Calculate KDE Scores 
    score_data = calculate_kde_scores(
        ligand_coords, 
        receptor_coords, 
        ligand_name=request.ligand_name, 
        receptor_name=request.receptor_name
    )

    avg_scores_response = {
        "co_occurrence": score_data.get("avg_co_occurrence_score"),
        "ligand_density": score_data.get("avg_ligand_density_score"),
        "receptor_density": score_data.get("avg_receptor_density_score"),
    }

    # 4. Determine Visualization Parameters based on Request Type
    score_type_for_calc = None
    custom_plot_desc = None
    custom_legend_label = None
    viz_type = request.visualization_type

    if viz_type == "Co-occurrence Density":
        score_type_for_calc = "co_occurrence"
        # Uses default labels
    elif viz_type == "Ligand Density": # User interpretation: Density of ligands around ligands
        score_type_for_calc = "ligand_density"
        # Uses default labels
    elif viz_type == "Receptor Density": # User interpretation: Density of receptors around receptors
        score_type_for_calc = "receptor_density"
        # Uses default labels
    elif viz_type == "Ligand Density at Receptors": # User interpretation
        score_type_for_calc = "ligand_density" # Calculation uses ligand density
        custom_plot_desc = "Ligand Density" # Title reflects the calculation
        custom_legend_label = "Ligand Density (near Receptors)" # Legend clarifies context
    elif viz_type == "Receptor Density at Ligands": # User interpretation
        score_type_for_calc = "receptor_density" # Calculation uses receptor density
        custom_plot_desc = "Receptor Density" # Title reflects the calculation
        custom_legend_label = "Receptor Density (near Ligands)" # Legend clarifies context
    else:
        logger.error(f"Invalid visualization_type received: {viz_type}")
        raise HTTPException(status_code=400, detail=f"Invalid visualization type: {viz_type}")

    # Define output directory (relative to a static path accessible by frontend)
    # !!! Placeholder: Define properly using settings !!!
    static_base = Path("backend/app/static") # Example
    output_sub_dir = f"spatial_viz/{request.scope_name}" # Organize by scope
    output_dir = static_base / output_sub_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # 5. Call Visualization Function
    image_result_path = visualize_spatial_density(
        ligand_coords=ligand_coords,
        receptor_coords=receptor_coords,
        score_data=score_data,
        ligand_name=request.ligand_name,
        receptor_name=request.receptor_name,
        scope_name=request.scope_name,
        output_dir=output_dir,
        score_type=score_type_for_calc, # Use the calculation type
        custom_plot_description=custom_plot_desc, # Pass custom labels
        custom_legend_label=custom_legend_label
    )

    # 6. Prepare Response
    if image_result_path:
        # Return path relative to static base for frontend access
        relative_image_path = str(image_result_path.relative_to(static_base))
        logger.info(f"Visualization successful. Image path: {relative_image_path}")
        return SpatialVisualizationResponse(
            image_path=relative_image_path,
            average_scores=avg_scores_response
        )
    else:
        logger.error(f"Visualization failed for {request.ligand_name}-{request.receptor_name} in {request.scope_name}.")
        return SpatialVisualizationResponse(
            error="Visualization generation failed.",
            average_scores=avg_scores_response # Still return scores if calculated
        ) 