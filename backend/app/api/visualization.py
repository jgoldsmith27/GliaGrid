from fastapi import APIRouter, HTTPException, Query, Depends
from typing import Dict, List
import pandas as pd
import numpy as np
import logging

# Import JobService and its dependency function
from app.services.job_service import JobService, get_job_service 
# Import AnalysisService only for its _load_and_standardize method (consider moving this)
from app.services.analysis_service import AnalysisService, STANDARDIZED_SPATIAL_COLS, get_analysis_service
from app.services.file_service import FileService
from app.models.analysis_models import AnalysisMapping # Import for type hinting/casting

logger = logging.getLogger(__name__)

router = APIRouter()

# TODO: Consider moving _load_and_standardize to a shared utility module 
#       or FileService to avoid depending on AnalysisService here.

@router.get("/visualization/{job_id}")
async def get_visualization_data(
    job_id: str,
    ligand: str = Query(..., description="Ligand gene name"),
    receptor: str = Query(..., description="Receptor gene name"),
    layer: str = Query(..., description="Layer name (required)"),
    job_service: JobService = Depends(get_job_service),
    # Inject AnalysisService mainly for _load_and_standardize
    analysis_service: AnalysisService = Depends(get_analysis_service) 
):
    logger.info(f"Received visualization request for job {job_id}, L: {ligand}, R: {receptor}, Layer: {layer}")
    try:
        # 1. Look up job context using JobService
        job_context = job_service.get_job_context(job_id)
        if not job_context:
            logger.warning(f"Job context not found for job ID: {job_id} in visualization endpoint.")
            raise HTTPException(status_code=404, detail=f"Job context for job ID {job_id} not found. Has the analysis run?")

        # Safely access context data
        spatial_file_id = job_context.get("spatialFileId")
        spatial_mapping_dict = job_context.get("spatialMapping")
        if not spatial_file_id or not spatial_mapping_dict:
             logger.error(f"Spatial file ID or mapping missing in job context for {job_id}")
             raise HTTPException(status_code=400, detail="Job context is incomplete. Missing spatial file ID or mapping.")
        
        # Convert mapping dict back to model if needed by loader function
        try:
            spatial_mapping = AnalysisMapping(**spatial_mapping_dict)
        except Exception as e:
             logger.error(f"Error converting spatial mapping from context for {job_id}: {e}")
             raise HTTPException(status_code=500, detail="Error processing job context mapping.")

        # 2. Load and standardize the spatial file
        logger.info(f"Loading spatial file {spatial_file_id} for job {job_id}")
        spatial_path = FileService.get_file_path(spatial_file_id)
        
        # Use the standardization method (currently from AnalysisService)
        df = analysis_service._load_and_standardize(
            spatial_path,
            spatial_mapping,
            STANDARDIZED_SPATIAL_COLS
        )
        logger.info(f"Spatial data loaded and standardized for job {job_id}. Shape: {df.shape}")

        # 3. Filter for the selected layer
        if "layer" not in df.columns:
            logger.error(f"Standardized spatial data for job {job_id} is missing 'layer' column.")
            raise HTTPException(status_code=400, detail="Layer column missing after standardization.")
        
        layer_df = df[df["layer"] == layer]
        if layer_df.empty:
            logger.warning(f"No data found for layer '{layer}' in job {job_id}.")
            # Return empty data instead of 404/400, frontend can handle this
            return {"ligand": [], "receptor": []}
            # raise HTTPException(status_code=404, detail=f"No data found for layer '{layer}'.")
        logger.info(f"Filtered data for layer '{layer}' in job {job_id}. Shape: {layer_df.shape}")

        # 4. For ligand and receptor, extract x/y for each
        if "gene" not in layer_df.columns or "x" not in layer_df.columns or "y" not in layer_df.columns:
            logger.error(f"Layer data for job {job_id} is missing required columns (gene, x, y). Columns: {layer_df.columns}")
            raise HTTPException(status_code=400, detail="Required columns (gene, x, y) missing after standardization.")
            
        ligand_points = layer_df[layer_df["gene"] == ligand][["x", "y"]]
        receptor_points = layer_df[layer_df["gene"] == receptor][["x", "y"]]
        
        # Handle cases where ligand/receptor might not be in the layer gracefully
        ligand_data = ligand_points.to_dict("records") if not ligand_points.empty else []
        receptor_data = receptor_points.to_dict("records") if not receptor_points.empty else []

        if not ligand_data:
             logger.warning(f"No data found for ligand '{ligand}' in layer '{layer}' for job {job_id}.")
        if not receptor_data:
             logger.warning(f"No data found for receptor '{receptor}' in layer '{layer}' for job {job_id}.")

        # 5. Return in frontend format
        logger.info(f"Returning visualization data for job {job_id}. Ligand points: {len(ligand_data)}, Receptor points: {len(receptor_data)}")
        return {
            "ligand": ligand_data,
            "receptor": receptor_data
        }
        
    except HTTPException as he:
        # Log HTTPExceptions raised explicitly
        logger.error(f"HTTP Exception in visualization for job {job_id}: {he.status_code} - {he.detail}")
        raise # Re-raise the exception
    except Exception as e:
        # Log unexpected errors
        logger.exception(f"Unexpected error generating visualization for job {job_id}, L:{ligand}, R:{receptor}, Layer:{layer}: {e}")
        raise HTTPException(status_code=500, detail=f"An internal error occurred while generating the visualization: {str(e)}") 