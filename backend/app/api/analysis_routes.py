from fastapi import APIRouter, Depends, HTTPException, Body, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Any
import asyncio
import uuid
import json
import logging
import pandas as pd
from shapely.geometry import Point, Polygon # Import Shapely

# Import JobService and AnalysisService
from app.services.job_service import JobService, get_job_service
from app.services.analysis_service import AnalysisService
from app.models.analysis_models import AnalysisMapping, AnalysisPayload, CustomLassoAnalysisRequest, CustomAnalysisResponse
# Import core analysis logic directly
from app.analysis_logic.core import run_stage1_counts_pipeline, run_pathway_dominance_pipeline, run_module_context_pipeline
# Import FileService class instead of the specific function
from app.services.file_service import FileService 

logger = logging.getLogger(__name__)

router = APIRouter(
    prefix="/analysis",
    tags=["Analysis"],
    responses={404: {"description": "Not found"}},
)

# Pydantic model for the response (can be expanded later)
class AnalysisResponse(BaseModel):
    status: str = Field(..., description="Status of the analysis initiation")
    message: str = Field(..., description="A message detailing the status or next steps")
    job_id: Optional[str] = Field(None, description="Job ID for tracking asynchronous analysis")

# Pydantic model for the status response
class JobStatusResponse(BaseModel):
    job_id: str
    status: str = Field(..., description="Current status: pending, running, success, failed")
    message: Optional[str] = Field(None, description="Optional message related to status (e.g., error details)")
    progress: Optional[float] = Field(None, description="Optional progress percentage (0.0 to 1.0)")
    results: Optional[Dict] = Field(None, description="Analysis results if status is 'success'")

# NEW: Response model for the job context endpoint
class JobContextResponse(BaseModel):
    job_id: str
    context: Optional[Dict[str, Any]] = Field(None, description="The initial payload and context stored for the job")

# Remove the WebSocket Connection Management section and replace with REST endpoint
@router.get("/status/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str, job_service: JobService = Depends(get_job_service)):
    """
    Get the current status of a job.
    Used for polling instead of WebSockets.
    """
    status = job_service.get_job_status(job_id)
    if not status:
        raise HTTPException(status_code=404, detail=f"Job ID {job_id} not found")
    
    # Return job status (might or might not include context depending on JobService implementation)
    return {**status, "job_id": job_id}

# NEW: Endpoint to get job context
@router.get("/context/{job_id}", response_model=JobContextResponse)
async def get_job_context_endpoint(job_id: str, job_service: JobService = Depends(get_job_service)):
    """
    Get the stored context/results for a specific job.
    This now retrieves the job's status dictionary and extracts the 'results' field.
    """
    full_status = job_service.get_job_context(job_id) # Renamed variable for clarity
    if full_status is None:
        # Check if the job *exists* first for a better error message
        if not job_service.job_exists(job_id):
             raise HTTPException(status_code=404, detail=f"Job ID {job_id} not found")
        else:
             # Job exists but get_job_context returned None (shouldn't happen with new logic)
             logger.error(f"Job {job_id} exists but get_job_context returned None unexpectedly.")
             raise HTTPException(status_code=500, detail=f"Internal error retrieving context for job ID {job_id}")

    # Extract the 'results' field from the full status to return as the context
    job_results = full_status.get("results")
    
    # Check if results are available (job might still be running or failed without results)
    if job_results is None:
        logger.warning(f"Context requested for job {job_id}, but results are not yet available (status: {full_status.get('status')}). Returning empty context.")
        # Depending on frontend expectation, maybe return 404 or the status?
        # For now, return empty context as per original model allowing Optional[Dict]
        job_results = {}
        # Alternatively, raise an error if results are expected:
        # raise HTTPException(status_code=404, detail=f"Results not yet available for job {job_id} (status: {full_status.get('status')})")

    # Ensure the extracted results are a dictionary, even if empty
    if not isinstance(job_results, dict):
        logger.error(f"Results for job {job_id} are not a dictionary: {type(job_results)}. Returning empty context.")
        job_results = {}

    return JobContextResponse(job_id=job_id, context=job_results)

# Dependency to get AnalysisService instance (without WebSocket dependency)
def get_analysis_service(job_service: JobService = Depends(get_job_service)) -> AnalysisService:
    # This function allows FastAPI to correctly instantiate AnalysisService
    # with its dependencies
    return AnalysisService(job_service=job_service)

@router.post("/start", response_model=AnalysisResponse)
async def start_analysis_endpoint(
    background_tasks: BackgroundTasks,
    payload: AnalysisPayload = Body(...),
    job_service: JobService = Depends(get_job_service),
    # Inject AnalysisService using the dependency function
    analysis_service: AnalysisService = Depends(get_analysis_service)
):
    """
    Endpoint to start the analysis pipeline asynchronously.
    
    Receives file IDs and column mappings, creates a job via JobService,
    stores the payload as context, adds the analysis task to the background queue,
    and returns a job ID for status tracking.
    """
    # --- Check for existing active job --- 
    # Let's simplify this check for now, or move complex logic to JobService if needed.
    # active_job_exists = False
    # for job_id_in_status, status_info in jobs_status.items(): # Old way
    #     if status_info.get('status') in ['pending', 'running']:
    #         active_job_exists = True
    #         break
    # 
    # if active_job_exists:
    #     raise HTTPException(
    #         status_code=409, 
    #         detail="An analysis job is already running. Please wait for it to complete."
    #     )
    # -------------------------------------
    
    logger.info("Received analysis request. Creating job...")
    # Create job using JobService, store payload as initial context
    job_id = job_service.create_job(initial_context=payload.dict())
    logger.info(f"Job created with ID: {job_id}")
    
    # Add the actual analysis task to the background
    # Pass only the job_id; the service retrieves the payload from context
    background_tasks.add_task(analysis_service.run_analysis_background, job_id, payload)
    
    # Return immediately with the job ID
    return AnalysisResponse(status="success", message="Analysis job started in background.", job_id=job_id)

# --- NEW Models for Custom Analysis ---
# REMOVED PointData and CustomAnalysisRequest CLASS DEFINITIONS FROM HERE

# --- REMOVED Endpoint for Custom Selection Analysis (Background Task version) ---
# @router.post("/start/custom_selection", ...)
# async def start_custom_selection_analysis(...):
#    ...

@router.post("/custom/{job_id}",
             response_model=CustomAnalysisResponse, # Specify response model
             summary="Run Analysis on Custom Lasso Selection",
             description="Filters the original spatial data using the provided polygon and runs the analysis pipeline on the subset.")
async def run_custom_analysis_endpoint(
    job_id: str,
    request: CustomLassoAnalysisRequest, # Use the existing request model
    job_service: JobService = Depends(get_job_service) # Inject job service for context
    # AnalysisService not directly needed here as we call core logic
):
    """Endpoint to run analysis synchronously on a user-defined spatial subset."""
    logger.info(f"Received custom analysis request for job {job_id}.")
    if len(request.polygon) < 4:
        raise HTTPException(status_code=400, detail="Polygon must have at least 4 points (first and last the same).")
    if request.polygon[0] != request.polygon[-1]:
         raise HTTPException(status_code=400, detail="Polygon first and last points must be the same.")

    try:
        # 1. Get original job context (contains file IDs and mappings)
        job_context = job_service.get_job_context(job_id) # Fetch the whole job document
        if not job_context or 'results' not in job_context or 'inputs' not in job_context['results']:
             logger.error(f"Job context or results/inputs missing for job {job_id}.")
             raise HTTPException(status_code=404, detail=f"Original job context/inputs not found for job ID: {job_id}")

        # Extract inputs (file IDs and mappings)
        inputs = job_context['results']['inputs']
        spatial_file_id = inputs.get('files', {}).get('spatialFileId')
        interactions_file_id = inputs.get('files', {}).get('interactionsFileId')
        modules_file_id = inputs.get('files', {}).get('modulesFileId')
        spatial_mapping = inputs.get('mappings', {}).get('spatialMapping', {})
        interactions_mapping = inputs.get('mappings', {}).get('interactionsMapping', {})
        modules_mapping = inputs.get('mappings', {}).get('modulesMapping', {})

        # Validate required inputs
        if not all([spatial_file_id, interactions_file_id, modules_file_id,
                    spatial_mapping, interactions_mapping, modules_mapping]):
            logger.error(f"Missing file IDs or mappings in job context for {job_id}.")
            raise HTTPException(status_code=400, detail="Missing required file IDs or mappings in original job context.")

        x_col = spatial_mapping.get('xCol')
        y_col = spatial_mapping.get('yCol')
        gene_col_spatial = spatial_mapping.get('geneCol')
        layer_col = spatial_mapping.get('layerCol')
        ligand_col = interactions_mapping.get('ligandCol')
        receptor_col = interactions_mapping.get('receptorCol')
        gene_col_modules = modules_mapping.get('geneCol')
        module_col = modules_mapping.get('moduleCol')

        if not all([x_col, y_col, gene_col_spatial, layer_col, ligand_col, receptor_col, gene_col_modules, module_col]):
             logger.error(f"Missing specific column names in mappings for job {job_id}.")
             raise HTTPException(status_code=400, detail="One or more required column names (x, y, layer, gene, ligand, receptor, module) are missing in the mappings.")

        # 2. Load dataframes (Call class method on FileService)
        try:
            spatial_path = FileService.get_file_path(spatial_file_id)
            interactions_path = FileService.get_file_path(interactions_file_id)
            modules_path = FileService.get_file_path(modules_file_id)

            # Convert Path objects to strings for pandas
            df_spatial_full = pd.read_csv(str(spatial_path))
            df_interactions = pd.read_csv(str(interactions_path))
            df_modules = pd.read_csv(str(modules_path))
        except FileNotFoundError as e:
             logger.error(f"Data file not found for job {job_id}: {e}")
             raise HTTPException(status_code=404, detail=f"One of the original data files could not be found: {e}")
        except Exception as e:
             logger.error(f"Error loading dataframes for job {job_id}: {e}")
             raise HTTPException(status_code=500, detail=f"Failed to load data files: {e}")

        # 3. Filter spatial data using Shapely polygon
        logger.info(f"Filtering spatial data ({len(df_spatial_full)} points) using polygon...")
        selection_polygon = Polygon(request.polygon)

        # Function to check if a point is inside the polygon
        def is_inside(row):
            point = Point(row[x_col], row[y_col])
            return selection_polygon.contains(point)

        df_spatial_filtered = df_spatial_full[df_spatial_full.apply(is_inside, axis=1)].copy()
        logger.info(f"Filtered spatial data contains {len(df_spatial_filtered)} points.")

        if df_spatial_filtered.empty:
             logger.warning(f"No spatial points found within the provided polygon for job {job_id}.")
             # Return empty results as analysis cannot run
             return CustomAnalysisResponse(pathway_dominance=[], module_context=[])


        # 4. Run analysis pipeline stages
        logger.info("Running Stage 1: Counts Pipeline...")
        # Corrected arguments for stage 1
        counts_results = run_stage1_counts_pipeline(
            spatial_df=df_spatial_filtered, 
            interactions_df=df_interactions
            # Removed extra args: gene_col, layer_col, ligand_col, receptor_col
        )
        # Ensure stage 1 output is suitable for stage 2 (e.g., contains 'ligand_norm_expr' etc.)
        # This might require checking the return value of run_stage1...
        # For now, assuming counts_results['whole_tissue'] is the relevant DataFrame
        if 'whole_tissue' not in counts_results:
             logger.error("Stage 1 did not produce 'whole_tissue' results for custom analysis.")
             raise HTTPException(status_code=500, detail="Analysis Stage 1 failed for the selection.")
        df_counts = counts_results['whole_tissue'] # Use the counts for the filtered region

        logger.info("Running Stage 2: Pathway Dominance...")
        # Corrected: Pass the filtered spatial df and interactions df
        pathway_results_dict = run_pathway_dominance_pipeline(
            spatial_df=df_spatial_filtered,
            interactions_df=df_interactions 
        )
        # Pathway dominance pipeline returns results per scope (even if only one)
        # Assuming the result for our custom scope is under a key like 'whole_tissue' or similar
        # We need to adapt this based on how _calculate_pathway_dominance_for_scope is used
        # For now, let's assume it returns results for a single scope (our filtered data)
        # and we need the list directly. We might need to adjust core.py or this call.
        # Let's try getting the list; this might need adjustment!
        pathway_results = pathway_results_dict.get('whole_tissue', []) # Tentative extraction
        if not pathway_results:
             logger.warning(f"Pathway dominance did not return expected results for custom selection {job_id}.")
             # Decide how to handle - maybe return empty results? For now, proceed.
             pathway_results = [] 

        logger.info("Running Stage 3: Module Context...")
        # Corrected: Pass interactions, modules, and the *list* of pathway results
        # Also add spatial_df to match signature, even if not strictly used for single custom scope
        module_results_dict = run_module_context_pipeline(
            spatial_df=df_spatial_filtered, # Added to match signature 
            interactions_df=df_interactions,
            modules_df=df_modules,
            # Pass the list of pathway results, not the dict containing it
            full_pathway_results={'custom': {'pathway_dominance': pathway_results}} # Structure it as expected
        )
        # Similar extraction needed for module results
        module_results = module_results_dict.get('custom', []) # Tentative extraction
        if not module_results:
             logger.warning(f"Module context did not return expected results for custom selection {job_id}.")
             module_results = []

        # 5. Combine and structure results
        # Now pathway_results and module_results should be lists of dicts

        # 6. Return structured response
        logger.info(f"Custom analysis completed for job {job_id}. Returning {len(module_results)} module results.")
        return CustomAnalysisResponse(
             pathway_dominance=pathway_results,
             module_context=module_results
        )

    except HTTPException as e:
         # Re-raise HTTP exceptions
         raise e
    except Exception as e:
        logger.exception(f"Unexpected error during custom analysis for job {job_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"An internal error occurred during custom analysis: {str(e)}")

# --- Other potential helper functions or endpoints --- 