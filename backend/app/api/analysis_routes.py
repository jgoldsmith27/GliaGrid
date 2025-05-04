from fastapi import APIRouter, Depends, HTTPException, Body, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Any, Union
import asyncio
import uuid
import json
import logging
import pandas as pd
from shapely.geometry import Point, Polygon # Import Shapely
import time

# Import JobService and AnalysisService
from app.services.job_service import JobService, get_job_service
from app.services.analysis_service import AnalysisService
from app.models.analysis_models import (
    AnalysisMapping, AnalysisPayload, 
    CustomLassoAnalysisRequest, # Use updated request model
    CustomAnalysisResultsBundle # Use NEW bundle response model
)
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

@router.post("/custom/{job_id}",
             response_model=CustomAnalysisResultsBundle, # MODIFIED response model
             summary="Run Analysis on Custom Lasso Selection",
             description="Filters the original spatial data using the provided polygon and runs the analysis pipeline on the subset, returning both whole-selection and layer-specific results.")
async def run_custom_analysis_endpoint(
    job_id: str,
    request: CustomLassoAnalysisRequest, # Uses updated model (no aggregation field)
    job_service: JobService = Depends(get_job_service), # Inject job service for context
    analysis_service: AnalysisService = Depends(get_analysis_service) # Inject analysis service
) -> CustomAnalysisResultsBundle: # MODIFIED return type hint
    """Endpoint to run analysis synchronously on a user-defined spatial subset.
       Uses the AnalysisService to perform the loading, filtering, and analysis.
    """
    request_start_time = time.time() # DEBUG: Record start time
    logger.debug(f"[Custom Analysis {job_id}] Received request with {len(request.polygon)} polygon points.") # DEBUG

    logger.info(f"Received custom analysis request for job {job_id}.") # Removed aggregation log
    if len(request.polygon) < 4:
        raise HTTPException(status_code=400, detail="Polygon must have at least 4 points (first and last the same).")
    if request.polygon[0] != request.polygon[-1]:
         raise HTTPException(status_code=400, detail="Polygon first and last points must be the same.")
    # REMOVED: Validation for request.aggregation

    try:
        # 1. Get original job context (contains file IDs and mappings)
        # The AnalysisService will handle retrieving and validating context
        # We need the *full* job status document which contains the results/inputs
        full_job_status = job_service.get_job_context(job_id)
        if not full_job_status:
             raise HTTPException(status_code=404, detail=f"Original job context not found for job ID: {job_id}")
        
        # Extract the results section which contains the inputs
        job_results_context = full_job_status.get('results')
        if not job_results_context:
            logger.error(f"Job results context missing for job {job_id}. Status: {full_job_status.get('status')}")
            raise HTTPException(status_code=404, detail=f"Original job results context not found for job ID: {job_id}. Job might not be complete or failed.")

        # 2. Call the AnalysisService method to run the custom analysis
        # Pass the original job ID, polygon coordinates, and the *results* context (which holds inputs)
        # REMOVED: aggregation_level argument
        custom_analysis_results: CustomAnalysisResultsBundle = await analysis_service.run_custom_analysis(
            original_job_id=job_id,
            polygon_coords=request.polygon,
            initial_context=job_results_context # Pass the results context containing inputs
            # REMOVED: aggregation_level=request.aggregation 
        )

        logger.info(f"Custom analysis completed for job {job_id}. Returning results.") # Removed aggregation log
        logger.debug(f"[Custom Analysis {job_id}] Endpoint duration: {time.time() - request_start_time:.4f} seconds.") # DEBUG
        return custom_analysis_results

    except HTTPException as http_exc:
        # Re-raise HTTPExceptions directly
        raise http_exc
    except ValueError as val_err:
        # Catch specific validation errors from service/core logic
        logger.error(f"Validation error during custom analysis for {job_id}: {val_err}")
        raise HTTPException(status_code=400, detail=str(val_err))
    except FileNotFoundError as fnf_err:
        logger.error(f"File not found during custom analysis for {job_id}: {fnf_err}")
        raise HTTPException(status_code=404, detail=str(fnf_err))
    except Exception as e:
        # Catch-all for other unexpected errors
        logger.exception(f"Unexpected error during custom analysis endpoint for job {job_id}: {e}")
        logger.debug(f"[Custom Analysis {job_id}] Endpoint failed after {time.time() - request_start_time:.4f} seconds.") # DEBUG
        raise HTTPException(status_code=500, detail="An unexpected error occurred during custom analysis.")

# --- Other potential helper functions or endpoints --- 