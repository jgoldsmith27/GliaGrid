from fastapi import APIRouter, Depends, HTTPException, Body, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Any
import asyncio
import uuid
import json
import logging
import pandas as pd # Import pandas for DataFrame creation

# Import JobService and AnalysisService
from app.services.job_service import JobService, get_job_service
# Only import the class, not the non-existent dependency function from the service file
from app.services.analysis_service import AnalysisService
# Import models from the correct location
from app.models.analysis_models import AnalysisMapping, AnalysisPayload, CustomLassoAnalysisRequest, CustomAnalysisResponse

# Assuming core analysis logic is importable
# from ..analysis_logic.core import run_analysis_pipeline # Placeholder - core logic needs adaptation

# REMOVING unused imports that were causing ModuleNotFoundError
# from ..utils.task_manager import TaskManager, TaskStatus 
# from ..utils.file_manager import FileManager, get_file_manager, get_task_manager

# REMOVING import from non-existent data_models module
# from .data_models import AnalysisRequest, FileMapping, JobStatusResponse, AnalysisResult # Assuming existing data_models

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
    request: CustomLassoAnalysisRequest, # Use the new request model
    analysis_service: AnalysisService = Depends(get_analysis_service), # Inject service
    job_service: JobService = Depends(get_job_service) # Inject job service for context?
):
    """Endpoint to run analysis on a user-defined spatial subset."""
    try:
        print(f"Received custom analysis request for job {job_id} with polygon vertex count: {len(request.polygon)}")
        # We might need the original job context (file IDs, mappings) 
        # Option 1: Fetch context using job_service
        initial_context = await job_service.get_job_context(job_id)
        if not initial_context:
            raise HTTPException(status_code=404, detail=f"Original job context not found for job ID: {job_id}")

        # Ensure we have necessary info (e.g., file IDs from context)
        # This depends on what `run_custom_analysis` needs
        # Example: 
        # spatial_file_id = initial_context.get('inputs', {}).get('files', {}).get('spatialFileId')
        # if not spatial_file_id:
        #     raise HTTPException(status_code=400, detail="Spatial file ID missing from original job context.")

        results = await analysis_service.run_custom_analysis(
            original_job_id=job_id, 
            polygon_coords=request.polygon,
            initial_context=initial_context # Pass context if needed by service
        )
        
        print(f"Custom analysis completed for job {job_id}.")
        return results # Should match CustomAnalysisResponse structure

    except FileNotFoundError as e:
        print(f"Error during custom analysis for job {job_id}: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        print(f"Value error during custom analysis for job {job_id}: {e}")
        raise HTTPException(status_code=400, detail=str(e)) 
    except Exception as e:
        print(f"Unexpected error during custom analysis for job {job_id}: {e}")
        # Log the full traceback for unexpected errors
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"An unexpected error occurred: {str(e)}") 