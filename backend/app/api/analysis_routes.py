from fastapi import APIRouter, Depends, HTTPException, Body, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Dict, Optional, List
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
from app.models.analysis_models import AnalysisMapping, AnalysisPayload, PointData, CustomAnalysisRequest

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
    
    # Add job_id for consistency
    return {**status, "job_id": job_id}

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

# --- NEW Endpoint for Custom Analysis ---
@router.post("/start/custom_selection", response_model=JobStatusResponse, status_code=202) 
async def start_custom_selection_analysis(
    # Reorder parameters: non-defaults first
    background_tasks: BackgroundTasks,
    request: CustomAnalysisRequest = Body(...),
    job_service: JobService = Depends(get_job_service),
    analysis_service: AnalysisService = Depends(get_analysis_service)
):
    """
    Starts a background analysis task using a custom-defined set of points 
    (e.g., from a user's lasso selection).
    """
    # Create a job ID using JobService (provides central tracking)
    # Pass the request payload as initial context for reference if needed later
    job_id = job_service.create_job(initial_context=request.dict())
    logger.info(f"Received custom selection analysis request. Assigning Job ID: {job_id}")
    
    if not request.ligands or not request.receptors:
        logger.error(f"Custom analysis request for job {job_id} missing ligand or receptor points.")
        # Update job status to failed
        asyncio.create_task(job_service.update_job_status(job_id, status='failed', message="Ligand and receptor point lists cannot be empty."))
        raise HTTPException(status_code=400, detail="Ligand and receptor point lists cannot be empty.")

    # --- Delegate to AnalysisService --- 
    # Add a new method to AnalysisService to handle this specific workflow.
    # This service method will be responsible for: 
    #   - Creating DataFrames from the point lists.
    #   - Handling column naming/validation.
    #   - Loading necessary auxiliary data (interactions, modules - requires thought on how to locate these).
    #   - Calling an adapted version of the core analysis logic.
    #   - Updating job status via JobService.
    
    # *** TODO: Implement `run_custom_analysis_background` in `AnalysisService` ***
    background_tasks.add_task(analysis_service.run_custom_analysis_background, job_id, request)

    # Initial status update via JobService
    initial_status_message = "Custom selection analysis job accepted."
    asyncio.create_task(job_service.update_job_status(
        job_id,
        status="pending",
        message=initial_status_message,
        progress=0.0
    ))
    
    # Return immediate response
    return JobStatusResponse(
        job_id=job_id,
        status="pending",
        message=initial_status_message
    ) 