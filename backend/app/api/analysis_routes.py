from fastapi import APIRouter, Depends, HTTPException, Body, BackgroundTasks, WebSocket, WebSocketDisconnect
from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Set
import asyncio
import uuid
import json
import logging

# Import JobService and AnalysisService
from app.services.job_service import JobService, get_job_service
from app.services.analysis_service import AnalysisService 
from app.models.analysis_models import AnalysisMapping, AnalysisPayload # Import models

logger = logging.getLogger(__name__)

router = APIRouter(
    prefix="/analysis",
    tags=["Analysis"],
    responses={404: {"description": "Not found"}},
)

# Separate router for WebSocket endpoint without prefix interference
ws_router = APIRouter(
    tags=["Analysis WebSocket"]
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

# --- WebSocket Connection Management ---
class ConnectionManager:
    def __init__(self):
        # Dictionary mapping job_id to a set of active WebSocket connections
        self.active_connections: Dict[str, Set[WebSocket]] = {}

    async def connect(self, websocket: WebSocket, job_id: str):
        await websocket.accept()
        if job_id not in self.active_connections:
            self.active_connections[job_id] = set()
        self.active_connections[job_id].add(websocket)
        print(f"WebSocket connected for job {job_id}. Total connections for job: {len(self.active_connections[job_id])}")

    def disconnect(self, websocket: WebSocket, job_id: str):
        if job_id in self.active_connections:
            self.active_connections[job_id].remove(websocket)
            print(f"WebSocket disconnected for job {job_id}. Remaining connections: {len(self.active_connections[job_id])}")
            # Clean up job_id entry if no more connections
            if not self.active_connections[job_id]:
                del self.active_connections[job_id]

    async def send_update(self, job_id: str, message: dict):
        if job_id in self.active_connections:
            # Create a task for each send operation to run them concurrently
            tasks = [conn.send_text(json.dumps(message)) for conn in self.active_connections[job_id]]
            await asyncio.gather(*tasks)
            # print(f"Sent update to {len(tasks)} connections for job {job_id}") # Optional: verbose logging

manager = ConnectionManager()
# Pass the manager instance to the AnalysisService or make it accessible
# Option 1: Make it globally accessible (simpler for now)
# This instantiation needs to happen where FastAPI can manage dependencies properly,
# typically not directly at the module level like this if services depend on each other.
# We will rely on FastAPI's dependency injection instead.
# analysis_service_instance = AnalysisService(connection_manager=manager)
# -------------------------------------

# Dependency to get AnalysisService instance
# Note: AnalysisService now depends on JobService, FastAPI handles this chain.
def get_analysis_service(manager_instance: ConnectionManager = Depends(lambda: manager), # Inject manager
                         job_service: JobService = Depends(get_job_service)) -> AnalysisService:
    # This function allows FastAPI to correctly instantiate AnalysisService
    # with its dependencies (manager and job_service)
    return AnalysisService(connection_manager=manager_instance, job_service=job_service)

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

@router.get("/status/{job_id}", response_model=JobStatusResponse)
async def get_analysis_status(
    job_id: str, 
    job_service: JobService = Depends(get_job_service)
):
    """
    Endpoint to check the status and potentially retrieve results of an analysis job.
    Uses JobService to retrieve status.
    """
    logger.info(f"Checking status for job ID: {job_id}")
    status_info = job_service.get_job_status(job_id)

    if not status_info:
        logger.warning(f"Job ID {job_id} not found for status check.")
        raise HTTPException(status_code=404, detail=f"Job ID {job_id} not found.")

    logger.info(f"Status found for job ID {job_id}: {status_info}")
    # Ensure all potential fields from the status_info dict are included
    return JobStatusResponse(
        job_id=job_id,
        status=status_info.get("status", "unknown"),
        message=status_info.get("message"),
        progress=status_info.get("progress"),
        results=status_info.get("results")
    )

@ws_router.websocket("/ws/analysis/status/{job_id}")
async def websocket_status_endpoint(
    websocket: WebSocket, 
    job_id: str,
    # Inject JobService directly into the WebSocket endpoint
    job_service: JobService = Depends(get_job_service) 
):
    """
    WebSocket endpoint for real-time analysis status updates.
    Uses JobService to check initial status.
    """
    await manager.connect(websocket, job_id)
    try:
        # Send the current status immediately upon connection using JobService
        current_status = job_service.get_job_status(job_id)
        if current_status:
            # Add job_id for consistency, although it's in the URL
            await manager.send_update(job_id, {**current_status, "job_id": job_id}) 
        else:
            # If job_id doesn't exist when WS connects
            logger.warning(f"WebSocket connection attempt for non-existent job ID: {job_id}")
            await manager.send_update(job_id, {"status": "error", "message": f"Job ID {job_id} not found.", "job_id": job_id})
            # Consider closing the connection if job doesn't exist after sending error
            # await websocket.close(code=1008) # Policy Violation or similar
            # return # Exit the handler

        # Keep the connection alive
        while True:
             # Keep connection open to push server updates
             # Wait for disconnect instead of client messages
             await websocket.receive_text() # Needed to detect disconnect properly

    except WebSocketDisconnect:
        manager.disconnect(websocket, job_id)
        logger.info(f"Client disconnected from job {job_id}")
    except Exception as e:
        manager.disconnect(websocket, job_id)
        logger.exception(f"Error in WebSocket connection for job {job_id}: {e}") 