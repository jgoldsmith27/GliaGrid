from fastapi import APIRouter, Depends, HTTPException, Body, BackgroundTasks, WebSocket, WebSocketDisconnect
from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Set
import asyncio
import uuid
import json

# Assume AnalysisService will be created in app/services/analysis_service.py
from app.services.analysis_service import AnalysisService, jobs_status # Import service and status dict
from app.models.analysis_models import AnalysisMapping, AnalysisPayload # Import models

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
analysis_service_instance = AnalysisService(connection_manager=manager)
# -------------------------------------

# Dependency to get AnalysisService instance (if needed, or use class methods)
# def get_analysis_service():
#     return AnalysisService() 

@router.post("/start", response_model=AnalysisResponse)
async def start_analysis_endpoint(
    background_tasks: BackgroundTasks,
    payload: AnalysisPayload = Body(...)
):
    """
    Endpoint to start the analysis pipeline asynchronously.
    
    Receives file IDs and column mappings, adds the analysis job to the 
    background queue, and returns a job ID for status tracking.
    """
    # --- Check for existing active job --- 
    active_job_exists = False
    for job_id_in_status, status_info in jobs_status.items():
        if status_info.get('status') in ['pending', 'running']:
            active_job_exists = True
            break
    
    if active_job_exists:
        raise HTTPException(
            status_code=409, 
            detail="An analysis job is already running. Please wait for it to complete."
        )
    # -------------------------------------
    
    job_id = str(uuid.uuid4())
    print(f"Received analysis request, assigning job ID: {job_id}")
    
    # Add job to status dictionary (initial state)
    initial_status = {"status": "pending", "message": "Analysis job queued.", "progress": 0.0}
    jobs_status[job_id] = initial_status

    # Instantiate service and add the actual analysis task to the background
    # Use the globally created instance that has the connection manager
    background_tasks.add_task(analysis_service_instance.run_analysis_background, job_id, payload)
    
    # Return immediately with the job ID
    return AnalysisResponse(status="success", message="Analysis job started in background.", job_id=job_id)

@router.get("/status/{job_id}", response_model=JobStatusResponse)
async def get_analysis_status(job_id: str):
    """
    Endpoint to check the status and potentially retrieve results of an analysis job.
    (This is still useful for non-WebSocket clients or quick checks)
    """
    print(f"Checking status for job ID: {job_id}")
    status_info = jobs_status.get(job_id)

    if not status_info:
        raise HTTPException(status_code=404, detail=f"Job ID {job_id} not found.")

    print(f"Status found for job ID {job_id}: {status_info}")
    # Ensure all potential fields from the status_info dict are included
    return JobStatusResponse(
        job_id=job_id,
        status=status_info.get("status", "unknown"),
        message=status_info.get("message"),
        progress=status_info.get("progress"),
        results=status_info.get("results")
    ) 

@ws_router.websocket("/ws/analysis/status/{job_id}")
async def websocket_status_endpoint(websocket: WebSocket, job_id: str):
    """
    WebSocket endpoint for real-time analysis status updates.
    """
    await manager.connect(websocket, job_id)
    try:
        # Send the current status immediately upon connection
        current_status = jobs_status.get(job_id)
        if current_status:
            await manager.send_update(job_id, {**current_status, "job_id": job_id}) # Include job_id in message
        else:
            # If job_id doesn't exist when WS connects (e.g., race condition or invalid ID)
            await manager.send_update(job_id, {"status": "error", "message": f"Job ID {job_id} not found.", "job_id": job_id})

        # Keep the connection alive, listening for client messages (optional)
        while True:
            # You might receive messages from the client here if needed
            # For now, we just keep the connection open to push server updates
            await websocket.receive_text() # Waits for a message, keeps connection open
            # If you don't expect client messages, you could use asyncio.sleep() 
            # await asyncio.sleep(60) # Keep alive by sleeping, less reactive to disconnects

    except WebSocketDisconnect:
        manager.disconnect(websocket, job_id)
        print(f"Client disconnected from job {job_id}")
    except Exception as e:
        manager.disconnect(websocket, job_id)
        print(f"Error in WebSocket connection for job {job_id}: {e}") 