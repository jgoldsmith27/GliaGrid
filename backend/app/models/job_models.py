from typing import Optional, Dict, Any
from pydantic import BaseModel

class JobStatus(BaseModel):
    progress: float = 0.0
    message: str = "Pending"
    results: Optional[Dict[str, Any]] = None
    stage_id: Optional[str] = None
    current_scope: Optional[str] = None

class JobContext(BaseModel):
    # Define any fields needed for job context here, or leave empty if none yet
    pass

    current_scope: Optional[str] = None # Added for detailed progress