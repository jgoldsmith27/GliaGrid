"""File handling API routes."""
from fastapi import APIRouter, UploadFile, File, HTTPException, Depends
from pathlib import Path
import os
import tempfile
from contextlib import asynccontextmanager

from app.services.file_service import FileService
from app.models.file_data import FilePreviewResult, ErrorResponse

# Temporary directory for file uploads (could be configured)
TEMP_DIR = Path("temp") 
TEMP_DIR.mkdir(exist_ok=True) # Ensure directory exists

# Create API router
router = APIRouter(
    prefix="/file",
    tags=["File Handling"],
    responses={404: {"description": "Not found"}}
)

# Register cleanup on application shutdown
@asynccontextmanager
async def lifespan(app):
    # Setup: temp directory is created on startup
    yield
    # Cleanup: remove temp directory on shutdown
    if TEMP_DIR.exists():
        import shutil
        shutil.rmtree(TEMP_DIR)

def get_temp_dir() -> Path:
    """Dependency to provide the temp directory.
    
    Returns a Path to a secure temporary directory that's created at app startup
    and cleaned up on app shutdown.
    """
    return TEMP_DIR

@router.post("/preview", response_model=FilePreviewResult)
async def upload_and_preview_file(file: UploadFile = File(...)):
    """
    Uploads a file, saves it temporarily, generates a preview, 
    and returns the preview data along with a unique file ID.
    
    Handles both CSV and H5AD files.
    """
    try:
        # Use the refactored service method
        return await FileService.save_and_preview_file(file)
    except HTTPException as e:
        # Re-raise known HTTP exceptions
        raise e
    except Exception as e:
        # Catch unexpected errors during processing
        raise HTTPException(status_code=500, detail=f"Error processing file: {str(e)}") 