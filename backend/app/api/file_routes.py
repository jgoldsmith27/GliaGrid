"""File handling API routes."""
from fastapi import APIRouter, UploadFile, File, HTTPException, Depends
from pathlib import Path
import os
import tempfile
from contextlib import asynccontextmanager

from app.services.file_service import FileService
from app.models.file_data import FilePreviewResult, ErrorResponse

# Create router
router = APIRouter(prefix="/api/file", tags=["file"])

# Global temporary directory that persists during application lifetime
TEMP_DIR = Path(tempfile.mkdtemp(prefix="gliagrid_"))

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


@router.post("/preview", response_model=FilePreviewResult, responses={400: {"model": ErrorResponse}, 500: {"model": ErrorResponse}})
async def preview_file(
    file: UploadFile = File(...),
    temp_dir: Path = Depends(get_temp_dir)
):
    """Preview file contents (headers and sample rows)
    
    Args:
        file: The uploaded file
        temp_dir: Directory to store temporary files
        
    Returns:
        FilePreviewResult: Preview data including headers and rows
        
    Raises:
        HTTPException: If file processing fails
    """
    return await FileService.process_file(file, temp_dir) 