"""File handling API routes."""
from fastapi import APIRouter, UploadFile, File, HTTPException, Depends, Response, status
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

@router.delete("/{file_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_temporary_file(file_id: str):
    """
    Deletes a specific temporary file from the backend storage based on its ID.
    Returns 204 No Content on success.
    Returns 404 Not Found if the file doesn't exist.
    Returns 500 Internal Server Error on other deletion errors.
    """
    try:
        # Use FileService to find the path associated with the ID
        file_path = FileService.get_file_path(file_id)
        
        # Check if file actually exists at the path found
        if not file_path.is_file():
             # This case might occur if get_file_path logic has issues or file was already deleted
             raise HTTPException(status_code=404, detail=f"File with ID {file_id} not found at expected path.")

        # Attempt to delete the file
        os.unlink(file_path)
        print(f"Deleted temporary file: {file_path}") # Optional logging
        # No need to return anything on successful deletion with 204 status code
        return Response(status_code=status.HTTP_204_NO_CONTENT)
        
    except HTTPException as e:
        # Re-raise specific HTTP exceptions (like 404 from get_file_path)
        raise e
    except OSError as e:
        # Handle specific OS errors during deletion (e.g., permissions)
        print(f"Error deleting file {file_path}: {e}") # Optional logging
        raise HTTPException(status_code=500, detail=f"Failed to delete file {file_id}: OS Error")
    except Exception as e:
        # Catch any other unexpected errors
        print(f"Unexpected error deleting file {file_id}: {e}") # Optional logging
        raise HTTPException(status_code=500, detail=f"Failed to delete file {file_id}: Server Error") 