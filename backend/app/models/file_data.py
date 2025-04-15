"""Data models for file handling."""
from typing import Dict, List, Optional, Any, Union
from pydantic import BaseModel, Field


class FileInfo(BaseModel):
    """Metadata information for files, especially H5AD files."""
    shape: str = Field(..., description="Data shape information")
    obs_keys: List[str] = Field(default_factory=list, description="Observation keys")
    var_keys: List[str] = Field(default_factory=list, description="Variable keys")
    obsm_keys: List[str] = Field(default_factory=list, description="Observation matrix keys")


class FilePreviewResult(BaseModel):
    """Result of file preview operation."""
    fileId: str = Field(..., description="Unique identifier for the uploaded file")
    headers: List[str] = Field(..., description="Column headers")
    previewRows: List[Dict[str, Any]] = Field(..., description="Preview rows")
    fileInfo: Optional[FileInfo] = Field(None, description="Additional file metadata")


class ErrorResponse(BaseModel):
    """Error response model."""
    detail: str = Field(..., description="Error details") 