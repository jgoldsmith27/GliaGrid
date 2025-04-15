"""Service for handling different file formats."""
import os
import tempfile
import shutil
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
import pandas as pd
import numpy as np
import anndata as ad
from fastapi import UploadFile, HTTPException
from contextlib import contextmanager

from app.models.file_data import FileInfo, FilePreviewResult


class FileService:
    """Service for handling different file formats and converting them to a unified data structure.
    
    This service provides methods to:
    1. Process various file formats (CSV, H5AD)
    2. Extract preview data for UI display
    3. Convert files to pandas DataFrames for uniform handling
    4. Store uploaded files temporarily with a unique ID
    5. Retrieve file paths based on their ID
    """
    
    # Supported file extensions
    SUPPORTED_FORMATS = {
        '.csv': 'process_csv',
        '.h5ad': 'process_h5ad',
    }
    
    # Directory for storing uploaded files (relative to backend script location)
    # Go up one level from backend/app/services to backend/, then up to project root, then into .temp_uploads
    TEMP_STORAGE_DIR = Path(__file__).resolve().parent.parent.parent / ".temp_uploads"
    
    @classmethod
    def initialize_storage(cls):
        """Create the temporary storage directory if it doesn't exist."""
        cls.TEMP_STORAGE_DIR.mkdir(parents=True, exist_ok=True)

    @classmethod
    async def save_and_preview_file(cls, file: UploadFile) -> FilePreviewResult:
        """Saves the uploaded file to temporary storage, generates a preview, 
        and returns the preview data along with a unique file ID.
        
        Args:
            file: The uploaded file
            
        Returns:
            FilePreviewResult: Preview data including headers, sample rows, and fileId.
            
        Raises:
            HTTPException: If file processing fails or format is unsupported.
        """
        cls.initialize_storage() # Ensure storage directory exists

        if not file.filename:
            raise HTTPException(status_code=400, detail="No filename provided")
            
        # Get file extension and check if supported
        file_ext = os.path.splitext(file.filename.lower())[1]
        if file_ext not in cls.SUPPORTED_FORMATS:
            raise HTTPException(
                status_code=400, 
                detail=f"Unsupported file format: {file_ext}. Supported formats: {', '.join(cls.SUPPORTED_FORMATS.keys())}"
            )
        
        # Generate a unique filename to avoid collisions
        file_id = str(uuid.uuid4())
        unique_filename = f"{file_id}{file_ext}"
        file_path = cls.TEMP_STORAGE_DIR / unique_filename

        try:
            # Save the uploaded file content to the persistent temp location
            with open(file_path, "wb") as buffer:
                shutil.copyfileobj(file.file, buffer)
                
            # Process the file to get preview data components
            process_method = getattr(cls, cls.SUPPORTED_FORMATS[file_ext])
            
            headers: List[str]
            preview_rows: List[Dict[str, Any]]
            file_info: Optional[FileInfo] = None
            
            if file_ext == '.csv':
                headers, preview_rows = process_method(file_path)
            elif file_ext == '.h5ad':
                headers, preview_rows, file_info = process_method(file_path)
            else:
                # Should not happen due to earlier check, but handle defensively
                raise HTTPException(status_code=500, detail="Internal error: Unexpected file type during processing.")
            
            # Construct the final result including fileId
            return FilePreviewResult(
                fileId=file_id,
                headers=headers,
                previewRows=preview_rows,
                fileInfo=file_info
            )

        except Exception as e:
            # Clean up the saved file if processing fails
            if os.path.exists(file_path):
                os.unlink(file_path)
            raise HTTPException(status_code=500, detail=f"Failed to process {file_ext} file: {str(e)}")
        finally:
            if file.file and not file.file.closed:
                 file.file.close()

    @classmethod
    def get_file_path(cls, file_id: str) -> Path:
        """Retrieves the full path of a stored file based on its ID.
        
        Args:
            file_id: The unique identifier of the file.
            
        Returns:
            Path: The full path to the stored file.
            
        Raises:
            HTTPException: If the file corresponding to the ID is not found or has an unexpected extension.
        """
        # Search for the file with the given ID in the temp directory
        potential_files = list(cls.TEMP_STORAGE_DIR.glob(f"{file_id}.*"))
        
        if not potential_files:
            raise HTTPException(status_code=404, detail=f"File with ID {file_id} not found.")
            
        if len(potential_files) > 1:
             # This shouldn't happen with UUIDs, but handle defensively
             raise HTTPException(status_code=500, detail=f"Multiple files found for ID {file_id}. Storage inconsistency.")
             
        file_path = potential_files[0]
        
        # Optional: Verify the extension matches supported formats if needed
        # file_ext = file_path.suffix.lower()
        # if file_ext not in cls.SUPPORTED_FORMATS:
        #     raise HTTPException(status_code=400, detail=f"File {file_id} has an unsupported extension: {file_ext}")
            
        return file_path

    @staticmethod
    def process_csv(file_path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
        """Process CSV file and return headers and preview rows.
        
        Args:
            file_path: Path to the CSV file
            
        Returns:
            Tuple[List[str], List[Dict[str, Any]]]: Headers and preview rows
            
        Raises:
            Exception: If CSV processing fails
        """
        # Read only the first 6 rows (header + 5 data rows) for preview efficiency
        try:
            df = pd.read_csv(file_path, nrows=6)
        except pd.errors.EmptyDataError:
            # Handle case where file has fewer than 6 rows or is empty
            df = pd.read_csv(file_path) # Read whatever is there
        except Exception as e:
             raise HTTPException(status_code=400, detail=f"Error reading CSV headers/preview: {str(e)}")

        # Get headers and preview rows
        headers = df.columns.tolist()
        preview_rows = df.head(5).to_dict('records')
        
        return headers, preview_rows
    
    @staticmethod
    def process_h5ad(file_path: Path) -> Tuple[List[str], List[Dict[str, Any]], Optional[FileInfo]]:
        """Process H5AD file and return headers, preview rows, and file info.
        
        Args:
            file_path: Path to the H5AD file
            
        Returns:
            Tuple containing:
            - List of headers
            - List of preview rows
            - Optional FileInfo metadata
            
        Raises:
            Exception: If H5AD processing fails
        """
        # Read H5AD file
        adata = ad.read_h5ad(file_path)
        
        # Convert to a dataframe-like structure for consistency with CSV
        headers, preview_rows, file_info = FileService._extract_h5ad_preview(adata)
        
        if not headers:
            raise Exception("No usable data found in H5AD file")
        
        return headers, preview_rows, file_info
    
    @staticmethod
    def _extract_h5ad_preview(adata: ad.AnnData) -> Tuple[List[str], List[Dict[str, Any]], FileInfo]:
        """Extract preview data from AnnData object.
        
        Args:
            adata: The AnnData object
            
        Returns:
            Tuple containing:
            - List of headers
            - List of preview rows
            - FileInfo metadata
        """
        headers = []
        preview_rows = []
        
        # Create file info
        file_info = FileInfo(
            shape=f"{adata.shape[0]} cells x {adata.shape[1]} genes",
            obs_keys=list(adata.obs.columns) if adata.obs is not None and hasattr(adata.obs, 'columns') else [],
            var_keys=list(adata.var.columns) if adata.var is not None and hasattr(adata.var, 'columns') else [],
            obsm_keys=list(adata.obsm.keys()) if hasattr(adata, 'obsm') and adata.obsm is not None else []
        )
        
        # Default to obs if available (cell metadata)
        if adata.obs is not None and not adata.obs.empty:
            # Include index as a column
            df = adata.obs.reset_index()
            headers = df.columns.tolist()
            preview_rows = df.head(5).to_dict('records')
        
        # Also extract X matrix data and convert it to columns
        if len(preview_rows) == 0 and adata.X is not None:
            # Convert sparse matrix to dense if needed
            if hasattr(adata.X, 'toarray'):
                X_dense = adata.X.toarray()
            else:
                X_dense = adata.X
            
            # Create dataframe with gene/feature names
            if X_dense.shape[1] <= 50:  # Only use all columns if reasonable number
                gene_cols = [f"gene_{i}" for i in range(X_dense.shape[1])]
                if adata.var is not None and not adata.var.empty:
                    gene_cols = adata.var.index.tolist()[:X_dense.shape[1]]
                
                df = pd.DataFrame(X_dense[:5], columns=gene_cols)
                # Add cell IDs if available
                if adata.obs is not None and not adata.obs.empty:
                    df.insert(0, 'cell_id', adata.obs.index.tolist()[:5])
                
                headers = df.columns.tolist()
                preview_rows = df.head(5).to_dict('records')
            else:
                # Just sample a few genes/features to avoid overwhelming response
                sampled_indices = np.random.choice(X_dense.shape[1], 20, replace=False)
                sampled_cols = [f"gene_{i}" for i in sampled_indices]
                if adata.var is not None and not adata.var.empty:
                    var_index = adata.var.index.tolist()
                    sampled_cols = [var_index[i] for i in sampled_indices]
                
                df = pd.DataFrame(X_dense[:5, sampled_indices], columns=sampled_cols)
                # Add cell IDs if available
                if adata.obs is not None and not adata.obs.empty:
                    df.insert(0, 'cell_id', adata.obs.index.tolist()[:5])
                
                headers = df.columns.tolist()
                preview_rows = df.head(5).to_dict('records')
        
        # Add spatial coordinates if available
        if 'spatial' in adata.obsm:
            spatial_df = pd.DataFrame(
                adata.obsm['spatial'][:5], 
                columns=['x', 'y'] + [f'dim_{i}' for i in range(2, adata.obsm['spatial'].shape[1])]
            )
            
            # Add cell IDs if available
            if adata.obs is not None and not adata.obs.empty:
                spatial_df.insert(0, 'cell_id', adata.obs.index.tolist()[:5])
            
            # If we already have data, just return both sets of headers
            if len(headers) > 0:
                added_headers = spatial_df.columns.tolist()
                headers += [h for h in added_headers if h not in headers]
            else:
                headers = spatial_df.columns.tolist()
                preview_rows = spatial_df.head(5).to_dict('records')
                
        return headers, preview_rows, file_info
        
    @staticmethod
    def dataframe_from_file(file_path: Path) -> pd.DataFrame:
        """Convert any supported file to a pandas DataFrame.
        
        This method serves as a unified interface for getting a DataFrame
        from any supported file format.
        
        Args:
            file_path: Path to the file
            
        Returns:
            pd.DataFrame: The data as a pandas DataFrame
            
        Raises:
            ValueError: If file format is not supported
        """
        file_ext = os.path.splitext(file_path.name.lower())[1]
        
        if file_ext == '.csv':
            return pd.read_csv(file_path)
        elif file_ext == '.h5ad':
            adata = ad.read_h5ad(file_path)
            
            # Combine obs data with spatial data if available
            result_df = adata.obs.reset_index()
            
            if 'spatial' in adata.obsm:
                spatial_cols = ['x', 'y'] + [f'dim_{i}' for i in range(2, adata.obsm['spatial'].shape[1])]
                spatial_df = pd.DataFrame(adata.obsm['spatial'], columns=spatial_cols)
                # Ensure it has the same index as obs
                spatial_df.index = adata.obs.index
                # Join the dataframes
                result_df = pd.concat([result_df, spatial_df], axis=1)
                
            return result_df
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
            
    @classmethod
    def add_file_format(cls, extension: str, method_name: str) -> None:
        """Register a new file format handler.
        
        This method can be used to extend the service with new file formats.
        
        Args:
            extension: File extension (e.g., '.txt')
            method_name: Name of the method to handle this file format
            
        Raises:
            ValueError: If the method doesn't exist
        """
        if not hasattr(cls, method_name):
            raise ValueError(f"Method {method_name} does not exist in FileService")
            
        cls.SUPPORTED_FORMATS[extension] = method_name 