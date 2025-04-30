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
import logging
import datetime
from rtree import index
import json

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

    # ==== ENHANCED DATA ACCESS METHODS ====
    
    @classmethod
    def get_spatial_data_chunked(cls, file_id: str, options: Dict[str, Any] = None) -> pd.DataFrame:
        """Get spatial data with support for chunking, sampling, and region filtering.
        
        Args:
            file_id: The file ID to load from
            options: Dict containing options for data loading:
                - resolution: float (0.0-1.0) for random sampling rate
                - offset: int for starting row 
                - limit: int for maximum rows to return
                - region: Dict with x_min, y_min, x_max, y_max for spatial filtering
                
        Returns:
            DataFrame of the requested data chunk/sample
            
        Raises:
            HTTPException: If file not found or processing fails
        """
        options = options or {}
        
        # Get logging context
        log_ctx = {
            'file_id': file_id,
            'resolution': options.get('resolution', 1.0),
            'offset': options.get('offset'),
            'limit': options.get('limit'),
            'has_region': 'region' in options
        }
        
        logger.info(f"Loading chunked spatial data", extra=log_ctx)
        
        try:
            # Get the file path 
            file_path = cls.get_file_path(file_id)
            
            # Determine file type and call appropriate chunked reader
            if file_path.suffix.lower() == '.csv':
                return cls._get_csv_chunked(file_path, options)
            elif file_path.suffix.lower() == '.h5ad':
                return cls._get_h5ad_chunked(file_path, options)
            else:
                raise HTTPException(
                    status_code=400, 
                    detail=f"Unsupported file type for chunked access: {file_path.suffix.lower()}"
                )
        
        except Exception as e:
            logger.error(f"Error getting chunked data: {str(e)}", extra=log_ctx)
            if isinstance(e, HTTPException):
                raise e
            raise HTTPException(status_code=500, detail=f"Failed to get chunked data: {str(e)}")
    
    @staticmethod
    def _get_csv_chunked(file_path: Path, options: Dict[str, Any] = None) -> pd.DataFrame:
        """Get chunked data from CSV file with resolution, limit and region support.
        
        Args:
            file_path: Path to the CSV file
            options: Dict of options (see get_spatial_data_chunked)
            
        Returns:
            DataFrame containing the requested chunk/sample
        """
        options = options or {}
        resolution = float(options.get('resolution', 1.0))
        offset = int(options.get('offset', 0)) if options.get('offset') is not None else None
        limit = int(options.get('limit', 0)) if options.get('limit') is not None else None
        region = options.get('region')
        
        logger.debug(f"Loading CSV chunked: path={file_path}, resolution={resolution}, offset={offset}, limit={limit}")
        
        # Get total line count (only if needed for sampling)
        total_rows = None
        if resolution < 1.0:
            with open(file_path, 'r') as f:
                # Count total lines (subtract 1 for header)
                total_rows = sum(1 for _ in f) - 1
                sample_size = max(1, int(total_rows * resolution))
                logger.debug(f"File has {total_rows} rows, sampling {sample_size} rows at resolution {resolution}")
        
        # Handle different loading strategies based on options
        if resolution >= 1.0 and not region and offset is not None and limit:
            # Simple case: just load the specific chunk with skiprows and nrows
            logger.debug(f"Using direct chunk loading with skiprows={offset+1}, nrows={limit}")
            return pd.read_csv(file_path, skiprows=range(1, offset+1), nrows=limit)
        else:
            # Complex case: need to handle sampling and/or region filtering
            result_df = None
            
            # Handle sampling with or without region filtering
            if resolution < 1.0:
                # Calculate which rows to sample
                assert total_rows is not None, "Total rows not calculated"
                sample_size = max(1, int(total_rows * resolution))
                
                # Random sampling of indices
                import random
                random.seed(42)  # For reproducibility
                
                # Generate random row indices (accounting for header row)
                rows_to_sample = sorted(random.sample(range(1, total_rows+1), sample_size))
                logger.debug(f"Sampling {len(rows_to_sample)} rows")
                
                # Use skiprows with an index list to skip everything except our samples
                # Create a list of all indices not in rows_to_sample
                skip_rows = [i for i in range(total_rows+1) if i != 0 and i not in rows_to_sample]
                result_df = pd.read_csv(file_path, skiprows=skip_rows)
            else:
                # No sampling, but may have region filter
                if offset is not None and limit:
                    result_df = pd.read_csv(file_path, skiprows=range(1, offset+1), nrows=limit)
                else:
                    result_df = pd.read_csv(file_path)

            # Apply region filtering if specified
            if region and result_df is not None:
                x_col = options.get('x_col', 'x')
                y_col = options.get('y_col', 'y')
                
                # Validate the region dict
                try:
                    x_min = float(region['x_min'])
                    y_min = float(region['y_min'])
                    x_max = float(region['x_max'])
                    y_max = float(region['y_max'])
                    
                    logger.debug(f"Filtering by region: x=[{x_min},{x_max}], y=[{y_min},{y_max}]")
                    
                    # Filter by region
                    result_df = result_df[
                        (result_df[x_col] >= x_min) & (result_df[x_col] <= x_max) &
                        (result_df[y_col] >= y_min) & (result_df[y_col] <= y_max)
                    ]
                except (KeyError, ValueError) as e:
                    logger.warning(f"Invalid region specification: {str(e)}")
                    raise HTTPException(status_code=400, detail=f"Invalid region specification: {str(e)}")
                
            # Apply limit after filtering if needed
            if result_df is not None and offset is None and limit and len(result_df) > limit:
                logger.debug(f"Applying post-filter limit: {limit}")
                result_df = result_df.iloc[:limit]
                
            return result_df
    
    @staticmethod
    def _get_h5ad_chunked(file_path: Path, options: Dict[str, Any] = None) -> pd.DataFrame:
        """Get chunked data from H5AD file with resolution, limit and region support.
        
        Args:
            file_path: Path to the H5AD file
            options: Dict of options (see get_spatial_data_chunked)
            
        Returns:
            DataFrame containing the requested chunk/sample
        """
        options = options or {}
        resolution = float(options.get('resolution', 1.0))
        offset = int(options.get('offset', 0)) if options.get('offset') is not None else None
        limit = int(options.get('limit', 0)) if options.get('limit') is not None else None
        region = options.get('region')
        
        logger.debug(f"Loading H5AD chunked: path={file_path}, resolution={resolution}, offset={offset}, limit={limit}")
        
        # Read the H5AD file
        try:
            adata = ad.read_h5ad(file_path)
            
            # Convert to a DataFrame for consistent handling
            logger.debug(f"Loaded H5AD with shape {adata.shape}")
            
            # Extract the data we need into a DataFrame for consistent processing
            result_df = pd.DataFrame(adata.X, columns=adata.var_names)
            
            # Add observation metadata
            for col in adata.obs.columns:
                result_df[col] = adata.obs[col].values
                
            # Apply the same filtering as for CSV
            if resolution < 1.0:
                # Random sampling
                import random
                random.seed(42)  # For reproducibility
                sample_size = max(1, int(len(result_df) * resolution))
                sample_indices = sorted(random.sample(range(len(result_df)), sample_size))
                result_df = result_df.iloc[sample_indices]
                logger.debug(f"Sampled to {len(result_df)} rows at resolution {resolution}")
                
            # Apply offset/limit
            if offset is not None and limit:
                end = min(offset + limit, len(result_df))
                result_df = result_df.iloc[offset:end]
                logger.debug(f"Applied offset/limit: {offset}:{end}, resulting in {len(result_df)} rows")
                
            # Apply region filtering if specified
            if region:
                x_col = options.get('x_col', 'x')
                y_col = options.get('y_col', 'y')
                
                # Validate that columns exist
                if x_col not in result_df.columns or y_col not in result_df.columns:
                    missing = []
                    if x_col not in result_df.columns:
                        missing.append(x_col)
                    if y_col not in result_df.columns:
                        missing.append(y_col)
                    raise HTTPException(
                        status_code=400, 
                        detail=f"Spatial columns not found in H5AD: {', '.join(missing)}"
                    )
                
                # Validate the region dict
                try:
                    x_min = float(region['x_min'])
                    y_min = float(region['y_min'])
                    x_max = float(region['x_max'])
                    y_max = float(region['y_max'])
                    
                    logger.debug(f"Filtering by region: x=[{x_min},{x_max}], y=[{y_min},{y_max}]")
                    
                    # Filter by region
                    result_df = result_df[
                        (result_df[x_col] >= x_min) & (result_df[x_col] <= x_max) &
                        (result_df[y_col] >= y_min) & (result_df[y_col] <= y_max)
                    ]
                except (KeyError, ValueError) as e:
                    logger.warning(f"Invalid region specification: {str(e)}")
                    raise HTTPException(
                        status_code=400, 
                        detail=f"Invalid region specification: {str(e)}"
                    )
                    
            return result_df
            
        except Exception as e:
            logger.error(f"Error reading H5AD file: {str(e)}")
            if isinstance(e, HTTPException):
                raise e
            raise HTTPException(
                status_code=500, 
                detail=f"Failed to read H5AD file: {str(e)}"
            )
            
    @classmethod
    def create_spatial_index(cls, file_id: str, x_col: str = 'x', y_col: str = 'y') -> None:
        """Create a spatial index for the file to enable efficient spatial queries.
        
        Args:
            file_id: The ID of the file to index
            x_col: The column name for x coordinates
            y_col: The column name for y coordinates
            
        Returns:
            None
            
        Raises:
            HTTPException: If file not found or indexing fails
        """
        try:
            # Get the file path
            file_path = cls.get_file_path(file_id)
            
            # Create the index directory if it doesn't exist
            index_dir = cls.TEMP_STORAGE_DIR / "spatial_indices"
            index_dir.mkdir(parents=True, exist_ok=True)
            
            # Index file path
            index_filename = f"{file_id}_spatial.idx"
            index_path = index_dir / index_filename
            
            logger.info(f"Creating spatial index for file {file_id} at {index_path}")
            
            # Load the data
            if file_path.suffix.lower() == '.csv':
                df = pd.read_csv(file_path)
            elif file_path.suffix.lower() == '.h5ad':
                adata = ad.read_h5ad(file_path)
                df = pd.DataFrame()
                df[x_col] = adata.obs[x_col].values if x_col in adata.obs.columns else None
                df[y_col] = adata.obs[y_col].values if y_col in adata.obs.columns else None
            else:
                raise HTTPException(
                    status_code=400, 
                    detail=f"Unsupported file type for spatial indexing: {file_path.suffix.lower()}"
                )
                
            # Check if spatial columns exist
            if x_col not in df.columns or y_col not in df.columns:
                missing = []
                if x_col not in df.columns:
                    missing.append(x_col)
                if y_col not in df.columns:
                    missing.append(y_col)
                raise HTTPException(
                    status_code=400, 
                    detail=f"Spatial columns not found: {', '.join(missing)}"
                )
                
            # Create spatial index using rtree
            from rtree import index
            
            # Function to generate index entries
            def generate_entries():
                for i, (_, row) in enumerate(df.iterrows()):
                    x, y = row[x_col], row[y_col]
                    if pd.notnull(x) and pd.notnull(y):
                        yield (i, (x, y, x, y), i)
            
            # Create and save the index
            idx = index.Index(str(index_path), generate_entries())
            
            # Store index metadata
            import json
            metadata = {
                'file_id': file_id,
                'x_col': x_col,
                'y_col': y_col,
                'row_count': len(df),
                'bounds': {
                    'x_min': float(df[x_col].min()),
                    'y_min': float(df[y_col].min()),
                    'x_max': float(df[x_col].max()),
                    'y_max': float(df[y_col].max()),
                },
                'created_at': datetime.datetime.now().isoformat(),
            }
            
            with open(str(index_path) + '_meta.json', 'w') as f:
                json.dump(metadata, f)
                
            logger.info(f"Spatial index created successfully with {len(df)} points", extra={
                'file_id': file_id,
                'row_count': len(df),
                'bounds': metadata['bounds']
            })
            
            return metadata
            
        except Exception as e:
            logger.error(f"Error creating spatial index: {str(e)}")
            if isinstance(e, HTTPException):
                raise e
            raise HTTPException(
                status_code=500, 
                detail=f"Failed to create spatial index: {str(e)}"
            )
    
    @classmethod
    def query_spatial_index(cls, file_id: str, region: Dict[str, float], 
                           x_col: str = 'x', y_col: str = 'y') -> List[int]:
        """Query the spatial index to get indices of points within a bounding box.
        
        Args:
            file_id: The ID of the file to query
            region: Dict with x_min, y_min, x_max, y_max for spatial filtering
            x_col: The column name for x coordinates
            y_col: The column name for y coordinates
            
        Returns:
            List of row indices that match the spatial query
            
        Raises:
            HTTPException: If index not found or query fails
        """
        try:
            # Index file path
            index_dir = cls.TEMP_STORAGE_DIR / "spatial_indices"
            index_filename = f"{file_id}_spatial.idx"
            index_path = index_dir / index_filename
            
            # Check if index exists
            if not index_path.with_suffix('.dat').exists():
                logger.warning(f"Spatial index not found for file {file_id}")
                
                # Try to create the index
                cls.create_spatial_index(file_id, x_col, y_col)
                
                # Check again
                if not index_path.with_suffix('.dat').exists():
                    raise HTTPException(
                        status_code=404, 
                        detail=f"Spatial index not found for file {file_id}"
                    )
            
            logger.info(f"Querying spatial index for file {file_id}", extra={
                'file_id': file_id,
                'region': region
            })
            
            # Validate the region dict
            try:
                x_min = float(region['x_min'])
                y_min = float(region['y_min'])
                x_max = float(region['x_max'])
                y_max = float(region['y_max'])
            except (KeyError, ValueError) as e:
                logger.warning(f"Invalid region specification: {str(e)}")
                raise HTTPException(
                    status_code=400, 
                    detail=f"Invalid region specification: {str(e)}"
                )
            
            # Open the index
            from rtree import index
            idx = index.Index(str(index_path))
            
            # Query the index
            matching_indices = list(idx.intersection((x_min, y_min, x_max, y_max)))
            
            logger.info(f"Spatial query returned {len(matching_indices)} points", extra={
                'file_id': file_id,
                'region': region,
                'point_count': len(matching_indices)
            })
            
            return matching_indices
            
        except Exception as e:
            logger.error(f"Error querying spatial index: {str(e)}")
            if isinstance(e, HTTPException):
                raise e
            raise HTTPException(
                status_code=500, 
                detail=f"Failed to query spatial index: {str(e)}"
            ) 