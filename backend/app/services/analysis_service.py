import os
from pathlib import Path
import pandas as pd
import anndata as ad
from fastapi import HTTPException
import time # For simulating work
from typing import Dict, Any 
import logging # Import logging

from app.services.file_service import FileService
from app.models.analysis_models import AnalysisPayload, AnalysisMapping # Use the new models file

# Get a logger instance
logger = logging.getLogger(__name__)

# --- Placeholder for importing actual analysis functions ---
# Assume these functions exist and are structured correctly
# You will need to replace these with your actual imports and function calls
try:
    # Attempt to import from the layer_analysis directory
    # Ensure layer_analysis is in Python's path or adjust import structure
    from layer_analysis.analysis_logic import ( 
        run_stage1_counts,
        run_pathway_dominance,
        run_module_context
    )
    print("Successfully imported analysis functions from layer_analysis")
except ImportError as e:
    logger.warning(f"Could not import from layer_analysis: {e}. Using placeholder functions.")
    # Define placeholder functions if import fails, so the service still runs
    def run_stage1_counts(spatial_df, interactions_df):
        logger.warning("Using placeholder run_stage1_counts")
        time.sleep(2)
        # Simulate counts for whole tissue and layers
        layers = list(spatial_df['layer'].unique())
        results = {
            'whole_tissue': {'unique_ligands': 100, 'unique_receptors': 120}
        }
        for layer in layers:
             results[layer] = {'unique_ligands': 50, 'unique_receptors': 60}
        return results
        
    def run_pathway_dominance(spatial_df, interactions_df):
        logger.warning("Using placeholder run_pathway_dominance")
        time.sleep(3)
        # Simulate pathway results for whole tissue and layers
        layers = list(spatial_df['layer'].unique())
        results = {
            'whole_tissue': pd.DataFrame([{'ligand': 'L1', 'receptor': 'R1', 'score': 0.95}, {'ligand': 'L2', 'receptor': 'R2', 'score': 0.91}])
        }
        for layer in layers:
             results[layer] = pd.DataFrame([{'ligand': 'L1', 'receptor': 'R1', 'score': 0.88}, {'ligand': 'L3', 'receptor': 'R3', 'score': 0.85}])
        # Convert DataFrames to dicts for JSON serialization
        return {scope: df.to_dict('records') for scope, df in results.items()}
        
    def run_module_context(spatial_df, interactions_df, modules_df, pathway_results):
        logger.warning("Using placeholder run_module_context")
        time.sleep(3)
        # Simulate module context results for whole tissue and layers
        layers = list(spatial_df['layer'].unique())
        results = {
            'whole_tissue': pd.DataFrame([{'ligand': 'L1', 'receptor': 'R1', 'type': 'intra-module'}, {'ligand': 'L2', 'receptor': 'R2', 'type': 'inter-module'}])
        }
        for layer in layers:
            results[layer] = pd.DataFrame([{'ligand': 'L1', 'receptor': 'R1', 'type': 'intra-module'}, {'ligand': 'L3', 'receptor': 'R3', 'type': 'inter-module'}])
        # Convert DataFrames to dicts for JSON serialization
        return {scope: df.to_dict('records') for scope, df in results.items()}
# ---------------------------------------------------------

# --- Job Status Store (In-Memory) ---
# Simple dictionary to store job status. 
# NOTE: This is in-memory and will be lost on server restart.
# For persistent jobs, use a database or a proper task queue backend (like Redis with Celery).
jobs_status: Dict[str, Dict[str, Any]] = {}
# -------------------------------------

# Define standardized column names expected by the analysis functions
STANDARDIZED_SPATIAL_COLS = {
    'geneCol': 'gene',
    'xCol': 'x',
    'yCol': 'y',
    'layerCol': 'layer'
}
STANDARDIZED_INTERACTION_COLS = {
    'ligandCol': 'ligand',
    'receptorCol': 'receptor'
}
STANDARDIZED_MODULE_COLS = {
    'geneCol': 'gene',
    'moduleCol': 'module'
}

class AnalysisService:
    """Service to handle the analysis pipeline logic."""

    def __init__(self, connection_manager: Any):
        """Initialize the service with a WebSocket connection manager."""
        self.manager = connection_manager
        logger.info("AnalysisService initialized with ConnectionManager.")

    async def _update_job_status(self, job_id: str, status_update: Dict[str, Any]):
        """Helper function to update in-memory status and notify WebSocket clients."""
        if job_id in jobs_status:
            jobs_status[job_id].update(status_update)
            # Add job_id to the message being sent
            message_to_send = {**jobs_status[job_id], "job_id": job_id}
            logger.info(f"Updating job {job_id} status: {status_update}. Sending via WebSocket.")
            await self.manager.send_update(job_id, message_to_send)
        else:
            logger.warning(f"Attempted to update status for non-existent job ID: {job_id}")

    async def run_analysis_background(self, job_id: str, payload: AnalysisPayload):
        """The actual analysis function intended to be run in the background.
        
        Loads data, runs analysis stages, and updates job status via WebSocket.
        Args:
            job_id: The unique identifier for this analysis job.
            payload: The analysis request payload containing file IDs and mappings.
        """
        logger.info(f"Background task started for job ID: {job_id}")
        final_results = {}
        try:
            # 0. Update status to running
            await self._update_job_status(job_id, {"status": "running", "message": "Loading data...", "progress": 0.1})
            
            # 1. Get file paths and load/standardize data
            spatial_path = FileService.get_file_path(payload.spatialFileId)
            interactions_path = FileService.get_file_path(payload.interactionsFileId)
            modules_path = FileService.get_file_path(payload.modulesFileId)
            
            await self._update_job_status(job_id, {"message": "Standardizing data...", "progress": 0.2})
            spatial_df = self._load_and_standardize(spatial_path, payload.spatialMapping, STANDARDIZED_SPATIAL_COLS)
            interactions_df = self._load_and_standardize(interactions_path, payload.interactionsMapping, STANDARDIZED_INTERACTION_COLS)
            modules_df = self._load_and_standardize(modules_path, payload.modulesMapping, STANDARDIZED_MODULE_COLS)
            logger.info(f"Job {job_id}: Data loaded successfully.")

            # 2. Run Stage 1: Initial Counts
            await self._update_job_status(job_id, {"message": "Calculating ligand/receptor counts...", "progress": 0.4})
            logger.info(f"Job {job_id}: Running Stage 1 Counts...")
            count_results = run_stage1_counts(spatial_df, interactions_df)
            # Structure results: {scope: {analysis_type: data}}
            for scope, counts in count_results.items():
                if scope not in final_results: final_results[scope] = {}
                final_results[scope]['ligand_receptor_counts'] = counts
            logger.info(f"Job {job_id}: Stage 1 Counts Complete.")
            await self._update_job_status(job_id, {"progress": 0.5})

            # 3. Run Stage 2 (Concurrent - simulated sequentially here)
            
            # 3a. Pathway Dominance
            await self._update_job_status(job_id, {"message": "Calculating pathway dominance...", "progress": 0.6})
            logger.info(f"Job {job_id}: Running Pathway Dominance...")
            pathway_results = run_pathway_dominance(spatial_df, interactions_df)
            for scope, pathways in pathway_results.items():
                 if scope not in final_results: final_results[scope] = {}
                 final_results[scope]['pathway_dominance'] = pathways # pathways assumed to be list of dicts
            logger.info(f"Job {job_id}: Pathway Dominance Complete.")
            await self._update_job_status(job_id, {"progress": 0.8})
            
            # 3b. Module Context
            await self._update_job_status(job_id, {"message": "Calculating module context...", "progress": 0.85})
            logger.info(f"Job {job_id}: Running Module Context...")
            # Module context might depend on pathway results (e.g., significant pairs)
            # Pass pathway_results if needed by the actual function
            module_context_results = run_module_context(spatial_df, interactions_df, modules_df, pathway_results)
            for scope, modules in module_context_results.items():
                 if scope not in final_results: final_results[scope] = {}
                 final_results[scope]['module_context'] = modules # modules assumed to be list of dicts
            logger.info(f"Job {job_id}: Module Context Complete.")

            # 4. Update status to success
            await self._update_job_status(job_id, {
                "status": "success",
                "message": "Analysis completed successfully.",
                "progress": 1.0,
                "results": final_results
            })
            logger.info(f"Job {job_id}: Analysis complete. Final status sent via WebSocket.")

        except HTTPException as e:
            # Handle exceptions raised during loading/standardization or analysis
            error_message = f"Error during analysis: {e.detail}"
            logger.error(f"Job {job_id}: Failed with HTTPException - {error_message}")
            await self._update_job_status(job_id, {"status": "failed", "message": error_message})
        except Exception as e:
            # Catch any other unexpected errors
            error_message = f"An unexpected error occurred during analysis: {str(e)}"
            logger.exception(f"Job {job_id}: Failed with Exception - {error_message}")
            await self._update_job_status(job_id, {"status": "failed", "message": error_message})
        finally:
            # --- Cleanup uploaded files --- 
            logger.info(f"Job {job_id}: Cleaning up temporary files...")
            files_to_delete = []
            if hasattr(payload, 'spatialFileId') and payload.spatialFileId:
                try:
                    files_to_delete.append(FileService.get_file_path(payload.spatialFileId))
                except HTTPException as e:
                     logger.warning(f"Job {job_id}: Could not find spatial file {payload.spatialFileId} for cleanup: {e.detail}")
            if hasattr(payload, 'interactionsFileId') and payload.interactionsFileId:
                try:
                    files_to_delete.append(FileService.get_file_path(payload.interactionsFileId))
                except HTTPException as e:
                     logger.warning(f"Job {job_id}: Could not find interactions file {payload.interactionsFileId} for cleanup: {e.detail}")
            if hasattr(payload, 'modulesFileId') and payload.modulesFileId:
                 try:
                    files_to_delete.append(FileService.get_file_path(payload.modulesFileId))
                 except HTTPException as e:
                     logger.warning(f"Job {job_id}: Could not find modules file {payload.modulesFileId} for cleanup: {e.detail}")
            
            for file_path in files_to_delete:
                try:
                    if file_path.exists():
                        file_path.unlink()
                        logger.info(f"Job {job_id}: Deleted temporary file: {file_path}")
                except OSError as e:
                    logger.error(f"Job {job_id}: Error deleting temporary file {file_path}: {e}")
            # -----------------------------

    def _load_and_standardize(self, file_path: Path, mapping: AnalysisMapping, standard_cols_map: Dict[str, str]) -> pd.DataFrame:
        """Loads data from CSV or H5AD and standardizes columns based on mapping."""
        file_ext = file_path.suffix.lower()
        user_col_map = {k: v for k, v in mapping.dict().items() if v is not None} # User's selection {frontendKey: userColumnName}
        required_user_cols = list(user_col_map.values())
        rename_map = {user_col_name: standard_cols_map[frontend_key] 
                      for frontend_key, user_col_name in user_col_map.items()}

        if file_ext == '.csv':
            try:
                df = pd.read_csv(file_path, usecols=required_user_cols)
                df.rename(columns=rename_map, inplace=True)
                return df
            except ValueError as e:
                 # More specific error for missing columns
                 raise HTTPException(status_code=400, detail=f"Error loading CSV {file_path.name}: Missing required columns ({e}). Selected: {required_user_cols}")
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"Error reading CSV file {file_path.name}: {str(e)}")

        elif file_ext == '.h5ad':
            try:
                adata = ad.read_h5ad(file_path)
                
                # --- H5AD Data Extraction Logic --- 
                # This needs careful handling based on where data might be (obs, var, obsm, index)
                data_dict = {}
                
                for frontend_key, user_col_name in user_col_map.items():
                    standard_col_name = standard_cols_map[frontend_key]
                    
                    # Case 1: Data in obs
                    if user_col_name in adata.obs.columns:
                        data_dict[standard_col_name] = adata.obs[user_col_name].values
                    # Case 2: Data is the index (usually gene or cell ID)
                    elif user_col_name == 'index' or user_col_name == adata.obs.index.name:
                         data_dict[standard_col_name] = adata.obs.index.values
                    elif user_col_name == adata.var.index.name: # Check var index too (e.g., for gene names)
                         # This usually doesn't align row-wise with obs, handle carefully. 
                         # For modules/interactions linking gene names, this might be okay if read separately.
                         # If spatial data needs gene names, it's more complex.
                         # Let's assume for modules/interactions, mapping 'index' or 'var_names' is intended.
                         if standard_col_name == 'gene': # Special case for gene names
                            # We might need the full var index if we're just mapping genes to modules
                            # This part needs refinement based on how modules/interactions data uses gene names
                            logger.warning(f"Mapping '{user_col_name}' (var index) directly. Ensure this aligns with expected data structure.")
                            # For simplicity now, let's assume it's okay for non-spatial gene mapping
                            if frontend_key in ['geneCol']: # Check if it's for modules/interactions gene mapping 
                                data_dict[standard_col_name] = adata.var.index.values
                            else:
                                raise ValueError(f"Cannot directly map var index '{user_col_name}' to '{standard_col_name}' in this context.")
                         else:
                            raise ValueError(f"Mapping '{user_col_name}' (var index) to '{standard_col_name}' is not supported yet.")
                            
                    # Case 3: Data in obsm (e.g., spatial coordinates)
                    elif '.' in user_col_name: # Heuristic: check for format like 'obsm_key.column_index_or_name'
                        try:
                            obsm_key, col_ref = user_col_name.split('.', 1)
                            if obsm_key in adata.obsm:
                                obsm_data = adata.obsm[obsm_key]
                                try:
                                    col_idx = int(col_ref) # Try interpreting as index
                                    if col_idx < obsm_data.shape[1]:
                                        data_dict[standard_col_name] = obsm_data[:, col_idx]
                                    else:
                                        raise ValueError(f"Column index {col_idx} out of bounds for obsm key '{obsm_key}'")
                                except ValueError: 
                                    # Try interpreting as column name if obsm is a DataFrame (less common)
                                    if isinstance(obsm_data, pd.DataFrame) and col_ref in obsm_data.columns:
                                         data_dict[standard_col_name] = obsm_data[col_ref].values
                                    else:
                                         # Add heuristic for common spatial keys like 'spatial'
                                         if obsm_key == 'spatial':
                                             if col_ref.lower() == 'x' and obsm_data.shape[1] > 0:
                                                 data_dict[standard_col_name] = obsm_data[:, 0]
                                             elif col_ref.lower() == 'y' and obsm_data.shape[1] > 1:
                                                  data_dict[standard_col_name] = obsm_data[:, 1]
                                             else:
                                                  raise ValueError(f"Cannot find column reference '{col_ref}' in obsm key '{obsm_key}'")
                                         else:
                                             raise ValueError(f"Cannot find column reference '{col_ref}' in obsm key '{obsm_key}'")
                            else:
                                raise ValueError(f"obsm key '{obsm_key}' not found in H5AD file.")
                        except Exception as e:
                            raise ValueError(f"Error parsing obsm mapping '{user_col_name}': {e}")
                    else:
                        raise ValueError(f"Column '{user_col_name}' not found in H5AD .obs, index, or supported .obsm formats.")
                
                # Check if all required data was extracted
                extracted_cols = list(data_dict.keys())
                required_standard_cols = list(standard_cols_map.values())
                if not all(col in extracted_cols for col in required_standard_cols):
                    missing = [col for col in required_standard_cols if col not in extracted_cols]
                    raise ValueError(f"Could not extract all required data from H5AD. Missing: {missing}")
                    
                # Ensure all arrays have the same length (relevant for obs/obsm based data)
                lengths = {len(v) for k, v in data_dict.items() if k != 'gene' or frontend_key not in ['geneCol']} # Exclude var-based gene list length check
                if len(lengths) > 1:
                    raise ValueError(f"Extracted columns from H5AD have inconsistent lengths: {lengths}")

                df = pd.DataFrame(data_dict)
                return df
                # --- End H5AD Extraction --- 
                
            except ValueError as e:
                raise HTTPException(status_code=400, detail=f"Error mapping H5AD {file_path.name}: {str(e)}")
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"Error reading H5AD file {file_path.name}: {str(e)}")

        else:
            raise HTTPException(status_code=400, detail=f"Unsupported file type for analysis: {file_ext}")

