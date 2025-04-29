import os
from pathlib import Path
import pandas as pd
import anndata as ad
from fastapi import HTTPException, Depends
import time # For simulating work
from typing import Dict, Any
import logging # Import logging

# Import the new JobService and its dependency function
from app.services.job_service import JobService, get_job_service
from app.services.file_service import FileService
from app.models.analysis_models import AnalysisPayload, AnalysisMapping # Use the new models file

# Get a logger instance
logger = logging.getLogger(__name__)

# --- Placeholder for importing actual analysis functions ---
# Assume these functions exist and are structured correctly
# You will need to replace these with your actual imports and function calls
try:
    # Attempt to import from the correct location within the app
    from app.analysis_logic.core import ( 
        run_stage1_counts_pipeline,
        run_pathway_dominance_pipeline,
        run_module_context_pipeline
    )
    logger.info("Successfully imported analysis functions from app.analysis_logic.core")
except ImportError as e:
    logger.warning(f"Could not import from app.analysis_logic.core: {e}. Using placeholder functions.")
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
            'whole_tissue': pd.DataFrame([{'ligand': 'L1', 'receptor': 'R1', 'interaction_type': 'intra-module'}, {'ligand': 'L2', 'receptor': 'R2', 'interaction_type': 'inter-module'}])
        }
        for layer in layers:
            results[layer] = pd.DataFrame([{'ligand': 'L1', 'receptor': 'R1', 'interaction_type': 'intra-module'}, {'ligand': 'L3', 'receptor': 'R3', 'interaction_type': 'inter-module'}])
        # Convert DataFrames to dicts for JSON serialization
        return {scope: df.to_dict('records') for scope, df in results.items()}
# ---------------------------------------------------------

# --- Job Status and Context Stores REMOVED ---
# These are now managed by JobService
# jobs_status: Dict[str, Dict[str, Any]] = {}
# job_contexts: Dict[str, Any] = {}
# ---------------------------------------------

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

    # Inject JobService using FastAPI's Depends
    # Make connection_manager optional
    def __init__(self, job_service: JobService = Depends(get_job_service), connection_manager: Any = None):
        """Initialize the service with JobService and optional WebSocket connection manager."""
        self.manager = connection_manager
        self.job_service = job_service # Store the injected JobService instance
        logger.info(f"AnalysisService initialized. Manager: {'Present' if self.manager else 'Absent'}, JobService: Present")

    # REMOVE _update_job_status - use self.job_service.update_job_status directly
    # async def _update_job_status(self, job_id: str, status_update: Dict[str, Any]):
    #     """Helper function to update in-memory status and notify WebSocket clients."""
    #     # ... old implementation ...

    async def run_analysis_background(self, job_id: str, payload: AnalysisPayload):
        """The actual analysis function intended to be run in the background.

        Loads data, runs analysis stages, and updates job status via JobService and WebSocket.
        Args:
            job_id: The unique identifier for this analysis job (already created).
            payload: The analysis request payload containing file IDs and mappings (already stored in context).
        """
        # Payload should already be in context via job_service.store_job_context in start_analysis
        logger.info(f"Background task started for job ID: {job_id}")
        final_results = {}
        try:
            # 0. Update status to running via JobService
            status_update = {"status": "running", "message": "Loading data...", "progress": 0.1}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})

            # 1. Get file paths and load/standardize data
            # Retrieve payload from context
            job_context = self.job_service.get_job_context(job_id)
            if not job_context:
                 raise Exception(f"Job context not found for job ID {job_id}")
            # Recreate payload object if needed, or directly use dict keys
            # Assuming job_context IS the payload dict here
            spatial_path = FileService.get_file_path(job_context['spatialFileId'])
            interactions_path = FileService.get_file_path(job_context['interactionsFileId'])
            modules_path = FileService.get_file_path(job_context['modulesFileId'])
            spatial_mapping = AnalysisMapping(**job_context['spatialMapping'])
            interactions_mapping = AnalysisMapping(**job_context['interactionsMapping'])
            modules_mapping = AnalysisMapping(**job_context['modulesMapping'])

            status_update = {"message": "Standardizing data...", "progress": 0.2}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
            
            spatial_df = self._load_and_standardize(spatial_path, spatial_mapping, STANDARDIZED_SPATIAL_COLS)
            interactions_df = self._load_and_standardize(interactions_path, interactions_mapping, STANDARDIZED_INTERACTION_COLS)
            modules_df = self._load_and_standardize(modules_path, modules_mapping, STANDARDIZED_MODULE_COLS)
            logger.info(f"Job {job_id}: Data loaded successfully.")

            # 2. Run Stage 1: Initial Counts
            status_update = {"message": "Calculating ligand/receptor counts...", "progress": 0.4}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
            
            logger.info(f"Job {job_id}: Running Stage 1 Counts...")
            count_results = run_stage1_counts_pipeline(spatial_df, interactions_df)
            # Structure results: {scope: {analysis_type: data}}
            for scope, counts in count_results.items():
                if scope not in final_results: final_results[scope] = {}
                final_results[scope]['ligand_receptor_counts'] = counts
            logger.info(f"Job {job_id}: Stage 1 Counts Complete.")
            
            status_update = {"progress": 0.5}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})

            # 3. Run Stage 2 (Concurrent - simulated sequentially here)

            # 3a. Pathway Dominance
            status_update = {"message": "Calculating pathway dominance...", "progress": 0.6}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
            
            logger.info(f"Job {job_id}: Running Pathway Dominance...")
            pathway_results = run_pathway_dominance_pipeline(spatial_df, interactions_df)
            for scope, pathways in pathway_results.items():
                 if scope not in final_results: final_results[scope] = {}
                 final_results[scope]['pathway_dominance'] = pathways # pathways assumed to be list of dicts
            logger.info(f"Job {job_id}: Pathway Dominance Complete.")
            
            status_update = {"progress": 0.8}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})

            # 3b. Module Context
            status_update = {"message": "Calculating module context...", "progress": 0.85}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
            
            logger.info(f"Job {job_id}: Running Module Context...")
            # Module context might depend on pathway results (e.g., significant pairs)
            # Pass pathway_results if needed by the actual function
            module_context_results = run_module_context_pipeline(spatial_df, interactions_df, modules_df, pathway_results)
            for scope, modules in module_context_results.items():
                 if scope not in final_results: final_results[scope] = {}
                 final_results[scope]['module_context'] = modules # modules assumed to be list of dicts
            logger.info(f"Job {job_id}: Module Context Complete.")

            # 4. Update status to success via JobService
            status_update = {
                "status": "success",
                "message": "Analysis completed successfully.",
                "progress": 1.0,
                "results": final_results
            }
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
            logger.info(f"Job {job_id}: Analysis complete. Final status sent via WebSocket.")

        except HTTPException as e:
            # Handle exceptions raised during loading/standardization or analysis
            error_message = f"Error during analysis: {e.detail}"
            logger.error(f"Job {job_id}: Failed with HTTPException - {error_message}")
            status_update = {"status": "failed", "message": error_message}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
        except Exception as e:
            # Catch any other unexpected errors
            error_message = f"An unexpected error occurred during analysis: {str(e)}"
            logger.exception(f"Job {job_id}: Failed with Exception - {error_message}")
            status_update = {"status": "failed", "message": error_message}
            self.job_service.update_job_status(job_id, **status_update)
            await self.manager.send_update(job_id, {**self.job_service.get_job_status(job_id), "job_id": job_id})
        finally:
            # Cleanup logic remains the same for now
            # TODO: Implement delayed or on-demand cleanup
            pass

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
                                # If not mapping gene names specifically, unclear how to handle var index
                                raise HTTPException(status_code=400, detail=f"Cannot map variable index '{user_col_name}' to '{standard_col_name}' in H5AD {file_path.name}. Expected column in .obs.")
                         else:
                            # If not mapping gene names specifically, unclear how to handle var index
                            raise HTTPException(status_code=400, detail=f"Cannot map variable index '{user_col_name}' to '{standard_col_name}' in H5AD {file_path.name}. Expected column in .obs.")
                    # Case 3: Data in obsm (e.g., spatial coordinates)
                    elif user_col_name in adata.obsm:
                        # Assuming obsm[user_col_name] is structured appropriately (e.g., 2D array for coords)
                        # Need to handle potential multi-dimensional data
                        # Example for spatial coordinates often in 'spatial' key as array
                         if standard_col_name in ['x', 'y']:
                             if 'spatial' in adata.obsm and isinstance(adata.obsm['spatial'], pd.DataFrame):
                                coords_df = adata.obsm['spatial']
                                if 'x' in coords_df.columns and 'y' in coords_df.columns:
                                    data_dict['x'] = coords_df['x'].values
                                    data_dict['y'] = coords_df['y'].values
                                else: 
                                    # Attempt default coord names if 'x', 'y' not present
                                    # Assuming first two columns are x, y
                                    if coords_df.shape[1] >= 2:
                                        data_dict['x'] = coords_df.iloc[:, 0].values
                                        data_dict['y'] = coords_df.iloc[:, 1].values
                                    else:
                                         raise HTTPException(status_code=400, detail=f"Could not find x, y coordinates in .obsm['spatial'] of H5AD {file_path.name}.")
                             elif 'spatial' in adata.obsm and isinstance(adata.obsm['spatial'], np.ndarray):
                                 coords_array = adata.obsm['spatial']
                                 if coords_array.shape[1] >= 2:
                                     data_dict['x'] = coords_array[:, 0]
                                     data_dict['y'] = coords_array[:, 1]
                                 else:
                                     raise HTTPException(status_code=400, detail=f"Could not find x, y coordinates in .obsm['spatial'] numpy array of H5AD {file_path.name}.")
                             else:
                                # If user specified a different obsm key containing coords, 
                                # requires more specific handling or user input for column indices/names
                                 raise HTTPException(status_code=400, detail=f"Spatial coordinates key '{user_col_name}' not found or incorrect format in .obsm of H5AD {file_path.name}. Expected 'spatial' key with DataFrame or NumPy array.")
                         else:
                             # Handle other potential obsm data mappings if needed
                             raise HTTPException(status_code=400, detail=f"Mapping from .obsm key '{user_col_name}' to '{standard_col_name}' not implemented yet for H5AD {file_path.name}.")

                    # Case 4: Data in var (less common for row-wise alignment with obs)
                    elif user_col_name in adata.var.columns:
                         raise HTTPException(status_code=400, detail=f"Cannot map from .var column '{user_col_name}' in H5AD {file_path.name}. Standard columns should generally map from .obs or index.")
                    
                    else:
                        raise HTTPException(status_code=400, detail=f"Specified column/key '{user_col_name}' not found in H5AD file {file_path.name} (.obs, .obs.index, .var.index, .obsm['spatial'] checked). Available obs columns: {list(adata.obs.columns)}")

                # Check if all required standard columns were successfully populated
                missing_standard_cols = set(standard_cols_map.values()) - set(data_dict.keys())
                if missing_standard_cols:
                    # This might indicate an issue with the loop logic or mapping assumptions
                    raise HTTPException(status_code=500, detail=f"Internal error: Failed to extract required standard columns ({missing_standard_cols}) from H5AD {file_path.name} after mapping.")
                
                # Ensure all extracted arrays/series have the same length (number of observations)
                lengths = [len(v) for v in data_dict.values()]
                if len(set(lengths)) > 1:
                    # This could happen if mapping var index with obs data - needs careful thought
                     raise HTTPException(status_code=400, detail=f"Inconsistent data lengths extracted from H5AD {file_path.name}. Check mappings. Lengths found: {lengths}")
                
                df = pd.DataFrame(data_dict)
                return df

            except FileNotFoundError:
                 raise HTTPException(status_code=404, detail=f"H5AD file not found: {file_path.name}")
            except Exception as e:
                logger.exception(f"Error reading or processing H5AD file {file_path.name}: {str(e)}")
                raise HTTPException(status_code=500, detail=f"Error reading or processing H5AD file {file_path.name}: {str(e)}")

        else:
            raise HTTPException(status_code=400, detail=f"Unsupported file type: {file_ext}")
            
    # REMOVED: start_analysis - This logic moves to the API route
    # def start_analysis(self, payload: AnalysisPayload, background_tasks: BackgroundTasks) -> Dict[str, Any]:
    #     """Initiates the analysis in a background task."""
    #     # ... old implementation ...

    # REMOVED: get_job_status - Use job_service.get_job_status directly from API/WebSocket handler
    # def get_job_status(self, job_id: str) -> Dict[str, Any]:
    #     """Retrieve the status of a specific analysis job."""
    #     # ... old implementation ...
    
    # REMOVED: store_job_context - Functionality moved to JobService
    # @staticmethod
    # def store_job_context(job_id: str, payload: AnalysisPayload):
    #    """Stores the job context (payload) for later use."""
    #    # ... old implementation ...
    
    # REMOVED: get_job_context - Use job_service.get_job_context directly from API endpoints
    # @staticmethod
    # def get_job_context(job_id: str) -> dict:
    #     """Retrieves the job context using the job ID."""
    #     # ... old implementation ...

    # --- MOVED the dependency function OUTSIDE the class --- 
    # def get_analysis_service_dependency(job_service: JobService = Depends(get_job_service)) -> "AnalysisService":
    #    # ...

# --- Dependency function MOVED HERE (outside the class) ---
def get_analysis_service_dependency(job_service: JobService = Depends(get_job_service)) -> "AnalysisService":
    """FastAPI dependency provider for AnalysisService (without ConnectionManager)."""
    # Note: This returns an instance without the manager, suitable for non-WebSocket endpoints
    # It also needs to resolve the forward reference to AnalysisService correctly
    return AnalysisService(job_service=job_service)

