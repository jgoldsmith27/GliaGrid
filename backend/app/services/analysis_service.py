import os
from pathlib import Path
import pandas as pd
import anndata as ad
from fastapi import HTTPException, Depends
import time # For simulating work
from typing import Dict, Any
import logging # Import logging
import asyncio
from scipy.spatial import ConvexHull # ADDED: Import ConvexHull
import numpy as np # ADDED: Import numpy for array handling

# Import the new JobService and its dependency function
from app.services.job_service import JobService, get_job_service
from app.services.file_service import FileService
# Import models from the models file
from app.models.analysis_models import AnalysisPayload, AnalysisMapping, CustomAnalysisRequest, PointData

# REMOVED OLD IMPORT: from ..api.analysis_routes import CustomAnalysisRequest, PointData 
# Assume settings are available for standard file paths
from ..core.config import settings 

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
    # Remove connection_manager parameter
    def __init__(self, job_service: JobService = Depends(get_job_service)):
        """Initialize the service with JobService."""
        self.job_service = job_service # Store the injected JobService instance
        logger.info(f"AnalysisService initialized with JobService")

    # Remove _update_job_status - use self.job_service.update_job_status directly

    async def run_analysis_background(self, job_id: str, payload: AnalysisPayload):
        logger.info(f"Background task started for job ID: {job_id}")
        final_results = {}
        try:
            # 0. Update status to running
            await self.job_service.update_job_status(job_id, status="running", message="Loading data...", progress=0.1)
            logger.info(f"[Job {job_id}] Status: Running, Progress: 10% (Loading data...)")

            # 1. Get file paths and mappings directly from the payload argument
            # REMOVED: job_context = self.job_service.get_job_context(job_id)
            # REMOVED: if not job_context: ...
            
            # Use payload attributes directly
            spatial_path = FileService.get_file_path(payload.spatialFileId)
            interactions_path = FileService.get_file_path(payload.interactionsFileId)
            modules_path = FileService.get_file_path(payload.modulesFileId)
            # Mappings are already AnalysisMapping objects within the payload
            spatial_mapping = payload.spatialMapping
            interactions_mapping = payload.interactionsMapping
            modules_mapping = payload.modulesMapping

            await self.job_service.update_job_status(job_id, message="Standardizing data...", progress=0.2)
            logger.info(f"[Job {job_id}] Status: Running, Progress: 20% (Standardizing data...)")
            
            spatial_df = self._load_and_standardize(spatial_path, spatial_mapping, STANDARDIZED_SPATIAL_COLS)
            interactions_df = self._load_and_standardize(interactions_path, interactions_mapping, STANDARDIZED_INTERACTION_COLS)
            modules_df = self._load_and_standardize(modules_path, modules_mapping, STANDARDIZED_MODULE_COLS)
            logger.info(f"Job {job_id}: Data loaded successfully.")

            # ADDED: Calculate Layer Boundaries (Convex Hulls)
            layer_boundaries = {}
            if 'layer' in spatial_df.columns and 'x' in spatial_df.columns and 'y' in spatial_df.columns:
                logger.info(f"[Job {job_id}] Calculating layer boundaries...")
                try:
                    # Ensure layer names are strings for consistent dictionary keys
                    spatial_df['layer'] = spatial_df['layer'].astype(str) 
                    grouped = spatial_df.groupby('layer')
                    for name, group in grouped:
                        # Need at least 3 points to form a convex hull
                        if len(group) >= 3:
                            try:
                                # Extract points as numpy array
                                points = group[['x', 'y']].to_numpy()
                                hull = ConvexHull(points)
                                # Get vertices coordinates and convert to list of [x, y] pairs
                                boundary_coords = points[hull.vertices].tolist()
                                # Ensure the polygon is closed by adding the first point at the end if it's not already there
                                if boundary_coords[0] != boundary_coords[-1]:
                                    boundary_coords.append(boundary_coords[0])
                                layer_boundaries[name] = boundary_coords
                            except Exception as hull_error:
                                logger.warning(f"    Could not compute convex hull for layer '{name}': {hull_error}")
                                layer_boundaries[name] = None # Indicate failure for this layer
                        else:
                            logger.warning(f"    Skipping boundary calculation for layer '{name}': requires at least 3 points, found {len(group)}.\n")
                            layer_boundaries[name] = None # Indicate insufficient points
                except Exception as e:
                    logger.exception(f"[Job {job_id}] Error during layer boundary calculation: {e}")
                    logger.info(f"[Job {job_id}] Layer boundary calculation failed after {time.time() - boundary_calc_start_time:.2f} seconds.")
                    # Continue analysis even if boundary calculation fails
            else:
                logger.warning(f"[Job {job_id}] Skipping layer boundary calculation: Missing 'layer', 'x', or 'y' columns in spatial data.\n")
            # --- END Boundary Calculation ---

            # 2. Run Stage 1: Initial Counts (Sequential)
            await self.job_service.update_job_status(job_id, message="Calculating ligand/receptor counts...", progress=0.4)
            logger.info(f"[Job {job_id}] Status: Running, Progress: 40% (Calculating counts...)")
            
            logger.info(f"Job {job_id}: Running Stage 1 Counts...")
            count_results = run_stage1_counts_pipeline(spatial_df, interactions_df)
            for scope, counts in count_results.items():
                scope_str = str(scope)
                if scope_str not in final_results: final_results[scope_str] = {}
                final_results[scope_str]['ligand_receptor_counts'] = counts
            logger.info(f"Job {job_id}: Stage 1 Counts Complete.")
            
            await self.job_service.update_job_status(job_id, progress=0.5)
            logger.info(f"[Job {job_id}] Status: Running, Progress: 50%")

            # 3. Run Stage 2a: Pathway Dominance (Sequential)
            await self.job_service.update_job_status(job_id, message="Calculating pathway dominance...", progress=0.6)
            logger.info(f"[Job {job_id}] Status: Running, Progress: 60% (Calculating pathway dominance...)")
            
            logger.info(f"Job {job_id}: Running Pathway Dominance...")
            # This now returns scores only, no points
            pathway_results = run_pathway_dominance_pipeline(spatial_df, interactions_df)
            
            # Store the pathway dominance results directly
            for scope, pathways in pathway_results.items():
                 scope_str = str(scope)
                 if scope_str not in final_results: final_results[scope_str] = {}
                 final_results[scope_str]['pathway_dominance'] = pathways 
            
            logger.info(f"Job {job_id}: Pathway Dominance Complete.")
            await self.job_service.update_job_status(job_id, progress=0.8)
            logger.info(f"[Job {job_id}] Status: Running, Progress: 80%")

            # 4. Run Stage 2b: Module Context (Sequential)
            await self.job_service.update_job_status(job_id, message="Calculating module context...", progress=0.85)
            logger.info(f"Job {job_id}: Running Module Context...")
            # Pass the results dict containing pathway_dominance per scope
            module_context_results = run_module_context_pipeline(spatial_df, interactions_df, modules_df, final_results)
            for scope, modules in module_context_results.items():
                 scope_str = str(scope)
                 if scope_str not in final_results: final_results[scope_str] = {}
                 final_results[scope_str]['module_context'] = modules
            logger.info(f"Job {job_id}: Module Context Complete.")

            # --- DEBUG: Log final_results before structuring ---
            # logger.debug(f"[Job {job_id}] Final analysis outputs before structuring: {final_results}")
            # -------------------------------------------------

            # ADDED: Include boundaries in final results
            if layer_boundaries:
                final_results['layer_boundaries'] = layer_boundaries

            # 5. Update status to success
            # Structure results to include inputs (payload) and outputs (final_results)
            structured_results = {
                "inputs": {
                    "files": {
                        "spatialFileId": payload.spatialFileId,
                        "interactionsFileId": payload.interactionsFileId,
                        "modulesFileId": payload.modulesFileId
                    },
                    "mappings": {
                        "spatialMapping": payload.spatialMapping.dict(), # Convert Pydantic model to dict
                        "interactionsMapping": payload.interactionsMapping.dict(),
                        "modulesMapping": payload.modulesMapping.dict()
                    }
                },
                "outputs": final_results # Contains results grouped by scope
            }
            
            # --- DEBUG: Log the final structured results to be saved ---
            # logger.debug(f"[Job {job_id}] Structured results being saved: {structured_results}")
            # ---------------------------------------------------------

            status_update = {
                "status": "success",
                "message": "Analysis completed successfully.",
                "progress": 1.0,
                "results": structured_results 
            }
            await self.job_service.update_job_status(job_id, **status_update)
            logger.info(f"[Job {job_id}] Status: Success, Progress: 100% (Analysis complete)")

        except Exception as e: # Catch exceptions from any stage
            error_message = f"Analysis failed: {str(e)}"
            logger.exception(f"Job {job_id}: Failed - {error_message}")
            status_update = {"status": "failed", "message": error_message}
            await self.job_service.update_job_status(job_id, **status_update)
            logger.info(f"[Job {job_id}] Status: Failed")
        finally:
            pass # Cleanup if needed

    # --- NEW Method for Custom Selection Analysis --- 
    async def run_custom_analysis_background(self, job_id: str, request: CustomAnalysisRequest):
        """Background task for running analysis on custom-selected points.
        
        Args:
            job_id: The ID of the analysis job.
            request: The request data containing selected points and analysis parameters.
        """
        logger.info(f"Starting custom selection analysis for job: {job_id}")
        
        try:
            # Update job status to running
            await self.job_service.update_job_status(job_id, 
                                              status="running", 
                                              message="Processing custom selection...",
                                              progress=0.1)
            
            # Create DataFrames from the request points
            logger.info(f"Creating DataFrames from {len(request.ligands)} ligand and {len(request.receptors)} receptor points")
            
            # Convert PointData lists to DataFrames
            ligand_df = pd.DataFrame([p.dict() for p in request.ligands])
            receptor_df = pd.DataFrame([p.dict() for p in request.receptors])
            
            # Prepare "spatial" DataFrame by combining points with standard column names
            spatial_data = []
            
            # Add ligand genes
            for ligand_name in set(ligand_df['gene']):
                # Filter points for this ligand
                ligand_points = ligand_df[ligand_df['gene'] == ligand_name]
                for _, point in ligand_points.iterrows():
                    spatial_data.append({
                        'gene': ligand_name, 
                        'x': point['x'], 
                        'y': point['y'],
                        'layer': 'custom_selection' # All points in same "layer" for custom selection
                    })
            
            # Add receptor genes
            for receptor_name in set(receptor_df['gene']):
                # Filter points for this receptor
                receptor_points = receptor_df[receptor_df['gene'] == receptor_name]
                for _, point in receptor_points.iterrows():
                    spatial_data.append({
                        'gene': receptor_name, 
                        'x': point['x'], 
                        'y': point['y'],
                        'layer': 'custom_selection'
                    })
                    
            # Create combined spatial DataFrame
            spatial_df = pd.DataFrame(spatial_data)
            
            # If original job ID was provided, try to get interaction data from it
            interactions_df = None
            if request.original_job_id:
                logger.info(f"Attempting to retrieve interactions data from original job: {request.original_job_id}")
                try:
                    original_job_context = self.job_service.get_job_context(request.original_job_id)
                    if original_job_context:
                        interactions_path = FileService.get_file_path(original_job_context['interactionsFileId'])
                        interactions_mapping = AnalysisMapping(**original_job_context['interactionsMapping'])
                        interactions_df = self._load_and_standardize(interactions_path, interactions_mapping, STANDARDIZED_INTERACTION_COLS)
                        logger.info(f"Successfully loaded interactions data from original job")
                    else:
                        logger.warning(f"Could not find context for original job: {request.original_job_id}")
                except Exception as e:
                    logger.error(f"Error loading interactions data from original job: {e}")
            
            # If we couldn't get interactions from original job, create a minimal one
            if interactions_df is None:
                logger.info("Creating minimal interactions DataFrame from point selection")
                # Create interaction pairs for all ligand-receptor combinations
                interaction_data = []
                for ligand_name in set(ligand_df['gene']):
                    for receptor_name in set(receptor_df['gene']):
                        interaction_data.append({
                            'ligand': ligand_name,
                            'receptor': receptor_name
                        })
                        
                interactions_df = pd.DataFrame(interaction_data)
                
            # Update progress
            await self.job_service.update_job_status(job_id, 
                                                 progress=0.3, 
                                                 message="Running analysis on selection...")
                
            # Run a simplified analysis pipeline on the custom data
            # For now we'll just do pathway dominance
            logger.info(f"Running pathway dominance on custom selection")
            # Use a dummy value for 'custom_selection' in place of multiple tissue layers
            results = {
                'custom_selection': {} # Dict to store results for this selection
            }
            
            # Actual analysis - depending on the points provided
            try:
                pathway_results = run_pathway_dominance_pipeline(spatial_df, interactions_df)
                results['custom_selection']['pathway_dominance'] = pathway_results.get('custom_selection', [])
                
                # Optional: You could also run module_context if the data supports it
                # But that likely needs module information
                
                # Update with completed results
                await self.job_service.update_job_status(job_id,
                                                    status="success",
                                                    progress=1.0,
                                                    message="Custom analysis complete",
                                                    results=results)
                                                    
                logger.info(f"Custom analysis for job {job_id} completed successfully")
                
            except Exception as analysis_error:
                logger.exception(f"Error in custom analysis calculations: {analysis_error}")
                await self.job_service.update_job_status(job_id,
                                                    status="failed",
                                                    message=f"Analysis calculation error: {str(analysis_error)}")
        
        except Exception as e:
            logger.exception(f"Error in custom selection analysis: {e}")
            await self.job_service.update_job_status(job_id,
                                                status="failed",
                                                message=f"Analysis error: {str(e)}")

    # Helper method removed - we now directly call self.job_service.update_job_status
            
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

