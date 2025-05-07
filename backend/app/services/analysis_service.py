import os
from pathlib import Path
import pandas as pd
import anndata as ad
from fastapi import HTTPException, Depends
import time # For simulating work
from typing import Dict, Any, List, Union, Tuple, Optional, Set
import logging # Import logging
import asyncio
from scipy.spatial import ConvexHull # ADDED: Import ConvexHull
import numpy as np # ADDED: Import numpy for array handling
# Need geopandas and shapely for polygon filtering
import geopandas as gpd
from shapely.geometry import Point, Polygon
# MODIFIED: Import new bundle response model and scope result, remove deprecated ones
from ..models.analysis_models import (
    AnalysisResultItem, AnalysisMapping, 
    CustomAnalysisScopeResult, CustomAnalysisResultsBundle,
    CustomLassoAnalysisRequest, 
    ComparisonRequest, ComparisonResponse, ComparisonResults, DifferentialExpressionResult, SelectionData, # Added Comparison models
    ComparisonJobResponse # Keep
)
# Import core logic functions
from ..analysis_logic import core

# Import the new JobService and its dependency function
from app.services.job_service import JobService, get_job_service
from app.services.file_service import FileService
# Import models from the models file
from app.models.analysis_models import AnalysisPayload # REMOVED CustomAnalysisRequest, PointData

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

# Statistics imports
from scipy.stats import fisher_exact, chi2_contingency # For comparing counts, ADDED chi2_contingency
from statsmodels.stats.multitest import multipletests # For FDR correction
from collections import Counter # For counting molecules
import math # For log2
from collections import defaultdict # For Reverse Lookup Maps

# ADDED: Helper function from core.py (or make it a shared utility)
def _split_receptor_complex(receptor_name: str) -> List[str]:
    """Split a receptor name into components if it's a complex (contains '_')."""
    return receptor_name.split('_')

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

            # 2. Run the Full Analysis Pipeline using the standardized function
            await self.job_service.update_job_status(job_id, message="Running core analysis pipeline...", progress=0.4, stage_id='analysis')
            logger.info(f"[Job {job_id}] Status: Running, Progress: 40% (Running core analysis...)")
            
            # Define the progress update callback function
            async def _update_progress_callback(stage_id: str, progress: float, current_scope: str = None, detail: str = None):
                progress_scaled = 0.4 + (progress * 0.5) # Scale pipeline progress (0-1) to overall progress (0.4-0.9)
                message = f"Processing {stage_id}..."
                if current_scope:
                    message += f" (Scope: {current_scope})"
                if detail:
                     message += f" - {detail}"
                await self.job_service.update_job_status(job_id, message=message, progress=progress_scaled, stage_id=stage_id, current_scope=current_scope)
                logger.info(f"[Job {job_id}] Progress: {progress_scaled*100:.1f}% ({message})")

            # Call the single, standardized pipeline function from core logic
            core_results = await core.run_full_analysis_pipeline(
                spatial_df=spatial_df,
                interactions_df=interactions_df,
                modules_df=modules_df,
                update_progress=_update_progress_callback
            )
            
            # final_results variable will hold the direct output from the core function
            final_results = core_results
            logger.info(f"Job {job_id}: Core analysis pipeline complete.")
            await self.job_service.update_job_status(job_id, message="Finalizing results...", progress=0.95, stage_id='cleanup')
            logger.info(f"[Job {job_id}] Status: Running, Progress: 95% (Finalizing results...)")

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

    # --- REMOVED Method for Custom Selection Analysis --- 
    # async def run_custom_analysis_background(self, job_id: str, request: CustomAnalysisRequest):
    #    ...
    #    ...

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
            
    # ADDED: Helper method to load and standardize interactions data
    def _load_standardized_interactions(self, file_id: str, mapping: AnalysisMapping) -> Tuple[Optional[pd.DataFrame], List[str]]:
        """Loads and standardizes interactions data."""
        errors = []
        if not file_id or not mapping:
            errors.append("Missing interactionsFileId or interactionsMapping.")
            return None, errors
        
        # Ensure ligandCol and receptorCol are present in the mapping
        if not mapping.ligandCol or not mapping.receptorCol:
            errors.append("Missing ligandCol or receptorCol in interactionsMapping.")
            return None, errors

        standard_interaction_cols = {
            'ligandCol': 'ligand',
            'receptorCol': 'receptor'
            # Add other interaction-specific columns if needed in the future
        }
        
        try:
            file_path = FileService.get_file_path(file_id)
            logger.info(f"Loading interactions data from: {file_path}")
            
            # Use a simplified mapping for interactions, only taking required ones
            interaction_mapping_dict = {'ligandCol': mapping.ligandCol, 'receptorCol': mapping.receptorCol}
            # We need to pass an AnalysisMapping object to _load_and_standardize
            # Create a minimal one for this purpose.
            # We can't directly create AnalysisMapping with a subset of fields if some are non-optional
            # in the model definition. Let's pass the original mapping but _load_and_standardize
            # should only use the relevant keys based on standard_interaction_cols.
            # This assumes _load_and_standardize correctly handles mapping if only a subset of its 
            # standard_cols_map keys are present in the provided AnalysisMapping's fields.
            # A safer way is to construct a new AnalysisMapping with only relevant fields if possible,
            # or adjust _load_and_standardize.
            # For now, we rely on _load_and_standardize to pick the correct fields.

            df_interactions = self._load_and_standardize(
                file_path,
                mapping, # Pass the original full mapping
                standard_interaction_cols # But guide it to only use these standard cols
            )
            
            if not all(col in df_interactions.columns for col in ['ligand', 'receptor']):
                raise ValueError("Standardized interactions data is missing required columns: ligand, receptor.")
            
            # Ensure receptor column is string for splitting
            df_interactions['receptor'] = df_interactions['receptor'].astype(str)
            logger.info(f"Successfully loaded and standardized interactions data with {len(df_interactions)} pairs.")
            return df_interactions, errors
        except FileNotFoundError:
            errors.append(f"Interactions file (ID: {file_id}) not found.")
            logger.error(f"Interactions file (ID: {file_id}) not found.")
            return None, errors
        except ValueError as ve:
            errors.append(f"Interactions data loading/validation error: {ve}")
            logger.error(f"Interactions data loading/validation error: {ve}")
            return None, errors
        except Exception as e:
            errors.append(f"Unexpected error loading interactions data: {e}")
            logger.exception(f"Unexpected error loading interactions data: {e}")
            return None, errors

    # ADDED: Helper method to get interaction data and classify molecules
    async def _get_interactions_and_molecule_types(
        self, 
        selection_for_interactions: SelectionData, 
        all_spatial_molecules: Set[str]
    ) -> Tuple[Optional[pd.DataFrame], Dict[str, str], List[str]]:
        """
        Loads interactions data and classifies molecules found in spatial data based on their roles.
        """
        errors = []
        molecule_classifications: Dict[str, str] = {mol: 'unknown' for mol in all_spatial_molecules} # Default to unknown
        interactions_df = None

        if not selection_for_interactions.files.interactionsFileId or \
           not selection_for_interactions.mappings.interactionsMapping:
            errors.append("Interactions file ID or mapping missing in selection data used for interactions.")
            # Still return default classifications for all_spatial_molecules
            return None, molecule_classifications, errors

        interactions_df, load_errors = self._load_standardized_interactions(
            selection_for_interactions.files.interactionsFileId,
            selection_for_interactions.mappings.interactionsMapping
        )
        errors.extend(load_errors)

        if interactions_df is None:
            logger.warning(f"Could not load interactions_df. Skipping molecule classification based on interactions. Errors: {errors}")
            # Return default classifications
            return None, molecule_classifications, errors
        
        try:
            db_ligands = set(interactions_df['ligand'].unique())
            db_receptor_components = set()
            for receptor_str in interactions_df['receptor'].astype(str).unique():
                for component in _split_receptor_complex(receptor_str):
                    db_receptor_components.add(component)

            logger.info(f"[Job classification] db_ligands size: {len(db_ligands)}, db_receptor_components size: {len(db_receptor_components)}")
            # Log first few elements, ensure list conversion for slicing
            log_db_ligands = list(db_ligands)
            log_db_receptor_components = list(db_receptor_components)
            logger.info(f"[Job classification] First 10 db_ligands: {log_db_ligands[:10]}")
            logger.info(f"[Job classification] First 10 db_receptor_components: {log_db_receptor_components[:10]}")

            example_molecules_to_trace = ['HFM1', 'NPY', 'SLCO5A1'] # Add more if needed from user output or expected data

            for molecule in all_spatial_molecules:
                is_ligand = molecule in db_ligands
                is_receptor_component = molecule in db_receptor_components

                # Updated classification logic
                classification = 'unknown' # Default if not in interactions db
                if is_ligand and not is_receptor_component:
                    classification = 'single_ligand'
                elif not is_ligand and is_receptor_component:
                    classification = 'single_receptor'
                elif is_ligand and is_receptor_component:
                    classification = 'ligand_and_receptor' # Molecule plays both roles
                # else: remains 'unknown' if neither ligand nor receptor component
                
                molecule_classifications[molecule] = classification

                if molecule in example_molecules_to_trace:
                    logger.info(f"[Job classification] Tracing molecule: {molecule} -> In db_ligands: {is_ligand}, In db_receptors: {is_receptor_component}, Classified as: {classification}")
            
            logger.info("Molecule classification based on interactions complete.")

        except Exception as e:
            error_msg = f"Error during molecule classification based on interactions: {e}"
            logger.exception(error_msg)
            errors.append(error_msg)
            # In case of error here, molecule_classifications might be partially updated or default
            # We return the current state of classifications and the error. interactions_df might be valid.

        return interactions_df, molecule_classifications, errors
            
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

    async def run_custom_analysis(
        self,
        original_job_id: str,
        polygon_coords: List[List[float]],
        initial_context: Dict[str, Any]
    ) -> CustomAnalysisResultsBundle:
        """Runs the analysis pipeline on a subset of spatial data defined by a polygon,
           returning results aggregated across the whole selection AND by layer within the selection.
        """
        
        service_start_time = time.time() # DEBUG
        logger.debug(f"[Custom Analysis {original_job_id}] Service layer started (calculating whole + layers).") # DEBUG

        try:
            # 1. Extract context (file IDs, mappings)
            inputs = initial_context.get('inputs', {})
            files = inputs.get('files', {})
            mappings = inputs.get('mappings', {})

            spatial_file_id = files.get('spatialFileId')
            interactions_file_id = files.get('interactionsFileId')
            modules_file_id = files.get('modulesFileId')

            spatial_mapping = mappings.get('spatialMapping', {})
            interactions_mapping = mappings.get('interactionsMapping', {})
            modules_mapping = mappings.get('modulesMapping', {})

            if not all([spatial_file_id, interactions_file_id, modules_file_id, 
                        spatial_mapping, interactions_mapping, modules_mapping]):
                raise ValueError("Missing required file IDs or mappings in the original job context.")

            # Convert mapping dicts back to AnalysisMapping objects for _load_and_standardize
            spatial_mapping_obj = AnalysisMapping(**spatial_mapping)
            interactions_mapping_obj = AnalysisMapping(**interactions_mapping)
            modules_mapping_obj = AnalysisMapping(**modules_mapping)

            if not all([spatial_file_id, interactions_file_id, modules_file_id, 
                        spatial_mapping_obj, interactions_mapping_obj, modules_mapping_obj]): # Check objects
                raise ValueError("Missing required file IDs or mappings in the original job context.")

            # 2. Load original data using _load_and_standardize
            load_start_time = time.time() # DEBUG
            logger.debug(f"[Custom Analysis {original_job_id}] Loading original data files...") # DEBUG
            spatial_path = FileService.get_file_path(spatial_file_id)
            interactions_path = FileService.get_file_path(interactions_file_id)
            modules_path = FileService.get_file_path(modules_file_id)
            
            spatial_df = self._load_and_standardize(spatial_path, spatial_mapping_obj, STANDARDIZED_SPATIAL_COLS)
            interactions_df = self._load_and_standardize(interactions_path, interactions_mapping_obj, STANDARDIZED_INTERACTION_COLS)
            modules_df = self._load_and_standardize(modules_path, modules_mapping_obj, STANDARDIZED_MODULE_COLS)
            logger.debug(f"[Custom Analysis {original_job_id}] Data loading duration: {time.time() - load_start_time:.4f} seconds.") # DEBUG
            logger.debug(f"[Custom Analysis {original_job_id}] Spatial data loaded with {len(spatial_df)} points initially.") # DEBUG

            # 3. Filter Spatial Data by Polygon
            filter_start_time = time.time() # DEBUG
            logger.debug(f"[Custom Analysis {original_job_id}] Filtering spatial data using lasso polygon...") # DEBUG
            if not polygon_coords or len(polygon_coords) < 4:
                 raise ValueError("Lasso polygon must have at least 4 vertices (including closing point).")
            if 'x' not in spatial_df.columns or 'y' not in spatial_df.columns:
                 raise ValueError("Spatial data must contain 'x' and 'y' columns for polygon filtering.")
            
            # Convert DataFrame to GeoDataFrame
            try:
                # OPTIMIZATION: Use efficient points_from_xy to create geometry
                gdf_creation_start_time = time.time() # DEBUG
                spatial_gdf = gpd.GeoDataFrame(
                    spatial_df, geometry=gpd.points_from_xy(spatial_df.x, spatial_df.y)
                )
                logger.debug(f"[Custom Analysis {original_job_id}] GeoDataFrame creation duration: {time.time() - gdf_creation_start_time:.4f} seconds.") # DEBUG
                
                # OPTIMIZATION: Create a spatial index
                index_start_time = time.time() # DEBUG
                spatial_gdf.sindex
                logger.debug(f"[Custom Analysis {original_job_id}] Spatial index creation duration: {time.time() - index_start_time:.4f} seconds.") # DEBUG
                
                # Create Shapely Polygon from coordinates
                lasso_polygon = Polygon(polygon_coords)
                
                # OPTIMIZATION: Pre-filter using bounding box intersection
                bbox_filter_start_time = time.time() # DEBUG
                possible_matches_index = list(spatial_gdf.sindex.intersection(lasso_polygon.bounds))
                possible_matches = spatial_gdf.iloc[possible_matches_index]
                logger.debug(f"[Custom Analysis {original_job_id}] BBox pre-filter reduced points from {len(spatial_gdf)} to {len(possible_matches)}. Duration: {time.time() - bbox_filter_start_time:.4f} seconds.") # DEBUG
                
                # Perform precise spatial query only on the pre-filtered subset
                precise_filter_start_time = time.time() # DEBUG
                points_within_lasso = possible_matches[possible_matches.within(lasso_polygon)]
                logger.debug(f"[Custom Analysis {original_job_id}] Precise filter duration: {time.time() - precise_filter_start_time:.4f} seconds.") # DEBUG
                
                # Drop the temporary geometry column if needed, keeping the filtered DataFrame
                filtered_spatial_df = pd.DataFrame(points_within_lasso.drop(columns='geometry'))
            except Exception as geo_error:
                 logger.exception(f"Error during geometric filtering: {geo_error}")
                 raise HTTPException(status_code=500, detail=f"Error during geometric filtering: {geo_error}")

            logger.debug(f"[Custom Analysis {original_job_id}] Filtering duration: {time.time() - filter_start_time:.4f} seconds.") # DEBUG
            logger.debug(f"[Custom Analysis {original_job_id}] Filtered spatial data contains {len(filtered_spatial_df)} points.") # DEBUG

            # Initialize results containers
            whole_results_data: CustomAnalysisScopeResult = CustomAnalysisScopeResult(pathway_dominance=[], module_context=[])
            layered_results_data: Dict[str, CustomAnalysisScopeResult] = {}

            # Check if any points remain after filtering
            if filtered_spatial_df.empty:
                 logger.warning(f"[Custom Analysis {original_job_id}] No spatial points found within the lasso selection. Returning empty results bundle.")
                 # Return empty bundle
                 return CustomAnalysisResultsBundle(whole_results=whole_results_data, layered_results=layered_results_data)

            # 4. Run Analysis Pipeline for Whole Selection
            whole_call_start_time = time.time() # DEBUG
            logger.debug(f"[Custom Analysis {original_job_id}] Calling core single-scope analysis pipeline for whole selection...") # DEBUG
            
            whole_custom_results_dict = await core.run_analysis_pipeline_from_dataframes(
                all_spatial_df=filtered_spatial_df, # Pass the polygon-filtered spatial data
                interactions_df=interactions_df,   # Pass the full interactions data
                modules_df=modules_df            # Pass the full modules data
            )
            
            # Extract results (assuming the function returns a dict with a single key like 'custom')
            whole_analysis_output = next(iter(whole_custom_results_dict.values()), {}) if whole_custom_results_dict else {}
            whole_results_data = CustomAnalysisScopeResult(
                pathway_dominance=whole_analysis_output.get('pathway_dominance', []),
                module_context=whole_analysis_output.get('module_context', [])
            )
            logger.debug(f"[Custom Analysis {original_job_id}] Core analysis duration (whole_custom): {time.time() - whole_call_start_time:.4f} seconds.") # DEBUG

            # 5. Run Analysis Pipeline Per Layer (if layer column exists)
            layers_call_start_time = time.time() # DEBUG
            if 'layer' in filtered_spatial_df.columns:
                logger.debug(f"[Custom Analysis {original_job_id}] Running layer-specific analysis within the custom selection...") # DEBUG
                
                # Ensure layer column is string type for consistent grouping
                filtered_spatial_df['layer'] = filtered_spatial_df['layer'].astype(str)
                unique_layers = filtered_spatial_df['layer'].unique()
                logger.debug(f"[Custom Analysis {original_job_id}] Found layers in selection: {unique_layers}")

                for layer_name in unique_layers:
                    layer_start_time = time.time() # DEBUG
                    logger.debug(f"[Custom Analysis {original_job_id}] Processing layer: {layer_name}")
                    layer_df = filtered_spatial_df[filtered_spatial_df['layer'] == layer_name]
                    
                    if layer_df.empty:
                        logger.warning(f"[Custom Analysis {original_job_id}] Skipping empty layer: {layer_name}")
                        layered_results_data[layer_name] = CustomAnalysisScopeResult(pathway_dominance=[], module_context=[]) # Add empty result
                        continue
                        
                    logger.debug(f"[Custom Analysis {original_job_id}] Calling core analysis for layer '{layer_name}' ({len(layer_df)} points)... ")
                    # Run analysis for this layer's subset of points
                    layer_results_dict = await core.run_analysis_pipeline_from_dataframes(
                        all_spatial_df=layer_df, # Pass the layer-specific spatial data
                        interactions_df=interactions_df, # Pass full interactions
                        modules_df=modules_df # Pass full modules
                    )
                    
                    # Extract results for this layer
                    layer_output = next(iter(layer_results_dict.values()), {}) if layer_results_dict else {}
                    
                    # Store results for this layer using the ScopeResult model
                    layered_results_data[layer_name] = CustomAnalysisScopeResult(
                        pathway_dominance=layer_output.get('pathway_dominance', []),
                        module_context=layer_output.get('module_context', [])
                    )
                    logger.debug(f"[Custom Analysis {original_job_id}] Layer '{layer_name}' processing time: {time.time() - layer_start_time:.4f} seconds.") # DEBUG
                
                logger.debug(f"[Custom Analysis {original_job_id}] Total core analysis duration (custom_by_layer): {time.time() - layers_call_start_time:.4f} seconds.") # DEBUG
            else:
                logger.warning(f"[Custom Analysis {original_job_id}] 'layer' column not found in filtered spatial data. Skipping layer-specific analysis.")

            # 6. Return the results bundle
            logger.debug(f"[Custom Analysis {original_job_id}] Service layer duration: {time.time() - service_start_time:.4f} seconds.") # DEBUG
            return CustomAnalysisResultsBundle(
                whole_results=whole_results_data,
                layered_results=layered_results_data
            )

        except Exception as e:
            logger.exception(f"Error during custom analysis service execution for job {original_job_id}: {e}")
            logger.debug(f"[Custom Analysis {original_job_id}] Service layer failed after {time.time() - service_start_time:.4f} seconds.") # DEBUG
            # Re-raise as a generic internal server error if not already an HTTPException or specific ValueError
            if not isinstance(e, (HTTPException, ValueError, FileNotFoundError)):
                raise HTTPException(status_code=500, detail=f"Internal server error during custom analysis: {str(e)}")
            else:
                raise e # Re-raise specific handled exceptions

    # --- Method for Comparison Analysis (now runs in background) --- #
    async def run_comparison_background(self, job_id: str, request: ComparisonRequest):
        """Performs the comparison logic asynchronously and updates job status."""
        logger.info(f"[Job {job_id}] Starting comparison analysis background task...")
        final_comparison_results = ComparisonResults(differential_expression=[])
        errors = []
        start_time = time.time()
        
        interactions_df: Optional[pd.DataFrame] = None
        molecule_classifications: Dict[str, str] = {}

        try:
            await self.job_service.update_job_status(job_id, status="running", message="Loading comparison data (Selection 1)...", progress=0.1)
            
            # 1. Process Selection 1
            df1, errors1 = await self._load_and_filter_selection(request.selection1)
            errors.extend(errors1)
            if df1 is None or df1.empty:
                errors.append("Could not load or filter data for Selection 1.")
                logger.error(f"[Job {job_id}] Failed to get data for Selection 1. Errors: {errors1}")
                raise ValueError("Failed to get data for Selection 1.")
            logger.info(f"[Job {job_id}] Selection 1 data loaded with {len(df1)} points.")

            await self.job_service.update_job_status(job_id, message="Loading comparison data (Selection 2)...", progress=0.2)
            # 2. Process Selection 2
            df2, errors2 = await self._load_and_filter_selection(request.selection2)
            errors.extend(errors2)
            if df2 is None or df2.empty:
                errors.append("Could not load or filter data for Selection 2.")
                logger.error(f"[Job {job_id}] Failed to get data for Selection 2. Errors: {errors2}")
                raise ValueError("Failed to get data for Selection 2.")
            logger.info(f"[Job {job_id}] Selection 2 data loaded with {len(df2)} points.")

            # Get all unique molecules from spatial data to guide classification
            all_spatial_molecules = set(df1['gene'].unique()) | set(df2['gene'].unique())
            if not all_spatial_molecules:
                errors.append("No molecules found in spatial data of either selection.")
                logger.error(f"[Job {job_id}] No molecules found in spatial data.")
                raise ValueError("No molecules to compare after loading spatial data.")

            await self.job_service.update_job_status(job_id, message="Loading interactions and classifying molecules...", progress=0.3)
            # 3. Get Interactions and Classify Molecules
            # Assuming selection1's interaction file is the reference. This could be a parameter.
            # Or, ensure both selections use the same job_id if "same project" comparison.
            # For now, using selection1 as the source for interaction definitions.
            interactions_df, molecule_classifications, typing_errors = await self._get_interactions_and_molecule_types(
                request.selection1, 
                all_spatial_molecules
            )
            errors.extend(typing_errors)
            if typing_errors:
                 logger.warning(f"[Job {job_id}] Errors during molecule typing/interaction loading: {typing_errors}")
            if interactions_df is None:
                logger.warning(f"[Job {job_id}] Interactions DataFrame could not be loaded. L-R pair analysis will be skipped.")
            else:
                logger.info(f"[Job {job_id}] Interactions data loaded with {len(interactions_df)} pairs. Molecules classified.")


            await self.job_service.update_job_status(job_id, message="Comparing individual molecule counts...", progress=0.4)
            
            # --- 4. INDIVIDUAL MOLECULE COMPARISON ---
            logger.info(f"[Job {job_id}] Comparing individual molecules (Total unique: {len(all_spatial_molecules)})...")
            
            mol_comp_overall_start_time = time.time() # Overall start for this section

            counter_start_time = time.time()
            counts1 = Counter(df1['gene'])
            logger.info(f"[Job {job_id}] Counter for Selection 1 (size {len(df1)}) created in {time.time() - counter_start_time:.4f}s.")

            counter_start_time = time.time()
            counts2 = Counter(df2['gene'])
            logger.info(f"[Job {job_id}] Counter for Selection 2 (size {len(df2)}) created in {time.time() - counter_start_time:.4f}s.")

            N1_mol = len(df1) # Total observations in selection 1 (spots/cells)
            N2_mol = len(df2) # Total observations in selection 2

            if N1_mol == 0 or N2_mol == 0:
                 errors.append("One or both selections have zero spatial entities for molecule counting.")
                 logger.error(f"[Job {job_id}] Zero spatial entities in one/both selections for molecule counts (N1_mol={N1_mol}, N2_mol={N2_mol}).")
            
            p_values_mol = []
            test_details_mol = []

            logger.info(f"[Job {job_id}] Starting individual molecule comparison loop for {len(all_spatial_molecules)} molecules...")
            loop_start_time = time.time()

            for mol in all_spatial_molecules:
                count1 = counts1.get(mol, 0)
                count2 = counts2.get(mol, 0)
                
                # Contingency table for this molecule vs all others
                #        Sel 1 | Sel 2
                # Gene X   a     |   b
                # Other  c     |   d
                # a = count1 (gene X in Sel 1)
                # b = count2 (gene X in Sel 2)
                # c = N1_mol - count1 (other genes in Sel 1)
                # d = N2_mol - count2 (other genes in Sel 2)
                table = [[count1, N1_mol - count1], [count2, N2_mol - count2]]
                
                p_value = 1.0 # Default p-value
                try:
                    # Use chi2_contingency for individual molecule comparison due to large N
                    # It's generally faster and a good approximation when expected frequencies are not too small.
                    if (N1_mol > 0 and N2_mol > 0) and (count1 + count2 > 0): # Ensure some data to test
                        # Check for non-negative values in table, which chi2_contingency expects
                        if table[0][1] < 0 or table[1][1] < 0:
                            # This can happen if count > N_mol, which indicates an issue with N_mol calculation or counts
                            logger.warning(f"[Job {job_id}] Negative value in contingency table for molecule '{mol}'. Table: {table}. N1_mol={N1_mol}, N2_mol={N2_mol}. Assigning p=1.")
                            # This case should ideally not happen if N1_mol and N2_mol are total counts including 'mol'
                        else:
                            chi2, p_value, dof, expected = chi2_contingency(table, correction=False) # correction=False is often recommended for 2x2
                    else:
                        # If no counts or no totals, no basis for test, assign non-significant p-value
                        pass # p_value remains 1.0
                                        
                    p_values_mol.append(p_value)
                    test_details_mol.append({
                        'id': mol, 
                        'p_value': p_value, 
                        'count1': count1, 
                        'count2': count2,
                        'n1': N1_mol, 
                        'n2': N2_mol, 
                        'type': 'molecule'
                    })                    
                except ValueError as stat_err: # e.g., if a row/col sum is zero for chi2 or other issues
                    logger.warning(f"[Job {job_id}] Chi-squared test failed for molecule '{mol}'. Table: {table}. Error: {stat_err}. Assigning p=1.")
                    p_values_mol.append(1.0) # Ensure p_value is 1.0 before appending
                    test_details_mol.append({
                        'id': mol, 'p_value': 1.0, 'count1': count1, 'count2': count2, 
                        'n1':N1_mol, 'n2':N2_mol, 'type': 'molecule'
                    })
            
            logger.info(f"[Job {job_id}] Individual molecule comparison loop finished in {time.time() - loop_start_time:.4f}s.")
            logger.info(f"[Job {job_id}] Individual molecule comparison section (including counters and loop) took {time.time() - mol_comp_overall_start_time:.4f}s.")
            logger.info(f"[Job {job_id}] Individual molecule comparison processed {len(test_details_mol)} molecules.") # This is an existing log
            # --- END INDIVIDUAL MOLECULE COMPARISON ---

            # --- 5. L-R PAIR CO-EXPRESSION COMPARISON ---
            test_details_lr = []
            p_values_lr = []
            if interactions_df is not None and not interactions_df.empty and N1_mol > 0 and N2_mol > 0 : # Check N1_mol, N2_mol as well
                await self.job_service.update_job_status(job_id, message="Comparing L-R pair co-expression...", progress=0.6)
                logger.info(f"[Job {job_id}] Comparing L-R pair co-expression using {len(interactions_df)} pairs from DB.")
                lr_comparison_start_time = time.time()

                # --- Pre-computation using Reverse Lookup Maps ---
                logger.info(f"[Job {job_id}] Calculating unique spot counts...")
                calc_spots_start_time = time.time()
                total_spots_s1 = df1[['x', 'y']].drop_duplicates().shape[0]
                total_spots_s2 = df2[['x', 'y']].drop_duplicates().shape[0]
                logger.info(f"[Job {job_id}] Unique spots: S1={total_spots_s1}, S2={total_spots_s2}. Calculation took {time.time() - calc_spots_start_time:.4f}s.")

                logger.info(f"[Job {job_id}] Building gene-to-spots map for Selection 1...")
                map_s1_start_time = time.time()
                gene_to_spots_map_s1 = defaultdict(set)
                # Use itertuples for efficiency
                for row in df1[['gene', 'x', 'y']].itertuples(index=False, name='SpatialPoint'):
                    gene_to_spots_map_s1[row.gene].add((row.x, row.y))
                logger.info(f"[Job {job_id}] Gene map S1 built in {time.time() - map_s1_start_time:.4f}s. Unique genes mapped: {len(gene_to_spots_map_s1)}")

                logger.info(f"[Job {job_id}] Building gene-to-spots map for Selection 2...")
                map_s2_start_time = time.time()
                gene_to_spots_map_s2 = defaultdict(set)
                for row in df2[['gene', 'x', 'y']].itertuples(index=False, name='SpatialPoint'):
                    gene_to_spots_map_s2[row.gene].add((row.x, row.y))
                logger.info(f"[Job {job_id}] Gene map S2 built in {time.time() - map_s2_start_time:.4f}s. Unique genes mapped: {len(gene_to_spots_map_s2)}")
                # --- End Pre-computation ---

                # Check if any spots exist before proceeding
                if total_spots_s1 > 0 or total_spots_s2 > 0: 
                    logger.info(f"[Job {job_id}] Starting L-R pair processing loop for {len(interactions_df)} pairs...")
                    loop_start_time = time.time()
                    processed_pairs_count = 0

                    # Pre-calculate receptor components cache (already added in previous step)
                    receptor_components_cache = {
                        r_str: _split_receptor_complex(r_str) 
                        for r_str in interactions_df['receptor'].astype(str).unique()
                    }

                    for _, interaction_row in interactions_df.iterrows():
                        pair_processing_start_time = time.time()
                        ligand = interaction_row['ligand']
                        receptor_complex_str = str(interaction_row['receptor'])
                        receptor_components = receptor_components_cache.get(receptor_complex_str, [])
                        if not receptor_components: 
                           receptor_components = _split_receptor_complex(receptor_complex_str)

                        # --- Optimized Co-expression Counting using Maps --- 
                        count_s1_start_time = time.time()
                        coexp_s1 = 0
                        if total_spots_s1 > 0:
                            ligand_spots_s1 = gene_to_spots_map_s1.get(ligand, set())
                            if ligand_spots_s1: # Only proceed if ligand is present
                                rc_spot_sets_s1 = [gene_to_spots_map_s1.get(rc, set()) for rc in receptor_components]
                                # Check if any receptor component is completely absent
                                if all(rc_spot_sets_s1):
                                    coexpressing_spots_s1 = ligand_spots_s1.intersection(*rc_spot_sets_s1)
                                    coexp_s1 = len(coexpressing_spots_s1)
                                # else: coexp_s1 remains 0 if a component is missing
                        count_s1_duration = time.time() - count_s1_start_time
                        
                        count_s2_start_time = time.time()
                        coexp_s2 = 0
                        if total_spots_s2 > 0:
                            ligand_spots_s2 = gene_to_spots_map_s2.get(ligand, set())
                            if ligand_spots_s2:
                                rc_spot_sets_s2 = [gene_to_spots_map_s2.get(rc, set()) for rc in receptor_components]
                                if all(rc_spot_sets_s2):
                                    coexpressing_spots_s2 = ligand_spots_s2.intersection(*rc_spot_sets_s2)
                                    coexp_s2 = len(coexpressing_spots_s2)
                        count_s2_duration = time.time() - count_s2_start_time
                        # --- End Optimized Counting --- #

                        # Fisher's Exact Test for co-expression counts
                        fisher_start_time = time.time()
                        table_lr = [[coexp_s1, total_spots_s1 - coexp_s1], [coexp_s2, total_spots_s2 - coexp_s2]]
                        
                        try:
                            if (total_spots_s1 > 0 and total_spots_s2 > 0) or (coexp_s1 + coexp_s2 > 0):
                                _, p_value_lr = fisher_exact(table_lr, alternative='two-sided')
                            else:
                                p_value_lr = 1.0
                            p_values_lr.append(p_value_lr)
                            test_details_lr.append({
                                'id': f"{ligand}<->{receptor_complex_str}", 
                                'p_value': p_value_lr, 
                                'count1': coexp_s1, 
                                'count2': coexp_s2,
                                'n1': total_spots_s1,
                                'n2': total_spots_s2,
                                'type': 'lr_pair',
                                'ligand_id': ligand,
                                'receptor_id': receptor_complex_str
                            })
                        except ValueError as fisher_err_lr:
                            logger.warning(f"[Job {job_id}] Fisher test failed for L-R pair {ligand}<->{receptor_complex_str}. Table: {table_lr}. Error: {fisher_err_lr}. Assigning p=1.")
                            p_values_lr.append(1.0)
                            test_details_lr.append({
                                'id': f"{ligand}<->{receptor_complex_str}", 'p_value': 1.0, 'count1': coexp_s1, 'count2': coexp_s2,
                                'n1':total_spots_s1, 'n2':total_spots_s2, 'type': 'lr_pair', 
                                'ligand_id': ligand, 'receptor_id': receptor_complex_str
                            })
                        fisher_duration = time.time() - fisher_start_time
                        pair_processing_duration = time.time() - pair_processing_start_time
                        processed_pairs_count +=1
                        if processed_pairs_count % 100 == 0: # Log progress every 100 pairs
                            logger.info(f"[Job {job_id}] Processed {processed_pairs_count}/{len(interactions_df)} L-R pairs. Last pair ({ligand}<->{receptor_complex_str} [Components: {receptor_components}]) took {pair_processing_duration:.4f}s (S1_count_t:{count_s1_duration:.4f}s, S2_count_t:{count_s2_duration:.4f}s, fisher_t:{fisher_duration:.4f}s).")

                    logger.info(f"[Job {job_id}] L-R pair processing loop finished in {time.time() - loop_start_time:.4f}s.")
                else:
                    logger.warning(f"[Job {job_id}] Skipping L-R pair co-expression calculation as total_spots_s1 and total_spots_s2 are both zero.")

                logger.info(f"[Job {job_id}] Full L-R pair co-expression comparison section took {time.time() - lr_comparison_start_time:.4f}s.")
            
            logger.info(f"[Job {job_id}] L-R pair co-expression comparison processed {len(test_details_lr)} pairs.") # This log might be redundant now or need adjustment
            # --- END L-R PAIR CO-EXPRESSION COMPARISON ---

            await self.job_service.update_job_status(job_id, message="Performing FDR correction...", progress=0.8)

            all_test_details = test_details_mol + test_details_lr
            all_p_values = p_values_mol + p_values_lr
            
            if not all_p_values:
                 # errors.append("No p-values generated from tests for FDR correction.") # Not a fatal error for the job
                 logger.warning(f"[Job {job_id}] No p-values were generated from any tests. Skipping FDR correction and results formatting.")
                 # Proceed to set success status with empty results if everything else was fine.
            else:
                # 6. Apply FDR Correction
                reject, q_values, _, _ = multipletests(all_p_values, alpha=request.fdr_threshold, method='fdr_bh')
                logger.info(f"[Job {job_id}] FDR correction done. Total tests: {len(all_p_values)}, Significant results before filtering: {sum(reject)}.")

                # 7. Format Significant Results
                significant_results_list = [] # Changed from significant_results to avoid type hint conflict
                pseudo_count = 1e-9 # For log2fc calculation to avoid division by zero

                for i, detail in enumerate(all_test_details):
                    if not reject[i]: 
                        continue 

                    item_id = detail['id']
                    item_type_str = detail['type']
                    
                    final_type = 'unknown' 
                    if item_type_str == 'molecule':
                        final_type = molecule_classifications.get(item_id, 'unknown')
                    elif item_type_str == 'lr_pair':
                        final_type = 'ligand_receptor_pair'
                    
                    # --- EXPLICIT FILTERING & FORMATTING --- 
                    # Only proceed if the type is one we want to show
                    if final_type in ['single_ligand', 'single_receptor', 'ligand_receptor_pair']:
                        
                        # Common calculations
                        count1 = detail['count1']
                        count2 = detail['count2']
                        n1_total = detail['n1'] 
                        n2_total = detail['n2'] 
                        mean1 = count1 / n1_total if n1_total > 0 else 0
                        mean2 = count2 / n2_total if n2_total > 0 else 0
                        
                        try:
                            log2fc = math.log2((mean2 + pseudo_count) / (mean1 + pseudo_count))
                        except ValueError: 
                             log2fc = None
                             logger.warning(f"[Job {job_id}] Log2FC calculation failed for item {item_id}. mean1={mean1}, mean2={mean2}")
                        
                        result_item_args = {
                            "molecule_id": item_id, 
                            "type": final_type, 
                            "log2_fold_change": log2fc,
                            "p_value": detail['p_value'],
                            "q_value": q_values[i],
                            "mean_selection1": mean1, 
                            "mean_selection2": mean2
                        }
                        # Add L-R specific fields if needed
                        if final_type == 'ligand_receptor_pair':
                            result_item_args["ligand_id"] = detail.get('ligand_id') # Use .get for safety
                            result_item_args["receptor_id"] = detail.get('receptor_id')

                        result_item = DifferentialExpressionResult(**result_item_args)
                        significant_results_list.append(result_item)
                    # --- END EXPLICIT FILTERING & FORMATTING --- 
                        
                significant_results_list.sort(key=lambda x: x.q_value if x.q_value is not None else 1.0)
                final_comparison_results.differential_expression = significant_results_list
                logger.info(f"[Job {job_id}] Formatted {len(significant_results_list)} significant results after FDR.")


            # 8. Update final job status to success
            final_message = f"Comparison finished in {time.time() - start_time:.2f}s. Found {len(final_comparison_results.differential_expression)} significant differences."
            if errors: # Append non-fatal errors to message
                final_message += f" Encountered non-fatal issues: {'; '.join(errors[:3])}{'...' if len(errors) > 3 else ''}"
            
            await self.job_service.update_job_status(
                job_id,
                status="success",
                message=final_message,
                progress=1.0,
                # Store full response in results, including any non-fatal errors collected
                results=ComparisonResponse(results=final_comparison_results, errors=errors).dict() 
            )
            logger.info(f"[Job {job_id}] Comparison analysis successful. {final_message}")

        except ValueError as ve: # Catch specific ValueErrors raised for critical failures
            logger.error(f"[Job {job_id}] Comparison analysis failed critically due to ValueError: {ve}")
            errors.append(f"Critical error: {str(ve)}")
            await self.job_service.update_job_status(
                job_id,
                status="failed",
                message=f"Comparison critical failure: {str(ve)}",
                errors=errors # Include all errors collected so far
            )
        except Exception as e:
            logger.exception(f"[Job {job_id}] Comparison analysis failed with unexpected exception: {e}")
            errors.append(f"Unexpected runtime error: {str(e)}")
            await self.job_service.update_job_status(
                job_id,
                status="failed",
                message=f"Comparison failed unexpectedly: {str(e)}",
                errors=errors # Include all errors collected so far
            )
            
    async def _load_and_filter_selection(self, selection: SelectionData) -> Tuple[pd.DataFrame | None, List[str]]:
        """Helper function to load, standardize, and filter spatial data for one selection."""
        errors = []
        df_spatial = None
        standard_spatial_cols = {"geneCol": "gene", "xCol": "x", "yCol": "y", "layerCol": "layer"}

        if not selection.files.spatialFileId or not selection.mappings.spatialMapping:
            errors.append(f"Missing spatialFileId or spatialMapping for selection (Source Job: {selection.source_job_id or 'N/A'})")
            return None, errors
        if not selection.mappings.spatialMapping.geneCol or not selection.mappings.spatialMapping.xCol or not selection.mappings.spatialMapping.yCol:
             errors.append(f"Missing required spatial mapping (geneCol, xCol, yCol) for selection (Source Job: {selection.source_job_id or 'N/A'})")
             return None, errors

        try:
            spatial_path = FileService.get_file_path(selection.files.spatialFileId)
            logger.info(f"Loading spatial data for selection from: {spatial_path}")
            df_spatial = self._load_and_standardize(
                spatial_path,
                selection.mappings.spatialMapping, 
                standard_spatial_cols
            )
            logger.info(f"Loaded spatial data with {len(df_spatial)} rows.")

            if not all(col in df_spatial.columns for col in ['gene', 'x', 'y']):
                 raise ValueError("Standardized spatial data is missing required columns: gene, x, y")

            if selection.type == 'layer':
                if not selection.definition.layer_name:
                    errors.append("Layer name missing in definition for layer-type selection.")
                    return None, errors
                layer_col_standard = standard_spatial_cols.get('layerCol', 'layer')
                if layer_col_standard not in df_spatial.columns:
                     errors.append(f"Standardized layer column '{layer_col_standard}' not found in spatial data. Cannot filter by layer.")
                     return None, errors
                layer_name = selection.definition.layer_name
                original_count = len(df_spatial)
                df_spatial = df_spatial[df_spatial[layer_col_standard] == layer_name]
                logger.info(f"Filtered for layer '{layer_name}'...")
                if df_spatial.empty:
                     errors.append(f"No data for layer '{layer_name}'.")

            elif selection.type == 'lasso':
                if not selection.definition.polygon_coords or len(selection.definition.polygon_coords) < 4:
                    errors.append("Polygon coordinates missing or insufficient for lasso-type selection.")
                    return None, errors
                original_count = len(df_spatial)
                polygon = Polygon(selection.definition.polygon_coords)
                try:
                    geometry = [Point(xy) for xy in zip(df_spatial['x'], df_spatial['y'])]
                    gdf = gpd.GeoDataFrame(df_spatial, geometry=geometry)
                    gdf_filtered = gdf[gdf.geometry.within(polygon)]
                    df_spatial = pd.DataFrame(gdf_filtered.drop(columns='geometry'))
                    logger.info(f"Filtered by lasso polygon...")
                    if df_spatial.empty:
                         errors.append(f"No data points in lasso.")
                except Exception as geo_err:
                    logger.exception(f"Error during GeoPandas filtering: {geo_err}")
                    errors.append(f"Failed to apply lasso filter: {geo_err}")
                    return None, errors

            elif selection.type == 'whole_tissue':
                logger.info("Using whole tissue...")
                pass 

            return df_spatial, errors # Success return

        except FileNotFoundError:
            errors.append(f"Spatial file not found...")
            return None, errors
        except ValueError as ve:
            errors.append(f"Data loading/validation error...: {ve}")
            return None, errors
        except Exception as e:
            errors.append(f"Unexpected error loading/filtering...: {e}")
            return None, errors
            
    # --- END: Method for Comparison Analysis --- #

    # ... (Keep run_custom_analysis, etc.) ...

# ... (Keep get_analysis_service_dependency) ...

