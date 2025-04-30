from pydantic import BaseModel, Field
from typing import Dict, Optional, List

# Pydantic model for the mapping part of the payload
class AnalysisMapping(BaseModel):
    geneCol: Optional[str] = None
    xCol: Optional[str] = None
    yCol: Optional[str] = None
    layerCol: Optional[str] = None
    ligandCol: Optional[str] = None
    receptorCol: Optional[str] = None
    moduleCol: Optional[str] = None

# Pydantic model for the main analysis request payload
class AnalysisPayload(BaseModel):
    spatialFileId: str = Field(..., description="Unique ID of the uploaded spatial data file")
    spatialMapping: AnalysisMapping = Field(..., description="Column mappings for the spatial data")
    interactionsFileId: str = Field(..., description="Unique ID of the uploaded interactions data file")
    interactionsMapping: AnalysisMapping = Field(..., description="Column mappings for the interactions data")
    modulesFileId: str = Field(..., description="Unique ID of the uploaded modules data file")
    modulesMapping: AnalysisMapping = Field(..., description="Column mappings for the modules data")

# --- Models moved from analysis_routes.py ---
# REMOVING old PointData and CustomAnalysisRequest definitions
# class PointData(BaseModel):
#     x: float
#     y: float
#     gene: str
#     # layer: Optional[str] = None # Layer might not be relevant from lasso
# 
# class CustomAnalysisRequest(BaseModel):
#     ligands: List[PointData]
#     receptors: List[PointData]
#     # Potentially add original jobId if needed to access other files like interactions/modules?
#     # original_job_id: Optional[str] = None 

# --- NEW Models for Custom Lasso Analysis ---

class CustomLassoAnalysisRequest(BaseModel):
    polygon: List[List[float]] = Field(..., description="List of [x, y] coordinates defining the lasso polygon")
    # We might need original file IDs/mappings here if they aren't easily retrievable from jobId context
    # spatialFileId: Optional[str] = None
    # interactionsFileId: Optional[str] = None
    # modulesFileId: Optional[str] = None
    # spatialMapping: Optional[AnalysisMapping] = None
    # interactionsMapping: Optional[AnalysisMapping] = None
    # modulesMapping: Optional[AnalysisMapping] = None

# Response model for analysis results (can be reused or adapted)
# Assuming the structure returned by the core logic is consistent
class AnalysisResultItem(BaseModel):
    # Define fields based on pathway_dominance and module_context structure
    # Example fields (adjust based on actual core.py output):
    ligand: str
    receptor: str
    score: Optional[float] = None
    ligand_norm_expr: Optional[float] = None
    receptor_avg_norm_expr: Optional[float] = None
    interaction_type: Optional[str] = None
    ligand_module: Optional[str] = None
    receptor_modules: Optional[List[str]] = None
    is_same_module: Optional[bool] = None

class CustomAnalysisResponse(BaseModel):
    pathway_dominance: List[AnalysisResultItem] = Field(..., description="Pathway dominance results for the lasso region")
    module_context: List[AnalysisResultItem] = Field(..., description="Module context results for the lasso region") 