from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Any

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
    polygon: List[List[float]] = Field(..., description="List of [x, y] coordinates defining the lasso polygon. First and last points must be the same.")

# Represents a single item in the results list for pathway/module analysis
class AnalysisResultItem(BaseModel):
    # Common fields likely present in both pathway_dominance and module_context results
    # Using Optional[...] and default=None allows flexibility if fields aren't always present
    ligand: Optional[str] = None
    receptor: Optional[str] = None
    # Pathway specific fields
    score: Optional[float] = None 
    ligand_norm_expr: Optional[float] = None
    receptor_avg_norm_expr: Optional[float] = None
    # Module context specific fields
    ligand_module: Optional[Any] = None # Can be str, int, float, null
    receptor_modules: Optional[Any] = None # Can be list, str, int, float, null
    interaction_type: Optional[str] = None

    # Allow extra fields to accommodate variations in analysis results
    class Config:
        extra = 'allow'

# --- REPLACED CustomAnalysisResponse and LayeredCustomAnalysisResponse ---

# Re-purposed: Represents results for a single scope (whole custom OR one layer)
class CustomAnalysisScopeResult(BaseModel):
    pathway_dominance: List[AnalysisResultItem] = Field(default_factory=list)
    module_context: List[AnalysisResultItem] = Field(default_factory=list)

# NEW: Response bundle containing both whole and layered results
class CustomAnalysisResultsBundle(BaseModel):
    # Results aggregated for the whole custom selection
    whole_results: CustomAnalysisScopeResult = Field(...)
    # Dictionary where keys are layer names (str)
    # and values are the results for that layer within the lasso.
    # Optional: Will be empty if no layers found or layer analysis fails.
    layered_results: Dict[str, CustomAnalysisScopeResult] = Field(default_factory=dict) 