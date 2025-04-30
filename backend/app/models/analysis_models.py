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
class PointData(BaseModel):
    x: float
    y: float
    gene: str
    # layer: Optional[str] = None # Layer might not be relevant from lasso

class CustomAnalysisRequest(BaseModel):
    ligands: List[PointData]
    receptors: List[PointData]
    # Potentially add original jobId if needed to access other files like interactions/modules?
    # original_job_id: Optional[str] = None 