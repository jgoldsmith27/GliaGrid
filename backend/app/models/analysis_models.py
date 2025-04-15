from pydantic import BaseModel, Field
from typing import Optional

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