from pydantic import BaseModel, Field
from typing import Dict, Optional, List, Any, Literal
import uuid

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

# --- ADDED: Models for Comparison Analysis --- #

class FileSet(BaseModel):
    spatialFileId: Optional[str] = None
    interactionsFileId: Optional[str] = None
    modulesFileId: Optional[str] = None

class MappingSet(BaseModel):
    # Use the existing AnalysisMapping model for structure consistency
    spatialMapping: Optional[AnalysisMapping] = None 
    interactionsMapping: Optional[AnalysisMapping] = None
    modulesMapping: Optional[AnalysisMapping] = None

class SelectionDefinition(BaseModel):
    layer_name: Optional[str] = None
    polygon_coords: Optional[List[List[float]]] = None # List of [x, y] coordinates

class SelectionData(BaseModel):
    source_job_id: Optional[str] = None # Original job ID, if relevant
    files: FileSet
    type: Literal['whole_tissue', 'layer', 'lasso'] 
    definition: SelectionDefinition
    mappings: MappingSet

class ComparisonRequest(BaseModel):
    comparison_name: Optional[str] = Field(None, description="Optional user-defined name for the comparison")
    selection1: SelectionData
    selection2: SelectionData
    # Define parameters for the analysis - simplified to FDR threshold as per user request
    fdr_threshold: float = Field(0.05, ge=0, le=1, description="False Discovery Rate threshold for significance")
    # analyses_to_perform: List[Dict[str, Any]] # Removed - using simplified approach

# Define structures for comparison results (adjust as needed)
class DifferentialExpressionResult(BaseModel):
    molecule_id: str
    # UPDATED: More specific molecule types, including new ones
    type: Literal[
        'single_ligand', 
        'single_receptor', 
        'complex_ligand', 
        'complex_receptor', 
        'ligand_receptor_pair', 
        'ligand_and_receptor', # NEW
        'unknown' # NEW (replaces 'other' for individual molecules)
    ]
    log2_fold_change: Optional[float] = None # Log2 Fold Change (Selection2 / Selection1)
    p_value: Optional[float] = None # Raw p-value from test
    q_value: Optional[float] = None # FDR adjusted p-value
    mean_selection1: Optional[float] = None # Mean normalized value/count in selection 1
    mean_selection2: Optional[float] = None # Mean normalized value/count in selection 2
    # Add flags for complex? depends on implementation
    # ADDED: Fields for L-R pair components
    ligand_id: Optional[str] = None
    receptor_id: Optional[str] = None

class ComparisonResults(BaseModel):
    differential_expression: List[DifferentialExpressionResult] = []
    # Add other result types here if needed later

class ComparisonResponse(BaseModel):
    comparison_id: str = Field(default_factory=lambda: f"comp_{uuid.uuid4()}", description="Unique ID for this comparison run")
    results: ComparisonResults
    errors: List[str] = []

# --- END: Models for Comparison Analysis --- #

# --- ADDED: Model for initiating comparison job --- #
class ComparisonJobResponse(BaseModel):
    status: str = Field(..., description="Status of the comparison job initiation")
    message: str = Field(..., description="A message detailing the status or next steps")
    job_id: Optional[str] = Field(None, description="Job ID for tracking the asynchronous comparison task")
# --- END: Model for initiating comparison job --- # 