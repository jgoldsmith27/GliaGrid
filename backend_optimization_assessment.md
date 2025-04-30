# GliaGrid Backend Architecture Assessment

## Current Architecture Overview

The backend appears well-organized with a clean separation of concerns:

1. **API Layer** (`backend/app/api/`):
   - `analysis_routes.py`: Handles analysis job submission and status
   - `visualization.py`: Serves visualization data
   - `file_routes.py`: Manages file uploads and processing

2. **Service Layer** (`backend/app/services/`):
   - `analysis_service.py`: Coordinates the analysis workflow
   - `job_service.py`: Manages job state (in-memory)
   - `file_service.py`: Handles file storage and processing

3. **Core Analysis** (`backend/app/analysis_logic/`):
   - `core.py`: Contains the computational pipeline for data analysis

4. **Data Models** (`backend/app/models/`):
   - `analysis_models.py`: Defines data structures
   - `file_data.py`: Defines file-related data structures

## Key Bottlenecks Identified

1. **Memory Management**:
   - The entire spatial dataset (potentially 9M rows) is loaded into memory
   - No chunking/streaming for large files
   - No spatial indexing for efficient area selection

2. **Redundant Processing**:
   - No caching mechanism for processed data
   - Reloading and reprocessing files on each request
   - Lack of pre-computed views for visualization

3. **State Management**:
   - In-memory job state (lost on server restart)
   - No client-side persistence

4. **Visualization Data Flow**:
   - Inefficient fetching of all points for visualization
   - No progressive loading or level-of-detail

## Architectural Recommendations

### 1. Data Processing Architecture

#### Current Flow:
```
Frontend -> API -> AnalysisService -> Load full files -> Process all data -> Return results
```

#### Recommended Flow:
```
Frontend -> Check local cache -> If missing -> API request with specific needs -> 
Backend serves chunked/indexed data -> Frontend caches and processes
```

### 2. Modular Service Improvements

#### A. Enhanced FileService

```python
class EnhancedFileService:
    """Handles file operations with performance optimizations."""
    
    @classmethod
    def get_chunked_data(cls, file_id: str, offset: int, limit: int) -> pd.DataFrame:
        """Reads a specific chunk of data from a file."""
        pass
        
    @classmethod
    def get_spatial_points_in_region(cls, file_id: str, bounding_box: dict) -> pd.DataFrame:
        """Efficiently retrieves only points within a specific region."""
        pass
        
    @classmethod
    def precompute_visualization_data(cls, file_id: str) -> Dict[str, Any]:
        """Precomputes data needed for visualization during initial analysis."""
        pass
        
    @classmethod
    def create_spatial_index(cls, file_id: str) -> None:
        """Creates an R-tree spatial index for the file to enable efficient spatial queries."""
        pass
```

#### B. Persistent JobService

```python
class PersistentJobService:
    """Manages job state with persistence."""
    
    def __init__(self, storage_path: str):
        """Initialize with a path for persistent storage."""
        self.storage_path = storage_path
        self._load_state()
        
    def _load_state(self) -> None:
        """Load job state from disk."""
        pass
        
    def _save_state(self) -> None:
        """Save job state to disk."""
        pass
    
    def create_precomputed_assets(self, job_id: str) -> Dict[str, Any]:
        """Generate assets for client-side storage."""
        pass
```

#### C. DataProcessingService

```python
class DataProcessingService:
    """Handles data transformation and processing."""
    
    @staticmethod
    def create_multi_resolution_data(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """Creates multiple resolution samples of the dataset."""
        pass
        
    @staticmethod
    def run_algorithm_on_subset(points: List[Dict], algorithm: str) -> Any:
        """Runs a specific algorithm on a subset of points."""
        pass
```

### 3. API Endpoint Improvements

#### A. Enhanced Points Endpoint

```python
@router.get("/points/{job_id}")
async def get_points(
    job_id: str,
    resolution: float = Query(1.0, description="Data resolution (0.1 to 1.0)"),
    region: str = Query(None, description="JSON bounding box for spatial filtering"),
    limit: int = Query(1000, description="Maximum points to return")
):
    """Get spatial points with resolution control and region filtering."""
```

#### B. Progressive Loading API

```python
@router.get("/points/{job_id}/progressive")
async def get_progressive_points(
    job_id: str,
    zoom_level: int = Query(..., description="Zoom level (1-10)"),
    viewport: str = Query(..., description="JSON viewport coordinates")
):
    """Returns appropriate detail level based on zoom."""
```

#### C. Custom Selection Analysis API

```python
@router.post("/analysis/custom/{job_id}")
async def analyze_custom_selection(
    job_id: str,
    selection: CustomSelectionRequest,
    background_tasks: BackgroundTasks,
    compute_location: str = Query("server", description="Where to perform computation (server/client)")
):
    """Analyze a custom selection with location preference."""
```

### 4. Client-Side Processing with Electron

#### A. Direct File Access

```typescript
// electron/main.js
ipcMain.handle('read-spatial-file-raw', async (event, filePath, options) => {
  // Read data directly from disk in Node.js process
  // Return only what's needed based on options
})
```

#### B. Background Processing

```typescript
// electron/background.js
function processDataInBackground(data, algorithm) {
  // Heavy computation in Node.js process
  // Return results via IPC
}
```

#### C. Client-Side Caching

```typescript
// frontend/src/services/localStorageService.ts
export async function cacheProcessedData(jobId: string, data: any): Promise<void> {
  // Store data in IndexedDB or local file system
}

export async function getProcessedData(jobId: string): Promise<any> {
  // Retrieve cached data if available
}
```

## Implementation Strategy

### Phase 1: Backend Optimization

1. **Immediate Performance Improvements**:
   - Implement chunked file reading in `FileService._load_and_standardize`
   - Add resolution parameter to `/points/{jobId}/all` endpoint
   - Create spatial indexing for point data

2. **Persistent State**:
   - Extend `JobService` with file-based persistence
   - Add pre-computed results storage

### Phase 2: Client-Side Integration

1. **Electron IPC Channels**:
   - Create direct file access methods
   - Implement background processing workers

2. **Frontend Data Management**:
   - Build caching layer for processed data
   - Add state persistence between sessions

### Phase 3: Advanced Visualization

1. **Multi-Resolution Data**:
   - Generate tile pyramids for spatial data
   - Implement level-of-detail rendering

2. **Progressive Loading**:
   - Create APIs for viewport-based data fetching
   - Implement streaming data processing

## Specific Technical Implementations

### 1. Memory-Efficient Point Loading

```python
def load_points_efficient(file_path, chunk_size=10000):
    """Process large spatial files in chunks."""
    for chunk in pd.read_csv(file_path, chunksize=chunk_size):
        # Process each chunk 
        yield chunk

# Usage in endpoints
@router.get("/points/{job_id}/all")
async def get_all_points_data(job_id: str, resolution: float = 1.0):
    """Get points data with controllable resolution."""
    # Implement sampling based on resolution
    if resolution < 1.0:
        sample_size = int(total_points * resolution)
        return random_sample(points, sample_size)
```

### 2. Spatial Indexing with PyRTree

```python
from rtree import index

class SpatialIndex:
    def __init__(self, points_df):
        self.idx = index.Index()
        for i, (_, row) in enumerate(points_df.iterrows()):
            self.idx.insert(i, (row.x, row.y, row.x, row.y))
        self.points_df = points_df
        
    def query_region(self, x_min, y_min, x_max, y_max):
        """Query points within a bounding box."""
        ids = list(self.idx.intersection((x_min, y_min, x_max, y_max)))
        return self.points_df.iloc[ids]
```

### 3. Electron Direct File Access

```typescript
// Add to preload.js
contextBridge.exposeInMainWorld('electron', {
  readSpatialFile: (filePath, options) => ipcRenderer.invoke('read-spatial-file', filePath, options),
  saveCachedResults: (jobId, data) => ipcRenderer.invoke('save-cached-results', jobId, data),
  loadCachedResults: (jobId) => ipcRenderer.invoke('load-cached-results', jobId)
})

// Main process implementation
ipcMain.handle('read-spatial-file', async (_, filePath, options) => {
  const { resolution, region } = options
  // Use efficient Node.js file operations
  // Return appropriate data based on options
})
```

## Priority Queue Approach for Data Loading

One optimal strategy for balancing initial data loading with responsive user interactions is implementing a priority queue for data loading:

### 1. Initial Load Strategy

- Load a low-resolution sample (~5-10%) of all points immediately
- Begin background loading of all points in chunks
- Store loaded chunk IDs to prevent duplication

### 2. Request Prioritization

- When user selects an area or needs specific points, interrupt background loading
- Load requested region at high priority
- Resume background loading where it left off

### 3. Implementation Pattern

```typescript
class DataLoadManager {
  private loadQueue = new PriorityQueue();
  private loadedChunkIds = new Set();
  private isLoading = false;

  // Add region to load queue with priority
  requestRegion(region, priority) {
    this.loadQueue.add({ region, priority });
    if (!this.isLoading) this.processQueue();
  }

  async processQueue() {
    this.isLoading = true;
    while (!this.loadQueue.isEmpty()) {
      const { region, priority } = this.loadQueue.pop();
      const chunks = this.getChunksForRegion(region);
      
      // Skip already loaded chunks
      const chunksToLoad = chunks.filter(c => !this.loadedChunkIds.has(c.id));
      await this.loadChunks(chunksToLoad);
      
      // Mark as loaded
      chunksToLoad.forEach(c => this.loadedChunkIds.add(c.id));
    }
    this.isLoading = false;
  }
}
```

This approach enables responsive UI while progressively building the complete dataset. The initial pipeline could run with a small sample, then continue loading in the background while allowing immediate responses to user interactions.

## Optimizations for Current Implementation

After reviewing the core visualization components (`SpatialOverviewVisualization` and `InteractionVisualization`), I can provide more targeted optimization recommendations:

### 1. Data Storage Architecture

The current implementation has two primary visualization components:
- **SpatialOverviewVisualization**: Shows all spatial points with layer filtering, handles lasso selection
- **InteractionVisualization**: Shows ligand-receptor interactions based on selection

These components need different views of the same underlying spatial data, which creates an opportunity for shared data storage:

```
┌─────────────────┐         ┌───────────────────┐         ┌───────────────────┐
│                 │         │                   │         │                   │
│ Backend (9M pts)│ ───────▶│ Shared Data Store │────────▶│ Custom Selection  │
│                 │         │                   │         │                   │
└─────────────────┘         └───────────────────┘         └───────────────────┘
                                     │                              ▲
                                     │                              │
                                     ▼                              │
                              ┌─────────────┐                ┌─────────────┐
                              │             │                │             │
                              │ Spatial     │                │ Interaction │
                              │ Overview    │────────────────│ Viz         │
                              │             │                │             │
                              └─────────────┘                └─────────────┘
```

### 2. Priority-Based Data Loading for Shared Storage

The priority queue approach can work across both components:

```typescript
class SharedSpatialDataStore {
  // In-memory storage of loaded data
  private loadedPoints: Record<string, AllPointsData[]> = {};
  private loadedChunks: Set<string> = new Set();
  private priorityQueue = new PriorityQueue<DataLoadRequest>();
  
  // Key methods
  async requestAllPoints(jobId: string) {
    // Queue low-priority request for all points
    this.enqueue({
      jobId,
      type: 'all',
      priority: 'low'
    });
    
    // Start with a low-resolution sample
    return this.fetchSample(jobId, 0.1);
  }
  
  async requestInteractionPoints(jobId: string, ligand: string, receptor: string) {
    // High priority - needed for user-initiated interaction view
    this.enqueue({
      jobId,
      type: 'interaction',
      ligand,
      receptor,
      priority: 'high'
    });
    
    // Check if we already have these points cached
    return this.getCachedInteractionPoints(jobId, ligand, receptor) ||
           this.fetchInteractionPoints(jobId, ligand, receptor);
  }
  
  async requestLassoSelection(jobId: string, polygon: Point[]) {
    // Highest priority - immediate user action
    this.enqueue({
      jobId,
      type: 'lasso',
      region: polygon, 
      priority: 'highest'
    });
    
    // Use spatial index if available or fetch from backend
    return this.getSpatialPointsInPolygon(jobId, polygon);
  }
}
```

### 3. Avoiding Redundant Computation

For the custom selection feature, we can store already processed data:

```typescript
// After a lasso selection is made
async function handleLassoSelect(polygon: Point[]) {
  // 1. Get points in selection
  const selectedPoints = await dataStore.requestLassoSelection(jobId, polygon);
  
  // 2. Store selection for future reference
  dataStore.storeSelection('current', selectedPoints);
  
  // 3. When interaction viz needs this data later, it can use:
  const cachedSelection = dataStore.getStoredSelection('current');
  
  // 4. Filter cached selection for specific ligand/receptor
  const { ligandPoints, receptorPoints } = dataStore.filterSelectionForInteraction(
    cachedSelection,
    currentLigand,
    currentReceptor
  );
}
```

### 4. Implementation in Current Components

The `SpatialOverviewVisualization` component would:
1. Request low-resolution sample initially (~10% of points)
2. Start background loading of all points in chunks
3. Immediately render what's available
4. Update as more chunks arrive

The `InteractionVisualization` component would:
1. First check if points for current ligand/receptor are in shared store
2. Fall back to API if not found
3. Store results back in shared store

## FastAPI vs Electron: Optimizing the Architecture Split

GliaGrid has a unique advantage with its Electron-based frontend, which allows for powerful client-side capabilities typically unavailable in browser-only applications. Here's how to optimally balance responsibilities:

### Recommended Architecture Division

| Processing Stage | FastAPI Backend | Electron Frontend |
|------------------|----------------|-------------------|
| Initial file parsing | ✅ (Python's data science ecosystem) | ⚠️ (Only for quick previews) |
| Heavy analysis algorithms | ✅ (Optimized with NumPy/pandas) | ❌ (Keep in Python) |
| Data serving | ✅ (Optimized, compressed) | ❌ (Request from backend) |
| File I/O after initial analysis | ❌ (Network bottleneck) | ✅ (Direct file access) |
| Visualization processing | ❌ (Send minimal data) | ✅ (Transform locally) |
| Selection & filtering | ❌ (Except complex queries) | ✅ (Fast on pre-loaded data) |
| Caching & persistence | ⚠️ (Job results only) | ✅ (All processed data) |

### 1. Direct File Access in Electron

Once initial processing is complete, Electron can directly access the processed result files:

```typescript
// In main process (main.js)
ipcMain.handle('read-processed-file', async (event, jobId, fileType) => {
  // Construct path to processed files based on job ID
  const filePath = path.join(app.getPath('userData'), 'processed', jobId, `${fileType}.parquet`);
  
  // Read directly using Node.js (e.g., with arrow-js for Parquet files)
  const fileData = await parquet.readFile(filePath);
  
  // Return only what's needed based on current view
  return optimizeForCurrentView(fileData);
});
```

### 2. Hybrid Analysis Pipeline

Optimize the analysis pipeline by splitting work:

```
┌───────────────┐   ┌───────────────┐   ┌───────────────────────┐
│ 1. FastAPI    │   │ 2. FastAPI    │   │ 3. Electron (Node.js) │
│ Initial Parse │──▶│ Core Analysis │──▶│ Post-processing &     │
│ & Validation  │   │ Algorithms    │   │ Visualization Prep    │
└───────────────┘   └───────────────┘   └───────────────────────┘
```

This approach:
1. Uses Python's strengths for data science
2. Minimizes data transfer over network
3. Leverages client-side computing power

### 3. Implementation Examples

#### Backend API Optimization

```python
@router.get("/analysis/{job_id}/metadata")
async def get_analysis_metadata(job_id: str):
    """Return lightweight metadata about analysis results"""
    # This includes file paths, summary statistics, and access info
    # Electron can use this to directly access result files
    return {
        "files": {
            "spatial": f"{job_id}/spatial_results.parquet",
            "interactions": f"{job_id}/interaction_results.parquet"
        },
        "summary": {
            "point_count": 9000000,
            "layer_count": 5,
            "interaction_count": 2500
        },
        "access_info": {
            "base_path": f"/path/to/results/{job_id}/",
            "checksum": "sha256:abc123..."  # For validating file integrity
        }
    }
```

#### Electron Direct Data Access

```typescript
// In renderer process (React app)
export async function loadSpatialData(jobId: string, viewportBounds?: BoundingBox) {
  // First check if we have cached this data already
  const cachedData = await localCache.get(`spatial:${jobId}`);
  if (cachedData) return filterForViewport(cachedData, viewportBounds);
  
  // Check if we're in Electron and can access files directly
  if (window.electron) {
    try {
      // Get metadata first (small request)
      const metadata = await api.getAnalysisMetadata(jobId);
      
      // Direct file access through Electron main process
      const data = await window.electron.readProcessedFile(
        jobId, 
        'spatial', 
        viewportBounds
      );
      
      // Cache for future use
      await localCache.set(`spatial:${jobId}`, data);
      return data;
    } catch (err) {
      console.error('Direct file access failed, falling back to API', err);
    }
  }
  
  // Fallback to API if direct access fails or in web mode
  return api.getSpatialData(jobId, viewportBounds);
}
```

### 4. Performance-Critical Functions to Move to Electron

1. **Point filtering & selection**: Move lasso selection processing to client
   ```typescript
   // Much faster in Electron than sending polygons to backend
   function getPointsInPolygon(points, polygon) {
     // Use point-in-polygon npm package
     return points.filter(pt => pointInPolygon([pt.x, pt.y], polygon));
   }
   ```

2. **Spatial data chunking**: Client-managed progressive loading
   ```typescript
   // Main process handles chunking based on current view
   function loadVisibleTilesForViewport(viewport, resolution) {
     const tiles = getTilesForViewport(viewport, resolution);
     return Promise.all(tiles.map(tile => 
       fs.promises.readFile(getTilePath(tile))
     ));
   }
   ```

3. **Visualization-specific transformations**: Transform once, reuse many times
   ```typescript
   // Process once, store results for different visualizations
   function preprocessForAllVisualizations(points) {
     return {
       spatialOverview: prepareSpatialPoints(points),
       layerBoundaries: extractLayerBoundaries(points),
       densityMap: computeDensityGrid(points)
     };
   }
   ```

### 5. Strategic Caching

Implement a two-level caching strategy:

1. **Backend Cache**: Store long-term results
   ```python
   # In FastAPI service
   class AnalysisService:
       def __init__(self):
           self.cache = Cache(ttl=3600*24*7)  # 1 week cache
           
       async def get_analysis_results(self, job_id):
           cached = self.cache.get(job_id)
           if cached:
               return cached
           # Otherwise compute and cache
   ```

2. **Electron Persistent Cache**: Store everything the user has seen
   ```typescript
   // In Electron app
   class PersistentCache {
     constructor() {
       this.dbPath = path.join(app.getPath('userData'), 'cache.db');
       this.db = new SQLite(this.dbPath);
     }
     
     async store(key, data) {
       await this.db.run(
         'INSERT OR REPLACE INTO cache VALUES (?, ?, ?)',
         [key, JSON.stringify(data), Date.now()]
       );
     }
     
     async retrieve(key) {
       const row = await this.db.get(
         'SELECT data FROM cache WHERE key = ?', 
         [key]
       );
       return row ? JSON.parse(row.data) : null;
     }
   }
   ```

### 6. Background Processing with IPC

Utilize Electron's background processes for heavy work:

```typescript
// In main.js
app.whenReady().then(() => {
  const workerPool = new WorkerPool(4); // 4 worker threads
  
  ipcMain.handle('process-data', async (event, jobId, operation, data) => {
    // Offload to worker thread
    return workerPool.submit({
      operation,
      data,
      jobId
    });
  });
});

// In renderer (React)
async function computeCustomSelection() {
  const result = await window.electron.processData(
    currentJobId,
    'customSelectionAnalysis', 
    {points: selectedPoints, parameters: analysisParams}
  );
  
  setAnalysisResults(result);
}
```

## Additional Frontend Optimizations

Some critical frontend-specific optimizations from the previous performance plan that deserve implementation:

### 1. WebGL Rendering Optimizations

```typescript
// In SpatialOverviewVisualization.tsx
const layers = useMemo(() => {
  return [
    new ScatterplotLayer({
      id: 'points-layer',
      data: pointsData,
      // Performance optimizations
      getPosition: d => [d.x, d.y],
      // Use instanced rendering for large point sets
      _instanced: true,
      // Adjust point rendering based on zoom
      radiusScale: Math.max(1, 10 / Math.pow(2, zoom)),
      // Only render points in current viewport
      _filterData: ({startRow, endRow}) => {
        // Filter visible
        return filterVisible(pointsData, startRow, endRow, viewport);
      },
      // Binary data for faster GPU transfer
      getColor: {
        type: 'accessor',
        value: new Uint8Array(colorData.buffer)  // Pre-allocated color array
      }
    })
  ];
}, [pointsData, zoom, viewport]);
```

### 2. Lazy Component Initialization

```typescript
// Only initialize heavy components when visible
function LazyVisualization({ isVisible, ...props }) {
  // Only create the component when visible
  if (!isVisible) return <div className="placeholder">Select data to view</div>;
  
  return <InteractionVisualization {...props} />;
}

// Use with React's Suspense for even better loading behavior
const LazyLoadedVisualization = React.lazy(() => 
  import('./InteractionVisualization')
);
```

### 3. Specialized Binary Formats 

For maximum rendering performance with large datasets, consider specialized data formats:

```typescript
// Convert points to binary format for WebGL
function pointsToBinary(points) {
  // Create a single ArrayBuffer for all points
  // 2 floats (x,y) + 1 uint32 (color) per point
  const buffer = new ArrayBuffer(points.length * 12);
  
  // Create typed array views
  const positions = new Float32Array(buffer, 0, points.length * 2);
  const colors = new Uint32Array(buffer, points.length * 8, points.length);
  
  // Fill arrays
  points.forEach((point, i) => {
    positions[i*2] = point.x;
    positions[i*2+1] = point.y;
    colors[i] = encodeColor(point.layer); // Convert color to packed RGBA
  });
  
  return {
    buffer,
    positions,
    colors
  };
}
```

### 4. Component-Specific Optimizations

#### SpatialOverviewVisualization

```typescript
// In SpatialOverviewVisualization.tsx

// Use multi-resolution approach
useEffect(() => {
  // Start with low resolution
  loadData(jobId, { resolution: 0.1 }).then(initialData => {
    setPoints(initialData);
    
    // Then load higher resolution in background
    loadData(jobId, { resolution: 0.5 }).then(mediumData => {
      if (mounted.current) setPoints(mediumData);
      
      // Finally load full resolution if needed
      if (zoom > 5) {
        loadData(jobId, { resolution: 1.0 }).then(fullData => {
          if (mounted.current) setPoints(fullData);
        });
      }
    });
  });
}, [jobId, zoom]);

// Implement layer boundary caching
const [layerBoundaries, setLayerBoundaries] = useState({});
useEffect(() => {
  if (points.length > 0 && !layerBoundaries[currentLayer]) {
    // Compute convex hull for layer boundaries
    const boundary = computeLayerBoundary(points, currentLayer);
    setLayerBoundaries(prev => ({...prev, [currentLayer]: boundary}));
  }
}, [points, currentLayer]);
```

#### InteractionVisualization

```typescript
// In InteractionVisualization.tsx

// Use viewport-based filtering
const visiblePoints = useMemo(() => {
  return filterPointsInViewport(allPoints, viewport);
}, [allPoints, viewport]);

// Pre-compute density maps
const densityMap = useMemo(() => {
  if (!visiblePoints.length) return null;
  // Generate grid-based density metrics
  return generateDensityGrid(visiblePoints, 50); // 50x50 grid
}, [visiblePoints]);
```

## Conclusion

By implementing these architectural changes, GliaGrid will be able to:

1. **Handle Massive Datasets**: Process 9M+ row datasets efficiently through chunking, indexing, and multi-resolution sampling.

2. **Accelerate Visualization**: Provide responsive visualizations through progressive loading and client-side caching.

3. **Enable Custom Analysis**: Support efficient lasso selection and analysis through spatial indexing and optimized data flow.

4. **Leverage Electron**: Utilize Electron's capabilities for client-side processing and direct file access.

The most critical first step is optimizing the `/points/{jobId}/all` endpoint to efficiently handle large datasets, followed by implementing client-side caching of processed results to avoid redundant work. 