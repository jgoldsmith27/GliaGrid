# GliaGrid Architecture & Performance Optimization Plan

## Table of Contents
1. [Current Architecture Assessment](#current-architecture-assessment)
2. [Performance Bottlenecks](#performance-bottlenecks)
3. [Optimization Strategy](#optimization-strategy)
4. [Implementation Roadmap](#implementation-roadmap)
   - [Phase 1: Data Access Optimization](#phase-1-data-access-optimization-1-2-weeks)
   - [Phase 2: Shared Data Management](#phase-2-shared-data-management-1-week)
   - [Phase 3: Visualization Enhancements](#phase-3-visualization-enhancements-1-week)
5. [Technical Implementation Details](#technical-implementation-details)
   - [Backend Optimizations](#backend-optimizations)
   - [Electron Integration](#electron-integration)
   - [Frontend Components](#frontend-components)

## Current Architecture Assessment

The backend is organized with a clean separation of concerns:

1. **API Layer** (`backend/app/api/`):
   - `analysis_routes.py`: Handles analysis job submission and status
   - `visualization.py`: Serves visualization data
   - `file_routes.py`: Manages file uploads and processing

2. **Service Layer** (`backend/app/services/`):
   - `analysis_service.py`: Coordinates the analysis workflow
   - `job_service.py`: Manages job state (in-memory)
   - `file_service.py`: Handles file storage and processing (in `.temp_uploads/`)

3. **Core Analysis** (`backend/app/analysis_logic/`):
   - `core.py`: Contains the computational pipeline for data analysis

4. **Frontend Visualization Components**:
   - `SpatialOverviewVisualization`: Shows all spatial points with layer filtering, handles lasso selection
   - `InteractionVisualization`: Shows ligand-receptor interactions based on selection

5. **Electron Integration**:
   - Projects are saved as JSON files in `app.getPath('userData')/projects/`
   - Projects store references to backend file IDs
   - Projects can be loaded/deleted through Electron IPC channels
   - No direct file access from Electron to backend storage

## Performance Bottlenecks

1. **Memory Management**:
   - The entire spatial dataset (9M rows) is loaded into memory
   - No chunking/streaming for large files
   - No spatial indexing for efficient area selection

2. **Redundant Processing**:
   - No caching mechanism for processed data
   - Files reloaded and reprocessed on each request
   - Lack of pre-computed views for visualization

3. **Data Transfer Inefficiency**:
   - All data flows through HTTP API, even in Electron
   - No progressive loading or resolution control
   - No direct file access despite Electron capabilities

4. **Visualization Rendering**:
   - No WebGL optimizations for large point datasets
   - All points rendered at once instead of level-of-detail approach
   - Heavy components not lazily initialized

## Optimization Strategy

The optimization strategy focuses on three key approaches:

1. **Priority Queue Data Loading**:
   - Start with low-resolution samples (~10%) for immediate feedback
   - Queue full dataset loading in background with cancellable requests
   - Prioritize user-requested regions/interactions
   - Track loaded chunks to prevent duplicate loading

2. **Balanced Workload Distribution**:

   | Processing Stage | FastAPI Backend | Electron Frontend |
   |------------------|----------------|-------------------|
   | Initial file parsing | ✅ (Python's data science ecosystem) | ⚠️ (Only for previews) |
   | Heavy analysis algorithms | ✅ (Optimized with NumPy/pandas) | ❌ (Keep in Python) |
   | Data serving | ✅ (Optimized, compressed) | ❌ (Request from backend) |
   | File I/O after initial analysis | ❌ (Network bottleneck) | ✅ (Direct file access) |
   | Visualization processing | ❌ (Send minimal data) | ✅ (Transform locally) |
   | Selection & filtering | ❌ (Except complex queries) | ✅ (Fast on pre-loaded data) |
   | Caching & persistence | ⚠️ (Job results only) | ✅ (All processed data) |

3. **Shared Data Architecture**:

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
                                 │ Spatial     │────────────────│ Interaction │
                                 │ Overview    │                │ Viz         │
                                 │             │                │             │
                                 └─────────────┘                └─────────────┘
   ```

## Implementation Roadmap

### Phase 1: Data Access Optimization (1-2 weeks)

1. **Backend Chunked Processing**:
   - Modify `FileService` to support chunked data retrieval
   - Add spatial indexing for region queries
   - Implement resolution parameter in visualization API endpoints

2. **Direct File Access in Electron**:
   - Create IPC handlers for reading backend files
   - Implement chunked CSV parsing
   - Add resolution/sampling parameters

3. **API Endpoints Enhancement**:
   - Add metadata endpoint to return file information
   - Update points endpoint to support resolution and regions
   - Implement job data file path endpoints

### Phase 2: Shared Data Management (1 week)

1. **Frontend Data Store**:
   - Create SharedDataStore singleton
   - Implement multi-level caching (memory and persistent)
   - Add priority queue for data requests

2. **Electron Background Processing**:
   - Setup worker threads for processing
   - Implement spatial operations in main process
   - Add IPC channels for communication

3. **Data Flow Optimization**:
   - Connect components to shared data store
   - Implement data transformation pipelines
   - Add selection caching

### Phase 3: Visualization Enhancements (1 week)

1. **WebGL Optimizations**:
   - Implement instanced rendering
   - Add level-of-detail based on zoom
   - Use binary formats for GPU transfer

2. **Progressive Loading**:
   - Update components to support resolution switching
   - Implement viewport-based filtering
   - Add loading indicators with progressive enhancement

3. **Component Optimizations**:
   - Add lazy initialization for heavy components
   - Implement shared layer boundary caching
   - Precompute visualization transformations

## Technical Implementation Details

### Backend Optimizations

#### Memory-Efficient Point Loading

```python
def load_points_efficient(file_path, chunk_size=10000, resolution=1.0):
    """Process large spatial files in chunks with resolution control."""
    # Calculate total rows
    total_rows = sum(1 for _ in open(file_path)) - 1  # Subtract header
    
    # Apply resolution sampling if needed
    if resolution < 1.0:
        sample_size = int(total_rows * resolution)
        sample_indices = set(random.sample(range(total_rows), sample_size))
        
    # Process in chunks
    result = []
    for i, chunk in enumerate(pd.read_csv(file_path, chunksize=chunk_size)):
        # Apply resolution filtering if needed
        if resolution < 1.0:
            start_idx = i * chunk_size
            chunk = chunk.iloc[[idx - start_idx for idx in sample_indices 
                               if start_idx <= idx < start_idx + len(chunk)]]
        result.append(chunk)
        
    return pd.concat(result) if result else pd.DataFrame()
```

#### Spatial Indexing

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

#### Enhanced API Endpoints

```python
@router.get("/points/{job_id}")
async def get_points(
    job_id: str,
    resolution: float = Query(1.0, description="Data resolution (0.1 to 1.0)"),
    region: str = Query(None, description="JSON bounding box for spatial filtering"),
    limit: int = Query(1000, description="Maximum points to return")
):
    """Get spatial points with resolution control and region filtering."""
    # Implementation using the optimized data loading methods
```

### Electron Integration

#### Direct File Access

```javascript
// In main.js
ipcMain.handle('read-backend-file', async (event, fileId, options = {}) => {
  // Construct path to temporary uploads directory
  const backendDir = path.join(app.getPath('userData'), '..', '.temp_uploads');
  
  try {
    // Find the file with this ID
    const files = await fsPromises.readdir(backendDir);
    const matchingFile = files.find(file => file.startsWith(fileId));
    
    if (!matchingFile) {
      throw new Error(`File with ID ${fileId} not found`);
    }
    
    const filePath = path.join(backendDir, matchingFile);
    
    // Handle different file types
    const ext = path.extname(filePath).toLowerCase();
    if (ext === '.csv') {
      // Read and parse CSV with chunking
      return readCsvChunked(filePath, options);
    } else if (ext === '.h5ad') {
      // For H5AD, we might need to use Python through child_process
      throw new Error('H5AD direct reading not implemented yet');
    }
    
    throw new Error(`Unsupported file type: ${ext}`);
  } catch (error) {
    console.error('Error reading backend file:', error);
    throw error;
  }
});

// Helper function to read CSV in chunks
async function readCsvChunked(filePath, options = {}) {
  const { offset = 0, limit = 1000, sampleRate = 1.0 } = options;
  
  return new Promise((resolve, reject) => {
    const results = [];
    let rowCount = 0;
    let skippedCount = 0;
    
    Papa.parse(fs.createReadStream(filePath), {
      header: true,
      skipEmptyLines: true,
      step: (row) => {
        rowCount++;
        
        // Apply sampling and offset/limit
        if (rowCount <= offset) return;
        if (limit && results.length >= limit) return;
        
        // Apply sampling
        if (sampleRate < 1.0 && Math.random() > sampleRate) {
          skippedCount++;
          return;
        }
        
        results.push(row.data);
      },
      complete: () => {
        resolve({
          data: results,
          totalRows: rowCount,
          sampledRows: results.length,
          skippedRows: skippedCount
        });
      },
      error: (error) => reject(error)
    });
  });
}
```

#### Preload Configuration

```javascript
// In preload.js
contextBridge.exposeInMainWorld('electron', {
  // Project management
  saveProject: (projectName, projectData) => ipcRenderer.invoke('save-project', projectName, projectData),
  listProjects: () => ipcRenderer.invoke('list-projects'),
  loadProject: (filePath) => ipcRenderer.invoke('load-project', filePath),
  deleteProject: (filePath) => ipcRenderer.invoke('delete-project', filePath),
  
  // Direct file access
  readBackendFile: (fileId, options) => ipcRenderer.invoke('read-backend-file', fileId, options),
  
  // Background processing
  processData: (jobId, operation, data) => ipcRenderer.invoke('process-data', jobId, operation, data)
});
```

### Frontend Components

#### SharedDataStore

```typescript
// In frontend/src/services/data/SharedDataStore.ts
export class SharedDataStore {
  private static instance: SharedDataStore;
  private cache: Map<string, any> = new Map();
  private loadingPromises: Map<string, Promise<any>> = new Map();
  private priorityQueue = new PriorityQueue<DataLoadRequest>();
  private isProcessing = false;
  
  // Singleton pattern
  static getInstance() {
    if (!this.instance) {
      this.instance = new SharedDataStore();
    }
    return this.instance;
  }
  
  // Priority queue methods
  enqueue(request: DataLoadRequest) {
    this.priorityQueue.add(request);
    if (!this.isProcessing) {
      this.processQueue();
    }
  }
  
  async processQueue() {
    this.isProcessing = true;
    while (!this.priorityQueue.isEmpty()) {
      const request = this.priorityQueue.pop();
      await this.handleRequest(request);
    }
    this.isProcessing = false;
  }
  
  // Data access methods
  async getChunk(jobId: string, fileType: string, options: { 
    offset?: number, 
    limit?: number, 
    resolution?: number,
    priority?: 'high' | 'medium' | 'low'
  } = {}) {
    const { priority = 'medium', ...fetchOptions } = options;
    const cacheKey = `${jobId}:${fileType}:${JSON.stringify(fetchOptions)}`;
    
    // Check cache first
    if (this.cache.has(cacheKey)) {
      return this.cache.get(cacheKey);
    }
    
    // Check if already loading
    if (this.loadingPromises.has(cacheKey)) {
      return this.loadingPromises.get(cacheKey);
    }
    
    // Add to priority queue
    const request = { 
      jobId, 
      fileType, 
      options: fetchOptions, 
      priority, 
      cacheKey 
    };
    
    this.enqueue(request);
    
    // Create a promise that will resolve when the request is processed
    let resolvePromise, rejectPromise;
    const promise = new Promise((resolve, reject) => {
      resolvePromise = resolve;
      rejectPromise = reject;
    });
    
    this.loadingPromises.set(cacheKey, promise);
    request.resolve = resolvePromise;
    request.reject = rejectPromise;
    
    return promise;
  }
  
  // Data loading implementation
  private async handleRequest(request: DataLoadRequest) {
    const { jobId, fileType, options, cacheKey, resolve, reject } = request;
    
    try {
      // Try direct file access if in Electron
      let result;
      if (window.electron?.readBackendFile) {
        try {
          // Get fileId from jobId mapping
          const fileId = await this.getFileIdFromJobId(jobId, fileType);
          if (fileId) {
            result = await window.electron.readBackendFile(fileId, options);
          }
        } catch (err) {
          console.warn('Direct file access failed, falling back to API', err);
        }
      }
      
      // Fallback to API if direct access failed or not available
      if (!result) {
        result = await api.getDataChunk(jobId, fileType, options);
      }
      
      // Cache the result
      this.cache.set(cacheKey, result);
      resolve(result);
      
    } catch (error) {
      reject(error);
    } finally {
      this.loadingPromises.delete(cacheKey);
    }
  }
  
  // Helper method to map job IDs to file IDs
  private async getFileIdFromJobId(jobId: string, fileType: string) {
    // This would need to be implemented based on your app's data structure
    // Could request from an API endpoint or store mapping in local storage
    // For now, returning a dummy implementation
    return jobId + '_' + fileType;
  }
}

// Priority queue implementation
class PriorityQueue<T extends { priority: string }> {
  private items: T[] = [];
  
  add(item: T) {
    // Convert priority string to numeric value
    const priorityValue = 
      item.priority === 'high' ? 3 :
      item.priority === 'medium' ? 2 : 1;
    
    // Add item with numeric priority
    this.items.push({...item, numericPriority: priorityValue});
    
    // Sort by priority (higher values first)
    this.items.sort((a, b) => b.numericPriority - a.numericPriority);
  }
  
  pop(): T | undefined {
    return this.items.shift();
  }
  
  isEmpty(): boolean {
    return this.items.length === 0;
  }
}
```

#### SpatialOverviewVisualization Component

```tsx
// In SpatialOverviewVisualization.tsx
function SpatialOverviewVisualization({ jobId }) {
  const [points, setPoints] = useState([]);
  const [resolution, setResolution] = useState(0.1); // Start with 10%
  const [isLoading, setIsLoading] = useState(true);
  const dataStore = useSharedDataStore();
  
  // Progressive loading based on zoom level
  useEffect(() => {
    async function loadData() {
      setIsLoading(true);
      try {
        // Initial low-res sample
        const initialData = await dataStore.getChunk(jobId, 'spatial', { 
          resolution,
          priority: 'high' 
        });
        setPoints(initialData);
        
        // If we're at a high zoom level, load more detailed data
        if (resolution < 0.5) {
          // After showing initial data, load more detailed view
          const detailedData = await dataStore.getChunk(jobId, 'spatial', { 
            resolution: 0.5,
            priority: 'medium'
          });
          setPoints(detailedData);
        }
      } catch (err) {
        console.error('Failed to load data:', err);
      } finally {
        setIsLoading(false);
      }
    }
    
    loadData();
  }, [jobId, resolution]);
  
  // WebGL optimized layers
  const layers = useMemo(() => {
    return [
      new ScatterplotLayer({
        id: 'points-layer',
        data: points,
        // Performance optimizations
        getPosition: d => [d.x, d.y],
        // Use instanced rendering for large point sets
        _instanced: true,
        // Adjust point rendering based on zoom
        radiusScale: Math.max(1, 10 / Math.pow(2, zoom)),
        // Only render points in current viewport
        _filterData: ({startRow, endRow}) => {
          // Filter visible
          return filterVisible(points, startRow, endRow, viewport);
        },
        // Binary data for faster GPU transfer
        getColor: {
          type: 'accessor',
          value: new Uint8Array(colorData.buffer)  // Pre-allocated color array
        }
      })
    ];
  }, [points, zoom, viewport]);
  
  // Rest of component...
}
```
