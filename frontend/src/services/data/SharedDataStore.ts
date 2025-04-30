/**
 * SharedDataStore.ts
 * 
 * A centralized data store for spatial data with:
 * - Caching
 * - Priority queue for loading
 * - Progressive resolution control
 * - Support for direct file access in Electron
 */

// Define data request priority levels
type Priority = 'high' | 'medium' | 'low';

// Options for data retrieval
interface DataRequestOptions {
  offset?: number;
  limit?: number;
  resolution?: number;
  region?: {
    x_min: number;
    y_min: number;
    x_max: number;
    y_max: number;
  };
  x_col?: string; // Needed if fetching specific columns
  y_col?: string; // Needed if fetching specific columns
  geneCol?: string; // Needed for filtering
  layerCol?: string; // Needed for potential layer filtering
  priority?: Priority;
  // Add properties specific to visualization requests
  ligand?: string; 
  receptor?: string;
  layer?: string; // Used by visualization logic
  // Filter/Columns options for readBackendFile
  filter?: { column: string; values: string[] };
  columns?: string[];
  sampleRate?: number;
}

// Data request type for the priority queue
interface DataLoadRequest {
  jobId: string;
  fileType: string;
  options: Omit<DataRequestOptions, 'priority'>;
  priority: Priority;
  cacheKey: string;
  resolve?: (value: any) => void;
  reject?: (reason: any) => void;
  numericPriority?: number;
}

// --- Define Job Status Event/Error types --- 
interface JobStatusEvent {
    jobId: string;
    job_id?: string; // Sometimes backend uses this key
    status: string;
    message?: string;
    progress?: number;
    results?: any;
    // Add other potential fields if the backend sends them
}

interface JobStatusError {
    jobId: string;
    error: string;
}

// Simple in-memory priority queue
class PriorityQueue<T extends { priority: Priority }> {
  private items: (T & { numericPriority: number })[] = [];
  
  add(item: T): void {
    const numericPriority = 
      item.priority === 'high' ? 3 :
      item.priority === 'medium' ? 2 : 1;
    this.items.push({...item, numericPriority});
    this.items.sort((a, b) => b.numericPriority - a.numericPriority);
    console.log(`[SharedDataStore] Added item to queue with priority ${item.priority} (${numericPriority}). Queue size: ${this.items.length}`);
  }
  
  pop(): T | undefined {
    return this.items.shift();
  }
  
  isEmpty(): boolean {
    return this.items.length === 0;
  }

  size(): number {
    return this.items.length;
  }
  
  peek(): T | undefined {
    return this.items[0];
  }
}

/**
 * Shared Data Store - Singleton that manages data loading and caching
 */
export class SharedDataStore {
  private static instance: SharedDataStore;
  
  private cache: Map<string, any> = new Map();
  private loadingPromises: Map<string, Promise<any>> = new Map();
  private priorityQueue = new PriorityQueue<DataLoadRequest>();
  private isProcessing = false;
  private jobMetadataCache: Map<string, any> = new Map();
  private jobFileMapping: Map<string, Record<string, string>> = new Map(); // Still needed for getFileIdFromJobId

  private constructor() {
    console.log('[SharedDataStore] Created new instance');
  }
  
  static getInstance(): SharedDataStore {
    if (!this.instance) {
      this.instance = new SharedDataStore();
    }
    return this.instance;
  }
  
  async getChunk(
    jobId: string, 
    fileType: string, 
    options: DataRequestOptions = {}
  ): Promise<any> {
    if (!window.electronAPI) {
      throw new Error('Electron API is not available for data fetching.');
    }

    const { priority = 'medium', ...fetchOptions } = options;
    // Use the full options (including priority for now) to create cache key if priority affects fetch result
    // Or ensure cache key only uses fetchOptions if priority doesn't change the data fetched
    const cacheKey = this.createCacheKey(jobId, fileType, fetchOptions);
    
    if (this.cache.has(cacheKey)) {
      console.log(`[SharedDataStore] Cache hit for ${cacheKey}`);
      return this.cache.get(cacheKey);
    }
    
    if (this.loadingPromises.has(cacheKey)) {
      console.log(`[SharedDataStore] Already loading ${cacheKey}, returning existing promise`);
      return this.loadingPromises.get(cacheKey);
    }
    
    console.log(`[SharedDataStore] Creating new request for ${cacheKey} with priority ${priority}`);
    
    let resolvePromise: (value: any) => void;
    let rejectPromise: (reason: any) => void;
    
    const promise = new Promise((resolve, reject) => {
      resolvePromise = resolve;
      rejectPromise = reject;
    });
    
    const request: DataLoadRequest = {
      jobId,
      fileType,
      options: fetchOptions, // Pass fetchOptions (without priority)
      priority,
      cacheKey,
      resolve: resolvePromise!,
      reject: rejectPromise!
    };
    
    this.loadingPromises.set(cacheKey, promise);
    this.enqueue(request);
    
    return promise;
  }
  
  private enqueue(request: DataLoadRequest): void {
    this.priorityQueue.add(request);
    if (!this.isProcessing) {
      console.log('[SharedDataStore] Starting queue processing');
      this.processQueue();
    }
  }
  
  private async processQueue(): Promise<void> {
    this.isProcessing = true;
    while (!this.priorityQueue.isEmpty()) {
      const request = this.priorityQueue.pop();
      if (!request) break;
      try {
        console.log(`[SharedDataStore] Processing request: ${request.cacheKey} (priority: ${request.priority})`);
        await this.handleRequest(request);
      } catch (error) {
        console.error(`[SharedDataStore] Error processing request ${request.cacheKey}:`, error);
        // Ensure promise is rejected if handleRequest throws outside its own try/catch
        if (request.reject) {
           request.reject(error); 
        }
      } finally {
         // Clean up promise only after handling or error
         this.loadingPromises.delete(request.cacheKey);
      }
    }
    this.isProcessing = false;
    console.log('[SharedDataStore] Queue processing complete');
  }
  
  private async handleRequest(request: DataLoadRequest): Promise<void> {
    const { jobId, fileType, options, cacheKey, resolve, reject } = request;
    
    if (!resolve || !reject) {
      console.error('[SharedDataStore] Internal error: Missing resolve/reject in request', request);
      return;
    }
    
    try {
      if (!window.electronAPI) {
        throw new Error('Electron API is not available.');
      }

      // Check if this is a visualization request based on presence of ligand/receptor options
      if (options.ligand && options.receptor) {
        console.log(`[SharedDataStore] Handling L/R visualization request for job ${jobId}`, options);

        // --- Visualization Logic (fetches SPATIAL data for L/R) --- 
        // 1. Get Job Metadata (needs spatial file info)
        // Metadata is now structured as { inputs: { files: ..., mappings: ... }, outputs: ... }
        const metadata = await this.getJobMetadata(jobId);
        
        // Extract inputs from the metadata structure
        const inputs = metadata?.inputs;
        const spatialFileId = inputs?.files?.spatialFileId;
        const mapping = inputs?.mappings?.spatialMapping;
        
        // Validate necessary metadata exists within inputs
        if (!spatialFileId || !mapping || !mapping.geneCol || !mapping.xCol || !mapping.yCol) {
          const missing = [];
          if (!inputs) missing.push('inputs object');
          else {
              if (!inputs.files) missing.push('inputs.files');
              else if (!spatialFileId) missing.push('inputs.files.spatialFileId');
              
              if (!inputs.mappings) missing.push('inputs.mappings');
              else if (!mapping) missing.push('inputs.mappings.spatialMapping');
              else {
                 if (!mapping.geneCol) missing.push('inputs.mappings.spatialMapping.geneCol');
                 if (!mapping.xCol) missing.push('inputs.mappings.spatialMapping.xCol');
                 if (!mapping.yCol) missing.push('inputs.mappings.spatialMapping.yCol');
              }
          }
          // Improved error message indicating where the data was expected
          throw new Error(`Missing required input metadata for L/R visualization for job ${jobId}: ${missing.join(', ')}. Received metadata structure: ${JSON.stringify(metadata)}`);
        }
        
        // 2. Determine Genes and Columns needed from options
        const ligandName = options.ligand; 
        const receptorName = options.receptor;
        // This check might be redundant due to the main 'if' condition, but kept for clarity
        if (!ligandName || !receptorName) { 
          throw new Error('Ligand and Receptor names are required for visualization request.');
        }
        const genesToFetch = Array.from(new Set([ligandName, receptorName]));
        const columnsToFetch = [mapping.xCol, mapping.yCol, mapping.geneCol];
        
        // Add layer column if layer filtering is requested
        const layerFilter = options.layer && options.layer !== 'whole_tissue' && mapping.layerCol;
        if (layerFilter && mapping.layerCol) { // Ensure mapping.layerCol exists before pushing
            columnsToFetch.push(mapping.layerCol);
        } else if (layerFilter && !mapping.layerCol) {
             console.warn(`[SharedDataStore] Layer filtering requested for '${options.layer}', but no layerCol defined in spatialMapping.`);
        }
        
        console.log(`[SharedDataStore] Fetching spatial data: File ID: ${spatialFileId}, Filter: ${mapping.geneCol}=[${genesToFetch.join(',')}, Columns: ${columnsToFetch.join(',')}`);

        // 3. Fetch Filtered SPATIAL Data using readBackendFile
        if (!window.electronAPI.readBackendFile) {
            throw new Error('Electron API function readBackendFile is not available.');
        }
        const spatialDataResult = await window.electronAPI.readBackendFile(spatialFileId, {
           filter: {
               column: mapping.geneCol,
               values: genesToFetch
           },
           columns: columnsToFetch
           // Consider adding limit/offset/sampleRate from original options if needed for large results?
           // limit: options.limit,
           // offset: options.offset,
           // sampleRate: options.sampleRate 
        });
        
        console.log(`[SharedDataStore] Received ${spatialDataResult?.returnedRows || 0} points from readBackendFile for L/R pair`);
        
        // 4. Process fetched spatial data into ligand/receptor points
        const ligandPoints: { x: number; y: number }[] = [];
        const receptorPoints: { x: number; y: number }[] = [];
        
        if (spatialDataResult?.data) {
            for (const point of spatialDataResult.data) {
                // Apply layer filtering post-fetch if requested and layerCol exists
                if (layerFilter && mapping.layerCol && point[mapping.layerCol] !== options.layer) {
                    continue;
                }
                
                const gene = point[mapping.geneCol];
                const x = parseFloat(point[mapping.xCol]);
                const y = parseFloat(point[mapping.yCol]);

                if (isNaN(x) || isNaN(y)) {
                     console.warn(`[SharedDataStore] Skipping point with invalid coordinates:`, point);
                     continue;
                }

                if (gene === ligandName) {
                    ligandPoints.push({ x, y });
                }
                if (gene === receptorName) { 
                    receptorPoints.push({ x, y });
                }
            }
        }
        
        const visualizationData = {
            ligand: ligandPoints,
            receptor: receptorPoints,
            // Include any warnings from spatialDataResult if applicable
            warnings: spatialDataResult?.warnings || [] 
        };
        
        console.log(`[SharedDataStore] Processed points for ${options.layer || 'all layers'} - Ligands: ${ligandPoints.length}, Receptors: ${receptorPoints.length}`);

        // 5. Cache and Resolve with the processed visualization data
        this.cache.set(cacheKey, visualizationData);
        resolve(visualizationData);

      } else { 
          // --- Logic for standard file types (spatial, interactions, modules, etc.) --- 
          console.log(`[SharedDataStore] Handling standard file request: ${fileType} for job ${jobId}`);
          
          if (!window.electronAPI.readBackendFile) {
              throw new Error('Electron API function readBackendFile is not available.');
          }
          
          // Get the specific fileId for the requested fileType from metadata
          const fileId = await this.getFileIdFromJobId(jobId, fileType);
          if (!fileId) {
            // Throw error if the specific fileId (e.g., interactionsFileId) is missing in metadata
            throw new Error(`No file ID found for file type '${fileType}' in metadata for job ${jobId}`);
          }
          
          console.log(`[SharedDataStore] Found fileId ${fileId} for job ${jobId}, type ${fileType}. Options:`, options);
          
          // Read the specific file using its fileId and pass options like filter/columns
          const result = await window.electronAPI.readBackendFile(fileId, options);
          console.log(`[SharedDataStore] Direct file access successful for ${fileType}, got ${result?.returnedRows || 0} rows`);
          
          // Cache and resolve with the raw data array from the file read result
          // Assuming the direct file read returns { data: [], returnedRows: ..., totalRows: ... }
          this.cache.set(cacheKey, result.data); 
          resolve(result.data);
      }

    } catch (error) {
      console.error(`[SharedDataStore] Error handling request ${cacheKey}:`, error);
      reject(error); // Reject the promise for the original caller
    } 
    // Removed finally block, cleanup is handled in processQueue
  }

  private async getJobMetadata(jobId: string): Promise<any> {
      if (this.jobMetadataCache.has(jobId)) {
          console.log(`[SharedDataStore] Metadata cache hit for job ${jobId}`);
          return this.jobMetadataCache.get(jobId);
      }
      
      if (!window.electronAPI?.getJobMetadata) {
         throw new Error('Electron API function getJobMetadata is not available.');
      }
      
      console.log(`[SharedDataStore] Fetching metadata for job ${jobId}...`);
      const result = await window.electronAPI.getJobMetadata(jobId);
      
      if (!result.success || !result.metadata) {
          throw new Error(result.error || `Failed to fetch metadata for job ${jobId}`);
      }
      
      // Ensure metadata.files and metadata.mappings exist, even if empty, for safety
      result.metadata.files = result.metadata.files || {};
      result.metadata.mappings = result.metadata.mappings || {};

      this.jobMetadataCache.set(jobId, result.metadata);
      console.log(`[SharedDataStore] Metadata fetched and cached for job ${jobId}`);
      return result.metadata;
  }
  
  // --- Helper method to map job IDs to file IDs --- 
  private async getFileIdFromJobId(jobId: string, fileType: string): Promise<string | null> {
      console.warn(`[SharedDataStore] getFileIdFromJobId: Looking for file ID for type '${fileType}'.`);
      try {
          const metadata = await this.getJobMetadata(jobId);
          
          // Determine the key name in metadata.files based on fileType
          let fileIdKey: string | null = null;
          if (fileType === 'spatial') {
              fileIdKey = 'spatialFileId';
          } else if (fileType === 'interactions') {
              fileIdKey = 'interactionsFileId';
          } else if (fileType === 'modules') {
              fileIdKey = 'modulesFileId';
          }
          // Add other file types as needed (e.g., 'annotations', 'celltypes')
          
          if (!fileIdKey) {
               console.error(`[SharedDataStore] Unknown fileType '${fileType}' requested.`);
               return null;
          }

          // Check if the key exists in metadata.files
          if (metadata?.files?.[fileIdKey]) {
              return metadata.files[fileIdKey];
          }
          
          // Log specific error if the expected key is missing
          console.error(`[SharedDataStore] Could not find fileId key '${fileIdKey}' for type '${fileType}' in metadata.files for job ${jobId}. Metadata files:`, metadata?.files);
          return null;
          
      } catch (error) {
          console.error(`[SharedDataStore] Error getting metadata to find fileId for job ${jobId}, type ${fileType}:`, error);
          return null;
      }
  }

  private createCacheKey(jobId: string, fileType: string, options: Omit<DataRequestOptions, 'priority'>): string {
    // Ensure consistent key generation by sorting option keys or using stable stringify
    const optionsString = JSON.stringify(options, Object.keys(options).sort());
    return `${jobId}:${fileType}:${optionsString}`;
  }

  clearCache(jobId?: string): void {
    if (jobId) {
      const prefix = `${jobId}:`;
      const keysToDelete = Array.from(this.cache.keys()).filter(k => k.startsWith(prefix));
      keysToDelete.forEach(k => this.cache.delete(k));
      this.jobMetadataCache.delete(jobId); // Clear metadata cache too
      this.jobFileMapping.delete(jobId); // Clear legacy mapping if used
      console.log(`[SharedDataStore] Cleared cache entries for job ${jobId}`);
    } else {
      this.cache.clear();
      this.jobMetadataCache.clear();
      this.jobFileMapping.clear();
      console.log('[SharedDataStore] Cleared all cache entries');
    }
  }

  getCacheStats(): { size: number, keys: string[] } {
    return {
      size: this.cache.size,
      keys: Array.from(this.cache.keys())
    };
  }
}

export function useSharedData(): SharedDataStore {
  return SharedDataStore.getInstance();
}

// --- REMOVED Duplicate Window interface definition --- 