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
  x_col?: string;
  y_col?: string;
  priority?: Priority;
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

// Simple in-memory priority queue
class PriorityQueue<T extends { priority: Priority }> {
  private items: (T & { numericPriority: number })[] = [];
  
  add(item: T): void {
    // Convert priority string to numeric value
    const numericPriority = 
      item.priority === 'high' ? 3 :
      item.priority === 'medium' ? 2 : 1;
    
    // Add item with numeric priority
    this.items.push({...item, numericPriority});
    
    // Sort by priority (higher values first)
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
  
  // Data storage
  private cache: Map<string, any> = new Map();
  private loadingPromises: Map<string, Promise<any>> = new Map();
  
  // Request management
  private priorityQueue = new PriorityQueue<DataLoadRequest>();
  private isProcessing = false;

  // Job to file mapping cache
  private jobFileMapping: Map<string, Record<string, string>> = new Map();
  
  // Singleton pattern
  private constructor() {
    console.log('[SharedDataStore] Created new instance');
  }
  
  static getInstance(): SharedDataStore {
    if (!this.instance) {
      this.instance = new SharedDataStore();
    }
    return this.instance;
  }
  
  /**
   * Get data with caching and priority queue support
   */
  async getChunk(
    jobId: string, 
    fileType: string, 
    options: DataRequestOptions = {}
  ): Promise<any> {
    // Check for Electron environment
    if (!window.electronAPI) {
      throw new Error('Direct file access is required but Electron API is not available');
    }

    // Separate priority from other options
    const { priority = 'medium', ...fetchOptions } = options;
    
    // Create a cache key from the parameters
    const cacheKey = this.createCacheKey(jobId, fileType, fetchOptions);
    
    // Check cache first
    if (this.cache.has(cacheKey)) {
      console.log(`[SharedDataStore] Cache hit for ${cacheKey}`);
      return this.cache.get(cacheKey);
    }
    
    // Check if already loading
    if (this.loadingPromises.has(cacheKey)) {
      console.log(`[SharedDataStore] Already loading ${cacheKey}, returning existing promise`);
      return this.loadingPromises.get(cacheKey);
    }
    
    // Create a new promise for this request
    console.log(`[SharedDataStore] Creating new request for ${cacheKey} with priority ${priority}`);
    
    let resolvePromise: (value: any) => void;
    let rejectPromise: (reason: any) => void;
    
    const promise = new Promise((resolve, reject) => {
      resolvePromise = resolve;
      rejectPromise = reject;
    });
    
    // Create a request object
    const request: DataLoadRequest = {
      jobId,
      fileType,
      options: fetchOptions,
      priority,
      cacheKey,
      resolve: resolvePromise!,
      reject: rejectPromise!
    };
    
    // Add to loading promises to track it
    this.loadingPromises.set(cacheKey, promise);
    
    // Add to priority queue
    this.enqueue(request);
    
    return promise;
  }
  
  /**
   * Add a request to the priority queue and start processing if needed
   */
  private enqueue(request: DataLoadRequest): void {
    this.priorityQueue.add(request);
    
    if (!this.isProcessing) {
      console.log('[SharedDataStore] Starting queue processing');
      this.processQueue();
    }
  }
  
  /**
   * Process the priority queue
   */
  private async processQueue(): Promise<void> {
    this.isProcessing = true;
    
    while (!this.priorityQueue.isEmpty()) {
      const request = this.priorityQueue.pop();
      if (!request) break;
      
      try {
        console.log(`[SharedDataStore] Processing request: ${request.cacheKey} (priority: ${request.priority})`);
        await this.handleRequest(request);
      } catch (error) {
        console.error(`[SharedDataStore] Error processing request: ${request.cacheKey}`, error);
      }
    }
    
    this.isProcessing = false;
    console.log('[SharedDataStore] Queue processing complete');
  }
  
  /**
   * Handle an individual data request
   */
  private async handleRequest(request: DataLoadRequest): Promise<void> {
    const { jobId, fileType, options, cacheKey, resolve, reject } = request;
    
    if (!resolve || !reject) {
      console.error('[SharedDataStore] Missing resolve/reject in request', request);
      return;
    }
    
    try {
      if (!window.electronAPI) {
        throw new Error('Direct file access is required but Electron API is not available');
      }

      // Special handling for visualization data
      if (fileType === 'visualization') {
        console.log(`[SharedDataStore] Fetching visualization data for job ${jobId}`, options);
        const result = await window.electronAPI.readVisualizationData(jobId, options);
        
        if (!result.success || result.error) {
          throw new Error(result.error || 'Failed to fetch visualization data');
        }
        
        // Cache the visualization data
        this.cache.set(cacheKey, result.data);
        
        // Resolve the promise with visualization data
        resolve(result.data);
        return;
      }

      // For regular file data
      // Get file ID from job context
      const fileId = await this.getFileIdFromJobId(jobId, fileType);
      
      if (!fileId) {
        throw new Error(`No file ID found for job ${jobId} and file type ${fileType}`);
      }
      
      console.log(`[SharedDataStore] Found fileId ${fileId} for job ${jobId}`);
      
      // Access file directly through Electron
      const result = await window.electronAPI.readBackendFile(fileId, options);
      console.log(`[SharedDataStore] Direct file access successful, got ${result?.data?.length || 0} points`);
      
      // Cache the result
      this.cache.set(cacheKey, result);
      
      // Resolve the promise
      resolve(result);
      
    } catch (error) {
      console.error('[SharedDataStore] Error handling request:', error);
      reject(error);
    } finally {
      this.loadingPromises.delete(cacheKey);
    }
  }
  
  /**
   * Create a cache key from request parameters
   */
  private createCacheKey(jobId: string, fileType: string, options: Omit<DataRequestOptions, 'priority'>): string {
    return `${jobId}:${fileType}:${JSON.stringify(options)}`;
  }
  
  /**
   * Get file ID from job ID and file type
   */
  private async getFileIdFromJobId(jobId: string, fileType: string): Promise<string | null> {
    // Check if mapping is already cached
    if (this.jobFileMapping.has(jobId)) {
      const mapping = this.jobFileMapping.get(jobId);
      if (mapping?.[fileType]) {
        return mapping[fileType];
      }
    }
    
    // Get job metadata through Electron
    try {
      // Get file mapping from job metadata
      const result = await window.electronAPI!.getJobMetadata(jobId);
      
      if (!result || !result.success) {
        throw new Error(`Failed to get job metadata: ${result?.error || 'Unknown error'}`);
      }
      
      const mapping: Record<string, string> = {};
      
      // Extract file IDs based on file type
      if (result.context) {
        // Map API fileType to context property
        const fileIdProperty = fileType === 'spatial' 
          ? 'spatialFileId' 
          : fileType === 'interactions' 
            ? 'interactionsFileId' 
            : fileType === 'modules' 
              ? 'modulesFileId' 
              : null;
        
        if (fileIdProperty && result.context[fileIdProperty]) {
          mapping[fileType] = result.context[fileIdProperty];
        }
      }
      
      // Cache the mapping
      this.jobFileMapping.set(jobId, mapping);
      
      // Return the file ID or null if not found
      return mapping[fileType] || null;
    } catch (error) {
      console.error(`[SharedDataStore] Error getting file ID for job ${jobId}:`, error);
      return null;
    }
  }
  
  /**
   * Clear cache entries
   */
  clearCache(jobId?: string): void {
    // If jobId provided, only clear entries for that job
    if (jobId) {
      const prefix = `${jobId}:`;
      // Filter keys that start with jobId prefix
      const keysToRemove: string[] = [];
      this.cache.forEach((_, key) => {
        if (key.startsWith(prefix)) {
          keysToRemove.push(key);
        }
      });
      
      // Remove the entries
      keysToRemove.forEach(key => {
        this.cache.delete(key);
      });
      
      // Remove job mapping
      this.jobFileMapping.delete(jobId);
      
      console.log(`[SharedDataStore] Cleared ${keysToRemove.length} cache entries for job ${jobId}`);
    } else {
      // Clear all cache and mapping
      this.cache.clear();
      this.jobFileMapping.clear();
      console.log('[SharedDataStore] Cleared entire cache');
    }
  }
  
  /**
   * Get cache statistics
   */
  getCacheStats(): { size: number, keys: string[] } {
    return {
      size: this.cache.size,
      keys: Array.from(this.cache.keys())
    };
  }
}

/**
 * Hook to use the SharedDataStore
 */
export function useSharedData(): SharedDataStore {
  return SharedDataStore.getInstance();
}

// Extend the Window interface to include our Electron API
declare global {
  interface Window {
    electronAPI?: {
      readBackendFile: (fileId: string, options: any) => Promise<any>;
      getJobMetadata: (jobId: string) => Promise<any>;
      readVisualizationData: (jobId: string, options: any) => Promise<{ success: boolean; data?: any; error?: string }>;
      saveProject: (projectName: string, projectData: any) => Promise<any>;
      listProjects: () => Promise<any>;
      loadProject: (filePath: string) => Promise<any>;
      deleteProject: (filePath: string, projectName: string) => Promise<any>;
    }
  }
} 