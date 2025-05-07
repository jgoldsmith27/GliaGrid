/**
 * SharedDataStore.ts
 * 
 * A centralized data store for spatial data with:
 * - Caching
 * - Priority queue for loading
 * - Progressive resolution control
 * - Support for direct file access in Electron
 */
import { useState, useEffect } from 'react';
// import PriorityQueue from 'ts-priority-queue'; // Assume installed or handle later
// Use a basic array as fallback if library is missing or causes issues
type PriorityQueue<T> = T[]; 

// Options for data requests
export interface DataRequestOptions {
  filter?: { column: string; values: any[] };
  columns?: string[];
  limit?: number;
  offset?: number;
  sampleRate?: number;
  layer?: string; // For layer-specific filtering in handleRequest
  ligand?: string; // For interaction points
  receptor?: string; // For interaction points
  polygon?: [number, number][]; // ADDED: Optional polygon for spatial filtering
  signal?: AbortSignal; // Add signal to options
}

// Data request type for the priority queue
interface DataLoadRequest {
  jobId: string;
  fileType: string;
  options: DataRequestOptions; // Options now include the signal
  cacheKey: string;
  resolve?: (value: any) => void;
  reject?: (reason: any) => void;
  // signal?: AbortSignal; // Signal is now within options
}

/**
 * Shared Data Store - Singleton that manages data loading and caching
 */
export class SharedDataStore {
  private static instance: SharedDataStore;
  
  private cache: Map<string, any> = new Map();
  private loadingPromises: Map<string, Promise<any>> = new Map();
  private jobMetadataCache: Map<string, any> = new Map();
  private jobFileMapping: Map<string, Record<string, string>> = new Map(); // Still needed for getFileIdFromJobId

  // --- State for Ongoing Spatial Data Stream - REMOVED --- 

  private constructor() {
    console.log('[SharedDataStore] Created new instance');
  }

  public static getInstance(): SharedDataStore {
    if (!SharedDataStore.instance) {
      SharedDataStore.instance = new SharedDataStore();
    }
    return SharedDataStore.instance;
  }

  // Public method to request data - Simplified
  public requestData(jobId: string, fileType: string, options: DataRequestOptions = {}): Promise<any> {
    const { signal } = options; // Extract signal from options
    // Create a cache key *without* the signal, as the signal is transient
    const cacheKeyOptions = { ...options }; 
    delete cacheKeyOptions.signal; 
    const cacheKey = `${jobId}_${fileType}_${JSON.stringify(cacheKeyOptions)}`;

    // 1. Check Cache
    if (this.cache.has(cacheKey)) {
      console.log(`[SharedDataStore] Cache hit for: ${cacheKey}`);
      return Promise.resolve(this.cache.get(cacheKey));
    }
    
    // 2. Check Loading Promises (prevent duplicate concurrent requests)
    const loadingPromise = this.loadingPromises.get(cacheKey);
    if (loadingPromise) { 
      console.log(`[SharedDataStore] Request already loading: ${cacheKey}`);
      return loadingPromise;
    }

    console.log(`[SharedDataStore] Requesting data (Cache miss): ${cacheKey}`);
    // 3. Create and store the promise, then immediately start handling
    const promise = new Promise<any>(async (resolve, reject) => {
       // Check if aborted *before* even creating the request object
       if (signal?.aborted) {
           console.log(`[SharedDataStore] Request aborted before initiation: ${cacheKey}`);
           // Need to reject the promise so the caller knows it was cancelled.
           // Use a specific error type or message if possible.
           reject(new DOMException('Request aborted', 'AbortError')); 
           // Don't store this promise in loadingPromises
           return; 
       }

       const request: DataLoadRequest = {
         jobId,
         fileType,
         options, // Pass the full options including the signal
         cacheKey,
         resolve, // Pass resolve/reject to handleRequest
         reject,
       };

       // Add signal listener to reject the promise if aborted *while waiting* 
       // (though direct call makes this less likely, good practice)
       const abortHandler = () => {
         console.log(`[SharedDataStore] Abort signal received for loading request: ${cacheKey}`);
         // Check if the promise is still pending (i.e., not already resolved/rejected)
         if (this.loadingPromises.has(cacheKey)) {
             reject(new DOMException('Request aborted', 'AbortError'));
             this.loadingPromises.delete(cacheKey); // Clean up on abort
         }
       };
       signal?.addEventListener('abort', abortHandler);

       try {
           // Directly call handleRequest now, no queueing
           await this.handleRequest(request); // Pass the full request object
       } catch (error) {
           console.error(`[SharedDataStore] Error directly handling request ${cacheKey}:`, error);
           // Ensure rejection propagates, unless it's an AbortError we already handled
           if (!(error instanceof DOMException && error.name === 'AbortError')) {
               if (reject) reject(error); 
           }
           this.loadingPromises.delete(cacheKey); // Clean up promise if handling fails here
       } finally {
           // Remove the listener once done
           signal?.removeEventListener('abort', abortHandler);
       }
    });

    this.loadingPromises.set(cacheKey, promise);
    return promise;
  }
  
  public async getJobMetadata(jobId: string): Promise<any> {
    // Check cache first
    if (this.jobMetadataCache.has(jobId)) {
      return this.jobMetadataCache.get(jobId);
    }
    
    // Check loading promises
    const loadingPromise = this.loadingPromises.get(jobId + '_metadata');
    if (loadingPromise) {
        return loadingPromise;
    }
    
    console.log(`[SharedDataStore] Fetching metadata for job: ${jobId}`);
    
    const promise = new Promise(async (resolve, reject) => {
        try {
            if (!window.electronAPI?.getJobMetadata) { // Check optional chaining
                throw new Error('Electron API getJobMetadata is not available.');
            }
            const result = await window.electronAPI.getJobMetadata(jobId);
            if (result.success) {
                console.log(`[SharedDataStore] Fetched metadata successfully for ${jobId}:`, result.metadata);
                this.jobMetadataCache.set(jobId, result.metadata); 
                // Extract and store file mappings if needed
                if (result.metadata?.inputs?.files) {
                    this.jobFileMapping.set(jobId, result.metadata.inputs.files);
                }
                resolve(result.metadata);
            } else {
                 console.error(`[SharedDataStore] Failed to fetch metadata for ${jobId}:`, result.error);
                reject(new Error(result.error || 'Failed to fetch job metadata'));
            }
        } catch (error) {
            console.error(`[SharedDataStore] Error fetching metadata for ${jobId}:`, error);
            reject(error);
        } finally {
            this.loadingPromises.delete(jobId + '_metadata');
        }
    });
    
    this.loadingPromises.set(jobId + '_metadata', promise);
    return promise;
  }

  // Helper to get a specific file ID from cached metadata
  public async getFileIdFromJobId(jobId: string, fileTypeKey: string): Promise<string | null> {
    // Ensure metadata is loaded first
    const metadata = await this.getJobMetadata(jobId);
    
    if (metadata?.inputs?.files && metadata.inputs.files[fileTypeKey]) {
        return metadata.inputs.files[fileTypeKey];
    }
    
    console.warn(`[SharedDataStore] Could not find file ID for type '${fileTypeKey}' in metadata for job ${jobId}`);
    // Attempt fallback to direct mapping cache if metadata structure was unexpected
    if (this.jobFileMapping.has(jobId)) {
         const mapping = this.jobFileMapping.get(jobId);
         if (mapping && mapping[fileTypeKey]) {
             console.warn(`[SharedDataStore] Using fallback file mapping for ${fileTypeKey} on job ${jobId}`);
             return mapping[fileTypeKey];
         }
    }
    
    return null;
  }
  
  public async getFileHeaders(fileId: string): Promise<string[]> {
    const cacheKey = `headers_${fileId}`;
    if (this.cache.has(cacheKey)) {
      return this.cache.get(cacheKey);
    }
    const loadingPromise = this.loadingPromises.get(cacheKey);
    if (loadingPromise) {
      return loadingPromise;
    }

    console.log(`[SharedDataStore] Fetching headers for file ID: ${fileId}`);
    const promise = new Promise<string[]>(async (resolve, reject) => {
      try {
        if (!window.electronAPI?.readCsvHeaders) { // Use optional chaining
          throw new Error('Electron API readCsvHeaders is not available.');
        }
        const result = await window.electronAPI.readCsvHeaders(fileId);
        if (result.success && result.headers) {
          this.cache.set(cacheKey, result.headers);
          resolve(result.headers);
        } else {
          reject(new Error(result.error || 'Failed to fetch headers or headers missing'));
        }
      } catch (error) {
        reject(error);
      } finally {
        this.loadingPromises.delete(cacheKey);
      }
    });

    this.loadingPromises.set(cacheKey, promise);
    return promise;
  }

  public async readCsvChunk(fileId: string, options: any): Promise<any> {
    // No caching for chunks by default, always fetch
    const cacheKey = `chunk_${fileId}_${JSON.stringify(options)}`; // Key includes options
    
    // Optional: Check loading promises to prevent duplicate *concurrent* requests for the exact same chunk
    const loadingPromise = this.loadingPromises.get(cacheKey);
    if (loadingPromise) {
       console.log(`[SharedDataStore] Duplicate request detected for chunk: ${cacheKey}. Returning existing promise.`);
       return loadingPromise;
    }
    
    console.log(`[SharedDataStore] Reading CSV chunk for file ID: ${fileId} with options:`, options);
    
    const promise = new Promise<any>(async (resolve, reject) => {
        try {
          if (!window.electronAPI?.readCsvChunk) { // Use optional chaining
            throw new Error('Electron API readCsvChunk is not available.');
          }
          const result = await window.electronAPI.readCsvChunk(fileId, options);
          if (result.success) {
            resolve(result); // Resolve with the full result object (includes headers, rows, stats)
          } else {
            reject(new Error(result.error || 'Failed to read CSV chunk'));
          }
        } catch (error) {
          reject(error);
        } finally {
            this.loadingPromises.delete(cacheKey); // Remove promise once completed (success or fail)
        }
    });
    
    this.loadingPromises.set(cacheKey, promise);
    return promise;
  }

  // Handles the actual data fetching for INTERACTION points (keep this logic)
  private async handleRequest(request: DataLoadRequest): Promise<void> {
    const { jobId, fileType, options, cacheKey, resolve, reject } = request;
    const { signal } = options; // Extract signal from options

    if (!resolve || !reject) {
      console.error('[SharedDataStore] handleRequest called without resolve/reject functions.');
      this.loadingPromises.delete(cacheKey);
      return; 
    }
    
    // Check if aborted before starting the main logic
    if (signal?.aborted) {
        console.log(`[SharedDataStore] handleRequest aborted at start for: ${cacheKey}`);
        reject(new DOMException('Request aborted', 'AbortError'));
        this.loadingPromises.delete(cacheKey);
        return;
    }
    
    try {
        // 1. Get Job Metadata (which includes file IDs and mappings)
        const metadata = await this.getJobMetadata(jobId); 
        if (signal?.aborted) throw new DOMException('Request aborted', 'AbortError'); 
        
        if (!metadata?.inputs?.files || !metadata?.inputs?.mappings) {
            throw new Error(`Incomplete metadata for job ${jobId}. Cannot determine files or mappings.`);
        }

        // 2. Determine File ID, Mapping, and Genes based on fileType
        let spatialFileId: string | null = null;
        let mapping: any = null; // Type this better if possible
        let ligandName: string | null = null;
        let receptorName: string | null = null;
        let receptorComponentNames: string[] = []; // NEW: Store component names
        let isComplexReceptor = false; // NEW: Flag for complex receptor

        if (fileType === 'interactionPoints' && options.ligand && options.receptor) {
            spatialFileId = metadata.inputs.files.spatialFileId;
            mapping = metadata.inputs.mappings.spatialMapping;
            ligandName = options.ligand;
            receptorName = options.receptor; // Keep the original complex name if provided

            // NEW: Check for complex receptor
            if (receptorName.includes('_')) {
                isComplexReceptor = true;
                receptorComponentNames = receptorName.split('_');
                console.log(`[SharedDataStore] Detected complex receptor: ${receptorName} -> Components: ${receptorComponentNames.join(', ')}`);
            } else {
                isComplexReceptor = false;
                receptorComponentNames = [receptorName]; // Treat single receptor as a component list of one
            }

            if (!spatialFileId || !mapping || !mapping.geneCol || !mapping.xCol || !mapping.yCol) {
                throw new Error('Missing required spatial file ID or mapping info (gene, x, y) for interaction points.');
            }
        } else {
            throw new Error(`Unsupported fileType requested: ${fileType}. Only 'interactionPoints' is handled.`);
        }

        // Extract genes and columns needed
        // MODIFIED: Use ligandName and receptorComponentNames
        const genesToFetch = [ligandName, ...receptorComponentNames].filter(Boolean) as string[];
        const columnsToFetch = [mapping.geneCol, mapping.xCol, mapping.yCol].filter(Boolean);
        const layerFilter = options.layer; // Check if a specific layer is requested
        if (layerFilter && mapping.layerCol) {
            columnsToFetch.push(mapping.layerCol);
        } else if (layerFilter && !mapping.layerCol) {
             console.warn(`[SharedDataStore] Layer filtering requested for '${options.layer}', but no layerCol defined in spatialMapping.`);
        }
        
        console.log(`[SharedDataStore] Fetching spatial data for interaction points: File ID: ${spatialFileId}, Genes: ${genesToFetch.join(',')}, Columns: ${columnsToFetch.join(',')}, Layer Filter: ${layerFilter || 'None'}`);

        // 3. Fetch Filtered SPATIAL Data using readBackendFile
        if (!window.electronAPI?.readBackendFile) { 
            throw new Error('Electron API function readBackendFile is not available.');
        }
        
        // *** Pass the signal and polygon to the backend API ***
        const backendOptions = {
           filter: {
               column: mapping.geneCol,
               values: genesToFetch 
           },
           columns: columnsToFetch,
           polygon: options.polygon, // ADDED: Pass polygon
           xCol: mapping.xCol, // ADDED: Pass X column name
           yCol: mapping.yCol  // ADDED: Pass Y column name
        };
        console.log("[SharedDataStore] Calling readBackendFile with options:", backendOptions);

        const spatialDataResult = await window.electronAPI.readBackendFile(
            spatialFileId, 
            backendOptions // Pass constructed options
        );
        
        // Check if aborted immediately after the potentially long operation
        if (signal?.aborted) throw new DOMException('Request aborted', 'AbortError');

        console.log(`[SharedDataStore] Received ${spatialDataResult?.returnedRows || 0} points from readBackendFile for L/R pair ${ligandName}-${receptorName}`);
        
        // 4. Process fetched spatial data into ligand/receptor points
        const ligandPoints: { x: number; y: number }[] = [];
        // MODIFIED: Collect receptor component points with their gene names
        const receptorComponentPoints: { x: number; y: number; gene: string }[] = [];

        if (spatialDataResult?.data) {
            for (const point of spatialDataResult.data) {
                // Apply layer filtering post-fetch only if a specific layer (not 'whole_tissue') is requested
                const isSpecificLayerFilter = layerFilter && layerFilter !== 'whole_tissue';
                if (isSpecificLayerFilter && mapping.layerCol && point[mapping.layerCol] !== layerFilter) {
                    continue; // Skip points not matching the specific layer request
                }
                
                const gene = point[mapping.geneCol];
                // MODIFIED: Check ligand/component names for safety
                if (!ligandName || !receptorComponentNames || receptorComponentNames.length === 0) { console.warn("Ligand/Receptor component names missing in comparison"); }
                // Log only the first few points to avoid flooding console
                const currentPointIndex = spatialDataResult.data.indexOf(point);
                if (currentPointIndex < 5) { // Log first 5 points
                    console.log(`[SharedDataStore DEBUG] Comparing point gene: '${gene}' (type: ${typeof gene}) with Ligand: '${ligandName}', Receptor: '${receptorName}'`);
                }

                const x = parseFloat(point[mapping.xCol]);
                const y = parseFloat(point[mapping.yCol]);

                if (isNaN(x) || isNaN(y)) {
                     console.warn(`[SharedDataStore] Skipping point with invalid coordinates:`, point);
                     continue;
                }

                if (gene === ligandName) {
                    // Ligand points only need x, y
                    ligandPoints.push({ x, y });
                }
                // MODIFIED: Check if gene matches any of the receptor components
                if (receptorComponentNames.includes(gene)) {
                    // Receptor points need x, y, and gene name for complex processing
                    receptorComponentPoints.push({ x, y, gene });
                }
            }
        }
        
        // MODIFIED: Create the new data structure to resolve with
        const visualizationData = {
            ligand: ligandPoints,
            receptor: receptorComponentPoints, // Use the array containing component points
            isComplex: isComplexReceptor,     // Include the flag
            receptorName: receptorName,       // Include the original (potentially complex) name
            warnings: spatialDataResult?.warnings || []
        };

        console.log(`[SharedDataStore] Processed points for L/R pair ${ligandName}-${receptorName} (Layer: ${layerFilter || 'All'}) - Ligands: ${ligandPoints.length}, Receptor Components: ${receptorComponentPoints.length}`);

        // 5. Cache and Resolve with the processed visualization data
        this.cache.set(cacheKey, visualizationData);
        resolve(visualizationData);
        this.loadingPromises.delete(cacheKey); // Remove promise on success

      } catch (error) {
          // Handle AbortError specifically
          if (error instanceof DOMException && error.name === 'AbortError') {
              console.log(`[SharedDataStore] Caught AbortError in handleRequest for ${cacheKey}`);
              // Ensure reject was called by the listener or throw
              if (!signal?.aborted) { 
                  reject(error); // Should have been rejected by listener, but belt-and-suspenders
              }
          } else {
              console.error(`[SharedDataStore] Error in handleRequest for ${cacheKey}:`, error);
              reject(error); // Reject with other errors
          }
          // Ensure cleanup happens regardless of error type
          this.loadingPromises.delete(cacheKey); 
      }
  }

  // --- Spatial Data Streaming Methods - REMOVED --- 
  // startSpatialStream(jobId: string, options: { resolution?: number, columns?: string[] } = {}) { ... }
  // stopSpatialStream() { ... } // Assuming a stop method might have existed
  // private _cleanupStreamListeners() { ... }
  // getSpatialStreamState(): SpatialStreamState { ... }
  // subscribeSpatialStream(callback: (state: SpatialStreamState) => void): () => void { ... }
  // private _notifySpatialStreamSubscribers() { ... }
  
}

// --- React Hook for Spatial Stream State - REMOVED --- 
// export function useSpatialStreamData(): SpatialStreamState { ... }

// Existing hook to get the singleton instance (can be kept)
export function useSharedData(): SharedDataStore {
  return SharedDataStore.getInstance();
}

// --- REMOVED Duplicate Window interface definition --- 