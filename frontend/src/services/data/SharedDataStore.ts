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
}

// Data request type for the priority queue
interface DataLoadRequest {
  jobId: string;
  fileType: string;
  options: DataRequestOptions;
  cacheKey: string;
  resolve?: (value: any) => void;
  reject?: (reason: any) => void;
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
    const cacheKey = `${jobId}_${fileType}_${JSON.stringify(options)}`;

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
       const request: DataLoadRequest = {
         jobId,
         fileType,
         options,
         cacheKey,
         resolve, // Pass resolve/reject to handleRequest
         reject,
       };
       try {
           // Directly call handleRequest now, no queueing
           await this.handleRequest(request);
       } catch (error) {
           console.error(`[SharedDataStore] Error directly handling request ${cacheKey}:`, error);
           if (reject) reject(error); // Ensure rejection propagates
           this.loadingPromises.delete(cacheKey); // Clean up promise if handling fails here
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

    if (!resolve || !reject) {
      console.error('[SharedDataStore] handleRequest called without resolve/reject functions.');
      this.loadingPromises.delete(cacheKey);
      return; 
    }

    try {
        // 1. Get Job Metadata (which includes file IDs and mappings)
        const metadata = await this.getJobMetadata(jobId);
        if (!metadata?.inputs?.files || !metadata?.inputs?.mappings) {
            throw new Error(`Incomplete metadata for job ${jobId}. Cannot determine files or mappings.`);
        }

        // 2. Determine File ID and Mapping based on fileType (expecting 'interactionPoints' for this path)
        let spatialFileId: string | null = null;
        let mapping: any = null; // Type this better if possible
        let ligandName: string | null = null;
        let receptorName: string | null = null;

        if (fileType === 'interactionPoints' && options.ligand && options.receptor) {
            spatialFileId = metadata.inputs.files.spatialFileId;
            mapping = metadata.inputs.mappings.spatialMapping;
            ligandName = options.ligand;
            receptorName = options.receptor;

            if (!spatialFileId || !mapping || !mapping.geneCol || !mapping.xCol || !mapping.yCol) {
                throw new Error('Missing required spatial file ID or mapping info (gene, x, y) for interaction points.');
            }
        } else {
            throw new Error(`Unsupported fileType requested: ${fileType}. Only 'interactionPoints' is handled.`);
        }

        // Extract genes and columns needed
        const genesToFetch = [ligandName, receptorName].filter(Boolean) as string[];
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
        
        const spatialDataResult = await window.electronAPI.readBackendFile(spatialFileId, {
           filter: {
               column: mapping.geneCol,
               values: genesToFetch 
           },
           columns: columnsToFetch 
        });
        
        console.log(`[SharedDataStore] Received ${spatialDataResult?.returnedRows || 0} points from readBackendFile for L/R pair ${ligandName}-${receptorName}`);
        
        // 4. Process fetched spatial data into ligand/receptor points
        const ligandPoints: { x: number; y: number }[] = [];
        const receptorPoints: { x: number; y: number }[] = [];
        
        if (spatialDataResult?.data) {
            for (const point of spatialDataResult.data) {
                // Apply layer filtering post-fetch only if a specific layer (not 'whole_tissue') is requested
                const isSpecificLayerFilter = layerFilter && layerFilter !== 'whole_tissue';
                if (isSpecificLayerFilter && mapping.layerCol && point[mapping.layerCol] !== layerFilter) {
                    continue; // Skip points not matching the specific layer request
                }
                
                const gene = point[mapping.geneCol];
                if (!ligandName || !receptorName) { console.warn("Ligand/Receptor name missing in comparison"); }
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
            warnings: spatialDataResult?.warnings || [] 
        };
        
        console.log(`[SharedDataStore] Processed points for L/R pair ${ligandName}-${receptorName} (Layer: ${layerFilter || 'All'}) - Ligands: ${ligandPoints.length}, Receptors: ${receptorPoints.length}`);

        // 5. Cache and Resolve with the processed visualization data
        this.cache.set(cacheKey, visualizationData);
        resolve(visualizationData);
        this.loadingPromises.delete(cacheKey); // Remove promise on success

      } catch (error) {
          console.error(`[SharedDataStore] Error in handleRequest for ${cacheKey}:`, error);
          reject(error);
          this.loadingPromises.delete(cacheKey); // Also remove promise on handleRequest error
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