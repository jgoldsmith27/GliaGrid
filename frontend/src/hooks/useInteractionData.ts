import { useState, useEffect, useRef, useCallback } from 'react';
import { SharedDataStore, useSharedData, DataRequestOptions } from '../services/data/SharedDataStore';
import { isEqual } from 'lodash'; // Add lodash for deep equality check

// Type for the data returned by the /api/visualization endpoint
export interface InteractionVisualizationData {
  ligand: { x: number; y: number; layer: string }[];
  receptor: { x: number; y: number; layer: string }[];
  warnings?: string[];
}

interface UseInteractionDataReturn {
    interactionVizData: InteractionVisualizationData | null;
    isLoading: boolean;
    error: string | null;
    warnings: string[];
    fetchInteractionData: () => void; // Function to manually trigger fetch if needed (or could be internal)
    cancelFetch: () => void; // Function to cancel ongoing fetch
}

// Static map to track active visualization fetches across hook instances
// Key is a hash of job ID + scope to prevent conflicts between different visualizations
const activeVisualizations = new Map<string, AbortController>();

/**
 * Custom hook to fetch and manage interaction visualization data.
 * 
 * @param jobId The ID of the analysis job.
 * @param selectedPair The currently selected [ligand, receptor] pair. Null if none selected.
 * @param apiScopeName The specific scope name for the API ('whole_tissue' or a layer name like 'Layer1'). 'custom' is handled separately (no fetch).
 * @param lassoCoords Lasso coordinates for polygon filtering
 * @returns State variables for data, loading, error, warnings, and control functions.
 */
const useInteractionData = (
    jobId: string | null, 
    selectedPair: [string, string] | null, 
    apiScopeName: string | null, // e.g., 'whole_tissue', 'Layer1', 'custom'
    lassoCoords: [number, number][] | null // ADDED: Lasso coordinates for polygon filtering
): UseInteractionDataReturn => {
    const [interactionVizData, setInteractionVizData] = useState<InteractionVisualizationData | null>(null);
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [warnings, setWarnings] = useState<string[]>([]);
    const abortControllerRef = useRef<AbortController | null>(null);
    
    // Add a ref to track if a fetch is currently in progress for the same parameters
    const fetchInProgressRef = useRef<boolean>(false);
    // Add a ref to track current fetch params to avoid duplicate fetches
    const currentFetchParamsRef = useRef<string | null>(null);
    // Add a unique ID for this hook instance to help with debugging
    const hookInstanceId = useRef<string>(`hook_${Math.random().toString(36).substring(2, 9)}`);
    // Add a debounce timer ref
    const debounceTimerRef = useRef<number | null>(null);
    // Keep track of mounted state
    const isMountedRef = useRef(true);
    // Keep a stable ref of selectedPair
    const selectedPairRef = useRef<[string, string] | null>(null);
    // Keep a stable ref of lassoCoords
    const lassoCoordsRef = useRef<[number, number][] | null>(null);
    
    // Get the shared data store instance
    const dataStore = useSharedData();

    // Update stable refs with current values
    useEffect(() => {
        selectedPairRef.current = selectedPair;
    }, [selectedPair]);

    useEffect(() => {
        lassoCoordsRef.current = lassoCoords;
    }, [lassoCoords]);

    // Ensure cleanup on unmount
    useEffect(() => {
        return () => {
            isMountedRef.current = false;
            // Clear any pending timers
            if (debounceTimerRef.current !== null) {
                window.clearTimeout(debounceTimerRef.current);
                debounceTimerRef.current = null;
            }
            // Abort any in-progress fetches
            if (abortControllerRef.current) {
                abortControllerRef.current.abort();
                abortControllerRef.current = null;
            }
            
            // Remove this instance from the active visualizations map
            if (jobId && selectedPair) {
                const visKey = `${jobId}_${selectedPair.join('_')}_${apiScopeName}`;
                if (activeVisualizations.get(visKey)?.signal === abortControllerRef.current?.signal) {
                    activeVisualizations.delete(visKey);
                }
            }
        };
    }, [jobId, selectedPair, apiScopeName]);

    const cancelFetch = useCallback(() => {
        if (abortControllerRef.current) {
            console.log(`[${hookInstanceId.current}] Cancelling fetch.`);
            abortControllerRef.current.abort();
            abortControllerRef.current = null;
            fetchInProgressRef.current = false;
            if (isMountedRef.current) {
                setIsLoading(false); // Ensure loading stops on manual cancel
            }
        }
        
        // Also clear any pending debounced fetches
        if (debounceTimerRef.current !== null) {
            window.clearTimeout(debounceTimerRef.current);
            debounceTimerRef.current = null;
        }
    }, []);

    // Create a hash of the current fetch params to detect duplicate fetches
    const getFetchParamsHash = useCallback((
        jobId: string | null, 
        selectedPair: [string, string] | null, 
        apiScopeName: string | null,
        lassoCoords: [number, number][] | null
    ): string => {
        const polygonKey = lassoCoords ? 'withPolygon' : 'noPolygon';
        return `${jobId}_${selectedPair?.[0]}_${selectedPair?.[1]}_${apiScopeName}_${polygonKey}`;
    }, []);

    const fetchInteractionData = useCallback(async () => {
        // Prevent executing if component is unmounted
        if (!isMountedRef.current) return;
        
        // Get latest values from stable refs
        const currentSelectedPair = selectedPairRef.current;
        const currentLassoCoords = lassoCoordsRef.current;
        
        // Guard conditions: Don't fetch if no job or no pair selected.
        if (!jobId || !currentSelectedPair) {
            if (isMountedRef.current) {
                setInteractionVizData(null);
                setError(null);
                setWarnings([]);
            }
            // If a fetch was in progress for a previous valid state, cancel it
            if (abortControllerRef.current) {
                 cancelFetch();
            }
            return; 
        }

        // Create a hash of the current fetch params
        const paramsHash = getFetchParamsHash(jobId, currentSelectedPair, apiScopeName, currentLassoCoords);
        
        // Don't refetch if we're already fetching for these exact parameters
        if (fetchInProgressRef.current && currentFetchParamsRef.current === paramsHash) {
            console.log(`[${hookInstanceId.current}] Skipping duplicate fetch for ${currentSelectedPair.join('-')}, scope: ${apiScopeName}`);
            return;
        }

        // Check if another visualization is already fetching this exact data
        const visKey = `${jobId}_${currentSelectedPair.join('_')}_${apiScopeName}`;
        const existingController = activeVisualizations.get(visKey);
        if (existingController && !existingController.signal.aborted) {
            console.log(`[${hookInstanceId.current}] Another visualization is already fetching this data, waiting for it to complete`);
            // Don't start a new fetch, let the existing one finish
            return;
        }

        // Abort previous fetch if any
        if (abortControllerRef.current) {
            abortControllerRef.current.abort();
            abortControllerRef.current = null;
        }
        
        const controller = new AbortController();
        abortControllerRef.current = controller;
        
        // Register this fetch in the active visualizations map
        activeVisualizations.set(visKey, controller);
        
        fetchInProgressRef.current = true;
        currentFetchParamsRef.current = paramsHash;

        console.log(`[${hookInstanceId.current}] Fetching for ${currentSelectedPair.join('-')}, scope: ${apiScopeName}`);
        
        if (isMountedRef.current) {
            setIsLoading(true);
            setError(null);
            // Don't clear existing data during fetch to avoid flashing
            // setInteractionVizData(null);
            setWarnings([]);
        }

        try {
            const [ligand, receptor] = currentSelectedPair;
            
            // Use the SharedDataStore to fetch data
            const options: DataRequestOptions = {
                ligand,
                receptor,
                layer: apiScopeName ?? undefined, // Pass undefined if apiScopeName is null
                polygon: currentLassoCoords || undefined, // MODIFIED: Pass undefined if lassoCoords is null
                signal: controller.signal, // Pass the abort signal
            };
            
            // Create a URL for logging purposes only
            // MODIFIED: Log layer/polygon info better
            const scopeInfo = apiScopeName ? `layer=${encodeURIComponent(apiScopeName)}` : 'scope=unknown';
            const polygonInfo = currentLassoCoords ? `polygon=present` : 'polygon=absent';
            const logUrl = `Using SharedDataStore for: ${jobId}?ligand=${encodeURIComponent(ligand)}&receptor=${encodeURIComponent(receptor)}&${scopeInfo}&${polygonInfo}`;
            console.log(`[${hookInstanceId.current}] Fetching:`, logUrl);

            // Track abort controller state
            const abortSignal = controller.signal;
            if (abortSignal.aborted) {
                console.log(`[${hookInstanceId.current}] Fetch aborted before starting.`);
                fetchInProgressRef.current = false;
                activeVisualizations.delete(visKey);
                return;
            }

            // Use data store to get the data
            try {
                const data = await dataStore.requestData(jobId, 'interactionPoints', options);
                
                // Check if aborted during fetch or if component unmounted
                if (abortSignal.aborted || !isMountedRef.current) {
                    console.log(`[${hookInstanceId.current}] Fetch aborted during execution.`);
                    fetchInProgressRef.current = false;
                    activeVisualizations.delete(visKey);
                    return;
                }
                
                // Check if this is still the current request
                if (abortControllerRef.current !== controller) {
                    console.log(`[${hookInstanceId.current}] Stale fetch response ignored.`);
                    activeVisualizations.delete(visKey);
                    return;
                }
                
                console.log(`[${hookInstanceId.current}] Fetched data:`, data);
                
                // Create the visualization data object from the returned data
                const vizData: InteractionVisualizationData = {
                    ligand: data.ligand || [],
                    receptor: data.receptor || [],
                    warnings: data.warnings || []
                };
                
                if (isMountedRef.current) {
                    setInteractionVizData(vizData);
                    setWarnings(vizData.warnings || []);
                }
            } catch (fetchErr) {
                // Only handle error if this is still the current request
                if (abortControllerRef.current === controller && !abortSignal.aborted && isMountedRef.current) {
                    throw fetchErr;
                }
            }
        } catch (err) {
             // Check again if the error corresponds to the *current* controller
             if (abortControllerRef.current === controller && !controller.signal.aborted && isMountedRef.current) {
                 console.error(`[${hookInstanceId.current}] Catch Error:`, err);
                 if ((err as Error).name === 'AbortError') {
                     console.log(`[${hookInstanceId.current}] Fetch aborted (caught).`);
                 } else {
                     const message = err instanceof Error ? err.message : 'An unknown error occurred';
                     setError(message);
                     // Don't clear visualization data on error to prevent flashing
                     // setInteractionVizData(null);
                     setWarnings([]);
                 }
             } else {
                  console.log(`[${hookInstanceId.current}] Error from stale fetch ignored.`, err);
             }
        } finally {
            // Only clear controller and set loading false if this is the *latest* request controller
            if (abortControllerRef.current === controller && isMountedRef.current) {
                 abortControllerRef.current = null;
                 fetchInProgressRef.current = false;
                 setIsLoading(false);
                 
                 // Remove from active visualizations map
                 const visKey = `${jobId}_${currentSelectedPair.join('_')}_${apiScopeName}`;
                 activeVisualizations.delete(visKey);
            }
        }
    }, [jobId, apiScopeName, cancelFetch, dataStore, getFetchParamsHash, hookInstanceId]);

    // Debounced fetch function to prevent rapid fetch requests
    const debouncedFetch = useCallback(() => {
        // Clear any existing timer
        if (debounceTimerRef.current !== null) {
            window.clearTimeout(debounceTimerRef.current);
            debounceTimerRef.current = null;
        }
        
        // Set a new timer to execute the fetch after 500ms (increased from 250ms)
        debounceTimerRef.current = window.setTimeout(() => {
            debounceTimerRef.current = null;
            fetchInteractionData();
        }, 500);
    }, [fetchInteractionData]);

    // Previous params for deep comparison
    const prevParamsRef = useRef<{
        jobId: string | null;
        selectedPair: [string, string] | null;
        apiScopeName: string | null;
        lassoCoords: [number, number][] | null;
    }>({
        jobId: null,
        selectedPair: null,
        apiScopeName: null,
        lassoCoords: null
    });

    // Effect to trigger fetch when dependencies change
    useEffect(() => {
        // Don't trigger fetch if the component is unmounted
        if (!isMountedRef.current) return;
        
        // Deep compare current params with previous params to prevent unnecessary fetches
        const currentParams = {
            jobId,
            selectedPair,
            apiScopeName,
            lassoCoords
        };
        
        // Skip if params haven't changed (using deep equality)
        if (isEqual(currentParams, prevParamsRef.current)) {
            return;
        }
        
        // Update previous params
        prevParamsRef.current = currentParams;
        
        // Create a hash of the current fetch params
        const paramsHash = getFetchParamsHash(jobId, selectedPair, apiScopeName, lassoCoords);
        
        // Only fetch if we have valid fetch parameters and are not already fetching these params
        if (jobId && selectedPair && (currentFetchParamsRef.current !== paramsHash || !fetchInProgressRef.current)) {
            console.log(`[${hookInstanceId.current}] Effect triggered fetch for ${selectedPair.join('-')}, scope: ${apiScopeName}, params hash: ${paramsHash}`);
            debouncedFetch();
        }

        // Cleanup function to abort request when dependencies change or component unmounts
        return () => {
             if (abortControllerRef.current) {
                  console.log(`[${hookInstanceId.current}] Cleanup: Aborting fetch.`);
                  abortControllerRef.current.abort();
                  abortControllerRef.current = null;
                  fetchInProgressRef.current = false;
             }
             
             // Also clear any pending debounced fetches
             if (debounceTimerRef.current !== null) {
                 window.clearTimeout(debounceTimerRef.current);
                 debounceTimerRef.current = null;
             }
        };
    }, [jobId, selectedPair, apiScopeName, lassoCoords, debouncedFetch, getFetchParamsHash, hookInstanceId]);

    return { interactionVizData, isLoading, error, warnings, fetchInteractionData, cancelFetch };
};

export default useInteractionData; 