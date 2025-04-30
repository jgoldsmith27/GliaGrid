import { useState, useEffect, useRef, useCallback } from 'react';
import { SharedDataStore, useSharedData, DataRequestOptions } from '../services/data/SharedDataStore';

// Type for the data returned by the /api/visualization endpoint
interface InteractionVisualizationData {
  ligand: { x: number; y: number }[];
  receptor: { x: number; y: number }[];
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

/**
 * Custom hook to fetch and manage interaction visualization data.
 * 
 * @param jobId The ID of the analysis job.
 * @param selectedPair The currently selected [ligand, receptor] pair. Null if none selected.
 * @param apiScopeName The specific scope name for the API ('whole_tissue' or a layer name like 'Layer1'). 'custom' is handled separately (no fetch).
 * @returns State variables for data, loading, error, warnings, and control functions.
 */
const useInteractionData = (
    jobId: string | null, 
    selectedPair: [string, string] | null, 
    apiScopeName: string | null // e.g., 'whole_tissue', 'Layer1', 'custom'
): UseInteractionDataReturn => {
    const [interactionVizData, setInteractionVizData] = useState<InteractionVisualizationData | null>(null);
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [warnings, setWarnings] = useState<string[]>([]);
    const abortControllerRef = useRef<AbortController | null>(null);
    
    // Get the shared data store instance
    const dataStore = useSharedData();

    const cancelFetch = useCallback(() => {
        if (abortControllerRef.current) {
            console.log('[useInteractionData] Cancelling fetch.');
            abortControllerRef.current.abort();
            abortControllerRef.current = null;
            setIsLoading(false); // Ensure loading stops on manual cancel
        }
    }, []);

    const fetchInteractionData = useCallback(async () => {
        // Guard conditions: Don't fetch if no job, no pair, or scope is 'custom'
        if (!jobId || !selectedPair || !apiScopeName || apiScopeName === 'custom') {
            setInteractionVizData(null);
            setError(null);
            setWarnings([]);
            setIsLoading(false);
             // If a fetch was in progress for a previous valid state, cancel it
            if (abortControllerRef.current) {
                 cancelFetch();
            }
            return; 
        }

        // Abort previous fetch if any
        if (abortControllerRef.current) {
            abortControllerRef.current.abort();
        }
        const controller = new AbortController();
        abortControllerRef.current = controller;

        console.log(`[useInteractionData] Fetching for ${selectedPair.join('-')}, scope: ${apiScopeName}`);
        
        setIsLoading(true);
        setError(null);
        setInteractionVizData(null);
        setWarnings([]);

        try {
            const [ligand, receptor] = selectedPair;
            
            // Use the SharedDataStore to fetch data
            const options: DataRequestOptions = {
                ligand,
                receptor,
                layer: apiScopeName,
            };
            
            // Create a URL for logging purposes only
            const logUrl = `Using SharedDataStore for: ${jobId}?ligand=${encodeURIComponent(ligand)}&receptor=${encodeURIComponent(receptor)}&layer=${encodeURIComponent(apiScopeName)}`;
            console.log("[useInteractionData] Fetching:", logUrl);

            // Track abort controller state
            const abortSignal = controller.signal;
            if (abortSignal.aborted) {
                console.log('[useInteractionData] Fetch aborted before starting.');
                return;
            }

            // Use data store to get the data
            try {
                const data = await dataStore.requestData(jobId, 'interactionPoints', options);
                
                // Check if aborted during fetch
                if (abortSignal.aborted) {
                    console.log('[useInteractionData] Fetch aborted during execution.');
                    return;
                }
                
                // Check if this is still the current request
                if (abortControllerRef.current !== controller) {
                    console.log('[useInteractionData] Stale fetch response ignored.');
                    return;
                }
                
                console.log("[useInteractionData] Fetched data:", data);
                
                // Create the visualization data object from the returned data
                const vizData: InteractionVisualizationData = {
                    ligand: data.ligand || [],
                    receptor: data.receptor || [],
                    warnings: data.warnings || []
                };
                
                setInteractionVizData(vizData);
                setWarnings(vizData.warnings || []);
            } catch (fetchErr) {
                // Only handle error if this is still the current request
                if (abortControllerRef.current === controller && !abortSignal.aborted) {
                    throw fetchErr;
                }
            }
        } catch (err) {
             // Check again if the error corresponds to the *current* controller
             if (abortControllerRef.current === controller || !controller.signal.aborted) {
                 console.error("[useInteractionData] Catch Error:", err);
                 if ((err as Error).name === 'AbortError') {
                     console.log('[useInteractionData] Fetch aborted (caught).');
                 } else {
                     const message = err instanceof Error ? err.message : 'An unknown error occurred';
                     setError(message);
                     setInteractionVizData(null); // Clear data on error
                     setWarnings([]);
                 }
             } else {
                  console.log('[useInteractionData] Error from stale fetch ignored.', err);
             }
        } finally {
            // Only clear controller and set loading false if this is the *latest* request controller
            if (abortControllerRef.current === controller) {
                 abortControllerRef.current = null;
                 setIsLoading(false);
            }
        }
    }, [jobId, selectedPair, apiScopeName, cancelFetch, dataStore]); // Added dataStore as dependency

    // Effect to trigger fetch when dependencies change
    useEffect(() => {
        fetchInteractionData();

        // Cleanup function to abort request if dependencies change or component unmounts
        return () => {
             if (abortControllerRef.current) {
                  console.log('[useInteractionData] Cleanup: Aborting fetch.');
                  abortControllerRef.current.abort();
                  abortControllerRef.current = null;
             }
        };
    }, [fetchInteractionData]); // Depend on the memoized fetch function

    return { interactionVizData, isLoading, error, warnings, fetchInteractionData, cancelFetch };
};

export default useInteractionData; 