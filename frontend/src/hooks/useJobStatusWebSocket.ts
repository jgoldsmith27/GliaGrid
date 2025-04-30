import { useState, useEffect, useRef } from 'react';

// Define the expected structure of the job status/result data
// (This should match what the WebSocket sends)
interface JobStatusData {
    status: 'processing' | 'success' | 'failed' | 'error';
    progress?: number;
    message?: string;
    results?: any; // Define more specifically if possible
    job_id?: string; // Can be useful for verification
}

interface UseJobStatusWebSocketReturn {
    jobStatus: JobStatusData | null;
    isConnected: boolean;
    isLoading: boolean;
    error: string | null;
}

/**
 * Custom hook to manage a WebSocket connection for monitoring analysis job status.
 * 
 * @param jobId The ID of the job to monitor. If null or undefined, no connection is made.
 * @returns An object containing the latest job status, connection state, loading state, and error messages.
 */
const useJobStatusWebSocket = (jobId: string | undefined | null): UseJobStatusWebSocketReturn => {
    const [jobStatus, setJobStatus] = useState<JobStatusData | null>(null);
    const [isConnected, setIsConnected] = useState<boolean>(false);
    const [isLoading, setIsLoading] = useState<boolean>(true); // Start loading until connected/final status
    const [error, setError] = useState<string | null>(null);
    const wsRef = useRef<WebSocket | null>(null);

    useEffect(() => {
        // Don't connect if jobId is invalid
        if (!jobId) {
            setIsConnected(false);
            setIsLoading(false);
            setJobStatus(null);
            setError('No Job ID provided for WebSocket.');
            return;
        }

        // Close any existing connection before opening a new one
        if (wsRef.current && wsRef.current.readyState !== WebSocket.CLOSED) {
            console.log('[useJobStatusWebSocket] Closing previous WebSocket connection.');
            wsRef.current.close(1000, "New connection requested");
        }

        const wsUrl = `ws://localhost:8000/ws/analysis/status/${jobId}`;
        console.log(`[useJobStatusWebSocket] Attempting to connect to ${wsUrl}`);
        const ws = new WebSocket(wsUrl);
        wsRef.current = ws; // Store the current WebSocket instance

        // Reset state on new connection attempt
        setIsConnected(false);
        setIsLoading(true); 
        setError(null);
        // Keep previous jobStatus until a new one arrives?
        // setJobStatus(null); // Optional: Decide if status should clear immediately

        ws.onopen = () => {
            console.log(`[useJobStatusWebSocket] Connected to ${wsUrl}`);
            setIsConnected(true);
            setIsLoading(true); // Remain loading until first message or final status
            setError(null);
        };

        ws.onmessage = (event) => {
            try {
                const data: JobStatusData = JSON.parse(event.data);
                console.log('[useJobStatusWebSocket] Message received:', data);
                setJobStatus(data); // Update with the latest status

                // Determine loading state based on message
                if (data.status === 'success' || data.status === 'failed' || data.status === 'error') {
                    setIsLoading(false);
                } else {
                    setIsLoading(true); // Still processing
                }

                // Handle specific error status from backend
                if (data.status === 'error') {
                    setError(data.message || 'An error occurred during analysis (via WebSocket).');
                }

            } catch (parseError) {
                console.error('[useJobStatusWebSocket] Error parsing message:', parseError, 'Data:', event.data);
                setError('Received invalid data from server.');
                setIsLoading(false);
            }
        };

        ws.onerror = (event) => {
            console.error('[useJobStatusWebSocket] WebSocket Error:', event);
            setError('WebSocket connection error.');
            setIsConnected(false);
            setIsLoading(false);
            // Optionally clear status on error?
            // setJobStatus(null);
        };

        ws.onclose = (event) => {
            console.log(`[useJobStatusWebSocket] Disconnected from ${wsUrl}. Code: ${event.code}, Reason: ${event.reason}`);
            setIsConnected(false);
            // Only set loading false if we didn't reach a final state
            // Use functional state update to avoid stale closure issues with jobStatus
            setJobStatus(currentStatus => {
                 if (currentStatus?.status !== 'success' && currentStatus?.status !== 'failed' && currentStatus?.status !== 'error') {
                    setIsLoading(false);
                     // Report disconnection error only if it wasn't a normal closure or explicitly requested
                     if (event.code !== 1000 && event.reason !== "New connection requested" && event.reason !== "Component unmounting") {
                        setError(prevError => prevError || 'WebSocket disconnected unexpectedly.');
                     }
                 }
                 return currentStatus;
            });
        };

        // Cleanup function: close WebSocket when jobId changes or component unmounts
        return () => {
            if (wsRef.current && wsRef.current.readyState !== WebSocket.CLOSED) {
                console.log(`[useJobStatusWebSocket] Closing WebSocket connection to ${wsUrl} on cleanup.`);
                wsRef.current.close(1000, "Component unmounting");
                wsRef.current = null;
            }
        };
    }, [jobId]); // Effect runs only when jobId changes

    return { jobStatus, isConnected, isLoading, error };
};

export default useJobStatusWebSocket; 