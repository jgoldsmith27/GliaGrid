import { useState, useEffect, useCallback, useRef } from 'react';

// Define the job status data structure
interface JobStatusData {
  status: 'pending' | 'running' | 'success' | 'failed' | 'error';
  progress?: number;
  message?: string;
  results?: any;
  job_id?: string;
}

interface UseJobStatusReturn {
  jobStatus: JobStatusData | null;
  isLoading: boolean;
  error: string | null;
  refreshStatus: () => Promise<void>;
}

/**
 * Custom hook to monitor job status using event-driven Electron IPC
 * @param jobId The ID of the job to monitor
 * @returns Object containing job status, loading state, error state, and refresh function
 */
const useJobStatus = (jobId: string | undefined | null): UseJobStatusReturn => {
  const [jobStatus, setJobStatus] = useState<JobStatusData | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  
  // Refs for cleanup functions
  const statusListenerCleanupRef = useRef<(() => void) | null>(null);
  const errorListenerCleanupRef = useRef<(() => void) | null>(null);

  // Handle job status update event
  const handleStatusUpdate = useCallback((data: JobStatusData) => {
    console.log('[useJobStatus] Received event update:', data);
    
    if (data.job_id !== jobId) {
      console.log(`[useJobStatus] Ignoring update for different job: ${data.job_id}`);
      return; // Ignore updates for other jobs
    }
    
    setJobStatus(data);
    
    // Determine loading state based on status
    if (data.status === 'success' || data.status === 'failed' || data.status === 'error') {
      setIsLoading(false);
    } else {
      setIsLoading(true);
    }
    
    // Handle error state
    if (data.status === 'error' || data.status === 'failed') {
      setError(data.message || 'An error occurred during analysis');
    } else {
      setError(null);
    }
  }, [jobId]);

  // Handle error event
  const handleError = useCallback((data: { jobId: string, error: string }) => {
    console.error('[useJobStatus] Error event:', data);
    
    if (data.jobId !== jobId) {
      console.log(`[useJobStatus] Ignoring error for different job: ${data.jobId}`);
      return; // Ignore updates for other jobs
    }
    
    setError(data.error);
    setIsLoading(false);
  }, [jobId]);

  // Manual refresh function
  const refreshStatus = useCallback(async () => {
    if (!jobId) {
      console.warn('[useJobStatus] Cannot refresh without jobId');
      return;
    }
    
    try {
      console.log(`[useJobStatus] Manually refreshing status for job ${jobId}`);
      const response = await window.electronAPI.checkJobStatus(jobId);
      if (response.success) {
        handleStatusUpdate(response.status);
      } else {
        handleError({ jobId, error: response.error || 'Unknown error checking job status' });
      }
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Unknown error';
      handleError({ jobId, error: message });
    }
  }, [jobId, handleStatusUpdate, handleError]);

  useEffect(() => {
    // Clean up any existing listeners
    if (statusListenerCleanupRef.current) {
      statusListenerCleanupRef.current();
      statusListenerCleanupRef.current = null;
    }
    
    if (errorListenerCleanupRef.current) {
      errorListenerCleanupRef.current();
      errorListenerCleanupRef.current = null;
    }
    
    if (!jobId) {
      setIsLoading(false);
      setJobStatus(null);
      setError('No Job ID provided');
      return;
    }
    
    // Set initial loading state
    setIsLoading(true);
    setError(null);
    
    console.log(`[useJobStatus] Setting up listeners for job ${jobId}`);
    
    // Set up event listeners using the correct function from preload.js
    // Note: There's only one exposed listener 'onJobUpdate'. Error statuses are handled within it.
    statusListenerCleanupRef.current = window.electronAPI.onJobUpdate(handleStatusUpdate);
    
    // Subscribe to updates for this job
    window.electronAPI.subscribeJobStatus(jobId)
      .then(response => {
        if (!response.success) {
          console.error(`[useJobStatus] Subscription error: ${response.error}`);
          setError(`Failed to subscribe to job status updates: ${response.error}`);
          setIsLoading(false);
        } else {
          console.log(`[useJobStatus] Successfully subscribed to job ${jobId}`);
        }
      })
      .catch(err => {
        console.error('[useJobStatus] Failed to subscribe:', err);
        setError('Failed to subscribe to job status updates');
        setIsLoading(false);
      });
    
    // Cleanup: unsubscribe and remove listeners
    return () => {
      console.log(`[useJobStatus] Cleaning up listeners for job ${jobId}`);
      
      // Unsubscribe from job status events
      window.electronAPI.unsubscribeJobStatus(jobId)
        .catch(err => console.error(`[useJobStatus] Error unsubscribing from job ${jobId}:`, err));
      
      // Clean up event listeners
      if (statusListenerCleanupRef.current) {
        statusListenerCleanupRef.current();
        statusListenerCleanupRef.current = null;
      }
    };
  }, [jobId, handleStatusUpdate, handleError]);
  
  return { jobStatus, isLoading, error, refreshStatus };
};

export default useJobStatus; 