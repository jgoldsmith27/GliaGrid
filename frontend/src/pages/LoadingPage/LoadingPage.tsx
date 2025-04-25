import React, { useEffect, useState, useRef } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import styles from './LoadingPage.module.css';
import LoadingAnimation from '../../components/LoadingPage/LoadingPage';
import ErrorDisplay from '../../components/ErrorDisplay/ErrorDisplay';
import { getWsUrl } from '../../config';

// Define types for the enhanced progress updates
interface ProgressUpdate {
  status: 'running' | 'success' | 'failed' | 'error' | 'pending' | 'PROCESSING' | 'complete';
  progress: number;
  message: string;
  stage_id?: string;
  stage_name?: string;
  stage_type?: string;
  current_scope?: string;
  detail?: string;
  error_type?: string;
  error_detail?: string;
  job_id: string;
  results?: any;
}

// WebSocket connection singleton
declare global {
  interface Window {
    wsConnections?: Map<string, WebSocket>;
  }
}

// Use the global wsConnections if it exists, or create a new one
if (!window.wsConnections) {
  window.wsConnections = new Map<string, WebSocket>();
}
const wsConnections = window.wsConnections;

const LoadingPage: React.FC = () => {
  const { jobId } = useParams<{ jobId: string }>();
  const navigate = useNavigate();
  
  const [progressData, setProgressData] = useState<ProgressUpdate | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [navigating, setNavigating] = useState<boolean>(false);
  const [reconnectAttempts, setReconnectAttempts] = useState<number>(0);
  const maxReconnectAttempts = 5;
  const [isConnected, setIsConnected] = useState<boolean>(false);
  const pingIntervalRef = useRef<number | null>(null);
  const [showReconnect, setShowReconnect] = useState<boolean>(false);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [showBypass, setShowBypass] = useState<boolean>(false);
  
  // Add component lifecycle logging
  useEffect(() => {
    console.log('[LoadingPage] Component mounted with jobId:', jobId);
    return () => {
      console.log('[LoadingPage] Component unmounting');
    };
  }, [jobId]);
  
  // Format stage-specific messages
  const getStageMessage = (stage: string): string => {
    switch(stage) {
      case 'data_loading': return 'Loading data files...';
      case 'spatial_loading': return 'Loading spatial data...';
      case 'interaction_loading': return 'Loading interaction data...';
      case 'module_loading': return 'Loading module data...';
      case 'counts_analysis': return 'Calculating ligand-receptor counts...';
      case 'pathway_dominance': return 'Analyzing pathway dominance...';
      case 'module_context': return 'Analyzing module context...';
      case 'cleanup': return 'Cleaning up temporary files...';
      case 'complete': return 'Analysis complete!';
      default: return `Processing ${stage}...`;
    }
  };
  
  // Format a detailed message from the progress data
  const getDetailedMessage = (data: ProgressUpdate): string => {
    // If we have a specific stage, use that for a more detailed message
    if (data.stage_id) {
      let message = getStageMessage(data.stage_id);
      
      // Add current scope information (e.g., which layer is being processed)
      if (data.current_scope) {
        message = `${message} Scope: ${data.current_scope}`;
      }
      
      return message;
    }
    
    // Otherwise use the default message logic
    let message = data.message || 'Processing...';
    
    // Add stage information if available
    if (data.stage_name) {
      message = `${data.stage_name}: ${message}`;
    }
    
    // Add current scope information (e.g., which layer is being processed)
    if (data.current_scope) {
      message = `${message} - Scope: ${data.current_scope}`;
    }
    
    // Add any additional details
    if (data.detail) {
      message = `${message} - ${data.detail}`;
    }
    
    return message;
  };
  
  // Handle navigation based on progress data
  useEffect(() => {
    console.log('[LoadingPage] Navigation effect running with:', {
      hasProgressData: !!progressData,
      progressDataStatus: progressData?.status,
      progressDataProgress: progressData?.progress,
      navigating
    });
    
    if (!progressData || navigating) return;
    
    // Check for success condition with completed status
    const isComplete = 
      (progressData.status === 'success' || progressData.status === 'complete') && 
      progressData.progress >= 0.99;
    
    console.log('[LoadingPage] isComplete check result:', isComplete);
    
    if (isComplete && !navigating) {
      console.log('%c[LoadingPage] Analysis complete, preparing for navigation', 'background: purple; color: white');
      setNavigating(true);
      
      // Small delay to ensure user sees the 100% completion state
      const delay = 1500;
      console.log(`[LoadingPage] Will navigate to results in ${delay}ms`);
      
      setTimeout(() => {
        console.log('%c[LoadingPage] Navigating to results page now: /analysis/${jobId}', 'background: purple; color: white');
        navigate(`/analysis/${jobId}`, { replace: true });
      }, delay);
    }
  }, [progressData, navigating, jobId, navigate]);

  // Function to send ping message to server
  const sendPing = () => {
    if (wsConnections.has(jobId || '')) {
      const ws = wsConnections.get(jobId || '')!;
      if (ws.readyState === WebSocket.OPEN) {
        console.log('[LoadingPage] Sending ping to server');
        ws.send(JSON.stringify({ type: 'ping' }));
      }
    }
  };

  // WebSocket connection setup
  useEffect(() => {
    if (!jobId) {
      setError('No job ID provided');
      return;
    }

    // Define WebSocket URL using the helper function from config
    const wsUrl = getWsUrl(`/ws/analysis/status/${jobId}`);
    let ws: WebSocket | null = null;
    
    // Check if we already have a connection for this job
    if (wsConnections.has(jobId)) {
      ws = wsConnections.get(jobId)!;
      // If connection is already open, just use it
      if (ws.readyState === WebSocket.OPEN) {
        console.log(`[LoadingPage] Using existing connection for ${jobId}`);
        return;
      } else if (ws.readyState === WebSocket.CONNECTING) {
        console.log(`[LoadingPage] Connection already being established for ${jobId}`);
        return;
      } else {
        // Connection is closing or closed, remove it from map
        wsConnections.delete(jobId);
      }
    }
    
    // Create a new WebSocket connection
    console.log(`[LoadingPage] Creating new WebSocket connection to ${wsUrl}`);
    ws = new WebSocket(wsUrl);
    wsConnections.set(jobId, ws);
    
    // Connection established
    ws.onopen = () => {
      console.log(`[LoadingPage] WebSocket connected: ${wsUrl}`);
      setIsConnected(true);
      setReconnectAttempts(0); // Reset reconnect attempts on successful connection
      
      // Start pinging the server every 5 seconds to keep connection alive
      if (pingIntervalRef.current) {
        window.clearInterval(pingIntervalRef.current);
      }
      pingIntervalRef.current = window.setInterval(sendPing, 5000);
      
      // Send an immediate ping to verify connection
      setTimeout(sendPing, 500);
    };
    
    // Message received
    ws.onmessage = (event) => {
      // Log the raw message before any processing
      console.log('%c[LoadingPage] RAW WebSocket message received:', 'background: purple; color: white', event.data);
      
      try {
        const data: ProgressUpdate = JSON.parse(event.data);
        console.log('%c[LoadingPage] WebSocket message received:', 'color: blue', data);
        console.log('[LoadingPage] Current status:', data.status, 'Progress:', data.progress);
        
        // Debug log critical state transitions more visibly
        if (data.status === 'success' || data.status === 'complete') {
          console.log('%c[LoadingPage] SUCCESS/COMPLETE STATUS RECEIVED', 'background: green; color: white', data);
          console.log('[LoadingPage] Will this trigger navigation?', 
            'Status:', data.status, 
            'Progress:', data.progress, 
            'IsComplete check:', (data.status === 'success' || data.status === 'complete') && data.progress >= 0.99);
        }
        
        setProgressData(data);
        
        // Handle errors
        if (data.status === 'error' || data.status === 'failed') {
          const errorMsg = data.error_detail || data.message || 'Analysis failed';
          console.error('[LoadingPage] Analysis error:', errorMsg);
          setError(errorMsg);
        }
      } catch (err) {
        console.error('[LoadingPage] Error parsing WebSocket message:', err);
        console.error('[LoadingPage] Message that failed to parse:', event.data);
        setError('Error communicating with server');
      }
    };
    
    // Error handling
    ws.onerror = (event) => {
      console.error(`[LoadingPage] WebSocket error:`, event);
      
      // Implement reconnection logic
      if (reconnectAttempts < maxReconnectAttempts && !navigating) {
        const attemptNumber = reconnectAttempts + 1;
        console.log(`[LoadingPage] Will attempt to reconnect (${attemptNumber}/${maxReconnectAttempts})`);
        setReconnectAttempts(attemptNumber);
        // Remove from connections map to allow reconnect
        wsConnections.delete(jobId);
      } else {
        setError('Connection error. Please try again.');
      }
    };
    
    // Connection closed
    ws.onclose = (event) => {
      console.log(`[LoadingPage] WebSocket closed. Code: ${event.code}, Reason: ${event.reason}`);
      
      // Only show error if we didn't navigate away due to success
      if (!event.wasClean && event.code !== 1000 && !navigating) {
        if (reconnectAttempts < maxReconnectAttempts) {
          // Remove from connections map to allow reconnect
          wsConnections.delete(jobId);
        } else {
          setError('Connection closed unexpectedly. Maximum reconnection attempts reached.');
        }
      }
    };
    
    // Cleanup function when component unmounts
    return () => {
      // We don't close the connection when component unmounts to preserve it for the results page
      // However, we need to make sure we don't respond to messages anymore
      console.log(`[LoadingPage] Component unmounting, disconnecting from WebSocket events`);
      
      if (pingIntervalRef.current) {
        window.clearInterval(pingIntervalRef.current);
        pingIntervalRef.current = null;
      }
      
      if (ws) {
        // Remove our event handlers but keep connection open
        ws.onmessage = null;
        ws.onerror = null;
        ws.onclose = null;
      }
    };
  }, [jobId, reconnectAttempts, navigating]);
  
  // Function to show reconnect button after a delay if we haven't received messages
  useEffect(() => {
    // If we're connected but haven't received any progress data after 10 seconds
    // show a reconnect button
    if (isConnected && !progressData && !error) {
      const timer = setTimeout(() => {
        console.log('[LoadingPage] No progress data received after timeout, showing reconnect button');
        setShowReconnect(true);
      }, 10000); // 10 second timeout
      
      return () => clearTimeout(timer);
    }
    
    // Hide reconnect button once we get data
    if (progressData) {
      setShowReconnect(false);
    }
  }, [isConnected, progressData, error]);

  // Function to force reconnect WebSocket
  const handleForceReconnect = () => {
    console.log('[LoadingPage] Force reconnecting WebSocket');
    // Remove current connection from the map to allow reconnect
    if (jobId && wsConnections.has(jobId)) {
      const oldWs = wsConnections.get(jobId)!;
      console.log('[LoadingPage] Closing existing connection');
      // Actually close it this time
      oldWs.close();
      wsConnections.delete(jobId);
    }
    
    // Reset state
    setIsConnected(false);
    setShowReconnect(false);
    setReconnectAttempts(0);
    
    // Force component re-render to recreate the WebSocket
    setIsLoading(true);
  };
  
  // Add a timer to show a bypass button after 20 seconds
  useEffect(() => {
    const timer = setTimeout(() => {
      console.log('[LoadingPage] Analysis has been running for 20 seconds, showing bypass button');
      setShowBypass(true);
    }, 20000); // 20 seconds
    
    return () => clearTimeout(timer);
  }, []);

  // Function to bypass waiting and go directly to results
  const handleBypassToResults = () => {
    console.log('[LoadingPage] User requested to bypass waiting and go to results');
    if (jobId) {
      navigate(`/analysis/${jobId}`, { replace: true });
    }
  };
  
  // Error display
  if (error) {
    return <ErrorDisplay message="Analysis Error" details={error} />;
  }
  
  // Calculate progress percentage
  const progressPercentage = progressData 
    ? Math.round(progressData.progress * 100) 
    : 0;
    
  // Determine if we should show manual navigation button
  const showManualNavigationButton = progressData && 
    (progressData.status === 'success' || progressData.status === 'complete') && 
    progressData.progress >= 0.99 && 
    !navigating;
    
  // Format the current stage name
  const stageLabel = progressData?.stage_id 
    ? getStageMessage(progressData.stage_id)
    : progressData?.stage_name || 'Processing...';
    
  // Function to manually navigate to results page
  const handleManualNavigation = () => {
    console.log('[LoadingPage] Manual navigation requested');
    navigate(`/analysis/${jobId}`, { replace: true });
  };
  
  // Loading display with enhanced progress information
  return (
    <div className={styles.loadingPageContainer}>
      <h1 className={styles.loadingTitle}>Analysis in Progress</h1>
      
      <LoadingAnimation 
        message={progressData ? getDetailedMessage(progressData) : 'Connecting to analysis server...'}
        jobId={jobId} 
        progress={progressData?.progress}
      />
      
      {showReconnect && (
        <div className={styles.reconnectContainer}>
          <p>No updates received. There might be a connection issue.</p>
          <button 
            className={styles.reconnectButton}
            onClick={handleForceReconnect}
          >
            Reconnect
          </button>
        </div>
      )}
      
      {showBypass && (
        <div className={styles.bypassContainer}>
          <p>Taking longer than expected? You can try viewing results directly.</p>
          <button 
            className={styles.bypassButton}
            onClick={handleBypassToResults}
          >
            View Results Now
          </button>
        </div>
      )}
      
      {progressData && (
        <div className={styles.progressInfoContainer}>
          {/* Progress bar */}
          <div className={styles.progressBarContainer}>
            <div 
              className={styles.progressBar} 
              style={{ width: `${progressPercentage}%` }}
            ></div>
            <span className={styles.progressPercentage}>
              {progressPercentage}%
            </span>
          </div>
          
          {/* Stage information */}
          <div className={styles.stageInfo}>
            <span className={styles.stageLabel}>Current Stage:</span>
            <span className={styles.stageValue}>{stageLabel}</span>
          </div>
          
          {/* Current scope (if available) */}
          {progressData.current_scope && (
            <div className={styles.scopeInfo}>
              <span className={styles.scopeLabel}>Processing:</span>
              <span className={styles.scopeValue}>{progressData.current_scope}</span>
            </div>
          )}
          
          {/* Manual navigation button */}
          {showManualNavigationButton && (
            <div className={styles.manualNavigationContainer}>
              <p>Analysis complete! Automatic navigation didn't trigger.</p>
              <button 
                className={styles.manualNavigationButton}
                onClick={handleManualNavigation}
              >
                View Results
              </button>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default LoadingPage; 