// Create this file if it doesn't exist, otherwise update it
// Configuration file for API endpoints and other global settings

// Base URL for API calls - without /api suffix
export const API_BASE_URL = 'http://localhost:8000';

// Specific API endpoints 
export const API_ENDPOINTS = {
  // Analysis API
  ANALYSIS: {
    START: `${API_BASE_URL}/analysis/start`,
    STATUS: (jobId: string) => `${API_BASE_URL}/analysis/status/${jobId}`,
  },
  // Jobs API (alternative path for results)
  JOBS: {
    STATUS: (jobId: string) => `${API_BASE_URL}/api/analysis/status/${jobId}`,
    RESULTS: (jobId: string) => `${API_BASE_URL}/api/analysis/${jobId}/results`,
  }
};

// Base URL for WebSocket connections
export const WS_BASE_URL = 'ws://localhost:8000';

// Use secure WebSocket protocol if we're on HTTPS
export const WS_PROTOCOL = window.location.protocol === 'https:' ? 'wss' : 'ws';

// Constructed WebSocket URL (allows for protocol switching)
export const getWsUrl = (path: string) => `${WS_PROTOCOL}://localhost:8000${path}`; 