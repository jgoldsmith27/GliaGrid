const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods that allow the renderer process (React app)
// to communicate with the main process (Electron backend)
contextBridge.exposeInMainWorld('electronAPI', {
  // File Upload/Management
  uploadFile: (fileType, filePath) => ipcRenderer.invoke('upload-file', fileType, filePath),
  getFilePath: (fileId) => ipcRenderer.invoke('get-file-path', fileId),
  readCsvHeaders: (fileId) => ipcRenderer.invoke('read-csv-headers', fileId),
  readCsvChunk: (fileId, options) => ipcRenderer.invoke('read-csv-chunk', fileId, options),
  readBackendFile: (fileId, options, signal) => {
    // Generate a unique ID for this request for cancellation tracking
    const requestId = `readReq-${Date.now()}-${Math.random().toString(36).substring(2, 9)}`;

    // Create a promise that resolves/rejects based on the IPC invoke
    const invokePromise = ipcRenderer.invoke('read-backend-file', fileId, options, requestId);

    // Handle the AbortSignal if provided
    let abortListener = null;
    if (signal) {
      abortListener = () => {
        console.log(`[preload] Abort signal received for request: ${requestId}. Sending cancel IPC.`);
        // Send a message to the main process to cancel this specific request
        ipcRenderer.send('cancel-read-request', requestId);
      };

      // If already aborted before we even start, trigger immediately
      if (signal.aborted) {
        abortListener();
      } else {
        signal.addEventListener('abort', abortListener);
      }
    }

    // Clean up the listener when the invoke completes or is aborted
    const cleanup = () => {
      if (signal && abortListener) {
        signal.removeEventListener('abort', abortListener);
      }
    };

    invokePromise.then(cleanup, cleanup); // Cleanup on resolve or reject

    return invokePromise; // Return the original promise from invoke
  },
  removeFile: (fileId) => ipcRenderer.invoke('remove-file', fileId),
  
  // Project Management
  saveProject: (projectName, projectData) => ipcRenderer.invoke('save-project', projectName, projectData),
  listProjects: () => ipcRenderer.invoke('list-projects'),
  loadProject: (filePath) => ipcRenderer.invoke('load-project', filePath),
  deleteProject: (filePath, projectName) => ipcRenderer.invoke('delete-project', filePath, projectName),
  
  // Job Management & Status
  subscribeJobStatus: (jobId) => ipcRenderer.invoke('subscribe-job-status', jobId),
  unsubscribeJobStatus: (jobId) => ipcRenderer.invoke('unsubscribe-job-status', jobId),
  checkJobStatus: (jobId) => ipcRenderer.invoke('check-job-status', jobId),
  // Listener for ZMQ/Backend job updates
  onJobUpdate: (callback) => {
    const subscription = (event, jobStatus) => callback(jobStatus);
    ipcRenderer.on('job-update', subscription);
    return () => ipcRenderer.removeListener('job-update', subscription); // Return cleanup function
  },
  
  // Metadata
  getJobMetadata: (jobId) => ipcRenderer.invoke('get-job-metadata', jobId),
  
  // Backend Process Management (Example)
  startBackend: () => ipcRenderer.invoke('start-backend'),
  stopBackend: () => ipcRenderer.invoke('stop-backend'),
  getBackendStatus: () => ipcRenderer.invoke('get-backend-status'),
  onBackendLog: (callback) => {
    const subscription = (event, message) => callback(message);
    ipcRenderer.on('backend-log', subscription);
    return () => ipcRenderer.removeListener('backend-log', subscription);
  },
  
  // Spatial Data Streaming - REMOVED
  // streamSpatialData: (jobId, options) => ipcRenderer.send('stream-spatial-data', jobId, options),
  // // Listener for stream chunks
  // onSpatialDataChunk: (callback) => { ... },
  // // Listener for stream end
  // onSpatialDataEnd: (callback) => { ... },
  // // Listener for stream errors
  // onSpatialDataError: (callback) => { ... }
});

console.log('preload.js executed - Electron API exposed (with project save/load, file access, and streaming)');
