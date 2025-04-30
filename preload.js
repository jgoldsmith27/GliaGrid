const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods that allow the renderer process (React app)
// to communicate with the main process (Electron backend)
contextBridge.exposeInMainWorld('electronAPI', {
  // File Upload/Management
  uploadFile: (fileType, filePath) => ipcRenderer.invoke('upload-file', fileType, filePath),
  getFilePath: (fileId) => ipcRenderer.invoke('get-file-path', fileId),
  readCsvHeaders: (fileId) => ipcRenderer.invoke('read-csv-headers', fileId),
  readCsvChunk: (fileId, options) => ipcRenderer.invoke('read-csv-chunk', fileId, options),
  readBackendFile: (fileId, options) => ipcRenderer.invoke('read-backend-file', fileId, options),
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
