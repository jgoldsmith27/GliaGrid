const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods that allow the renderer process (React app)
// to communicate with the main process (Electron backend)
contextBridge.exposeInMainWorld('electronAPI', {
  // Project management methods
  saveProject: (projectName, projectData) => ipcRenderer.invoke('save-project', projectName, projectData),
  listProjects: () => ipcRenderer.invoke('list-projects'),
  loadProject: (filePath) => ipcRenderer.invoke('load-project', filePath),
  deleteProject: (filePath, projectName) => ipcRenderer.invoke('delete-project', filePath, projectName),

  // Direct file access methods
  readBackendFile: (fileId, options) => ipcRenderer.invoke('read-backend-file', fileId, options),
  getJobMetadata: (jobId) => ipcRenderer.invoke('get-job-metadata', jobId),
  
  // Visualization data access
  readVisualizationData: (jobId, options) => ipcRenderer.invoke('read-visualization-data', jobId, options),

  // Event-driven job status API
  subscribeJobStatus: (jobId) => ipcRenderer.invoke('subscribe-job-status', jobId),
  unsubscribeJobStatus: (jobId) => ipcRenderer.invoke('unsubscribe-job-status', jobId),
  checkJobStatus: (jobId) => ipcRenderer.invoke('check-job-status', jobId),
  onJobStatusEvent: (callback) => {
    const subscription = (event, data) => callback(data);
    ipcRenderer.on('job-status-event', subscription);
    return () => ipcRenderer.removeListener('job-status-event', subscription);
  },
  onJobStatusError: (callback) => {
    const subscription = (event, data) => callback(data);
    ipcRenderer.on('job-status-error', subscription);
    return () => ipcRenderer.removeListener('job-status-error', subscription);
  }
});

console.log('preload.js executed - Electron API exposed (with project save/load and direct file access)');
