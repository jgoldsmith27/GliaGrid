const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods that allow the renderer process (React app)
// to communicate with the main process (Electron backend)
contextBridge.exposeInMainWorld('electronAPI', {
  // Project management methods
  saveProject: (projectName, projectData) => ipcRenderer.invoke('save-project', projectName, projectData),
  listProjects: () => ipcRenderer.invoke('list-projects'),
  loadProject: (filePath) => ipcRenderer.invoke('load-project', filePath),
  deleteProject: (filePath) => ipcRenderer.invoke('delete-project', filePath),

  // Other APIs if needed (e.g., readFilePreview if it exists)
  // readFilePreview: (fileContent, fileName) => ipcRenderer.invoke('read-file-preview', fileContent, fileName)
});

console.log('preload.js executed - Electron API exposed (with project save/load)');
