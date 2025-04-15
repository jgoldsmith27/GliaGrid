const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods that allow the renderer process (React app)
// to communicate with the main process (Electron backend)
contextBridge.exposeInMainWorld('electronAPI', {
  // Add any future Electron-specific APIs here
  // For now, we're using direct fetch calls to the backend API
});

console.log('preload.js executed - Backend communication now uses fetch API');
