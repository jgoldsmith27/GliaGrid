const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods that allow the renderer process (React app)
// to communicate with the main process (Electron backend)
contextBridge.exposeInMainWorld('electronAPI', {
  // Example: Function to send a message to main process and get a response
  // We will replace this with actual backend calls later
  // invoke: (channel, ...args) => ipcRenderer.invoke(channel, ...args),

  // Example: Sending a message one-way to the main process
  // send: (channel, data) => ipcRenderer.send(channel, data),

  // Example: Receiving messages from the main process
  // on: (channel, func) => {
  //   const subscription = (event, ...args) => func(...args);
  //   ipcRenderer.on(channel, subscription);
  //   // Return a function to remove the listener
  //   return () => ipcRenderer.removeListener(channel, subscription);
  // }

  // Expose the function to request file preview data
  readFilePreview: (fileContent, fileName) => ipcRenderer.invoke('read-file-preview', fileContent, fileName)
});

console.log('preload.js executed - Updated readFilePreview, removed ping');
