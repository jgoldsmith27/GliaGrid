const { app, BrowserWindow, ipcMain } = require('electron');
const path = require('path');
const isDev = require('electron-is-dev'); // To check if running in development
const fs = require('fs'); // Added for file system access
const Papa = require('papaparse'); // Added for CSV parsing

function createWindow() {
  // Create the browser window.
  const win = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true, // Protect against prototype pollution
      nodeIntegration: false, // Recommended security practice
    },
  });

  // Load the index.html of the app.
  // In development, load from the Vite dev server.
  // In production, load the built frontend index.html file.
  const startUrl = isDev
    ? 'http://localhost:5173' // Default Vite port, adjust if your frontend runs elsewhere
    : `file://${path.join(__dirname, 'frontend/dist/index.html')}`; // Path to production build

  win.loadURL(startUrl);

  // Open the DevTools automatically if in development
  if (isDev) {
    win.webContents.openDevTools({ mode: 'detach' });
  }
}

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.whenReady().then(createWindow);

// Quit when all windows are closed, except on macOS. There, it's common
// for applications and their menu bar to stay active until the user quits
// explicitly with Cmd + Q.
app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  // On macOS it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow();
  }
});

// --- IPC Handlers ---

// Handle request to read file header and preview rows
ipcMain.handle('read-file-preview', async (event, fileContent, fileName) => { // Accept content and name
  console.log('Received request to read preview for:', fileName);
  if (!fileContent || !fileName) {
    return { error: 'No file content or name provided' };
  }

  // For now, only handle CSV via content
  if (fileName.toLowerCase().endsWith('.csv')) {
    try {
      // Parse the content directly
      let headers = [];
      let previewRows = [];
      let rowCount = 0;

      Papa.parse(fileContent, {
        header: true,
        preview: 6, // Read header + 5 rows
        skipEmptyLines: true,
        step: function(result, parser) {
          if (rowCount === 0) {
            headers = result.meta.fields || [];
            console.log(`[Main Process] Found Headers (${fileName}):`, headers);
          }
          if (rowCount < 5 && result.data) {
            previewRows.push(result.data);
          }
          if(rowCount >= 5) {
              parser.abort();
              console.log(`[Main Process] Aborting parse early (${fileName})`);
          }
          rowCount++;
        },
        complete: function(results) {
          console.log(`[Main Process] CSV Parsing complete (${fileName}). Headers:`, headers.length, 'Rows:', previewRows.length);
        },
        error: function(error, file) {
            console.error(`[Main Process] CSV Parsing error (${fileName}):`, error);
            // Return value will be handled by the outer check
        }
      });

      if (headers.length === 0 && previewRows.length === 0) {
         console.error(`[Main Process] No headers or preview rows obtained after parsing (${fileName}).`);
         return { error: 'Failed to parse CSV or file is empty/malformed.' };
      }

      return { headers, previewRows };
    } catch (error) {
      console.error('Error parsing CSV content:', error);
      return { error: `Failed to parse CSV content: ${error.message}` };
    }
  } else if (fileName.toLowerCase().endsWith('.h5ad')) {
    console.warn('H5AD preview via content not implemented yet.');
    // Placeholder: We'll need to send this content differently or maybe save temporarily
    // and send the path to the Python backend.
    return { error: 'H5AD preview not supported yet.', headers: [], previewRows: [] };
  } else {
    return { error: 'Unsupported file type' };
  }
});

