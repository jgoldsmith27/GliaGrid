const { app, BrowserWindow, ipcMain, dialog } = require('electron');
const path = require('path');
const fs = require('fs'); // Added for file system access
const fsPromises = require('fs').promises;
const Papa = require('papaparse'); // Added for CSV parsing
const zmq = require('zeromq'); // Restored for ZeroMQ IPC
const readline = require('readline'); // Added for line-by-line reading

// --- Cancellation Tracking for File Reads ---
const activeReadRequests = new Map(); // Map<requestId, isCancelled>

// --- Project File Directory Setup ---
const projectsDir = path.join(app.getPath('userData'), 'projects');
if (!fs.existsSync(projectsDir)) {
  fs.mkdirSync(projectsDir, { recursive: true });
}
const projectFileExtension = '.gliaproj';

// --- Backend Integration Setup ---
// Path to the backend's temporary uploads directory
const backendTempDir = path.join(__dirname, 'backend', '.temp_uploads');

// Setup ZeroMQ subscriber for direct IPC
const ZMQ_IPC_ADDRESS = "ipc:///tmp/gliagrid-job-status";
let zmqSocket = null;

// Add to main.js - Event-driven job status tracking
const activeJobSubscriptions = new Map(); // Track which jobs have active subscriptions

// Function to set up ZeroMQ subscriber
async function setupZmqSubscriber() {
  try {
    console.log(`[main.js] Setting up ZeroMQ subscriber...`);
    zmqSocket = new zmq.Subscriber();
    
    // Subscribe to all job status updates
    zmqSocket.subscribe('job.');
    
    await zmqSocket.connect(ZMQ_IPC_ADDRESS);
    console.log(`[main.js] ZeroMQ subscriber connected to ${ZMQ_IPC_ADDRESS}`);
    
    // Set up message processing
    (async () => {
      console.log(`[main.js] Starting ZeroMQ message listener`);
      for await (const [topic, message] of zmqSocket) {
        try {
          const topicStr = topic.toString();
          const jobId = topicStr.substring(4); // Remove 'job.' prefix
          const statusData = JSON.parse(message.toString());
          
          console.log(`[main.js] Received ZeroMQ message for job ${jobId}: ${statusData.status}`);
          
          // Only process if we have active subscribers for this job
          if (activeJobSubscriptions.has(jobId)) {
            const subscribers = activeJobSubscriptions.get(jobId);
            
            // Send update to all subscribers
            for (const subscriber of subscribers) {
              if (!subscriber.isDestroyed()) {
                subscriber.send('job-update', { jobId, ...statusData });
              } else {
                subscribers.delete(subscriber);
              }
            }
            
            // If job is complete or failed, remove from active subscriptions
            if (statusData.status === 'success' || statusData.status === 'failed' || statusData.status === 'error') {
              console.log(`[main.js] Job ${jobId} completed with status: ${statusData.status}. Removing from subscription list.`);
              activeJobSubscriptions.delete(jobId);
            }
          }
        } catch (err) {
          console.error(`[main.js] Error processing ZeroMQ message: ${err}`);
        }
      }
    })().catch(err => console.error(`[main.js] ZeroMQ listener error: ${err}`));
    
    return true;
  } catch (err) {
    console.error(`[main.js] Failed to setup ZeroMQ subscriber: ${err}`);
    return false;
  }
}

async function createWindow() {
  // Dynamically import isDev
  const { default: isDev } = await import('electron-is-dev');
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

// --- Revised Project File Handling ---

// Handle saving a project with a user-provided name
ipcMain.handle('save-project', async (event, projectName, projectData) => {
  if (!projectName || typeof projectName !== 'string' || projectName.trim().length === 0) {
    return { success: false, error: 'Invalid project name provided.' };
  }
  // Sanitize projectName slightly to make a safe filename
  const safePart = projectName.trim().replace(/[^a-zA-Z0-9_-]/g, '_').substring(0, 50);
  const fileName = `${safePart}-${Date.now()}${projectFileExtension}`;
  const filePath = path.join(projectsDir, fileName);

  // Add projectName to the data being saved
  const dataToSave = {
    projectName: projectName.trim(),
    savedAt: new Date().toISOString(),
    ...projectData // Should contain version, files (with fileIds), mappings
  };

  // <<< Log before saving >>>
  console.log(`[main.js] Attempting to save project '${projectName}' to path: ${filePath}`);
  console.log('[main.js] Project data being saved:', JSON.stringify(dataToSave, null, 2));

  try {
    const fileContent = JSON.stringify(dataToSave, null, 2);
    await fsPromises.writeFile(filePath, fileContent, 'utf-8');
    console.log('[main.js] Project saved successfully:', filePath);
    // Return info needed to update the list immediately
    return { success: true, filePath, name: dataToSave.projectName, savedAt: dataToSave.savedAt }; 
  } catch (error) {
    console.error(`[main.js] Failed to save project to ${filePath}:`, error);
    return { success: false, error: `Failed to save project file: ${error.message}` };
  }
});

// Handle listing saved projects (parsing name/date from files)
ipcMain.handle('list-projects', async () => {
  let projectFilesData = [];
  try {
    const files = await fsPromises.readdir(projectsDir);
    for (const file of files) {
      if (path.extname(file).toLowerCase() === projectFileExtension) {
        const filePath = path.join(projectsDir, file);
        try {
          const fileContent = await fsPromises.readFile(filePath, 'utf-8');
          const projectData = JSON.parse(fileContent);
          // Use saved name/date if available, otherwise fallback
          const name = projectData?.projectName || path.basename(file, projectFileExtension);
          const savedAt = projectData?.savedAt || (await fsPromises.stat(filePath)).mtime.toISOString();
          
          projectFilesData.push({
            name,
            savedAt,
            filePath
          });
        } catch (error) {
          console.error(`[main.js] Error reading project file ${filePath}:`, error);
          // Still include file but with error info
          projectFilesData.push({
            name: path.basename(file, projectFileExtension),
            savedAt: (await fsPromises.stat(filePath)).mtime.toISOString(),
            filePath,
            error: `Could not parse project file: ${error.message}`
          });
        }
      }
    }
    
    // Sort by newest first
    projectFilesData.sort((a, b) => new Date(b.savedAt) - new Date(a.savedAt));
    
    return { success: true, projects: projectFilesData };
  } catch (error) {
    console.error('[main.js] Error listing projects:', error);
    return { success: false, error: error.message, projects: [] };
  }
});

// Handle loading a specific project file
ipcMain.handle('load-project', async (event, filePath) => {
  try {
    const content = await fsPromises.readFile(filePath, 'utf-8');
    const projectData = JSON.parse(content);
    return { success: true, projectState: projectData };
  } catch (error) {
    console.error(`[main.js] Failed to load project from ${filePath}:`, error);
    return { success: false, error: `Failed to load project: ${error.message}` };
  }
});

// Handle deleting a specific project file
ipcMain.handle('delete-project', async (event, filePath, projectName) => {
  try {
    await fsPromises.unlink(filePath);
    console.log(`[main.js] Deleted project file: ${filePath}`);
    return { success: true };
  } catch (error) {
    console.error(`[main.js] Failed to delete project ${filePath}:`, error);
    return { success: false, error: `Failed to delete project: ${error.message}` };
  }
});

// --- Direct File Access (Direct IPC) ---

// Listener for cancellation requests from preload
ipcMain.on('cancel-read-request', (event, requestId) => {
  if (activeReadRequests.has(requestId)) {
    console.log(`[main.js] Received cancel request for read ID: ${requestId}`);
    activeReadRequests.set(requestId, true); // Mark as cancelled
  } else {
    console.warn(`[main.js] Received cancel request for unknown/completed read ID: ${requestId}`);
  }
});

ipcMain.handle('read-backend-file', async (event, fileId, options = {}, requestId) => {
  // Check if a requestId was provided (should be from preload)
  if (!requestId) {
    console.error('[main.js] read-backend-file invoked without a requestId.');
    return { success: false, error: 'Internal error: Missing request ID for cancellation.' };
  }

  // Construct path to temporary uploads directory
  const backendDir = path.join(__dirname, 'backend', '.temp_uploads');
  
  // Register this request for cancellation tracking
  activeReadRequests.set(requestId, false);
  // Define the cancellation check function for this request
  const isCancelled = () => activeReadRequests.get(requestId) === true;
  
  try {
    // Find the file with this ID
    const files = await fsPromises.readdir(backendDir);
    const matchingFile = files.find(file => file.startsWith(fileId));
    
    if (!matchingFile) {
      throw new Error(`File with ID ${fileId} not found`);
    }
    
    const filePath = path.join(backendDir, matchingFile);
    
    // Check for cancellation before starting potentially long read
    if (isCancelled()) {
      console.log(`[main.js] Read request ${requestId} cancelled before starting file read.`);
      throw new Error('Operation cancelled'); // Throw cancellation error
    }

    // Handle different file types
    const ext = path.extname(filePath).toLowerCase();
    if (ext === '.csv') {
      // Read and parse CSV, passing the cancellation checker
      return await readCsvChunked(filePath, options, isCancelled);
    } else if (ext === '.h5ad') {
      // For H5AD, we might need to use Python through child_process
      // Cancellation for child_process would need specific handling (e.g., process.kill)
      throw new Error('H5AD direct reading cancellation not implemented yet');
    }
    
    throw new Error(`Unsupported file type: ${ext}`);
  } catch (error) {
    console.error(`[main.js] Error reading backend file (Req ID: ${requestId}):`, error);
    // If the error is our specific cancellation error, return it structured
    if (error.message === 'Operation cancelled') {
        return { success: false, error: error.message, cancelled: true };
    }
    // Otherwise, return a generic error structure
    return { success: false, error: error.message };
  } finally {
    // Always clean up the request ID from the tracking map
    activeReadRequests.delete(requestId);
    console.log(`[main.js] Cleaned up read request ID: ${requestId}`);
  }
});

// --- Helper function for reading CSV with filtering and cancellation ---
async function readCsvChunked(filePath, options = {}, isCancelled) {
  console.log(`[main.js] readCsvChunked called for ${filePath} with options:`, options);
  const { filter, columns } = options; // Destructure options
  const data = [];
  const warnings = [];
  let returnedRows = 0;
  // Note: totalRows calculation might be inefficient for large files without reading fully
  // Consider skipping totalRows or finding a faster way if performance is critical.

  return new Promise((resolve, reject) => {
    const fileStream = fs.createReadStream(filePath, { encoding: 'utf-8' });
    const lineReader = readline.createInterface({
      input: fileStream,
      crlfDelay: Infinity
    });

    let headers = [];
    let headerProcessed = false;
    let filterColIndex = -1;
    let columnIndices = [];

    lineReader.on('line', (line) => {
      // *** Check for cancellation frequently ***
      if (isCancelled()) {
        console.log(`[main.js] Cancellation detected during CSV read of ${filePath}`);
        lineReader.close(); // Close the reader
        fileStream.destroy(); // Destroy the stream
        // Reject the promise with a specific cancellation error
        reject(new Error('Operation cancelled'));
        return; 
      }
      
      if (!headerProcessed) {
        headers = line.split(',').map(h => h.trim()); // Basic CSV split, consider PapaParse for robustness
        headerProcessed = true;

        // Determine indices based on headers
        if (filter && filter.column) {
          filterColIndex = headers.indexOf(filter.column);
          if (filterColIndex === -1) {
            warnings.push(`Filter column "${filter.column}" not found.`);
          }
        }
        if (columns && columns.length > 0) {
          columnIndices = columns.map(col => {
            const index = headers.indexOf(col);
            if (index === -1) {
              warnings.push(`Requested column "${col}" not found.`);
            }
            return index;
          }).filter(index => index !== -1); // Only keep valid indices
        } else {
          // If no specific columns requested, use all columns
          columnIndices = headers.map((_, index) => index);
        }
        return; // Skip processing header row as data
      }

      const values = line.split(',').map(v => v.trim());
      
      // Apply row filter if specified and valid
      if (filter && filter.values && filterColIndex !== -1) {
        const cellValue = values[filterColIndex];
        // Simple includes check, adjust if type conversion or specific logic needed
        if (!filter.values.includes(cellValue)) {
          return; // Skip row if it doesn't match filter
        }
      }

      // Construct row object with requested columns
      const rowObject = {};
      for (const index of columnIndices) {
        if (index < values.length) { // Check bounds
          rowObject[headers[index]] = values[index];
        } else {
            rowObject[headers[index]] = undefined; // Handle case where row has fewer columns than header
        }
      }
      data.push(rowObject);
      returnedRows++;

      // Optional: Add limits or sampling logic here if needed
    });

    lineReader.on('close', () => {
       // Check cancellation one last time after closing (in case it happened between last line and close)
      if (isCancelled()) {
           console.log(`[main.js] Cancellation detected just before resolving CSV read of ${filePath}`);
           reject(new Error('Operation cancelled'));
           return;
      }
      console.log(`[main.js] Finished reading ${returnedRows} rows from ${filePath}`);
      resolve({ 
          success: true, 
          data: data, 
          returnedRows: returnedRows, 
          // totalRows: could be estimated or determined differently if needed
          warnings: warnings 
        });
    });

    lineReader.on('error', (err) => {
      console.error(`[main.js] Error reading CSV file ${filePath}:`, err);
      // Check if the error is due to cancellation (might be wrapped)
      if (err.message === 'Operation cancelled') {
           reject(err); // Propagate cancellation error
      } else {
           reject(new Error(`Failed to read CSV: ${err.message}`)); // Wrap other errors
      }
    });
    
    fileStream.on('error', (err) => {
       console.error(`[main.js] File stream error for ${filePath}:`, err);
       lineReader.close(); // Ensure linereader is closed on stream error
        if (err.message === 'Operation cancelled') {
           reject(err);
       } else {
           reject(new Error(`File stream error: ${err.message}`));
       }
    });

  });
}

// --- Direct Job Status IPC ---

// Add handlers for job status management
ipcMain.handle('subscribe-job-status', async (event, jobId) => {
  const sender = event.sender;
  
  if (!activeJobSubscriptions.has(jobId)) {
    activeJobSubscriptions.set(jobId, new Set());
  }
  
  activeJobSubscriptions.get(jobId).add(sender);
  console.log(`[main.js] Subscribed to job status for ${jobId}`);
  
  // Ensure ZMQ subscriber is setup
  if (!zmqSocket) {
    await setupZmqSubscriber();
  }
  
  // Immediately fetch current status and send it
  try {
    const initialStatus = await fetchJobStatus(jobId);
    sender.send('job-update', { jobId, ...initialStatus });
  } catch (error) {
    console.error(`[main.js] Error fetching initial job status: ${error}`);
    sender.send('job-status-error', { jobId, error: error.message });
  }
  
  // Return success
  return { success: true };
});

ipcMain.handle('unsubscribe-job-status', (event, jobId) => {
  const sender = event.sender;
  
  if (activeJobSubscriptions.has(jobId)) {
    activeJobSubscriptions.get(jobId).delete(sender);
    
    // Clean up if no more subscribers
    if (activeJobSubscriptions.get(jobId).size === 0) {
      activeJobSubscriptions.delete(jobId);
    }
    
    console.log(`[main.js] Unsubscribed from job status for ${jobId}`);
  }
  
  return { success: true };
});

// Direct job status check using backend API
ipcMain.handle('check-job-status', async (event, jobId) => {
  try {
    const status = await fetchJobStatus(jobId);
    return { success: true, status };
  } catch (error) {
    console.error(`[main.js] Error checking job status for ${jobId}:`, error);
    return { success: false, error: error.message };
  }
});

// Helper function to fetch job status from API
async function fetchJobStatus(jobId) {
  const response = await fetch(`http://localhost:8000/api/analysis/status/${jobId}`);
  
  if (!response.ok) {
    throw new Error(`Failed to fetch job status: ${response.statusText}`);
  }
  
  return await response.json();
}

// --- Reusable Internal Function for Metadata Fetching ---
async function getJobMetadataInternal(jobId) {
  console.log(`[main.js Internal] Fetching metadata for job ${jobId}`);
  const apiUrl = `http://localhost:8000/api/analysis/context/${jobId}`;
  try {
    const response = await fetch(apiUrl);
    if (!response.ok) {
      let errorDetail = `API returned status ${response.status}`;
      try {
         const errorBody = await response.json();
         errorDetail = errorBody.detail || errorDetail;
      } catch (_) { /* Ignore */ }
      throw new Error(errorDetail + ` fetching context for job ${jobId}`);
    }
    const contextResponse = await response.json();
    if (!contextResponse || contextResponse.context === undefined) {
      throw new Error(`Context field missing in response from ${apiUrl}`);
    }
    console.log(`[main.js Internal] Successfully fetched context for job ${jobId}`);
    // Check the structure before returning
    if (typeof contextResponse.context !== 'object' || contextResponse.context === null) {
        console.warn(`[main.js Internal] Context received for job ${jobId} is not a valid object:`, contextResponse.context);
        // Decide how to handle - return error or empty object?
        // Let's return an error for now, as the structure is expected.
        throw new Error(`Invalid context structure received from backend for job ${jobId}`);
    }
    return { success: true, metadata: contextResponse.context };
  } catch (error) {
    console.error(`[main.js Internal] Error fetching metadata for job ${jobId}:`, error);
    return { success: false, error: error.message };
  }
}

// --- Get Job Metadata/Context (IPC Handler) ---
ipcMain.handle('get-job-metadata', async (event, jobId) => {
  console.log(`[main.js Handler] Received request for metadata for job ${jobId}`);
  // Call the internal reusable function
  return await getJobMetadataInternal(jobId);
});

// Start ZMQ subscriber when app is ready
app.whenReady().then(async () => {
  await setupZmqSubscriber();
  createWindow();
});