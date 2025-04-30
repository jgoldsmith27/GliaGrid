const { app, BrowserWindow, ipcMain, dialog } = require('electron');
const path = require('path');
const fs = require('fs'); // Added for file system access
const fsPromises = require('fs').promises;
const Papa = require('papaparse'); // Added for CSV parsing
const zmq = require('zeromq'); // Restored for ZeroMQ IPC

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
                subscriber.send('job-status-event', { jobId, ...statusData });
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
ipcMain.handle('read-backend-file', async (event, fileId, options = {}) => {
  // Construct path to temporary uploads directory
  const backendDir = path.join(__dirname, 'backend', '.temp_uploads');
  
  try {
    // Find the file with this ID
    const files = await fsPromises.readdir(backendDir);
    const matchingFile = files.find(file => file.startsWith(fileId));
    
    if (!matchingFile) {
      throw new Error(`File with ID ${fileId} not found`);
    }
    
    const filePath = path.join(backendDir, matchingFile);
    
    // Handle different file types
    const ext = path.extname(filePath).toLowerCase();
    if (ext === '.csv') {
      // Read and parse CSV with chunking
      return readCsvChunked(filePath, options);
    } else if (ext === '.h5ad') {
      // For H5AD, we might need to use Python through child_process
      throw new Error('H5AD direct reading not implemented yet');
    }
    
    throw new Error(`Unsupported file type: ${ext}`);
  } catch (error) {
    console.error('Error reading backend file:', error);
    throw error;
  }
});

// Helper function to read CSV in chunks
async function readCsvChunked(filePath, options = {}) {
  const { offset = 0, limit = 1000, sampleRate = 1.0 } = options;
  
  return new Promise((resolve, reject) => {
    const results = [];
    let rowCount = 0;
    let skippedCount = 0;
    
    Papa.parse(fs.createReadStream(filePath), {
      header: true,
      skipEmptyLines: true,
      step: (row) => {
        rowCount++;
        
        // Apply sampling and offset/limit
        if (rowCount <= offset) return;
        if (limit && results.length >= limit) return;
        
        // Apply sampling
        if (sampleRate < 1.0 && Math.random() > sampleRate) {
          skippedCount++;
          return;
        }
        
        results.push(row.data);
      },
      complete: () => {
        resolve({
          data: results,
          totalRows: rowCount,
          sampledRows: results.length,
          skippedRows: skippedCount
        });
      },
      error: (error) => reject(error)
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
    sender.send('job-status-event', { jobId, ...initialStatus });
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

// Get visualization data directly from backend
ipcMain.handle('read-visualization-data', async (event, jobId, options = {}) => {
  try {
    // Build query parameters
    const queryParams = new URLSearchParams();
    if (options.resolution) queryParams.append('resolution', options.resolution);
    if (options.region) queryParams.append('region', JSON.stringify(options.region));
    if (options.limit) queryParams.append('limit', options.limit);
    
    const queryString = queryParams.toString() ? `?${queryParams.toString()}` : '';
    const response = await fetch(`http://localhost:8000/api/visualization/${jobId}/points${queryString}`);
    
    if (!response.ok) {
      throw new Error(`API returned status ${response.status}`);
    }
    
    return await response.json();
  } catch (error) {
    console.error(`[main.js] Error fetching visualization data for ${jobId}:`, error);
    throw error;
  }
});

// Start ZMQ subscriber when app is ready
app.whenReady().then(async () => {
  await setupZmqSubscriber();
  createWindow();
});