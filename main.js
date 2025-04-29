const { app, BrowserWindow, ipcMain, dialog } = require('electron');
const path = require('path');
const isDev = require('electron-is-dev'); // To check if running in development
const fs = require('fs'); // Added for file system access
const fsPromises = require('fs').promises;
const Papa = require('papaparse'); // Added for CSV parsing

// --- Project File Directory Setup ---
const projectsDir = path.join(app.getPath('userData'), 'projects');
if (!fs.existsSync(projectsDir)) {
  fs.mkdirSync(projectsDir, { recursive: true });
}
const projectFileExtension = '.gliaproj';

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
          projectFilesData.push({ filePath, name, savedAt });
        } catch (parseError) {
          console.warn(`Failed to read or parse project file ${file}:`, parseError);
          // Add with fallback data if needed
          try {
              const stat = await fsPromises.stat(filePath);
              projectFilesData.push({ filePath, name: path.basename(file, projectFileExtension), savedAt: stat.mtime.toISOString()});
          } catch (statError) {
               console.error(`Failed to stat file ${file} after parse error:`, statError);
          }
        }
      }
    }
    // Sort by save date, newest first
    projectFilesData.sort((a, b) => new Date(b.savedAt).getTime() - new Date(a.savedAt).getTime());
    return { success: true, projects: projectFilesData };
  } catch (error) {
    if (error.code === 'ENOENT') { /* ... handle not found ... */ return { success: true, projects: [] }; }
    console.error('Failed to list projects:', error);
    return { success: false, error: `Failed to list projects: ${error.message}`, projects: [] };
  }
});

// Handle loading project state from a specific file path (reads the structure saved by save-project)
ipcMain.handle('load-project', async (event, filePath) => {
  if (!filePath || !filePath.startsWith(projectsDir)) { /* ... invalid path ... */ }
  try {
    const fileContent = await fsPromises.readFile(filePath, 'utf-8');
    const projectState = JSON.parse(fileContent); // Should include projectName, savedAt, version, files, mappings
    // Basic validation
    if (!projectState || typeof projectState !== 'object' || !projectState.version || !projectState.files || !projectState.mappings) {
        throw new Error('Invalid project file format.');
    }
    console.log('Loaded project from:', filePath);
    return { success: true, projectState }; // Return the full saved object
  } catch (error) { /* ... error handling ... */ }
});

// Handle deleting a specific project file by path
ipcMain.handle('delete-project', async (event, filePath, projectName) => {
   // Basic validation: Ensure the path is within the expected directory
   if (!filePath || !filePath.startsWith(projectsDir)) {
       return { success: false, error: 'Invalid project file path provided.' };
   }

   // Use projectName in the message, fallback if not provided
   const displayProjectName = projectName || path.basename(filePath, projectFileExtension);

  // Confirmation dialog
  const { response } = await dialog.showMessageBox({
    type: 'warning',
    buttons: ['Cancel', 'Delete'], // Index 0: Cancel, Index 1: Delete
    defaultId: 0, 
    title: 'Confirm Delete',
    message: `Are you sure you want to delete project '${displayProjectName}'?`, // Use the name
    detail: `File: ${path.basename(filePath)}\nThis action cannot be undone.`, // Show filename in detail
  });

  if (response === 0) { // User clicked Cancel (Index 0)
    console.log('Delete cancelled for:', filePath);
    return { success: false, error: 'Delete cancelled by user.' }; 
  }

  // --- User clicked Delete --- 
  let projectFileDeleted = false;
  let backendFilesDeleted = true; // Assume success initially
  let backendError = null;

  try {
    // 1. Read the project file to get file IDs *before* deleting it
    let fileIdsToDelete = [];
    try {
      const fileContent = await fsPromises.readFile(filePath, 'utf-8');
      const projectState = JSON.parse(fileContent);
      // Extract non-null/undefined file IDs from the 'files' object
      if (projectState && projectState.files) {
        fileIdsToDelete = Object.values(projectState.files).filter(id => id != null);
      }
      // <<< Log extracted file IDs >>>
      console.log(`[main.js] Extracted file IDs from ${path.basename(filePath)} for deletion:`, fileIdsToDelete);
    } catch (readError) {
      console.error(`Error reading project file ${filePath} before deletion:`, readError);
      // Decide if you want to proceed with deleting the project file anyway
      // Or return an error here? Let's proceed for now, but log the warning.
      backendError = `Could not read project file to find associated data files: ${readError.message}`;
    }

    // 2. Attempt to delete associated files on the backend
    if (fileIdsToDelete.length > 0) {
        console.log(`Attempting to delete ${fileIdsToDelete.length} associated backend files...`);
        const backendUrlBase = 'http://localhost:8000'; // Make sure this is correct
        
        const deletePromises = fileIdsToDelete.map(fileId => 
            fetch(`${backendUrlBase}/api/files/${fileId}`, { method: 'DELETE' })
        );
        const results = await Promise.allSettled(deletePromises);
        
        // Use a for...of loop to allow awaiting inside
        let index = 0;
        for (const result of results) {
            const fileId = fileIdsToDelete[index];
            if (result.status === 'fulfilled') {
                // <<< Log backend response status >>>
                console.log(`[main.js] Backend DELETE response for fileId ${fileId}: Status ${result.value.status}`);
                if (result.value.ok || result.value.status === 204) { 
                    console.log(`Successfully deleted backend file for ID: ${fileId}`);
                } else if (result.value.status === 404) {
                     console.warn(`Backend file for ID ${fileId} not found (already deleted?).`);
                } else {
                    // Other backend error - Await the error text
                    backendFilesDeleted = false;
                    let errorDetail = 'Unknown backend error';
                    try {
                         // Try to get error detail from response body
                         errorDetail = await result.value.text(); 
                         // <<< Log backend error detail >>>
                         console.log(`[main.js] Backend DELETE error detail for fileId ${fileId}:`, errorDetail);
                    } catch (textError) {
                         console.error(`Could not get error text for failed delete of ${fileId}:`, textError);
                         errorDetail = `${result.value.status} ${result.value.statusText}`;
                    }
                    console.error(`[main.js] Failed backend DELETE for fileId ${fileId}: Status ${result.value.status}`);
                    backendError = backendError ? `${backendError}; ${errorDetail}` : `Backend delete failed for ${fileId}: ${errorDetail}`;
                }
            } else {
                // Fetch itself failed (network error, etc.)
                backendFilesDeleted = false;
                console.error(`[main.js] Network/fetch error deleting backend file for ID ${fileId}:`, result.reason);
                 backendError = backendError ? `${backendError}; Fetch error for ${fileId}` : `Fetch error for ${fileId}`;
            }
            index++; // Increment index for the next fileId
        }
    }

    // 3. Delete the project file itself
    await fsPromises.unlink(filePath);
    projectFileDeleted = true;
    console.log('[main.js] Deleted project configuration file:', filePath);

    // 4. Return overall status
    if (projectFileDeleted && backendFilesDeleted) {
        return { success: true, filePath };
    } else {
        // Construct a meaningful error message
        let finalError = 'Project deleted, but failed to delete some associated data files.';
        if (!projectFileDeleted) finalError = 'Failed to delete project configuration file.';
        if (backendError) finalError += ` Backend errors: ${backendError}`;
        return { success: false, error: finalError };
    }

  } catch (error) {
    // Handle errors during the main try block (e.g., project file unlink error)
    console.error(`Failed to delete project file ${filePath}:`, error);
    if (error.code === 'ENOENT') { /* ... */ }
    return { success: false, error: `Failed to delete project file: ${error.message}` };
  }
});

