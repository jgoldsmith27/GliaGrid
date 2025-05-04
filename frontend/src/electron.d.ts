// Assuming ColumnMappings is defined elsewhere (e.g., types/dataTypes or types/projectTypes)
// If not, define it here or import it.
// interface ColumnMappings { ... }

// Structure of the data *saved* in the .gliaproj file
interface SavedProjectData {
    projectName: string;
    savedAt: string; // ISO date string
    version: string;
    files: { // Contains fileIds
      expression?: string; 
      metadata?: string;   
      modules?: string;    
    };
    mappings: any; // Use specific ColumnMappings type if available
}

// Structure returned by listProjects
interface ProjectListing {
    filePath: string;
    name: string;
    savedAt: string; // ISO date string
}

// Define Job Status type (align with backend/ZMQ messages)
interface JobStatusUpdate { 
    jobId: string;
    status: 'pending' | 'running' | 'success' | 'failed' | 'error';
    message?: string;
    progress?: number;
    results?: any; 
    job_id?: string; // Allow for backend key variation
}

declare global {
    interface Window {
      electronAPI: {
        // --- File Upload/Management (Matches preload.js) ---
        uploadFile: (fileType: string, filePath: string) => Promise<any>;
        getFilePath: (fileId: string) => Promise<any>;
        readCsvHeaders: (fileId: string) => Promise<{ success: boolean; headers?: string[]; error?: string }>;
        readCsvChunk: (fileId: string, options: any) => Promise<{ success: boolean; data?: any[]; headers?: string[]; stats?: any; error?: string }>;
        readBackendFile: (fileId: string, options: any, signal?: AbortSignal) => Promise<{ success: boolean; data?: any[]; returnedRows?: number; totalRows?: number; warnings?: string[]; error?: string }>;
        removeFile: (fileId: string) => Promise<any>;
        
        // --- Job Management & Status (Matches preload.js) ---
        subscribeJobStatus: (jobId: string) => Promise<any>; // Returns { success: boolean }
        unsubscribeJobStatus: (jobId: string) => Promise<any>; // Returns { success: boolean }
        checkJobStatus: (jobId: string) => Promise<{ success: boolean; status?: any; error?: string }>; // Returns full status payload on success
        // Listener for ZMQ/Backend job updates (Matches preload.js)
        onJobUpdate: (callback: (jobStatus: JobStatusUpdate) => void) => () => void; // Returns cleanup function
        
        // --- Metadata (Matches preload.js) ---
        getJobMetadata: (jobId: string) => Promise<{ success: boolean; metadata?: any; error?: string }>;
        
        // --- Backend Process Management (Example - Matches preload.js) ---
        startBackend: () => Promise<any>;
        stopBackend: () => Promise<any>;
        getBackendStatus: () => Promise<any>;
        onBackendLog: (callback: (message: string) => void) => () => void; // Returns cleanup function

        // --- REMOVED Old/Unused Functions --- 
        // saveProject: ... (Now handled by backend API or specific IPC if needed differently)
        // listProjects: ... 
        // loadProject: ...
        // deleteProject: ...
        // readVisualizationData: ... (Likely replaced by direct file reads or metadata)
        // onJobStatusEvent: ... (Replaced by onJobUpdate via ZMQ)
        // onJobStatusError: ... (Errors handled via ZMQ or direct calls)
        // streamSpatialData: ... (Removed)
        // onSpatialDataChunk: ... (Removed)
        // onSpatialDataEnd: ... (Removed)
        // onSpatialDataError: ... (Removed)
      };
    }
}

export {}; // Ensure this file is treated as a module 