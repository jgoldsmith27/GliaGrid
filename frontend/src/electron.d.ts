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

declare global {
    interface Window {
      electronAPI: {
        // Project management methods
        saveProject: (projectName: string, projectData: any) => Promise<{ success: boolean; error?: string; filePath?: string; name?: string; savedAt?: string }>; // projectData should match structure expected by main process (version, files{ids}, mappings)
        listProjects: () => Promise<{ success: boolean; error?: string; projects?: ProjectListing[] }>;
        loadProject: (filePath: string) => Promise<{ success: boolean; error?: string; projectState?: SavedProjectData }>; // Returns the full saved data structure
        deleteProject: (filePath: string, projectName: string) => Promise<{ success: boolean; error?: string; filePath?: string }>;

        // Direct file access methods
        readBackendFile: (fileId: string, options: any) => Promise<any>;
        getJobMetadata: (jobId: string) => Promise<any>;
        
        // Visualization data access
        readVisualizationData: (jobId: string, options: any) => Promise<{ success: boolean; data?: any; error?: string }>;

        // Event-driven job status API
        subscribeJobStatus: (jobId: string) => Promise<{ success: boolean; error?: string }>;
        unsubscribeJobStatus: (jobId: string) => Promise<{ success: boolean; error?: string }>;
        checkJobStatus: (jobId: string) => Promise<{ success: boolean; status?: any; error?: string }>;
        onJobStatusEvent: (callback: (data: any) => void) => (() => void); // Returns cleanup function
        onJobStatusError: (callback: (data: any) => void) => (() => void); // Returns cleanup function
      };
    }
}

export {}; // Ensure this file is treated as a module 