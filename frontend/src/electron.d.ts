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

        // Add other methods exposed from preload.js
      };
    }
}

export {}; // Ensure this file is treated as a module 