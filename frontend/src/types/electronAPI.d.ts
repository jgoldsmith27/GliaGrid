interface JobStatusData {
  status: 'pending' | 'running' | 'success' | 'failed' | 'error';
  progress?: number;
  message?: string;
  results?: any;
  job_id?: string;
}

interface JobStatusError {
  jobId: string;
  error: string;
}

// Define the Electron API interface globally
declare interface Window {
  electronAPI: {
    // Project management methods
    saveProject: (projectName: string, projectData: any) => Promise<{ success: boolean; error?: string; filePath?: string; name?: string; savedAt?: string }>;
    listProjects: () => Promise<{ success: boolean; error?: string; projects?: any[] }>;
    loadProject: (filePath: string) => Promise<{ success: boolean; error?: string; projectState?: any }>;
    deleteProject: (filePath: string, projectName: string) => Promise<{ success: boolean; error?: string; filePath?: string }>;

    // Direct file access methods
    readBackendFile: (fileId: string, options: any) => Promise<any>;
    getJobMetadata: (jobId: string) => Promise<any>;
    
    // Visualization data access
    readVisualizationData: (jobId: string, options: any) => Promise<{ success: boolean; data?: any; error?: string }>;

    // Event-driven job status API
    subscribeJobStatus: (jobId: string) => Promise<{ success: boolean; error?: string }>;
    unsubscribeJobStatus: (jobId: string) => Promise<{ success: boolean; error?: string }>;
    checkJobStatus: (jobId: string) => Promise<{ success: boolean; status?: JobStatusData; error?: string }>;
    onJobStatusEvent: (callback: (data: JobStatusData) => void) => (() => void);
    onJobStatusError: (callback: (data: JobStatusError) => void) => (() => void);
  };
} 