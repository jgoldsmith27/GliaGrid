import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import styles from './DataInputPage.module.css';
import FileUploader from '../../components/FileUploader/FileUploader';
import ColumnMapper from '../../components/ColumnMapper/ColumnMapper';
import DataPreview from '../../components/DataPreview/DataPreview';
import { FileState, FileType, FilePreviewResult, requiredColumns, mappingFields, AnalysisPayload, AnalysisMapping } from '../../types/dataTypes';
import { SharedDataStore } from '../../services/data/SharedDataStore';

// --- Define Interfaces Locally ---
interface ColumnMappings { 
  expression?: { 
    geneCol?: string;
    xCol?: string;
    yCol?: string;
    layerCol?: string;
  };
  metadata?: {
    ligandCol?: string;
    receptorCol?: string;
  };
  modules?: {
    geneCol?: string;
    moduleCol?: string;
  };
}

interface SavedProjectData { 
    projectName: string;
    savedAt: string; 
    version: string;
    files: { 
      expression?: string; 
      metadata?: string;   
      modules?: string;    
    };
    mappings: ColumnMappings; // Use the locally defined type
}

interface ProjectListing { 
    filePath: string;
    name: string;
    savedAt: string; 
}
// --- End Local Interface Definitions ---

const projectFileExtension = '.gliaproj';

const DataInputPage: React.FC = () => {
  const navigate = useNavigate();
  const dataStore = SharedDataStore.getInstance();

  // Helper to create initial file state
  const createInitialFileState = (): FileState => ({
    file: null,
    fileId: null, // Added fileId
    headers: [],
    previewRows: [],
    error: null,
    isLoading: false,
    geneCol: undefined, xCol: undefined, yCol: undefined, layerCol: undefined,
    ligandCol: undefined, receptorCol: undefined, moduleCol: undefined
  });
  
  const [spatialFile, setSpatialFile] = useState<FileState>(createInitialFileState());
  const [interactionsFile, setInteractionsFile] = useState<FileState>(createInitialFileState());
  const [modulesFile, setModulesFile] = useState<FileState>(createInitialFileState());
  const [analysisStatus, setAnalysisStatus] = useState<string>('');
  const [isAnalyzing, setIsAnalyzing] = useState<boolean>(false);
  const [lastJobId, setLastJobId] = useState<string | null>(null);
  const [savedProjects, setSavedProjects] = useState<ProjectListing[]>([]); 
  const [projectListStatus, setProjectListStatus] = useState<string>('Loading projects...');

  // --- State for Project Name Input ---
  const [showProjectNameInput, setShowProjectNameInput] = useState<boolean>(false);
  const [projectName, setProjectName] = useState<string>('');

  // --- Fetch Saved Projects on Mount --- 
  useEffect(() => {
    const fetchProjects = async () => {
      setProjectListStatus('Loading project history...');
      try {
        // Call the modified listProjects
        const result = await window.electronAPI.listProjects(); 
        if (result.success && result.projects) {
          setSavedProjects(result.projects);
          setProjectListStatus(result.projects.length > 0 ? '' : 'No saved projects found.');
        } else {
          setSavedProjects([]);
          setProjectListStatus(`Error loading projects: ${result.error || 'Unknown error'}`);
          console.error("List Projects Error:", result.error);
        }
      } catch (error) {
        setSavedProjects([]);
        setProjectListStatus(`Error loading projects: ${error instanceof Error ? error.message : String(error)}`);
        console.error("List Projects IPC Error:", error);
      }
    };
    fetchProjects();
  }, []);

  const updateFileState = (
    setter: React.Dispatch<React.SetStateAction<FileState>>,
    update: Partial<FileState>
  ) => {
    setter(prev => ({ ...prev, ...update }));
  };

  const handleFileChange = async (event: React.ChangeEvent<HTMLInputElement>, type: FileType) => {
    const setter = type === 'spatial' ? setSpatialFile :
                   type === 'interactions' ? setInteractionsFile :
                   setModulesFile;

    updateFileState(setter, createInitialFileState()); // Reset completely on new selection

    if (event.target.files && event.target.files[0]) {
      const file = event.target.files[0];
      const fileName = file.name;

      if (!fileName.endsWith('.csv') && !fileName.endsWith('.h5ad')) {
        alert('Please select a .csv or .h5ad file.');
        event.target.value = '';
        return;
      }

      updateFileState(setter, { file, isLoading: true });
      console.log(`Processing file for ${type}:`, fileName);

      try {
        // Create a FormData object to send the file to the backend
        const formData = new FormData();
        formData.append('file', file);

        // Call the backend API to process the file
        const response = await fetch('http://localhost:8000/api/files/preview', {
          method: 'POST',
          body: formData,
        });

        if (!response.ok) {
          const errorText = await response.text();
          throw new Error(`Server error: ${response.status} - ${errorText}`);
        }

        const result: FilePreviewResult = await response.json(); // Ensure type
        console.log(`[React] Preview result received for ${type}:`, result);

        if (result.error) {
          updateFileState(setter, { error: result.error, isLoading: false });
        } else if (result.headers && result.previewRows && result.fileId) { // Check for fileId
          updateFileState(setter, {
            fileId: result.fileId, // Store fileId
            headers: result.headers,
            previewRows: result.previewRows,
            fileInfo: result.fileInfo, // Store file info if available (for H5AD)
            error: null,
            isLoading: false,
            // Reset mappings
            geneCol: undefined, xCol: undefined, yCol: undefined, layerCol: undefined,
            ligandCol: undefined, receptorCol: undefined, moduleCol: undefined
          });
          
          // Auto-detect and map columns for H5AD files if possible
          if (fileName.endsWith('.h5ad') && result.headers) {
            const headerMap: Record<string, string> = {};
            
            // Try to auto-map common column names
            result.headers.forEach((header: string) => {
              const lowerHeader = header.toLowerCase();
              
              // Match column names based on common naming conventions
              if (type === 'spatial') {
                if (lowerHeader.includes('gene') || lowerHeader === 'index') 
                  headerMap.geneCol = header;
                else if (lowerHeader === 'x' || lowerHeader.includes('coord_x')) 
                  headerMap.xCol = header;
                else if (lowerHeader === 'y' || lowerHeader.includes('coord_y')) 
                  headerMap.yCol = header;
                else if (lowerHeader.includes('layer') || lowerHeader.includes('cluster')) 
                  headerMap.layerCol = header;
              } else if (type === 'interactions') {
                if (lowerHeader.includes('ligand')) 
                  headerMap.ligandCol = header;
                else if (lowerHeader.includes('receptor')) 
                  headerMap.receptorCol = header;
              } else if (type === 'modules') {
                if (lowerHeader.includes('gene') || lowerHeader === 'index') 
                  headerMap.geneCol = header;
                else if (lowerHeader.includes('module') || lowerHeader.includes('cluster')) 
                  headerMap.moduleCol = header;
              }
            });
            
            // Apply the auto-mappings
            if (Object.keys(headerMap).length > 0) {
              console.log(`[React] Auto-mapped columns for ${type}:`, headerMap);
              updateFileState(setter, headerMap);
            }
          }
        } else {
          updateFileState(setter, { error: 'Invalid response from server.', isLoading: false });
        }
      } catch (error) {
        const errorMessage = error instanceof Error ? error.message : String(error);
        updateFileState(setter, { error: `Error processing file: ${errorMessage}`, isLoading: false });
      }

      // Clear the input value so the same file can be selected again if needed
      event.target.value = '';
    }
  };

  const handleMappingChange = (
    event: React.ChangeEvent<HTMLSelectElement>,
    type: FileType,
    columnKey: keyof FileState // Revert to expecting keyof FileState
  ) => {
    const setter = type === 'spatial' ? setSpatialFile :
                   type === 'interactions' ? setInteractionsFile :
                   setModulesFile;
    const selectedValue = event.target.value;
    updateFileState(setter, { [columnKey]: selectedValue }); // This is type-safe now
  };

  const handleRemoveFile = (type: FileType) => {
    const setter = type === 'spatial' ? setSpatialFile :
                   type === 'interactions' ? setInteractionsFile :
                   setModulesFile;
    // Reset the state for this file type using the helper
    updateFileState(setter, createInitialFileState());
  };

  // Check if all required columns for a file type are mapped
  const checkRequiredMapping = (state: FileState, type: FileType): boolean => {
    if (!state.file || !state.fileId) return false; // Ensure file is loaded AND has an ID
    const reqCols = requiredColumns[type] || [];
    return reqCols.every(key => !!state[key as keyof FileState]);
  };

  const canProcess = checkRequiredMapping(spatialFile, 'spatial') &&
                     checkRequiredMapping(interactionsFile, 'interactions') &&
                     checkRequiredMapping(modulesFile, 'modules');

  // Convert FileState and mapping config to ColumnMapper props format
  const getColumnMapperProps = (fileState: FileState, type: FileType) => {
    if (!fileState.file) return null;
    
    // Keep track of the actual keys alongside the IDs for the callback
    const fieldsWithKeys = mappingFields[type].map(({ key, label }) => {
      return {
        id: key, // The ID used for the component key/mapping
        key: key, // The actual keyof FileState
        label,
        required: requiredColumns[type].includes(key)
      };
    });

    // Get current mappings using the key
    const currentMappings: Record<string, string> = {};
    // Iterate only over the actual mapping keys defined for the current file type
    mappingFields[type].forEach(({ key }) => { 
      // Check if the key corresponds to a mapping property (like 'geneCol', 'xCol', etc.)
      if (key in fileState) {
        const value = fileState[key]; // Access state directly with the key
        // Ensure the value is a string before assigning (it should be for mapping keys)
        if (typeof value === 'string') { 
            currentMappings[key] = value; // Map using the key as ID
        }
      }
    });
    
    return {
      headers: fileState.headers,
      isLoading: fileState.isLoading,
      error: fileState.error,
      // Pass only id, label, required to ColumnMapper props if it expects that shape
      requiredFields: fieldsWithKeys.map(({ id, label, required }) => ({ id, label, required })), 
      mappings: currentMappings,
      onMappingChange: (fieldId: string, selectedColumn: string) => {
        // Find the corresponding field object to get the correct keyof FileState
        const field = fieldsWithKeys.find(f => f.id === fieldId);
        if (field) {
          handleMappingChange(
            { target: { value: selectedColumn } } as React.ChangeEvent<HTMLSelectElement>,
            type,
            field.key // Pass the actual keyof FileState
          );
        } else {
           console.error(`Could not find field definition for ID: ${fieldId}`);
        }
      }
    };
  };

  // --- Modified handleStartAnalysis: Now just shows the input field --- 
  const handleInitiateSave = () => {
    if (!canProcess) {
      setAnalysisStatus('Please ensure all files are uploaded and required columns are mapped.');
      return;
    }
    setProjectName(''); // Initialize with an empty string
    setShowProjectNameInput(true);
    setAnalysisStatus('Please enter a name for the project configuration.');
  };

  // --- NEW: Function to start the backend analysis job ---
  const startAnalysisJob = async (
      spatialFileId: string,
      spatialMappingData: AnalysisMapping,
      interactionsFileId: string,
      interactionsMappingData: AnalysisMapping,
      modulesFileId: string,
      modulesMappingData: AnalysisMapping
  ): Promise<{ success: boolean; error?: string }> => {
      
      setAnalysisStatus('Starting analysis job...');
      // Ensure mappings are at least empty objects if null/undefined passed
      const safeSpatialMapping = spatialMappingData || {};
      const safeInteractionsMapping = interactionsMappingData || {};
      const safeModulesMapping = modulesMappingData || {};
      
      const payload: AnalysisPayload = {
          spatialFileId: spatialFileId,
          spatialMapping: safeSpatialMapping,
          interactionsFileId: interactionsFileId,
          interactionsMapping: safeInteractionsMapping,
          modulesFileId: modulesFileId,
          modulesMapping: safeModulesMapping
      };

      console.log("[DataInputPage] Sending analysis payload:", JSON.stringify(payload));

      try {
          const response = await fetch('http://localhost:8000/api/analysis/start', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify(payload),
          });

          console.log("[DataInputPage] Raw analysis start response status:", response.status);
          const responseText = await response.text(); 
          console.log("[DataInputPage] Raw analysis start response body:", responseText);

          if (!response.ok) {
              let errorDetail = `Server error: ${response.status} ${response.statusText}`;
              try { errorDetail = JSON.parse(responseText).detail || errorDetail; } catch {} 
              throw new Error(errorDetail);
          }

          const result = JSON.parse(responseText); 
          console.log("[DataInputPage] Parsed analysis start result:", result);

          if (!result.job_id) {
              throw new Error("Backend started analysis but did not return a job_id.");
          }

          setAnalysisStatus(`Analysis job started successfully! Job ID: ${result.job_id}`);
          setLastJobId(result.job_id);
          console.log("[DataInputPage] Navigating to results page for job:", result.job_id);
          navigate(`/analysis/${result.job_id}`);

          return { success: true }; // Analysis started successfully

      } catch (error: any) {
          const errorMessage = error instanceof Error ? error.message : String(error);
          console.error("[DataInputPage] Error during analysis start:", errorMessage);
          setAnalysisStatus(`Error starting analysis: ${errorMessage}`);
          setIsAnalyzing(false); // Allow user to try again or modify
          return { success: false, error: errorMessage };
      }
      // Note: setIsAnalyzing(false) is handled within the error case or implicitly on navigation
  };

  // --- Modified handler: Saves config and starts backend job --- 
  const handleSaveAndStart = async () => {
    if (!projectName || projectName.trim().length === 0) {
      setAnalysisStatus('Project name cannot be empty.');
      return; 
    }

    setShowProjectNameInput(false);
    setIsAnalyzing(true); // Set analyzing state here
    setAnalysisStatus('Saving project configuration...');
    setLastJobId(null);

    // Prepare Project Data for saving (using current component state)
    const currentFileIds = { 
      expression: spatialFile.fileId ?? undefined,
      metadata: interactionsFile.fileId ?? undefined,
      modules: modulesFile.fileId ?? undefined
    }; 
    
    const currentMappings: ColumnMappings = {};
    if (spatialFile.fileId) {
        currentMappings.expression = {};
        if (spatialFile.geneCol) currentMappings.expression.geneCol = spatialFile.geneCol;
        if (spatialFile.xCol) currentMappings.expression.xCol = spatialFile.xCol;
        if (spatialFile.yCol) currentMappings.expression.yCol = spatialFile.yCol;
        if (spatialFile.layerCol) currentMappings.expression.layerCol = spatialFile.layerCol;
        if (Object.keys(currentMappings.expression).length === 0) delete currentMappings.expression;
    }
    if (interactionsFile.fileId) {
        currentMappings.metadata = {};
        if (interactionsFile.ligandCol) currentMappings.metadata.ligandCol = interactionsFile.ligandCol;
        if (interactionsFile.receptorCol) currentMappings.metadata.receptorCol = interactionsFile.receptorCol;
        if (Object.keys(currentMappings.metadata).length === 0) delete currentMappings.metadata;
    }
    if (modulesFile.fileId) {
        currentMappings.modules = {};
        if (modulesFile.geneCol) currentMappings.modules.geneCol = modulesFile.geneCol;
        if (modulesFile.moduleCol) currentMappings.modules.moduleCol = modulesFile.moduleCol;
        if (Object.keys(currentMappings.modules).length === 0) delete currentMappings.modules;
    }

    const projectData = { version: "1.0", files: currentFileIds, mappings: currentMappings };

    try {
        // 1. Save Project Config
        console.log("[DataInputPage] Attempting to save project config...");
        const saveResult = await window.electronAPI.saveProject(projectName.trim(), projectData);
        if (!saveResult.success || !saveResult.filePath) {
            throw new Error(saveResult.error || 'Failed to save project configuration.');
        }
        const savedProjectInfo = { 
            filePath: saveResult.filePath, 
            name: saveResult.name!, 
            savedAt: saveResult.savedAt! 
        };
        setSavedProjects(prev => [savedProjectInfo!, ...prev.filter(p => p.filePath !== savedProjectInfo!.filePath)]);
        if (projectListStatus.includes('No saved projects')) setProjectListStatus('');
        console.log('[DataInputPage] Project configuration saved:', savedProjectInfo);
        setAnalysisStatus('Configuration saved. Starting analysis job...');

        // 2. Start Backend Analysis Job (using current state)
        const spatialMapping: AnalysisMapping = {};
        if (spatialFile.geneCol) spatialMapping.geneCol = spatialFile.geneCol;
        if (spatialFile.xCol) spatialMapping.xCol = spatialFile.xCol;
        if (spatialFile.yCol) spatialMapping.yCol = spatialFile.yCol;
        if (spatialFile.layerCol) spatialMapping.layerCol = spatialFile.layerCol;
        
        const interactionsMapping: AnalysisMapping = {};
        if (interactionsFile.ligandCol) interactionsMapping.ligandCol = interactionsFile.ligandCol;
        if (interactionsFile.receptorCol) interactionsMapping.receptorCol = interactionsFile.receptorCol;
        
        const modulesMapping: AnalysisMapping = {};
        if (modulesFile.geneCol) modulesMapping.geneCol = modulesFile.geneCol;
        if (modulesFile.moduleCol) modulesMapping.moduleCol = modulesFile.moduleCol;

        // Call the extracted function
        const analysisResult = await startAnalysisJob(
            spatialFile.fileId!, 
            spatialMapping, 
            interactionsFile.fileId!, 
            interactionsMapping, 
            modulesFile.fileId!, 
            modulesMapping
        );

        if (!analysisResult.success) {
             // Error message is already set by startAnalysisJob
             // isAnalyzing state is already set to false by startAnalysisJob
             console.error("[DataInputPage] Analysis job failed to start after saving project.");
             // No need to throw again, just let the function finish
        }
        // If analysisResult.success is true, navigation happens in startAnalysisJob

    } catch (error) {
        // This catch block handles errors from saveProject OR any synchronous errors before/between awaits
        const errorMessage = error instanceof Error ? error.message : String(error);
        console.error("[DataInputPage] Error during project save or analysis initiation:", errorMessage);
        setAnalysisStatus(`Error: ${errorMessage}`);
        setIsAnalyzing(false); // Ensure analyzing state is reset on any error in this process
    } 
  };

  // --- Modified handleLoadProjectUI --- (Added explicit return type void for clarity)
  const handleLoadProjectUI = async (projectFilePath: string, projectName: string): Promise<void> => { 
    console.log('Attempting to load project configuration and start analysis:', projectFilePath);
    setProjectListStatus(`Loading project: ${projectName}...`);
    setIsAnalyzing(true); 
    setAnalysisStatus(''); 
    setLastJobId(null);

    try {
        // 1. Load Project Config
        const result = await window.electronAPI.loadProject(projectFilePath);
        const loadResult = result as { success: boolean; error?: string; projectState?: SavedProjectData };

        if (!loadResult.success || !loadResult.projectState) {
             throw new Error(loadResult.error || 'Unknown error loading project file.');
        }
        
        const loadedData = loadResult.projectState;
        console.log("Loaded project data:", loadedData);

        // --- Check and assign file IDs --- 
        const spatialFileId = loadedData.files.expression;
        const interactionsFileId = loadedData.files.metadata;
        const modulesFileId = loadedData.files.modules;

        // Check non-null/undefined after assignment
        if (!spatialFileId || !interactionsFileId || !modulesFileId) {
             throw new Error('Loaded project is missing one or more essential file references. Cannot start analysis.');
        }
        // TypeScript should now know these are strings

        // --- 2. Update the UI state --- 
        let newSpatial = createInitialFileState();
        let newInteractions = createInitialFileState();
        let newModules = createInitialFileState();
        // Assign checked IDs
        newSpatial.fileId = spatialFileId;
        newInteractions.fileId = interactionsFileId;
        newModules.fileId = modulesFileId;
        // Restore Mappings
        newSpatial = { ...newSpatial, ...(loadedData.mappings?.expression || {}) };
        newInteractions = { ...newInteractions, ...(loadedData.mappings?.metadata || {}) };
        newModules = { ...newModules, ...(loadedData.mappings?.modules || {}) };
        
        const updateStateWithPlaceholder = (state: FileState, fileId: string): FileState => ({ 
            ...state, 
            file: { name: `[Loaded: ${fileId.substring(0, 8)}...]` } as File,
            headers: [], previewRows: [], fileInfo: undefined, error: null, isLoading: false,
        });
        // Use the already checked variables here too for consistency
        setSpatialFile(updateStateWithPlaceholder(newSpatial, spatialFileId));
        setInteractionsFile(updateStateWithPlaceholder(newInteractions, interactionsFileId));
        setModulesFile(updateStateWithPlaceholder(newModules, modulesFileId));
        
        setProjectListStatus(''); 
        setAnalysisStatus(`Project '${loadedData.projectName}' loaded. Starting analysis...`);
        setShowProjectNameInput(false); 
        console.log("[DataInputPage] State updated after load. Attempting to start analysis job automatically.");

        // --- 3. Attempt to start the analysis job --- (Use checked variables)
        const analysisResult = await startAnalysisJob(
             spatialFileId,
             loadedData.mappings?.expression || {},
             interactionsFileId,
             loadedData.mappings?.metadata || {},
             modulesFileId,
             loadedData.mappings?.modules || {}
        );

        // Check result from startAnalysisJob
        if (!analysisResult.success) {
            // Error message is already set by startAnalysisJob
            // isAnalyzing state is already set to false by startAnalysisJob
            console.error(`[DataInputPage] Failed to auto-start analysis for loaded project '${projectName}'.`);
        }
         // Success case handled by startAnalysisJob

    } catch (error) {
        // This catch handles errors from loadProject OR the file ID checks OR synchronous errors
        const msg = error instanceof Error ? error.message : String(error);
        setProjectListStatus(''); 
        setAnalysisStatus(`Error loading project: ${msg}`); 
        console.error("Load Project / Auto-Analysis Error:", error);
        setIsAnalyzing(false); 
    }
  };

  // handleDeleteSpecificProject (Ensure it updates projectListStatus)
  const handleDeleteSpecificProject = async (filePath: string, projectName: string) => {
    // Proceed directly to calling the Electron main process, which will show its own dialog
    setProjectListStatus(`Deleting project: ${projectName}...`);
    try {
        // Pass projectName to the Electron API call
        const result = await window.electronAPI.deleteProject(filePath, projectName);
        
        // Handle the result (success or cancellation/error)
        if (result.success) {
             setProjectListStatus(`Project '${projectName}' deleted.`);
             setSavedProjects(prev => prev.filter(p => p.filePath !== filePath));
             setTimeout(() => setProjectListStatus(savedProjects.length -1 === 0 ? 'No saved projects found.' : ''), 3000);
        } else {
             // Handle cancellation message specifically
             if (result.error === 'Delete cancelled by user.') {
                 setProjectListStatus('Delete operation cancelled.');
             } else {
                 // Handle other errors
                 throw new Error(result.error || 'Unknown error during deletion.');
             }
        }
    } catch (error) {
      const msg = error instanceof Error ? error.message : String(error);
      setProjectListStatus(`Error deleting project '${projectName}': ${msg}`);
      console.error("Delete Project Error:", error);
      // Set timeout to clear error message later
      setTimeout(() => setProjectListStatus(savedProjects.length > 0 ? '' : 'No saved projects found.'), 5000);
    }
  };

  // formatProjectDate (Helper for display)
  const formatProjectDate = (isoDateString: string): string => {
    try {
      return new Date(isoDateString).toLocaleString();
    } catch { 
      return isoDateString; // Fallback
    } 
  };

  return (
    <div className={styles.container}>
      <h1 className={styles.pageTitle}>Data Input</h1>

      {/* --- Project History Table --- */}
      <section className={styles.section}>
        <h2 className={styles.sectionTitle}>Project History</h2>
        {projectListStatus && <p className={styles.statusMessage}><i>{projectListStatus}</i></p>}
        {savedProjects.length > 0 ? (
          <table className={styles.projectTable}> {/* Add CSS for projectTable */} 
            <thead>
              <tr>
                <th>Project Name</th>
                <th>Date Saved</th>
                <th>Actions</th>
              </tr>
            </thead>
            <tbody>
              {savedProjects.map((proj) => (
                <tr key={proj.filePath}>
                  <td>{proj.name}</td>
                  <td>{formatProjectDate(proj.savedAt)}</td>
                  <td className={styles.projectActions}> {/* Add CSS */} 
                    <button
                      onClick={() => handleLoadProjectUI(proj.filePath, proj.name)}
                      className={styles.projectActionButton} /* Add CSS */ 
                      disabled={isAnalyzing} 
                      title="Load this project configuration"
                    >
                      Load Configuration
                    </button>
                    <button
                      onClick={() => handleDeleteSpecificProject(proj.filePath, proj.name)}
                      className={`${styles.projectActionButton} ${styles.deleteButton}`} /* Add CSS */
                      disabled={isAnalyzing} 
                      title="Delete this project configuration"
                    >
                      Delete
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        ) : (
          !projectListStatus.includes('Loading') && <p>No saved project configurations found.</p>
        )}
      </section>

      {/* Spatial Data Section */}
      <section className={styles.section}>
        <h2 className={styles.sectionTitle}>1. Spatial Data</h2>
        
        <FileUploader
          fileState={spatialFile}
          fileType="spatial"
          onFileChange={handleFileChange}
          onRemoveFile={handleRemoveFile}
          title="Spatial Data Upload"
          description="Required columns: Gene ID, X coordinate, Y coordinate, Layer."
        />
        
        {spatialFile.file && !spatialFile.error && spatialFile.headers.length > 0 && (
          <>
            {/* Column Mapper */}
            <ColumnMapper {...getColumnMapperProps(spatialFile, 'spatial')!} />
            
            {/* Data Preview */}
            <DataPreview 
              headers={spatialFile.headers}
              rows={spatialFile.previewRows}
              isLoading={spatialFile.isLoading}
              error={spatialFile.error}
              fileInfo={spatialFile.fileInfo}
            />
          </>
        )}
      </section>

      {/* Interactions Data Section */}
      <section className={styles.section}>
        <h2 className={styles.sectionTitle}>2. Interaction Data</h2>
        
        <FileUploader
          fileState={interactionsFile}
          fileType="interactions"
          onFileChange={handleFileChange}
          onRemoveFile={handleRemoveFile}
          title="Interaction Data Upload"
          description="Required columns: Ligand, Receptor (use '_' for complexes)."
        />
        
        {interactionsFile.file && !interactionsFile.error && interactionsFile.headers.length > 0 && (
          <>
            {/* Column Mapper */}
            <ColumnMapper {...getColumnMapperProps(interactionsFile, 'interactions')!} />
            
            {/* Data Preview */}
            <DataPreview 
              headers={interactionsFile.headers}
              rows={interactionsFile.previewRows}
              isLoading={interactionsFile.isLoading}
              error={interactionsFile.error}
              fileInfo={interactionsFile.fileInfo}
            />
          </>
        )}
      </section>

      {/* Modules Data Section */}
      <section className={styles.section}>
        <h2 className={styles.sectionTitle}>3. Module Data</h2>
        
        <FileUploader
          fileState={modulesFile}
          fileType="modules"
          onFileChange={handleFileChange}
          onRemoveFile={handleRemoveFile}
          title="Module Data Upload"
          description="Required columns: Gene ID, Module ID."
        />
        
        {modulesFile.file && !modulesFile.error && modulesFile.headers.length > 0 && (
          <>
            {/* Column Mapper */}
            <ColumnMapper {...getColumnMapperProps(modulesFile, 'modules')!} />
            
            {/* Data Preview */}
            <DataPreview 
              headers={modulesFile.headers}
              rows={modulesFile.previewRows}
              isLoading={modulesFile.isLoading}
              error={modulesFile.error}
              fileInfo={modulesFile.fileInfo}
            />
          </>
        )}
      </section>

      {/* Analysis Trigger Section */}
      <section className={`${styles.section} ${styles.analysisSection}`}>
        <h2 className={styles.sectionTitle}>4. Run Analysis</h2>
        
        {/* Conditionally show input or button */} 
        {!showProjectNameInput ? (
          <>
            <p>Upload or modify data above, then start a new analysis run.</p>
            <button
              className={styles.startButton}
              onClick={handleInitiateSave}
              disabled={!canProcess || isAnalyzing}
            >
              {isAnalyzing ? 'Processing...' : 'Save & Start Analysis'}
            </button>
          </>
        ) : (
          <div className={styles.projectNameInputArea}>
            <label htmlFor="projectName">Project Name:</label>
            <input 
              type="text"
              id="projectName"
              value={projectName}
              onChange={(e) => setProjectName(e.target.value)}
              className={styles.projectNameInput}
              placeholder="Enter a name for this configuration"
            />
            <button 
              onClick={handleSaveAndStart}
              disabled={isAnalyzing || !projectName.trim()}
              className={styles.confirmButton}
            >
              {isAnalyzing ? 'Saving...' : 'Confirm & Start'}
            </button>
            <button 
              onClick={() => { setShowProjectNameInput(false); setAnalysisStatus(''); }}
              disabled={isAnalyzing}
              className={styles.cancelButton}
              title="Cancel saving configuration"
            >
               Cancel
            </button>
          </div>
        )}
        
        {analysisStatus && (
          <p 
            className={styles.statusMessage} 
            style={{ color: analysisStatus.includes('failed') ? 'red' : analysisStatus.includes('successfully') ? 'green' : 'inherit' }}
          >
            {analysisStatus}
          </p>
        )}
        {lastJobId && analysisStatus.includes('successfully') && (
            <p className={styles.statusMessage}> 
               <a href="#" onClick={(e) => { e.preventDefault(); alert(`Navigate to results for job: ${lastJobId}`)} }>
                   View Results (Job: {lastJobId})
               </a>
            </p>
        )}
      </section>

    </div>
  );
};

export default DataInputPage; 