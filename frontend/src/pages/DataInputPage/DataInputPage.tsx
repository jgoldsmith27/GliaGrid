import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import styles from './DataInputPage.module.css';
import FileUploader from '../../components/FileUploader/FileUploader';
import ColumnMapper from '../../components/ColumnMapper/ColumnMapper';
import DataPreview from '../../components/DataPreview/DataPreview';
import { FileState, FileType, FilePreviewResult, requiredColumns, mappingFields, AnalysisPayload } from '../../types/dataTypes';

// Augment the window interface to tell TypeScript about electronAPI
declare global {
  interface Window {
    electronAPI: {
      readFilePreview: (fileContent: string, fileName: string) => Promise<FilePreviewResult>;
    }
  }
}

const DataInputPage: React.FC = () => {
  const navigate = useNavigate();

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
  const [analysisStatus, setAnalysisStatus] = useState<string>(''); // To display analysis status/errors
  const [isAnalyzing, setIsAnalyzing] = useState<boolean>(false); // Loading state for analysis
  const [lastJobId, setLastJobId] = useState<string | null>(null); // Store the ID of the last started job

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

  // --- Start Analysis Handler ---
  const handleStartAnalysis = async () => {
    if (!canProcess) {
      setAnalysisStatus('Please ensure all files are uploaded and required columns are mapped.');
      return;
    }

    setIsAnalyzing(true);
    setAnalysisStatus('Starting analysis...');
    setLastJobId(null); // Reset last job ID

    // Construct the payload
    const payload: AnalysisPayload = {
      spatialFileId: spatialFile.fileId!,
      spatialMapping: {
        geneCol: spatialFile.geneCol!,
        xCol: spatialFile.xCol!,
        yCol: spatialFile.yCol!,
        layerCol: spatialFile.layerCol!
      },
      interactionsFileId: interactionsFile.fileId!,
      interactionsMapping: {
        ligandCol: interactionsFile.ligandCol!,
        receptorCol: interactionsFile.receptorCol!
      },
      modulesFileId: modulesFile.fileId!,
      modulesMapping: {
        geneCol: modulesFile.geneCol!,
        moduleCol: modulesFile.moduleCol!
      }
    };

    console.log("[React] Sending analysis payload:", payload);

    try {
      const response = await fetch('http://localhost:8000/api/analysis/start', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(payload),
      });

      const result = await response.json();

      if (!response.ok) {
        throw new Error(result.detail || `Server error: ${response.status}`);
      }

      console.log("[React] Analysis response received:", result);
      setAnalysisStatus(`Analysis job started successfully! Job ID: ${result.job_id}`);
      setLastJobId(result.job_id);
      
      // Navigate to the status page with the job ID
      if (result.job_id) {
        navigate(`/analysis/${result.job_id}`); 
      }

    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : String(error);
      console.error("[React] Error starting analysis:", errorMessage);
      setAnalysisStatus(`Error starting analysis: ${errorMessage}`);
    } finally {
      setIsAnalyzing(false);
    }
  };
  // --- End Analysis Handler ---

  return (
    <div className={styles.container}>
      <h1 className={styles.pageTitle}>Data Input</h1>

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
        <button 
          className={styles.startButton} 
          onClick={handleStartAnalysis} 
          disabled={!canProcess || isAnalyzing} // Disable if not ready or already analyzing
        >
          {isAnalyzing ? 'Processing...' : 'Start Analysis'}
        </button>
        {analysisStatus && (
          <p 
            className={styles.statusMessage} 
            style={{ color: analysisStatus.includes('failed') ? 'red' : analysisStatus.includes('successfully') ? 'green' : 'inherit' }}
          >
            {analysisStatus}
          </p>
        )}
        {/* Simple link placeholder - replace with proper routing/navigation later */}
        {lastJobId && analysisStatus.includes('successfully') && (
            <p className={styles.statusMessage}> 
               {/* In a real app, this would be a Link component from react-router-dom */}
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