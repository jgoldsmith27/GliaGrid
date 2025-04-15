import React, { useState } from 'react';
import styles from './DataInputPage.module.css';
import FileUploader from '../../components/FileUploader/FileUploader';
import ColumnMapper from '../../components/ColumnMapper/ColumnMapper';
import DataPreview from '../../components/DataPreview/DataPreview';
import { FileState, FileType, FilePreviewResult, requiredColumns, mappingFields } from '../../types/dataTypes';

// Augment the window interface to tell TypeScript about electronAPI
declare global {
  interface Window {
    electronAPI: {
      readFilePreview: (fileContent: string, fileName: string) => Promise<FilePreviewResult>;
    }
  }
}

const DataInputPage: React.FC = () => {
  const [spatialFile, setSpatialFile] = useState<FileState>({ 
    file: null, 
    headers: [], 
    previewRows: [], 
    error: null, 
    isLoading: false 
  });
  
  const [interactionsFile, setInteractionsFile] = useState<FileState>({ 
    file: null, 
    headers: [], 
    previewRows: [], 
    error: null, 
    isLoading: false 
  });
  
  const [modulesFile, setModulesFile] = useState<FileState>({ 
    file: null, 
    headers: [], 
    previewRows: [], 
    error: null, 
    isLoading: false 
  });

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

    updateFileState(setter, { file: null, headers: [], previewRows: [], error: null, isLoading: false }); // Reset on new selection

    if (event.target.files && event.target.files[0]) {
      const file = event.target.files[0];
      const fileName = file.name;

      if (!fileName.endsWith('.csv') && !fileName.endsWith('.h5ad')) {
        alert('Please select a .csv or .h5ad file.');
        event.target.value = '';
        return;
      }

      updateFileState(setter, { file, isLoading: true });
      console.log(`Reading file content slice for ${type}:`, fileName);

      // --- Read file content SLICE using FileReader --- 
      const reader = new FileReader();
      const CHUNK_SIZE = 1024 * 512; // Read first 512KB - should be plenty for headers + few rows
      const blobSlice = file.slice(0, CHUNK_SIZE);

      reader.onload = async (e) => {
        const fileContentChunk = e.target?.result as string;
        if (!fileContentChunk) {
            updateFileState(setter, { error: 'Could not read file content chunk.', isLoading: false });
            return;
        }

        try {
            // Pass the chunk and the original filename
            const result: FilePreviewResult = await window.electronAPI.readFilePreview(fileContentChunk, fileName);
            console.log(`[React] Preview result received for ${type}:`, result);

            if (result.error) {
              updateFileState(setter, { error: result.error, isLoading: false });
            } else if (result.headers && result.previewRows) {
              updateFileState(setter, {
                 headers: result.headers,
                 previewRows: result.previewRows,
                 error: null,
                 isLoading: false,
                 // Reset mappings
                 geneCol: undefined, xCol: undefined, yCol: undefined, layerCol: undefined,
                 ligandCol: undefined, receptorCol: undefined, moduleCol: undefined
              });
            } else {
              updateFileState(setter, { error: 'Invalid response from preview reader.', isLoading: false });
            }
        } catch (ipcError) {
            const errorMessage = ipcError instanceof Error ? ipcError.message : String(ipcError);
            updateFileState(setter, { error: `Error communicating with main process: ${errorMessage}`, isLoading: false });
        }
      };

      reader.onerror = (e) => {
        updateFileState(setter, { error: 'Failed to read file chunk.', isLoading: false });
      };

      // Read the BLOB SLICE as text
      reader.readAsText(blobSlice);

      // H5AD still needs separate handling (cannot read slice as text)
      if (fileName.endsWith('.h5ad')) {
          reader.abort(); // Stop the text reader if it was started
          try {
            const result = await window.electronAPI.readFilePreview("", fileName);
            if (result.error) {
                updateFileState(setter, { error: result.error, isLoading: false });
            } else {
                updateFileState(setter, { error: 'Unexpected success from H5AD preview?', isLoading: false });
            }
          } catch (ipcError) {
              const errorMessage = ipcError instanceof Error ? ipcError.message : String(ipcError);
              updateFileState(setter, { error: `Error communicating with main process: ${errorMessage}`, isLoading: false });
          }
          event.target.value = '';
          return; // Don't proceed further for H5AD for now
      }

      // Clear the input value so the same file can be selected again if needed
      event.target.value = '';
    }
  };

  const handleMappingChange = (
    event: React.ChangeEvent<HTMLSelectElement>,
    type: FileType,
    columnKey: keyof FileState
  ) => {
    const setter = type === 'spatial' ? setSpatialFile :
                   type === 'interactions' ? setInteractionsFile :
                   setModulesFile;
    const selectedValue = event.target.value;
    updateFileState(setter, { [columnKey]: selectedValue });
  };

  const handleRemoveFile = (type: FileType) => {
    const setter = type === 'spatial' ? setSpatialFile :
                   type === 'interactions' ? setInteractionsFile :
                   setModulesFile;
    // Reset the state for this file type
    updateFileState(setter, {
        file: null,
        headers: [],
        previewRows: [],
        error: null,
        isLoading: false,
        geneCol: undefined, xCol: undefined, yCol: undefined, layerCol: undefined,
        ligandCol: undefined, receptorCol: undefined, moduleCol: undefined
    });
  };

  // Check if all required columns for a file type are mapped
  const checkRequiredMapping = (state: FileState, type: FileType): boolean => {
    if (!state.file) return false; // No file loaded
    const reqCols = requiredColumns[type] || [];
    return reqCols.every(key => !!state[key as keyof FileState]);
  };

  const canProcess = checkRequiredMapping(spatialFile, 'spatial') &&
                     checkRequiredMapping(interactionsFile, 'interactions');
                     // Module file is optional, but if present, must be mapped
                     // (Add checkRequiredMapping(modulesFile, 'modules') if modulesFile.file)

  // Convert FileState and mapping config to ColumnMapper props format
  const getColumnMapperProps = (fileState: FileState, type: FileType) => {
    if (!fileState.file) return null;
    
    // Create mapping fields in the format required by ColumnMapper
    const fields = mappingFields[type].map(({ key, label }) => {
      return {
        id: key,
        label,
        required: requiredColumns[type].includes(key as string)
      };
    });
    
    // Get current mappings
    const currentMappings: Record<string, string> = {};
    mappingFields[type].forEach(({ key }) => {
      const value = fileState[key as keyof FileState];
      if (value) currentMappings[key] = value as string;
    });
    
    return {
      headers: fileState.headers,
      isLoading: fileState.isLoading,
      error: fileState.error,
      requiredFields: fields,
      mappings: currentMappings,
      onMappingChange: (fieldId: string, selectedColumn: string) => {
        handleMappingChange(
          { target: { value: selectedColumn } } as React.ChangeEvent<HTMLSelectElement>,
          type,
          fieldId as keyof FileState
        );
      }
    };
  };

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
          description="Required columns: Gene ID, X coordinate, Y coordinate. Optional: Layer."
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
            />
          </>
        )}
      </section>

      {/* Modules Data Section */}
      <section className={styles.section}>
        <h2 className={styles.sectionTitle}>3. Module Data (Optional)</h2>
        
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
            />
          </>
        )}
      </section>

      {/* Process Button */}
      <button
         className={styles.processButton}
         disabled={!canProcess}
         onClick={() => console.log('Process files:', { spatialFile, interactionsFile, modulesFile })} // Placeholder action
         title={canProcess ? "Proceed to analysis" : "Please select files and map required columns"}
      >
        Start Analysis
      </button>
    </div>
  );
};

export default DataInputPage; 