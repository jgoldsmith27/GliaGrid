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
      console.log(`Processing file for ${type}:`, fileName);

      try {
        // Create a FormData object to send the file to the backend
        const formData = new FormData();
        formData.append('file', file);

        // Call the backend API to process the file
        const response = await fetch('http://localhost:8000/api/file/preview', {
          method: 'POST',
          body: formData,
        });

        if (!response.ok) {
          const errorText = await response.text();
          throw new Error(`Server error: ${response.status} - ${errorText}`);
        }

        const result = await response.json();
        console.log(`[React] Preview result received for ${type}:`, result);

        if (result.error) {
          updateFileState(setter, { error: result.error, isLoading: false });
        } else if (result.headers && result.previewRows) {
          updateFileState(setter, {
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
                     checkRequiredMapping(interactionsFile, 'interactions') &&
                     checkRequiredMapping(modulesFile, 'modules');

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