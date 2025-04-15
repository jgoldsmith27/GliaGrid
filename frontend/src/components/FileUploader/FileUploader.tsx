import React from 'react';
import styles from './FileUploader.module.css';
import { FileState, FileType } from '../../types/dataTypes';

interface FileUploaderProps {
  fileState: FileState;
  fileType: FileType;
  onFileChange: (event: React.ChangeEvent<HTMLInputElement>, type: FileType) => void;
  onRemoveFile: (type: FileType) => void;
  title: string;
  description: string;
}

const FileUploader: React.FC<FileUploaderProps> = ({
  fileState,
  fileType,
  onFileChange,
  onRemoveFile,
  title,
  description
}) => {
  const getFileName = (state: FileState) => {
    return state.file?.name ?? 'No file selected';
  };

  return (
    <div className={styles.container}>
      <h3 className={styles.title}>{title}</h3>
      <p className={styles.description}>{description}</p>
      
      <div className={styles.inputContainer}>
        <label className={styles.inputLabel}>
          Choose {fileType.charAt(0).toUpperCase() + fileType.slice(1)} File (.csv, .h5ad)
          <input 
            type="file" 
            accept=".csv,.h5ad"
            onChange={(e) => onFileChange(e, fileType)}
            className={styles.input}
          />
        </label>
        <span className={styles.fileName}>{getFileName(fileState)}</span>
        {fileState.file && (
          <button 
            onClick={() => onRemoveFile(fileType)} 
            className={styles.removeButton}
            aria-label={`Remove ${fileType} file`}
          >
            ×
          </button>
        )}
      </div>
      
      {fileState.isLoading && <p className={styles.loadingText}>Loading file...</p>}
      {fileState.error && <p className={styles.errorText}>Error: {fileState.error}</p>}
    </div>
  );
};

export default FileUploader; 