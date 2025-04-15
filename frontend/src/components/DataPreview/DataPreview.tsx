import React from 'react';
import styles from './DataPreview.module.css';
import { FileState } from '../../types/dataTypes';

interface DataPreviewProps {
  headers: string[];
  rows: Record<string, any>[];
  isLoading: boolean;
  error: string | null;
  fileInfo?: {
    shape: string;
    obs_keys: string[];
    var_keys: string[];
    obsm_keys: string[];
  };
}

const DataPreview: React.FC<DataPreviewProps> = ({ 
  headers, 
  rows, 
  isLoading, 
  error,
  fileInfo
}) => {
  if (headers.length === 0 || rows.length === 0) {
    return null;
  }

  if (isLoading) {
    return <p className={styles.loadingText}>Loading preview data...</p>;
  }

  if (error) {
    return <p className={styles.errorText}>Error: {error}</p>;
  }

  return (
    <div className={styles.container}>
      <h4 className={styles.title}>File Preview (First 5 Rows)</h4>
      
      {/* Display H5AD file info if available */}
      {fileInfo && (
        <div className={styles.fileInfo}>
          <h5>H5AD File Information</h5>
          <p><strong>File size:</strong> {fileInfo.shape}</p>
          
          {fileInfo.obs_keys.length > 0 && (
            <div>
              <p><strong>Observation keys:</strong></p>
              <div className={styles.keyList}>
                {fileInfo.obs_keys.map(key => <span key={key} className={styles.keyItem}>{key}</span>)}
              </div>
            </div>
          )}
          
          {fileInfo.var_keys.length > 0 && (
            <div>
              <p><strong>Variable keys:</strong></p>
              <div className={styles.keyList}>
                {fileInfo.var_keys.map(key => <span key={key} className={styles.keyItem}>{key}</span>)}
              </div>
            </div>
          )}
          
          {fileInfo.obsm_keys.length > 0 && (
            <div>
              <p><strong>Observation matrices:</strong></p>
              <div className={styles.keyList}>
                {fileInfo.obsm_keys.map(key => <span key={key} className={styles.keyItem}>{key}</span>)}
              </div>
            </div>
          )}
        </div>
      )}
      
      <div className={styles.tableContainer}>
        <table className={styles.table}>
          <thead>
            <tr>
              {headers.map(header => (
                <th key={header} className={styles.headerCell}>{header}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((row, rowIndex) => (
              <tr key={rowIndex} className={styles.tableRow}>
                {headers.map(header => (
                  <td key={`${rowIndex}-${header}`} className={styles.tableCell}>
                    {row[header] !== undefined && row[header] !== null ? String(row[header]) : ''}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export default DataPreview; 