import React from 'react';
import styles from './DataPreview.module.css';
import { FileState } from '../../types/dataTypes';

interface DataPreviewProps {
  headers: string[];
  rows: Record<string, any>[];
  isLoading: boolean;
  error: string | null;
}

const DataPreview: React.FC<DataPreviewProps> = ({ 
  headers, 
  rows, 
  isLoading, 
  error 
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