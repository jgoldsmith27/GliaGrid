import React from 'react';
import styles from './AnalysisTable.module.css';

// Generic type for column definitions
export interface ColumnDefinition<T> {
  key: keyof T;
  header: string;
  format?: (value: T[keyof T]) => string; // Optional formatter
}

// Generic props interface
interface AnalysisTableProps<T> {
  data: T[] | null;
  columns: ColumnDefinition<T>[];
  loading?: boolean;
  onRowClick?: (row: T, index: number) => void;
  selectedRowIndex?: number | null;
}

// Helper to format cell values - uses provided formatter or default
const formatCellValue = <T, K extends keyof T>(
    value: T[K], 
    formatter?: (value: T[K]) => string
): string => {
    if (formatter) {
        return formatter(value);
    }
    // Default formatting (can be expanded)
    if (typeof value === 'number') {
        return value.toFixed(4); 
    }
    if (typeof value === 'boolean') {
        return value ? 'Yes' : 'No';
    }
     if (Array.isArray(value)) {
        return value.join(', ') || 'N/A'; 
    }
    if (value === null || value === undefined) {
        return 'N/A';
    }
    return String(value);
};


// Generic table component using React.FC with generics
function AnalysisTable<T extends Record<string, any>>({
  data,
  columns,
  loading = false,
  onRowClick,
  selectedRowIndex,
}: AnalysisTableProps<T>): React.ReactElement | null { // Return type specified

  if (!data || data.length === 0) {
    // Consider passing a 'noDataMessage' prop for customization
    return <p className={styles.noData}>No data available.</p>;
  }

  return (
    // Use the tableWrapper for positioning the overlay
    <div className={styles.tableWrapper}> 
      {/* Conditional loading overlay */}
      {loading && (
        <div className={styles.loadingOverlay}>
          <div className={styles.spinner}></div>
        </div>
      )}
      {/* Actual table */}
      <table className={styles.table}>
        <thead>
          <tr>
            {columns.map((col) => (
              <th key={String(col.key)} className={styles.th}>
                {col.header}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {data.map((row, rowIndex) => (
            <tr
              key={rowIndex}
              className={`${styles.tr} ${
                selectedRowIndex === rowIndex ? styles.selectedRow : ''
              }`}
              onClick={() => !loading && onRowClick && onRowClick(row, rowIndex)}
              style={{ cursor: loading || !onRowClick ? 'default' : 'pointer' }}
            >
              {columns.map((col) => (
                <td key={String(col.key)} className={styles.td}>
                  {formatCellValue(row[col.key], col.format)}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

export default AnalysisTable; 