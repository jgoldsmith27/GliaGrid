import React from 'react';
import styles from './ModuleContextTable.module.css'; // Assume CSS module exists or create it
import { ModuleContextResult } from '../../types/analysisResults';

interface ModuleContextTableProps {
  data: ModuleContextResult[] | null;
}

// Helper to format header
const formatHeader = (key: keyof ModuleContextResult): string => {
    switch (key) {
        case 'ligand': return 'Ligand';
        case 'receptor': return 'Receptor';
        case 'interaction_type': return 'Module Context';
        case 'ligand_module': return 'Ligand Module';
        case 'receptor_modules': return 'Receptor Module(s)';
        case 'is_same_module': return 'Same Module?';
        // Add cases for other potential fields
        default:
            const strKey = String(key);
            return strKey.replace(/_/g, ' ').replace(/\b\w/g, (l: string) => l.toUpperCase());
    }
};

// Helper to format cell values
const formatCellValue = (value: any): string => {
    if (typeof value === 'boolean') {
        return value ? 'Yes' : 'No'; // Format boolean
    }
    if (Array.isArray(value)) {
        return value.join(', ') || 'N/A'; // Format array
    }
    if (value === null || value === undefined) {
        return 'N/A';
    }
    return String(value);
};

const ModuleContextTable: React.FC<ModuleContextTableProps> = ({ data }) => {
  if (!data || data.length === 0) {
    return <p className={styles.noData}>No Module Context data available.</p>;
  }

  // Explicitly define columns and their order
  const columns: (keyof ModuleContextResult)[] = [
      'ligand',
      'receptor',
      'ligand_module',
      'receptor_modules',
      'is_same_module',
      'interaction_type'
  ];

  return (
    <div className={styles.tableContainer}>
      <table className={styles.table}>
        <thead>
          <tr>
            {columns.map(col => (
              <th key={col} className={styles.th}>{formatHeader(col)}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {data.map((row, rowIndex) => (
            <tr key={rowIndex} className={styles.tr}>
              {columns.map(col => (
                <td key={col} className={styles.td}>{formatCellValue(row[col])}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default ModuleContextTable; 