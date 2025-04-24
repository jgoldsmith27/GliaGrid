import React from 'react';
import styles from './PathwayDominanceTable.module.css'; // Assume CSS module exists or create it
import { PathwayDominanceResult } from '../../types/analysisResults';

interface PathwayDominanceTableProps {
  data: PathwayDominanceResult[] | null;
}

// Helper to format header
const formatHeader = (key: keyof PathwayDominanceResult): string => {
    switch (key) {
        case 'ligand': return 'Ligand';
        case 'receptor': return 'Receptor';
        case 'score': return 'Interaction Score';
        case 'pathway': return 'Pathway';
        case 'ligand_norm_expr': return 'Ligand Norm Expr';
        case 'receptor_avg_norm_expr': return 'Receptor Avg Norm Expr';
        default:
            const strKey = String(key);
            return strKey.replace(/_/g, ' ').replace(/\b\w/g, (l: string) => l.toUpperCase());
    }
};

// Helper to format cell values
const formatCellValue = (value: any): string => {
    if (typeof value === 'number') {
        return value.toFixed(4); // Format score
    }
    if (value === null || value === undefined) {
        return 'N/A';
    }
    return String(value);
};

const PathwayDominanceTable: React.FC<PathwayDominanceTableProps> = ({ data }) => {
  if (!data || data.length === 0) {
    return <p className={styles.noData}>No Pathway Dominance data available.</p>;
  }

  // Explicitly define columns and their order for this specific table
  const columns: (keyof PathwayDominanceResult)[] = [
      'ligand',
      'receptor',
      'pathway',
      'ligand_norm_expr',
      'receptor_avg_norm_expr',
      'score' 
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

export default PathwayDominanceTable; 