import React from 'react';
import styles from './InteractionTable.module.css';

interface InteractionTableProps {
  data: Array<{ [key: string]: any }> | null; // Expects an array of interaction objects
}

// Helper to format keys/headers (e.g., ligand_norm_expr -> Ligand Norm Expr)
const formatHeader = (key: string): string => {
    // Handle specific cases for clarity
    if (key === 'ligand_norm_expr') return 'Ligand Norm Expr';
    if (key === 'receptor_avg_norm_expr') return 'Receptor Avg Norm Expr';
    if (key === 'score') return 'Interaction Score';
    if (key === 'receptor_modules') return 'Receptor Module(s)';
    if (key === 'is_same_module') return 'Same Module?';
    if (key === 'interaction_type') return 'Module Context';
    
    return key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
};

// Helper to format cell values
const formatCellValue = (value: any): string => {
    if (typeof value === 'number') {
        // Format numbers to a reasonable precision
        return value.toFixed(4);
    }
    if (typeof value === 'boolean') {
        return value ? 'Yes' : 'No';
    }
    if (value === null || value === undefined) {
        return 'N/A';
    }
    if (Array.isArray(value)) {
        return value.join(', ') || 'N/A'; // Join array elements like receptor modules
    }
    return String(value);
};

const InteractionTable: React.FC<InteractionTableProps> = ({ data }) => {
  if (!data || data.length === 0) {
    return <p className={styles.noData}>No interaction data available.</p>;
  }

  // Determine columns dynamically from the keys of the first object
  // Exclude keys that might be less relevant for direct display or handle them specially
  const columns = Object.keys(data[0]).filter(key => 
      ![ 'num_receptor_components_in_modules' ].includes(key)
  );
  
  // Reorder columns for better readability (put score near the end)
  const preferredOrder = ['ligand', 'receptor', 'pathway', 'ligand_module', 'receptor_modules', 'is_same_module', 'interaction_type', 'ligand_norm_expr', 'receptor_avg_norm_expr', 'score'];
  const sortedColumns = columns.sort((a, b) => {
      const indexA = preferredOrder.indexOf(a);
      const indexB = preferredOrder.indexOf(b);
      if (indexA === -1 && indexB === -1) return 0; // Keep original order if both not preferred
      if (indexA === -1) return 1; // Put non-preferred after preferred
      if (indexB === -1) return -1; // Put preferred before non-preferred
      return indexA - indexB; // Sort based on preferred order
  });

  return (
    <div className={styles.tableContainer}>
      <table className={styles.table}>
        <thead>
          <tr>
            {sortedColumns.map(col => (
              <th key={col} className={styles.th}>{formatHeader(col)}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {data.map((row, rowIndex) => (
            <tr key={rowIndex} className={styles.tr}>
              {sortedColumns.map(col => (
                <td key={col} className={styles.td}>{formatCellValue(row[col])}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default InteractionTable; 