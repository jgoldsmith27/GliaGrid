import React, { useState } from 'react';
// Removed unused CSS module import
// import styles from './ModuleContextTable.module.css'; 
import { ModuleContextResult } from '../../types/analysisResults';
// Import the generic table and its types
import AnalysisTable, { ColumnDefinition } from '../AnalysisTable/AnalysisTable';

// Update props interface to match AnalysisTable's expected handler signature
interface ModuleContextTableProps {
  data: ModuleContextResult[] | null;
  // Expect the handler from parent to take row and index
  onRowClick?: (row: ModuleContextResult, index: number) => void; 
  loading?: boolean;
  selectedRowIndex?: number | null; // Accept selected index from parent
}

// Define columns specifically for ModuleContextResult
const moduleContextColumns: ColumnDefinition<ModuleContextResult>[] = [
  { key: 'ligand', header: 'Ligand' },
  { key: 'receptor', header: 'Receptor' },
  { key: 'ligand_module', header: 'Ligand Module' },
  { key: 'receptor_modules', header: 'Receptor Module(s)' }, // Array formatting handled by AnalysisTable
  { key: 'is_same_module', header: 'Same Module?' }, // Boolean formatting handled by AnalysisTable
  { key: 'interaction_type', header: 'Module Context' },
];

const ModuleContextTable: React.FC<ModuleContextTableProps> = ({
  data,
  onRowClick,
  loading = false,
  selectedRowIndex, // Receive selected index from props
}) => {

  // No internal state or handler needed anymore

  // Let AnalysisTable handle the no data case directly
  // if (!data || data.length === 0) { ... }

  return (
    // Removed the outer container div, assuming layout is handled by DisplayPanel
    <AnalysisTable<ModuleContextResult> // Specify the type for generics
      data={data}
      columns={moduleContextColumns}
      loading={loading}
      onRowClick={onRowClick} // Pass parent handler directly
      selectedRowIndex={selectedRowIndex} // Pass parent selected index directly
    />
  );
};

export default ModuleContextTable; 