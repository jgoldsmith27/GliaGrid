import React, { useState } from 'react';
import { PathwayDominanceResult } from '../../types/analysisResults';
import AnalysisTable, { ColumnDefinition } from '../AnalysisTable/AnalysisTable';

interface PathwayDominanceTableProps {
  data: PathwayDominanceResult[] | null;
  onRowClick?: (row: PathwayDominanceResult, index: number) => void;
  loading?: boolean;
  selectedRowIndex?: number | null;
}

// Define columns specifically for PathwayDominanceResult
const pathwayDominanceColumns: ColumnDefinition<PathwayDominanceResult>[] = [
  { key: 'ligand', header: 'Ligand' },
  { key: 'receptor', header: 'Receptor' },
  { key: 'pathway', header: 'Pathway' },
  { key: 'ligand_norm_expr', header: 'Ligand Norm Expr' }, // Formatting handled by AnalysisTable default
  { key: 'receptor_avg_norm_expr', header: 'Receptor Avg Norm Expr' }, // Formatting handled by AnalysisTable default
  { key: 'score', header: 'Interaction Score' }, // Formatting handled by AnalysisTable default
];

const PathwayDominanceTable: React.FC<PathwayDominanceTableProps> = ({
  data,
  onRowClick,
  loading = false,
  selectedRowIndex,
}) => {
  return (
    <AnalysisTable<PathwayDominanceResult>
      data={data}
      columns={pathwayDominanceColumns}
      loading={loading}
      onRowClick={onRowClick}
      selectedRowIndex={selectedRowIndex}
    />
  );
};

export default PathwayDominanceTable; 