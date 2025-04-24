import React from 'react';
import styles from './SummaryStats.module.css';
import { SummaryStatsData } from '../../types/analysisResults';

interface SummaryStatsProps {
  data: SummaryStatsData | null;
}

const SummaryStats: React.FC<SummaryStatsProps> = ({ data }) => {
  if (!data) {
    return <p className={styles.noData}>No summary statistics available.</p>;
  }

  // Format keys for display (e.g., unique_ligands -> Unique Ligands)
  const formatKey = (key: string): string => {
    return key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
  };

  return (
    <div className={styles.summaryStatsContainer}>
      {Object.entries(data).map(([key, value]) => (
        <div key={key} className={styles.statItem}>
          <span className={styles.statLabel}>{formatKey(key)}:</span>
          <span className={styles.statValue}>{value}</span>
        </div>
      ))}
    </div>
  );
};

export default SummaryStats; 