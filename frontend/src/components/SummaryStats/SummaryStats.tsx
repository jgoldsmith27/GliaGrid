import React from 'react';
import styles from './SummaryStats.module.css';

interface SummaryStatsProps {
  data: { [key: string]: number } | null; // Expects e.g., { unique_ligands: 100, unique_receptors: 120 }
}

const SummaryStats: React.FC<SummaryStatsProps> = ({ data }) => {
  if (!data || Object.keys(data).length === 0) {
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