import React from 'react';
import styles from './LoadingPage.module.css';

interface LoadingPageProps {
  message?: string;
  jobId?: string;
  progress?: number;
}

const LoadingPage: React.FC<LoadingPageProps> = ({ 
  message = 'Loading...', 
  jobId,
  progress 
}) => {
  return (
    <div className={styles.loadingContainer}>
      <div className={styles.loadingContent}>
        <div className={styles.spinner}>
          <div className={styles.dot}></div>
          <div className={styles.dot}></div>
          <div className={styles.dot}></div>
        </div>
        
        <h2 className={styles.loadingMessage}>{message}</h2>
        
        {progress !== undefined && (
          <div className={styles.simpleProgressContainer}>
            <div 
              className={styles.simpleProgressBar} 
              style={{ width: `${progress * 100}%` }}
            ></div>
          </div>
        )}
        
        {jobId && (
          <p className={styles.jobIdText}>
            Job ID: <span className={styles.jobId}>{jobId}</span>
          </p>
        )}
      </div>
    </div>
  );
};

export default LoadingPage; 