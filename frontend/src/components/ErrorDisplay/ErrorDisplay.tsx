import React from 'react';
import './ErrorDisplay.css';

interface ErrorDisplayProps {
  message: string;
  details?: string;
  onRetry?: () => void;
  retryLabel?: string;
}

const ErrorDisplay: React.FC<ErrorDisplayProps> = ({ 
  message, 
  details, 
  onRetry, 
  retryLabel = "Try Again" 
}) => {
  return (
    <div className="error-container">
      <div className="error-icon">❌</div>
      <h2 className="error-title">Error</h2>
      <p className="error-message">{message}</p>
      {details && <p className="error-details">{details}</p>}
      {onRetry && (
        <button className="error-retry-button" onClick={onRetry}>
          {retryLabel}
        </button>
      )}
    </div>
  );
};

export default ErrorDisplay; 