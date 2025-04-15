import React from 'react';
import styles from './ColumnMapper.module.css';

export interface ColumnMapperProps {
  headers: string[];
  isLoading: boolean;
  error: string | null;
  requiredFields: { id: string, label: string, required: boolean }[];
  mappings: Record<string, string>;
  onMappingChange: (fieldId: string, selectedColumn: string) => void;
}

const ColumnMapper: React.FC<ColumnMapperProps> = ({
  headers,
  isLoading,
  error,
  requiredFields,
  mappings,
  onMappingChange,
}) => {
  if (isLoading) {
    return <div className={styles.loadingText}>Loading column headers...</div>;
  }

  if (error) {
    return <div className={styles.errorText}>Error: {error}</div>;
  }

  if (!headers || headers.length === 0) {
    return <div className={styles.errorText}>No headers found in the uploaded file.</div>;
  }

  return (
    <div className={styles.container}>
      <h3 className={styles.title}>Column Mapping</h3>
      <p className={styles.description}>
        Select the column from your file that corresponds to each required field:
      </p>

      <div className={styles.grid}>
        {requiredFields.map((field) => (
          <div key={field.id} className={styles.item}>
            <label className={styles.label} htmlFor={`field-${field.id}`}>
              {field.label}
              {field.required && <span className={styles.requiredIndicator}>*</span>}
            </label>
            <select
              id={`field-${field.id}`}
              className={styles.select}
              value={mappings[field.id] || ''}
              onChange={(e) => onMappingChange(field.id, e.target.value)}
            >
              <option value="">-- Select a column --</option>
              {headers.map((header) => (
                <option key={header} value={header}>
                  {header}
                </option>
              ))}
            </select>
          </div>
        ))}
      </div>
    </div>
  );
};

export default ColumnMapper; 