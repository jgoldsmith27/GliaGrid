import React, { useState, useCallback } from 'react';
import axios from 'axios';
import styles from './Visualization.module.css'; // Reuse styles for now, or create a dedicated one

// --- Visualization Types (copied from Visualization.tsx) ---
const visualizationTypes = [
  "Co-occurrence Density",
  "Ligand Density",
  "Receptor Density",
  "Ligand Density at Receptors",
  "Receptor Density at Ligands",
] as const;
type VisualizationType = typeof visualizationTypes[number];

// --- Average Scores Type (copied from Visualization.tsx) ---
interface AverageScores {
  co_occurrence?: number | null;
  ligand_density?: number | null;
  receptor_density?: number | null;
}

// --- Component Props ---
interface DensityMapControlsProps {
  ligandName: string;
  receptorName: string;
  currentScope: string; // Scope (e.g., layer) is needed for the API call
}

const DensityMapControls: React.FC<DensityMapControlsProps> = ({ 
    ligandName, 
    receptorName, 
    currentScope 
}) => {
    
  // --- State for Density Maps (moved from Visualization.tsx) ---
  const [selectedVizType, setSelectedVizType] = useState<VisualizationType>(visualizationTypes[0]);
  const [densityMapUrl, setDensityMapUrl] = useState<string | null>(null);
  const [averageScores, setAverageScores] = useState<AverageScores | null>(null);
  const [isLoadingMap, setIsLoadingMap] = useState<boolean>(false);
  const [mapError, setMapError] = useState<string | null>(null);
  // -----------------------------------------------------------

  // --- Function to Generate Density Map (moved from Visualization.tsx) ---
  const generateDensityMap = useCallback(async () => {
    if (!ligandName || !receptorName || !currentScope) {
        setMapError("Missing ligand, receptor, or scope information.");
        return;
    }
    setIsLoadingMap(true);
    setMapError(null);
    setDensityMapUrl(null);
    setAverageScores(null);

    const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000'; 

    try {
        console.log(`Requesting density map: L=${ligandName}, R=${receptorName}, Scope=${currentScope}, Type=${selectedVizType}`);
        const response = await axios.post(`${API_BASE_URL}/api/spatial/visualize/`, {
            ligand_name: ligandName,
            receptor_name: receptorName,
            scope_name: currentScope,
            visualization_type: selectedVizType,
        });

        console.log("Density Map API Response:", response.data);

        if (response.data.error) {
            setMapError(response.data.error);
        } else if (response.data.image_path) {
            // Construct full URL for the image - ensure static path matches backend setup
            setDensityMapUrl(`${API_BASE_URL}/static/${response.data.image_path}`); 
            setAverageScores(response.data.average_scores || {});
        } else {
            setMapError("Received an unexpected response from the server.");
        }
    } catch (error: any) {
        console.error("Error generating density map:", error);
        let errorMsg = "Failed to generate density map.";
        if (axios.isAxiosError(error) && error.response) {
            errorMsg = error.response.data.detail || errorMsg;
        } else if (error instanceof Error) {
            errorMsg = error.message;
        }
        setMapError(errorMsg);
    } finally {
        setIsLoadingMap(false);
    }
  }, [ligandName, receptorName, currentScope, selectedVizType]); // Dependencies for useCallback
  // -----------------------------------------------------------------------

  return (
    <div className={styles.densityMapSection}> 
      <h4>Density Map Visualization</h4>
      <div className={styles.densityMapControls}>
        <label htmlFor="vizTypeSelect">Select Visualization:</label>
        <select 
          id="vizTypeSelect"
          value={selectedVizType}
          onChange={(e) => setSelectedVizType(e.target.value as VisualizationType)}
          disabled={isLoadingMap}
          className={styles.vizSelect}
        >
          {visualizationTypes.map(type => (
            <option key={type} value={type}>{type}</option>
          ))}
        </select>
        <button 
          onClick={generateDensityMap}
          disabled={isLoadingMap || !ligandName || !receptorName || !currentScope}
          className={styles.generateButton}
          title={!ligandName || !receptorName || !currentScope ? "Select Ligand, Receptor, and Scope first" : "Generate Density Map"}
        >
          {isLoadingMap ? 'Generating...' : 'Generate Map'}
        </button>
      </div>

      {mapError && <p className={styles.errorMessage}>Error: {mapError}</p>}

      {/* Only show display area if loading or if map exists */}
      {(isLoadingMap || densityMapUrl) && (
          <div className={styles.densityMapDisplay}>
              {isLoadingMap && <p>Loading map...</p>}
              {densityMapUrl && !isLoadingMap && (
                  <img src={densityMapUrl} alt={`${selectedVizType} for ${ligandName}-${receptorName}`} className={styles.densityMapImage}/>
              )}
              {averageScores && !isLoadingMap && (
                  <div className={styles.scoresDisplay}>
                      <h5>Average Scores:</h5>
                      <ul>
                          {Object.entries(averageScores).map(([key, value]) => (
                              value !== null && value !== undefined && (
                                  <li key={key}>{key.replace('_', ' ').replace(/\b\w/g, l => l.toUpperCase())}: {value.toFixed(4)}</li>
                              )
                          ))}
                      </ul>
                  </div>
              )}
          </div>
      )}
    </div>
  );
};

export default DensityMapControls; 