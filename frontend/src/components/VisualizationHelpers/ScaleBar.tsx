import React from 'react';
import styles from './ScaleBar.module.css';

interface ScaleBarProps {
  currentZoom: number;
  unitsPerMicron: number;
  targetPixelWidth?: number; // Approximate desired width of the bar on screen
}

const NICE_MICRON_VALUES = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];

const ScaleBar: React.FC<ScaleBarProps> = ({
  currentZoom,
  unitsPerMicron,
  targetPixelWidth = 100, // Aim for a bar around 100px wide
}) => {
  // Calculate how many world units fit into the targetPixelWidth at the current zoom
  // targetPixelWidth = worldUnits * Math.pow(2, currentZoom)
  // worldUnits = targetPixelWidth / Math.pow(2, currentZoom)
  const worldUnitsForTargetWidth = targetPixelWidth / Math.pow(2, currentZoom);
  const micronsForTargetWidth = worldUnitsForTargetWidth / unitsPerMicron;

  // Find the "nicest" micron value that is less than or equal to micronsForTargetWidth
  let bestMicronValue = NICE_MICRON_VALUES[0];
  for (const niceValue of NICE_MICRON_VALUES) {
    if (niceValue <= micronsForTargetWidth) {
      bestMicronValue = niceValue;
    } else {
      break; // Stop when we exceed the target
    }
  }
  // If micronsForTargetWidth is very small, we might need to adjust or show a smaller unit (not handled here)
  if (micronsForTargetWidth < NICE_MICRON_VALUES[0]) {
      // Potentially handle cases where even 1 micron is too large for the targetPixelWidth
      // For now, default to smallest, or consider alternative rendering
      bestMicronValue = NICE_MICRON_VALUES[0]; 
      // Or, if it's extremely zoomed in, a scale bar might not make sense or needs different logic
      // This can happen if targetPixelWidth is very small or zoom is very high.
      // We could also try to find the closest nice value instead of just <=.
      // For simplicity now, we stick to this.
  }


  const scaleBarMicronLength = bestMicronValue;
  const scaleBarCoordinateLength = scaleBarMicronLength * unitsPerMicron;
  const scaleBarPixelLength = scaleBarCoordinateLength * Math.pow(2, currentZoom);

  if (scaleBarPixelLength < 1) { // Avoid rendering a bar that's too small to see
    return null;
  }

  return (
    <div className={styles.scaleBarContainer}>
      <div className={styles.scaleBar} style={{ width: `${scaleBarPixelLength}px` }} />
      <div className={styles.scaleBarLabel}>{`${scaleBarMicronLength} µm`}</div>
    </div>
  );
};

export default ScaleBar; 