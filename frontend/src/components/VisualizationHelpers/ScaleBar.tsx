import React from 'react';
import styles from './ScaleBar.module.css';

interface ScaleBarProps {
  /** The current zoom level of the deck.gl OrthographicView (zoom = log2(pixels/world_unit)). */
  currentZoom: number;
  /** The number of world coordinate units that correspond to one micrometer (µm). 
   *  Based on analysis, this is typically 2 for this application (1 µm = 2 units). */
  unitsPerMicron: number;
  /** The approximate desired width of the scale bar in screen pixels. Defaults to 100. */
  targetPixelWidth?: number; 
}

/**
 * Renders a dynamic scale bar overlay for a deck.gl OrthographicView.
 * It calculates an appropriate length and label (in µm) based on the current zoom level
 * and the known relationship between world coordinate units and micrometers.
 * The displayed value is rounded to the nearest whole micron.
 */
const ScaleBar: React.FC<ScaleBarProps> = ({
  currentZoom,
  unitsPerMicron, // e.g., 2 units = 1 µm
  targetPixelWidth = 100, // Aim for a bar around 100px wide
}) => {
  // --- Calculate the Micron Value for the Scale Bar ---

  // 1. Determine the relationship between screen pixels and world units at the current zoom.
  const pixelsPerWorldUnit = Math.pow(2, currentZoom);

  // 2. Calculate how many world units would correspond to our target pixel width.
  const worldUnitsForTargetWidth = targetPixelWidth / pixelsPerWorldUnit;

  // 3. Convert this world unit width into micrometers using the provided conversion factor.
  const micronsForTargetWidth = worldUnitsForTargetWidth / unitsPerMicron;

  // 4. Round the calculated micron value to the nearest whole number (minimum 1).
  //    This will be the value displayed in the label.
  const scaleBarMicronLength = Math.max(1, Math.round(micronsForTargetWidth));

  // --- Calculate the Final Pixel Length of the Scale Bar ---
  
  // 5. Convert the *rounded* micron length (from step 4) back into world coordinate units.
  const scaleBarCoordinateLength = scaleBarMicronLength * unitsPerMicron;
  
  // 6. Convert this world coordinate length into screen pixels at the current zoom.
  //    This ensures the physical bar length matches the rounded label value.
  const scaleBarPixelLength = scaleBarCoordinateLength * pixelsPerWorldUnit;

  // Avoid rendering a bar that's too small to see or has non-positive length.
  if (scaleBarPixelLength < 1) { 
    return null;
  }

  return (
    <div className={styles.scaleBarContainer}>
      {/* The visible scale bar line, width set by calculated pixel length */}
      <div className={styles.scaleBar} style={{ width: `${scaleBarPixelLength}px` }} />
      {/* The label showing the corresponding micron length */}
      <div className={styles.scaleBarLabel}>{`${scaleBarMicronLength} µm`}</div>
    </div>
  );
};

export default ScaleBar; 