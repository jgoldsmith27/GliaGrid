# Custom Selection Feature Implementation Plan

This document outlines the plan to implement the "Custom Selection" feature, allowing users to visualize all spatial data points, select a region using a lasso tool, and run analysis on that custom region.

**Phase 1: Display All Points by Layer**

The first goal is to modify the visualization to display all spatial points colored by their assigned layer, with a legend and toggles to control layer visibility.

## Proposed Architecture

Leverage the existing components (`SummaryTabContent`, `Visualization`) with modifications to handle a new "custom" scope.

1.  **Backend:** Introduce a new API endpoint to serve all spatial points for a job ID, including layer information.
2.  **Frontend Data Fetching (`SummaryTabContent.tsx`):** When the scope is "custom", fetch data from the new endpoint. Pass this data to the `Visualization` component via a new prop.
3.  **Frontend Visualization (`Visualization.tsx`):** Modify the component to accept and render the "all points" data when the scope is "custom". Handle layer coloring, legend display, and visibility toggling.

## Implementation Steps

### 1. Backend API Endpoint

*   **Endpoint:** `GET /api/points/{jobId}/all`
*   **Functionality:** Retrieve all spatial points (`x`, `y`, `layer`) associated with the given `jobId`.
*   **Response Format (Example):**
    ```json
    [
      { "x": 10.5, "y": 25.1, "layer": "Layer 1" },
      { "x": 12.3, "y": 28.9, "layer": "Layer 2" },
      // ... other points
    ]
    ```

### 2. Frontend: `SummaryTabContent.tsx` Modifications

*   **State:** Add state for `allPointsData`, `loadingAllPoints`, `allPointsError`.
    ```typescript
    interface AllPointsData { x: number; y: number; layer: string; }
    const [allPointsData, setAllPointsData] = useState<AllPointsData[] | null>(null);
    // ... loading/error states
    ```
*   **Fetching Logic:**
    *   Create `fetchAllPoints` function targeting the new `/api/points/{jobId}/all` endpoint.
    *   Modify the main `useEffect` (watching `selectedPair`, `jobId`, `currentScope`) to call `fetchAllPoints` when `currentScope === 'custom'`, and the existing `fetchVisualizationData` otherwise. Handle state resets and cleanup.
*   **Rendering `<Visualization>`:** Conditionally render based on `currentScope`:
    ```typescript
    {currentScope === 'custom' ? (
        // Render based on loadingAllPoints, allPointsError, allPointsData
        <Visualization
            allPointsData={allPointsData} // New prop
            currentScope={currentScope}
            // No ligandName/receptorName needed
        />
    ) : (
        // Existing rendering logic using visualizationData
        <Visualization
            data={visualizationData}
            ligandName={selectedPair[0]}
            receptorName={selectedPair[1]}
            currentScope={currentScope} // e.g., 'whole_tissue', 'layers'
        />
    )}
    ```

### 3. Frontend: `Visualization.tsx` Modifications

*   **Props (`VisualizationProps`):**
    *   Add `allPointsData?: AllPointsData[] | null;`
    *   Make `data`, `ligandName`, `receptorName` optional (`?`).
    *   Update `currentScope` type: `type ScopeType = 'whole_tissue' | 'layers' | 'custom';`
*   **State:** Add state for layer colors and visibility.
    ```typescript
    const [visibleLayers, setVisibleLayers] = useState<Set<string>>(new Set());
    const [layerColors, setLayerColors] = useState<Record<string, import('@deck.gl/core').Color>>({});
    ```
*   **Effect (`useEffect`):** When `currentScope === 'custom'` and `allPointsData` is available:
    *   Extract unique layer names.
    *   Generate and store colors for each layer in `layerColors` (use a color generation library/function).
    *   Initialize `visibleLayers` (e.g., all visible).
    *   Reset state otherwise.
*   **Layer Toggling:** Implement `handleLayerToggle` function to update `visibleLayers` state.
*   **Layer Memoization (`useMemo`):** Modify the `layers` calculation:
    *   **If `currentScope === 'custom'`:**
        *   Filter `allPointsData` based on `visibleLayers`.
        *   Return a single `ScatterplotLayer` using the filtered data.
        *   Use `getFillColor: d => layerColors[d.layer]`
    *   **Else (if `currentScope !== 'custom'`):**
        *   Use the *existing* logic based on `props.data` (ligand/receptor points, density scoring, etc.). Ensure this logic handles optional `props.data`.
    *   Return `[]` otherwise.
*   **Render Legend & Toggles:** Inside the `.deckOverlayControls`, add a conditional block for `currentScope === 'custom'`:
    *   Render a list/group of checkboxes, one for each layer found in `layerColors`.
    *   Bind `checked` to `visibleLayers.has(layerName)`.
    *   Bind `onChange` to `handleLayerToggle(layerName)`.
    *   Display the layer name and its color swatch.
    *   Conditionally render the *existing* controls (ligand/receptor legend, density options) only when `currentScope !== 'custom'`.
*   **CSS:** Add necessary styles for the new layer legend and toggles.

## Benefits

*   **Modular:** Keeps concerns separated (data fetching vs. rendering).
*   **Decoupled:** Backend, data fetching parent, and visualization child have distinct roles.
*   **Reuses Existing Structure:** Builds upon the current components.

## Phase 2 (Future)

*   Integrate Deck.gl's `EditableGeoJsonLayer` or a similar mechanism for lasso selection when `currentScope === 'custom'`.
*   Capture selected points based on the drawn polygon.
*   Implement backend endpoint and frontend logic to trigger analysis on the selected points. 