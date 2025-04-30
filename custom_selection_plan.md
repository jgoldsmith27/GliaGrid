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

### 1. Backend API Endpoint (Phase 1 - All Points Data)

*   **Requirement:** A way to fetch all raw spatial points (`x`, `y`, `layer`) for the initial display in the custom selection mode.
*   **Confirmation:** Code review confirms the existing `/api/visualization/{job_id}` endpoint is designed for specific L-R pairs within a predefined scope (layer or whole tissue) and cannot serve this purpose.
*   **Action:** Implement the **new endpoint** as planned:
    *   **Endpoint:** `GET /api/points/{jobId}/all`
    *   **Functionality:** Retrieve all spatial points (`x`, `y`, `layer`) associated with the given `jobId` from the underlying data source.
    *   **Response Format (Example):**
        ```json
        [
          { \"x\": 10.5, \"y\": 25.1, \"layer\": \"Layer 1\" },
          { \"x\": 12.3, \"y\": 28.9, \"layer\": \"Layer 2\" },
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

*   **Lasso Selection:** Integrate Deck.gl's `EditableGeoJsonLayer` or a similar mechanism into `SpatialOverviewVisualization` to allow users to draw a polygon (lasso) when `currentScope === 'custom'`.
*   **Point Capture:** Capture the coordinates of the points from `allPointsData` that fall within the user-drawn polygon.
*   **Custom Analysis Trigger:** Implement a mechanism (e.g., a button appearing after selection) for the user to initiate the analysis pipeline on the subset of captured points. This will involve:
    *   **Backend Endpoint & Service:** 
        *   **Confirmation:** Code review reveals an existing, partially implemented endpoint `POST /api/analysis/start/custom_selection` designed for this purpose.
        *   **Action:** The primary backend task is to **fully implement the service method** called by this endpoint: `AnalysisService.run_custom_analysis_background(job_id, request)`.
        *   This service method must handle the received point list, potentially adapt/call the core analysis pipeline, and manage job status updates via `JobService`.
    *   **Core Analysis Pipeline:** The underlying analysis logic (pathway dominance, module context, etc.) may need **adaptation** to correctly process an arbitrary subset of points instead of just predefined scopes (whole tissue/layers).
    *   **Frontend Logic:** Calling the *existing* `POST /api/analysis/start/custom_selection` endpoint with the captured points and handling the asynchronous analysis status (polling/WebSocket).
*   **Displaying Custom Results:** Once the custom analysis is complete and results are available:
    *   **State Management:** Introduce state in `ResultsPage.tsx` or a dedicated hook to store the results (scores, interaction coordinates, etc.) specifically for the custom run.
    *   **`SummaryTabContent.tsx` Update:**
        *   Modify the conditional rendering logic. The `AnalysisTable` should reappear when custom results are available, displaying the scores calculated *for the custom region*.
        *   The `InteractionVisualization` should be rendered when a user clicks a row in the custom results table, visualizing the interactions *within the custom region*. The necessary data (coordinates, scope details for the custom run) needs to be passed down.
    *   **Data Flow & Result Storage:** 
        *   Ensure the data fetching logic (potentially within `useInteractionData` or a new hook) can request and handle visualization data based on the completed custom analysis run.
        *   Define how custom analysis results are stored and identified on the backend (e.g., associated with the original `jobId` but under a unique custom run identifier or scope).
*   **Scope within Custom Results (Optional):** Consider adding controls (similar to the main Scope Selector) that appear *after* custom analysis is done, allowing the user to view the custom results aggregated across the whole selection (\"Whole Selection\") or broken down by the original layers present within that selection (\"Layers\"). This would require the custom analysis results to be structured accordingly (e.g., results available per layer within the selection).
*   **UI Refinements:** Adjust titles and UI elements to clearly indicate that the user is viewing results from their custom selection. 