# Comparison Tool Feature Documentation

## 1. Introduction & Core Goal

The core goal of the Comparison Tool is to enable users to perform detailed molecular comparisons between different spatial datasets or regions within the same dataset. This includes comparing whole tissues, specific layers, or user-defined lassoed sections. The tool aims to identify statistically significant differences in molecular expression, providing insights into varying biological states or conditions.

## 2. Comparison Scenarios

The tool must support a flexible range of comparison scenarios:

*   **Project A vs. Project B:**
    *   Whole Tissue (A) vs. Whole Tissue (B)
    *   Layer (A) vs. Layer (B)
    *   Lassoed Section (A) vs. Lassoed Section (B)
    *   Whole Tissue (A) vs. Layer (B)
    *   Whole Tissue (A) vs. Lassoed Section (B)
    *   Layer (A) vs. Lassoed Section (B)
    *   And all permutations thereof.
*   **Within the Same Project (Project A vs. Project A):**
    *   Whole Tissue vs. Layer
    *   Whole Tissue vs. Lassoed Section
    *   Layer vs. Lassoed Section
    *   Lassoed Section 1 vs. Lassoed Section 2
    *   Layer 1 vs. Layer 2

## 3. Data Points & Metrics for Comparison

The following metrics will be calculated and compared:

### 3.1. Individual Ligand/Receptor Expression
*   **Description:** Identify specific ligands or receptors that are significantly more or less abundant in one selection compared to the other.
*   **Output:** A list of differentially expressed molecules, showing fold change (or similar metric) and statistical significance (e.g., p-value, FDR). Only results passing a defined significance threshold will be highlighted.
*   **Considerations:** Normalization of expression values is crucial.

### 3.2. Ligand-Receptor Pair Analysis
*   **Description:** Analyze the co-occurrence or potential interaction frequency of specific ligand-receptor pairs.
*   **Output:** Metrics indicating the strength or frequency of LR pairing in each selection, and a comparison of these metrics.
*   **Considerations:** This might involve proximity analysis or scoring based on co-localization.

## 4. Statistical Methods & Analysis Pipeline

### 4.1. Differential Expression of Individual Molecules
*   **Input:** Counts or normalized expression values of ligands/receptors per cell or per unit area for each selected region.
*   **Normalization:**
    *   Consider methods like Counts Per Million (CPM), TPM, or library size normalization.
    *   For spatial data, normalization by area might also be relevant.
    *   The existing backend pipeline uses MinMaxScaler on relative frequency for pathway dominance; evaluate if this or a similar per-selection normalization is appropriate before comparison, or if group-wise normalization is needed.
*   **Statistical Tests:**
    *   For comparing two groups (Selection 1 vs. Selection 2):
        *   If data resembles normalized expression values: t-tests (if distributions are normal) or non-parametric tests like Wilcoxon rank-sum test / Mann-Whitney U.
        *   If data is count-based (e.g., occurrences of `geneID` per region): methods like those used in DESeq2/edgeR (e.g., Negative Binomial GLM) if applicable, or chi-squared/Fisher's exact test for proportions if comparing presence/absence or high/low categories.
    *   Consider the number of data points (e.g., cells, spots, or synthetic points if aggregating) in each selected region.
*   **Significance:** p-values adjusted for multiple comparisons (e.g., Benjamini-Hochberg for FDR). A threshold (e.g., FDR < 0.05) and a fold-change cutoff will be used.

### 4.2. Ligand-Receptor Pair Analysis
*   **Proximity-based:** Count pairs within a certain distance threshold.
*   **Statistical Significance:** Compare observed pair counts against expected counts (e.g., permutation testing by shuffling molecule labels).

## 5. Technical Implementation Plan

### 5.1. Data Requirements
*   **For each selection (tissue, layer, lasso):**
    *   Identifier for the base dataset: `fileId` (referencing a file in `backend/.temp_uploads/`, typically a CSV).
    *   Optionally, `job_id` if the selection is tied to a specific completed analysis run (useful if comparing across different original uploads/projects).
    *   Definition of the region:
        *   Type: 'whole_tissue', 'layer', 'lasso'.
        *   Specifics: `layer_name` (if 'layer'), or `polygon_coordinates` (if 'lasso').
    *   From the data file (e.g., CSV referenced by `fileId`):
        *   Coordinates (`x`, `y`) for each detected molecule.
        *   Molecule type (`geneID`).
        *   Each row is assumed to represent a single molecule detection.

### 5.2. Data Availability & Fetching
*   **Current State & Data Columns:** 
    *   The primary data files (CSVs in `.temp_uploads/`) contain `geneID` (string), `x` (float), `y` (float), and `layer` (string), among others. This provides the necessary molecular identity, coordinates, and pre-defined layer information.
    *   Lassoed region definitions (polygon coordinates) will come from the client-side UI.
*   **Fetching Strategy for Comparison Tool:**
    *   The backend `FileService` is capable of loading data from `backend/.temp_uploads/` given a `fileId`.
    *   The backend `AnalysisService` (e.g., in `run_custom_analysis`) demonstrates a pattern of loading data for a specific `job_id` (which implies access to its `fileId`s and mappings via `JobService`) and then applying polygon filtering using GeoPandas/Shapely. This is a strong precedent.
    *   For the comparison tool, the backend will:
        1.  Receive definitions for Selection 1 and Selection 2 (including `fileId`s, type, and region specifics).
        2.  For each selection, use `FileService` to load the base data DataFrame(s).
        3.  Apply layer or polygon filtering to these DataFrames to get the precise data subsets for Selection 1 and Selection 2.

### 5.3. Computation Strategy (Client-side vs. Server-side via IPC)
*   **Context:** Electron app with a main process and a Python backend.
*   **Existing Architecture:** The current system performs light data filtering and preparation in `main.js` for direct visualization needs (e.g., `readCsvChunked` with polygon filter). However, all substantial analysis (initial pipeline, custom lasso analysis) is handled by the Python backend.
*   **Recommendations for Comparison Tool:**
    *   **Data Subsetting:** Backend (Python) - leveraging `FileService` and `GeoPandas`/`Shapely` for polygon filtering, similar to `AnalysisService.run_custom_analysis`.
    *   **Statistical Calculations:** Backend (Python) - All new statistical tests for comparing the two selections should be implemented in Python. This aligns with the existing pattern and keeps computationally intensive tasks off the UI thread.
    *   **Results Transmission:** Results will be sent from the backend to `main.js` and then to the renderer process via IPC.

### 5.4. IPC Design
*   **UI to Main Process (Electron):**
    *   A new UI component will be developed for defining the two selections and triggering the comparison.
    *   This UI will likely use an existing or new `window.electronAPI` function (exposed via `preload.js`) to send the comparison request to `main.js`.
    *   `main.js` will then make an HTTP request to the new backend API endpoint.
*   **New Backend API Endpoint:**
    *   **Endpoint:** `POST /api/analysis/compare` (or similar)
    *   **Request Payload (`ComparisonRequest`):**
        ```json
        {
          "comparison_name": "UserDefinedComparison1", // Optional name for the comparison
          "selection1": {
            "source_job_id": "job_id_A_optional", // If comparing data from a specific previous job context
            "file_id": "file_id_from_spatialFile_A",
            "type": "lasso", // "whole_tissue", "layer", "lasso"
            "definition": { // Contents depend on type
              // "layer_name": "Layer1" (if type is 'layer')
              // "polygon_coords": [[x1,y1], [x2,y2], ...] (if type is 'lasso')
            },
            "column_mappings": { // Original column mappings for this file_id
                 "geneCol": "gene", "xCol": "x", "yCol": "y", "layerCol": "layer" // etc.
            }
          },
          "selection2": {
            "source_job_id": "job_id_B_optional",
            "file_id": "file_id_from_spatialFile_B",
            "type": "layer",
            "definition": {
              "layer_name": "Layer2"
            },
            "column_mappings": { // Original column mappings for this file_id
                 "geneCol": "gene", "xCol": "x", "yCol": "y", "layerCol": "layer" // etc.
            }
          },
          "analyses_to_perform": [
            {"type": "differential_expression", "params": {"fdr_threshold": 0.05}},
            {"type": "lr_pair_hotspots", "params": {"target_pairs": [["L1","R1"], ["L2","R2"]]}}
          ]
        }
        ```
    *   **Response Payload (`ComparisonResponse`):**
        ```json
        {
          "comparison_id": "unique_id_for_this_comparison_run",
          "results": {
            "differential_expression": { /* ... results ... */ },
            "lr_pair_hotspots": { /* ... results ... */ }
          },
          "errors": []
        }
        ```
*   **Status Updates:** For potentially long-running comparisons, the backend could optionally use the existing `JobService` to provide status through ZeroMQ, or the client might poll a status endpoint if the comparison is run asynchronously.

## 7. Next Steps & Open Questions

*   Refine choices of statistical tests based on data characteristics.
*   Detail the UI/UX for selecting regions and displaying comparison results.