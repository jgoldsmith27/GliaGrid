# GliaGrid Explorer (Placeholder Name)

## Overview

GliaGrid Explorer is an open-source desktop application designed to provide a user-friendly interface for analyzing spatial transcriptomics data, focusing on ligand-receptor interactions, pathway dominance, and module context within tissue layers or the whole sample. It offers powerful visualization tools for exploring spatial patterns and comparing different regions or datasets.

## Goal

To create an intuitive, cross-platform application that empowers researchers to:
1. Easily load and map their spatial, interaction, and module data.
2. Run a comprehensive analysis pipeline based on the methods in the `layer_analysis` scripts.
3. Visualize spatial data, interaction scores, and module context interactively.
4. Perform advanced comparisons between tissue regions, layers, or multiple samples.

## Technology Stack

*   **Frontend:** Electron, React, Plotly.js
*   **Backend:** Python 3.x, FastAPI
*   **Core Python Libraries:** pandas, numpy, scipy, scikit-learn, anndata, matplotlib, seaborn, networkx

## Features

### 1. Data Input & Management
    *   Supports CSV and H5AD file formats.
    *   Intuitive interface for uploading:
        *   **Spatial Data:** Requires columns for Gene ID, X coordinate, Y coordinate, and Layer. Allows multiple file uploads for comparison.
        *   **Interaction Data:** Requires columns for Ligand, Receptor. (Note: Complex receptors should use underscore separation, e.g., `ReceptorA_ReceptorB`).
        *   **Module Data:** Requires columns for Gene ID, Module ID.
    *   **Column Mapping:** User-friendly interface to map columns from uploaded files to the required data fields (Gene, X, Y, Layer, Ligand, Receptor, Module). All mappings are required.
    *   **Data Preview:** Display headers and the first 5 rows of uploaded files for verification.
    *   Ability to revisit and modify data inputs.

### 2. Analysis Pipeline
    *   Adapts logic from the reference `layer_analysis` scripts.
    *   **Analysis Scope:** Performs calculations for **both the entire tissue and layer-by-layer**. Results for both scopes will be clearly presented and accessible in the UI.
    *   **Processing Stages:**
        *   **Stage 1 (Initial):** Calculates counts of unique ligands and receptors present in the spatial data (whole tissue and per layer).
        *   **Stage 2 (Concurrent):** Once Stage 1 is complete, Pathway Dominance and Module Context analyses can run concurrently.
            *   **Pathway Dominance:** Uses Min-Max scaling on expression/proximity scores to identify and rank probable ligand-receptor interactions (whole tissue and per layer).
            *   **Module Context:** Determines module assignments for ligands/receptors in significant pairs, classifying interactions as intra-module or inter-module (whole tissue and per layer).
    *   **Asynchronous Processing:** Backend tasks run asynchronously with clear progress indicators in the UI.
    *   **Optimization Note:** Intermediate results (e.g., gene counts, normalized expression, coordinate subsets) will be cached or stored efficiently (e.g., using Feather/Parquet files or in-memory structures where feasible) to avoid redundant calculations, especially for visualization updates and comparisons.

### 3. Visualization Tools

    *   **Basic Sample View:**
        *   Displays spatial coordinates of spots.
        *   Color-codes spots by layer.
    *   **Advanced Interaction Visualization:**
        *   Overlay pathway/interaction information onto the spatial plot.
        *   Visualize interaction scores spatially using different metrics:
            *   **Co-occurrence Score:** Based on Ligand KDE * Receptor KDE.
            *   **Ligand Density Score:** Ligand KDE evaluated at receptor locations.
            *   **Receptor Density Score:** Receptor KDE evaluated at receptor locations.
        *   Option to view and compare these different score visualizations simultaneously.
    *   **Interactivity:**
        *   Click on layers/sections to zoom/focus (with ability to revert).
        *   Lasso tool to select arbitrary regions of interest.
        *   Line tool to divide the sample into two sections for comparison.

### 4. Comparison Features
    *   **Flexible Comparisons:** The UI will allow users to define comparisons between various selections:
        *   A lassoed region vs. the entire tissue.
        *   A lassoed region vs. one or more specific layers.
        *   A lassoed region vs. one or more other lassoed regions.
        *   A lassoed region vs. a section defined by the line tool.
        *   Comparisons between layers, sections, etc., are also implicitly supported through the selection mechanism.
        *   Comparison of multiple uploaded spatial datasets.
    *   **Comparative Analysis:** Identifies and visualizes differentially prominent interactions, pathways, or module statistics between the selected groups.
    *   **Backend Design:** Core comparison logic will be designed to handle comparisons between 2 to N distinct subsets of data points, regardless of how those subsets were defined (lasso, layer, slice, etc.).

## Dependencies (Preliminary)

*   **Python:** pandas, numpy, matplotlib, scipy, scikit-learn, seaborn, anndata, openpyxl, FastAPI, uvicorn, python-dotenv, networkx, pyarrow (for Feather/Parquet)
*   **Node.js/Electron:** electron, react, react-dom, plotly.js, node package manager (npm or yarn)

## Project Structure (Initial Idea)
