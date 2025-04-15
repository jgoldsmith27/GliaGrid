/**
 * Common types for data handling in the application.
 */

export interface FileState {
  file: File | null;
  fileId: string | null; // Unique ID received from backend
  headers: string[];
  previewRows: Record<string, any>[];
  fileInfo?: FilePreviewResult['fileInfo']; // Optional H5AD specific info
  error: string | null;
  isLoading: boolean;
  // Column mappings
  geneCol?: string;
  xCol?: string;
  yCol?: string;
  layerCol?: string;
  ligandCol?: string;
  receptorCol?: string;
  moduleCol?: string;
}

export interface FilePreviewResult {
  headers: string[];
  previewRows: Record<string, any>[];
  fileId?: string; // Unique ID assigned by backend after upload
  fileInfo?: {
    shape: string;
    obs_keys: string[];
    var_keys: string[];
    obsm_keys: string[];
  };
  error?: string;
}

export type FileType = 'spatial' | 'interactions' | 'modules';

// Mapping requirements by file type
export const requiredColumns: Record<FileType, (keyof FileState)[]> = {
  spatial: ['geneCol', 'xCol', 'yCol', 'layerCol'],
  interactions: ['ligandCol', 'receptorCol'],
  modules: ['geneCol', 'moduleCol'],
};

// Mapping field definitions for UI
export const mappingFields: Record<FileType, Array<{key: keyof FileState, label: string}>> = {
  spatial: [
    { key: 'geneCol', label: "Map 'Gene ID *' to:" },
    { key: 'xCol', label: "Map 'X Coordinate *' to:" },
    { key: 'yCol', label: "Map 'Y Coordinate *' to:" },
    { key: 'layerCol', label: "Map 'Layer *' to:" }
  ],
  interactions: [
    { key: 'ligandCol', label: "Map 'Ligand *' to:" },
    { key: 'receptorCol', label: "Map 'Receptor *' to:" }
  ],
  modules: [
    { key: 'geneCol', label: "Map 'Gene ID *' to:" },
    { key: 'moduleCol', label: "Map 'Module ID *' to:" }
  ],
};

// Payload for the /api/analysis/start endpoint
export interface AnalysisMapping {
  [key: string]: string; // e.g., { geneCol: "gene_name", xCol: "X", ... }
}

export interface AnalysisPayload {
  spatialFileId: string;
  spatialMapping: AnalysisMapping;
  interactionsFileId: string;
  interactionsMapping: AnalysisMapping;
  modulesFileId: string;
  modulesMapping: AnalysisMapping;
} 