// Shared data types for input components

export interface FileState {
  file: File | null;
  headers: string[];
  previewRows: Record<string, any>[];
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
  // File info (for H5AD files)
  fileInfo?: {
    shape: string;
    obs_keys: string[];
    var_keys: string[];
    obsm_keys: string[];
  };
}

export interface FilePreviewResult {
  headers?: string[];
  previewRows?: Record<string, any>[];
  error?: string;
  fileInfo?: {
    shape: string;
    obs_keys: string[];
    var_keys: string[];
    obsm_keys: string[];
  };
}

export type FileType = 'spatial' | 'interactions' | 'modules';

// Mapping requirements by file type
export const requiredColumns: Record<FileType, string[]> = {
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