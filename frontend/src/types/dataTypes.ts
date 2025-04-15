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
}

export interface FilePreviewResult {
  headers?: string[];
  previewRows?: Record<string, any>[];
  error?: string;
}

export type FileType = 'spatial' | 'interactions' | 'modules';

// Mapping requirements by file type
export const requiredColumns: Record<FileType, string[]> = {
  spatial: ['geneCol', 'xCol', 'yCol'], // layerCol is optional
  interactions: ['ligandCol', 'receptorCol'],
  modules: ['geneCol', 'moduleCol'],
};

// Mapping field definitions for UI
export const mappingFields: Record<FileType, Array<{key: keyof FileState, label: string}>> = {
  spatial: [
    { key: 'geneCol', label: "Map 'Gene ID *' to:" },
    { key: 'xCol', label: "Map 'X Coordinate *' to:" },
    { key: 'yCol', label: "Map 'Y Coordinate *' to:" },
    { key: 'layerCol', label: "Map 'Layer (Optional)' to:" }
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