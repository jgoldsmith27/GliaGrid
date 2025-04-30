import React, { useState } from 'react';
import { Box, Button, IconButton, Tooltip, Typography, Badge, Menu, MenuItem } from '@mui/material';
import { LassoSelect, Cancel, PlayArrow, Save, Download } from '@mui/icons-material';
import { LassoState, Point } from '../../hooks/useCustomSelection';

interface CustomSelectionToolbarProps {
  lassoState: LassoState;
  enabled: boolean;
  toggleLassoMode: () => void;
  clearSelection: () => void;
  onRunAnalysis: (selectedPoints: Point[]) => void;
  onSaveSelection: (selectedPoints: Point[]) => void;
  onExportSelection: (selectedPoints: Point[]) => void;
}

const CustomSelectionToolbar: React.FC<CustomSelectionToolbarProps> = ({
  lassoState,
  enabled,
  toggleLassoMode,
  clearSelection,
  onRunAnalysis,
  onSaveSelection,
  onExportSelection
}) => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const menuOpen = Boolean(anchorEl);
  
  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  
  const handleClose = () => {
    setAnchorEl(null);
  };
  
  const handleSave = () => {
    onSaveSelection(lassoState.selectedPoints);
    handleClose();
  };
  
  const handleExport = () => {
    onExportSelection(lassoState.selectedPoints);
    handleClose();
  };
  
  const selectionCount = lassoState.selectedPoints.length;
  const hasSelection = selectionCount > 0;

  return (
    <Box 
      sx={{ 
        display: 'flex', 
        alignItems: 'center', 
        p: 1, 
        backgroundColor: 'background.paper',
        borderRadius: 1,
        boxShadow: 1
      }}
    >
      <Tooltip title={enabled ? "Disable Lasso Selection" : "Enable Lasso Selection"}>
        <IconButton 
          color={enabled ? "primary" : "default"}
          onClick={toggleLassoMode}
          sx={{ mr: 1 }}
        >
          <LassoSelect />
        </IconButton>
      </Tooltip>
      
      {enabled && (
        <Tooltip title="Clear Selection">
          <span>
            <IconButton 
              onClick={clearSelection}
              disabled={!hasSelection}
              sx={{ mr: 1 }}
            >
              <Cancel />
            </IconButton>
          </span>
        </Tooltip>
      )}
      
      {enabled && (
        <Badge 
          badgeContent={selectionCount} 
          color="primary"
          sx={{ mr: 2 }}
        >
          <Typography variant="body2">
            Selected Points
          </Typography>
        </Badge>
      )}
      
      {hasSelection && (
        <>
          <Tooltip title="Run Analysis on Selection">
            <Button
              variant="contained"
              color="primary"
              startIcon={<PlayArrow />}
              onClick={() => onRunAnalysis(lassoState.selectedPoints)}
              size="small"
              sx={{ mr: 1 }}
            >
              Run Analysis
            </Button>
          </Tooltip>
          
          <Tooltip title="More Options">
            <IconButton
              onClick={handleClick}
            >
              <Save />
            </IconButton>
          </Tooltip>
          
          <Menu
            anchorEl={anchorEl}
            open={menuOpen}
            onClose={handleClose}
          >
            <MenuItem onClick={handleSave}>
              <Save fontSize="small" sx={{ mr: 1 }} />
              Save Selection
            </MenuItem>
            <MenuItem onClick={handleExport}>
              <Download fontSize="small" sx={{ mr: 1 }} />
              Export Data
            </MenuItem>
          </Menu>
        </>
      )}
    </Box>
  );
};

export default CustomSelectionToolbar; 