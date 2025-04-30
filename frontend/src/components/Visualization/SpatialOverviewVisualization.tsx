import React, { useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { PolygonLayer, ScatterplotLayer, GeoJsonLayer } from '@deck.gl/layers'; // Correct PolygonLayer import
import { MapView } from '@deck.gl/core';
import { useJobStatus } from '../../hooks/useJobStatus';
import { LayerBoundary, JobResultOutput, PointFeature } from '../../types/analysis'; // Assuming PointFeature is defined here or elsewhere
import { Paper, CircularProgress, Typography, Box, Button, Switch, FormControlLabel, FormGroup } from '@mui/material';
import { Feature, Polygon, Point, FeatureCollection, GeoJsonProperties } from 'geojson'; // Import necessary GeoJSON types
import { booleanPointInPolygon } from '@turf/boolean-point-in-polygon'; // Correct import
import { point as turfPoint } from '@turf/helpers'; // Correct import for turf point helper

const SpatialOverviewVisualization = () => {
  const [lassoPoints, setLassoPoints] = useState<[number, number][]>([]);
  const [isDragging, setIsDragging] = useState(false);
  const [lassoPolygon, setLassoPolygon] = useState<Feature<Polygon> | null>(null);
  const [isLassoMode, setIsLassoMode] = useState(false);
  const [selectedPoints, setSelectedPoints] = useState<PointFeature[]>([]); // Assuming PointFeature is your point type

  const handleDragStart = useCallback((info: any) => {
    if (!isLassoMode) return;
    setIsDragging(true);
    setLassoPoints([info.coordinate]);
    setLassoPolygon(null); // Clear previous lasso polygon
    setSelectedPoints([]); // Clear previous selection
  }, [isLassoMode]);

  const handleDrag = useCallback((info: any) => {
    if (!isLassoMode || !isDragging) return;
    setLassoPoints(prev => [...prev, info.coordinate]);
  }, [isLassoMode, isDragging]);

  const handleDragEnd = useCallback((info: any) => {
    if (!isLassoMode || !isDragging) return;
    setIsDragging(false);
    if (lassoPoints.length > 2) {
      const finalLassoPoints = [...lassoPoints, lassoPoints[0]]; // Close the polygon
      const polygonFeature: Feature<Polygon> = {
        type: 'Feature',
        properties: {},
        geometry: {
          type: 'Polygon',
          coordinates: [finalLassoPoints]
        }
      };
      setLassoPolygon(polygonFeature);

      // TODO: Replace this with actual points data when available for selection
      // For now, it uses dummy points for demonstration.
      // This section will need to be adapted once point data fetching is implemented
      // alongside the boundary visualization.
      const dummyPointsForSelection: PointFeature[] = [
        { type: 'Feature', geometry: { type: 'Point', coordinates: [-122.4, 37.8] }, properties: { id: 'dummy1' } },
        { type: 'Feature', geometry: { type: 'Point', coordinates: [-122.45, 37.75] }, properties: { id: 'dummy2' } },
        // Add more points or fetch actual points based on view/analysis
      ];

      const pointsInLasso = dummyPointsForSelection.filter(p =>
          booleanPointInPolygon(turfPoint(p.geometry.coordinates), polygonFeature) // Use correct import
      );
      setSelectedPoints(pointsInLasso);
      console.log('Selected points:', pointsInLasso);


    } else {
      // Not enough points for a polygon
      setLassoPolygon(null);
    }
    setLassoPoints([]); // Clear drawn points after finishing drag
  }, [isLassoMode, isDragging, lassoPoints]);


  // Prepare layers for DeckGL
  const layers = useMemo(() => {
    const generatedLayers: (PolygonLayer<LayerBoundary> | GeoJsonLayer | ScatterplotLayer)[] = []; // Add explicit type

    // Layer Boundary Polygons
    if (layerBoundaries && !isLassoMode) {
       Object.keys(layerBoundaries).forEach(layerName => {
        if (layerVisibility[layerName]) {
          const boundary = layerBoundaries[layerName];
          generatedLayers.push(
            new PolygonLayer<LayerBoundary>({ // Add type parameter
              id: `polygon-${layerName}`,
              data: [boundary], // Data should be an array
              getPolygon: d => d.coordinates,
              getFillColor: layerColors[layerName] ? [...layerColors[layerName], 150] : [128, 128, 128, 150], // Use state color, alpha 150
              getLineColor: [255, 255, 255],
              getLineWidth: 1,
              pickable: true,
              autoHighlight: true,
              highlightColor: [255, 255, 255, 50],
            })
          );
        }
      });
    }

    // Lasso Polygon Visualization Layer
    if (lassoPolygon) {
        const lassoData: FeatureCollection<Polygon> = { // Wrap in FeatureCollection
            type: 'FeatureCollection',
            features: [lassoPolygon]
        };
        generatedLayers.push(
            new GeoJsonLayer({
                id: 'lasso-layer',
                data: lassoData, // Use FeatureCollection
                filled: false,
                stroked: true,
                getLineColor: [255, 0, 0, 255], // Red color for lasso
                getLineWidth: 2,
                lineWidthUnits: 'pixels',
            })
        );
    }

     // Layer for Selected Points (example using ScatterplotLayer)
     if (selectedPoints.length > 0) {
        generatedLayers.push(
            new ScatterplotLayer<PointFeature>({ // Add type parameter
                id: 'selected-points-layer',
                data: selectedPoints,
                getPosition: d => d.geometry.coordinates,
                getFillColor: [0, 255, 0, 255], // Green for selected points
                getRadius: 5, // Adjust size as needed
                radiusUnits: 'pixels',
            })
        );
     }


    return generatedLayers;
  }, [layerBoundaries, layerVisibility, layerColors, isLassoMode, lassoPolygon, selectedPoints]);

  // Tooltip handler
  const getTooltip = useCallback(({ object, layer }: { object: any, layer: any }): string | null => { // Add explicit types
    if (!object || !layer) {
      return null;
    }

    // Check if hovering over a polygon layer
    if (layer.id.startsWith('polygon-')) {
        // Assuming the data object 'd' passed to PolygonLayer was the LayerBoundary
        // And LayerBoundary has a 'name' property
       const boundaryData = object as LayerBoundary; // Type assertion
       if (boundaryData && boundaryData.name) {
         return `Layer: ${boundaryData.name}`;
       }
    }

    // Add tooltip logic for other layers if needed (e.g., selected points)

    return null; // Default no tooltip
  }, []);


  // Loading and Error States
  if (jobStatus.loading || !jobStatus.results) {
    return <CircularProgress />;
  }

  if (jobStatus.error) {
    return <Typography color="error">Error loading analysis results: {jobStatus.error}</Typography>;
  }

  // Render component
  return (
    <Paper elevation={3} style={{ height: '80vh', position: 'relative' }}>
      <Box sx={{ position: 'absolute', top: 10, left: 10, zIndex: 1, backgroundColor: 'rgba(255, 255, 255, 0.8)', padding: '5px', borderRadius: '4px' }}>
        <Typography variant="h6">Spatial Overview</Typography>
        <FormGroup>
             <FormControlLabel
                control={<Switch checked={isLassoMode} onChange={(e) => setIsLassoMode(e.target.checked)} />}
                label="Lasso Select Mode"
             />
          </FormGroup>
          <Typography variant="subtitle1">Layer Visibility:</Typography>
        <FormGroup>
          {layerNames.map(name => (
            <FormControlLabel
              key={name}
              control={
                <Box sx={{ display: 'flex', alignItems: 'center' }}>
                   <Box sx={{ width: 16, height: 16, backgroundColor: `rgb(${layerColors[name]?.join(',')})`, marginRight: 1, border: '1px solid #ccc' }} />
                  <Switch
                    checked={layerVisibility[name] ?? false}
                    onChange={() => toggleLayerVisibility(name)}
                    size="small"
                  />
                </Box>
              }
              label={name}
            />
          ))}
        </FormGroup>
         {selectedPoints.length > 0 && (
             <Typography variant="body2" sx={{ mt: 1 }}>Selected Points: {selectedPoints.length}</Typography>
         )}

      </Box>
      <DeckGL
        initialViewState={initialViewState}
        controller={true} // Enable map controls (zoom, pan)
        layers={layers}
        getTooltip={getTooltip} // Use the memoized tooltip handler
        views={new MapView({ repeat: true })} // Use MapView
        onDragStart={handleDragStart}
        onDrag={handleDrag}
        onDragEnd={handleDragEnd}
        pickingRadius={5} // Adjust picking radius if needed
      >
          {/* No base map needed for this visualization type */}
      </DeckGL>
    </Paper>
  );
};

export default SpatialOverviewVisualization; 