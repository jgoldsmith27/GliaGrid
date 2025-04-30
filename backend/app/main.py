#!/usr/bin/env python
# -*- coding: utf-8 -*-
# print("--- EXECUTING main.py ---") # <<< REMOVED

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
# Import all necessary routers
from .api import visualization # <<< RE-ENABLED IMPORT
from .api import file_routes, analysis_routes # <<< Add spatial_data_routes


# --- DEBUG: Check if router object exists after import ---
# print(f"DEBUG: visualization.router object: {visualization.router}") # <<< REMOVED
# print(f"DEBUG: file_routes.router object: {file_routes.router}") # <<< REMOVED
# print(f"DEBUG: analysis_routes.router object: {analysis_routes.router}") # <<< REMOVED
# print(f"DEBUG: analysis_routes.ws_router object: {analysis_routes.ws_router}") # <<< REMOVED
# ------------------------------------------------------

app = FastAPI()

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, replace with specific origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
# Include file routes (adjust prefix if needed, assuming /api/files)
app.include_router(file_routes.router, prefix="/api/files", tags=["Files"]) 
# Include analysis HTTP routes
app.include_router(analysis_routes.router, prefix="/api", tags=["Analysis"]) # <<< ADDED PREFIX
# Include analysis WebSocket router (no prefix needed as path is defined in the router)
app.include_router(analysis_routes.ws_router, tags=["Analysis WebSocket"])
# Include visualization routes (already present, but good to confirm)
app.include_router(visualization.router, prefix="/api", tags=["Visualization"]) # <<< RE-ENABLED ROUTER
# Remove spatial data router inclusion
# app.include_router(spatial_data_routes.router, prefix="/api", tags=["Spatial Data"])
# app.include_router(density_maps.router, prefix="/api", tags=["Spatial Visualization"])

# --- REMOVE DEBUG: Print registered routes before starting --- 
# for route in app.routes:
#     if hasattr(route, "path"):
#         print(f"Registered Route: Path={route.path}, Name={route.name}, Methods={getattr(route, 'methods', 'N/A')}")
#     elif hasattr(route, "path_format"):
#         # Handle WebSocket routes
#         print(f"Registered WS Route: Path={route.path_format}")
# print("--- Finished listing routes ---", flush=True)
# ------------------------------------------------------

@app.get("/")
async def root():
    return {"message": "GliaGrid Explorer API"} 