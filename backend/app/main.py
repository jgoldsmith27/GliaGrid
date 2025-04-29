#!/usr/bin/env python
# -*- coding: utf-8 -*-
# print("--- EXECUTING main.py ---") # <<< REMOVED

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
# Import all necessary routers
from .api import visualization # <<< RE-ENABLED IMPORT
from .api import file_routes, analysis_routes

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

@app.get("/")
async def root():
    return {"message": "GliaGrid Explorer API"} 