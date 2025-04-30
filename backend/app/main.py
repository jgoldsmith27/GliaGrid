#!/usr/bin/env python
# -*- coding: utf-8 -*-
# print("--- EXECUTING main.py ---") # <<< REMOVED

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
# Import necessary routers
from .api import file_routes, analysis_routes # Removed visualization import

app = FastAPI(
    title="GliaGrid API",
    description="""
    GliaGrid backend API.
    Note: Many data retrieval endpoints have been replaced with direct Electron file access.
    This API now primarily handles:
    1. Initial file uploads and processing
    2. Analysis job creation and management 
    3. REST-based status updates (polling)
    4. Custom analysis computation
    """,
    version="2.0.0"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, replace with specific origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
# File routes for initial upload and management
app.include_router(file_routes.router, prefix="/api/files", tags=["Files"]) 
# Analysis job creation and custom analysis computation
app.include_router(analysis_routes.router, prefix="/api", tags=["Analysis"])
# Visualization routes for retrieving visualization data - REMOVED
# app.include_router(visualization.router, prefix="/api/visualization", tags=["Visualization"])

@app.get("/")
async def root():
    return {
        "message": "GliaGrid Explorer API",
        "architecture": "Hybrid architecture with direct Electron file access for data retrieval",
        "version": "2.0.0"
    } 

# Import job service instance
from .services.job_service import job_service_instance

@app.on_event("shutdown")
def cleanup_zmq():
    """Clean up ZeroMQ resources on application shutdown"""
    if hasattr(job_service_instance, 'cleanup_zmq'):
        job_service_instance.cleanup_zmq() 