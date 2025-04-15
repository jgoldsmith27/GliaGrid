"""Main entry point for the GliaGrid backend application."""
import uvicorn
import os

from app.core.config import settings
from app.core.app_factory import create_app

# Create FastAPI app
app = create_app()

if __name__ == "__main__":
    print(f"Starting {settings.PROJECT_NAME} backend on {settings.HOST}:{settings.PORT}...")
    uvicorn.run(
        "main:app", 
        host="127.0.0.1",  # Explicitly use loopback IP to avoid hostname issues
        port=settings.PORT, 
        reload=settings.RELOAD
    )
