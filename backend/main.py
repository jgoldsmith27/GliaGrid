"""Main entry point for the GliaGrid backend application."""
import uvicorn
import os
import logging # Import logging

from app.core.config import settings
from app.core.app_factory import create_app

# --- Logging Configuration ---
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
# ---------------------------

# Create FastAPI app
app = create_app()

# The server will be run using hypercorn command via package.json script
# if __name__ == "__main__":
#     logger.info(f"Starting {settings.PROJECT_NAME} backend on {settings.HOST}:{settings.PORT}...") # Use logger
#     uvicorn.run(
#         "main:app", 
#         host="127.0.0.1",  # Explicitly use loopback IP to avoid hostname issues
#         port=settings.PORT, 
#         reload=settings.RELOAD,
#         ws="wsproto" # Try forcing wsproto for WebSocket handling
#     )
