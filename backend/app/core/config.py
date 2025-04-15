"""Configuration settings for the application."""
import os
from pathlib import Path
from pydantic_settings import BaseSettings
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

class Settings(BaseSettings):
    """Application settings.
    
    These settings can be configured via environment variables.
    """
    # API settings
    API_V1_STR: str = "/api"
    PROJECT_NAME: str = "GliaGrid Explorer"
    
    # Server settings
    HOST: str = "127.0.0.1"  # Using IP address instead of hostname
    PORT: int = int(os.getenv("BACKEND_PORT", "8000"))
    RELOAD: bool = True
    
    # CORS settings
    ALLOW_ORIGINS: list[str] = [
        "http://localhost:5173",  # Default Vite dev server
        "http://localhost:8000",  # Backend API
        "file://"                 # Electron local files
    ]
    
    # File handling
    TEMP_DIR: Path = Path("./temp")
    
    class Config:
        env_file = ".env"
        case_sensitive = True


# Create global settings instance
settings = Settings()

# Ensure temp directory exists
settings.TEMP_DIR.mkdir(exist_ok=True) 