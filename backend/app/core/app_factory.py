"""Application factory for creating FastAPI instances."""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.core.config import settings
from app.api import file_routes


def create_app() -> FastAPI:
    """Create and configure a FastAPI application.
    
    Returns:
        FastAPI: Configured FastAPI application
    """
    # Create FastAPI app
    app = FastAPI(
        title=settings.PROJECT_NAME,
        openapi_url=f"{settings.API_V1_STR}/openapi.json",
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=file_routes.lifespan  # Register lifespan for temp directory cleanup
    )
    
    # Add CORS middleware
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.ALLOW_ORIGINS,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    
    # Include routers
    app.include_router(file_routes.router)
    
    # Add health check endpoint
    @app.get("/ping")
    async def ping():
        """Basic health check endpoint."""
        return {"message": "Backend is alive!"}
    
    return app 