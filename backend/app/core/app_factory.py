"""Application factory for creating FastAPI instances."""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware import Middleware

from app.core.config import settings
from app.api import file_routes
from app.api import analysis_routes

# Define middleware list
middleware = [
    Middleware(
        CORSMiddleware,
        allow_origins=[str(origin) for origin in settings.ALLOW_ORIGINS],
        allow_origin_regex=r"^(http://localhost:5173|http://localhost:8000|file://.*)$",
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
        # trust_origins=True # Removed - not supported by current Starlette version
    )
]

def create_app() -> FastAPI:
    """Create and configure a FastAPI application.
    
    Returns:
        FastAPI: Configured FastAPI application
    """
    # Create FastAPI app, passing middleware directly
    app = FastAPI(
        title=settings.PROJECT_NAME,
        openapi_url=f"{settings.API_V1_STR}/openapi.json",
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=file_routes.lifespan,  # Register lifespan for temp directory cleanup
        middleware=middleware # Pass middleware here
    )
    
    # Include routers
    app.include_router(file_routes.router, prefix=settings.API_V1_STR)
    app.include_router(analysis_routes.router, prefix=settings.API_V1_STR)
    app.include_router(analysis_routes.ws_router)
    
    # Add health check endpoint
    @app.get("/ping")
    async def ping():
        """Basic health check endpoint."""
        return {"message": "Backend is alive!"}
    
    return app 