from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from .api import visualization

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
app.include_router(visualization.router, prefix="/api")

@app.get("/")
async def root():
    return {"message": "GliaGrid Explorer API"} 