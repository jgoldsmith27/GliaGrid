from fastapi import FastAPI
import uvicorn
import os
from dotenv import load_dotenv

# Load environment variables (optional but good practice)
load_dotenv()

app = FastAPI()

@app.get("/ping")
async def ping():
    """Basic health check endpoint."""
    return {"message": "Backend is alive!"}

if __name__ == "__main__":
    port = int(os.getenv("BACKEND_PORT", 8000)) # Default to 8000 if not set
    print(f"Starting backend server on port {port}...")
    uvicorn.run("main:app", host="127.0.0.1", port=port, reload=True)
