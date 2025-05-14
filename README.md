# GliaGrid

GliaGrid is a desktop application for analyzing and visualizing spatial molecular data, with a focus on ligand-receptor interactions and complex receptor modeling.

## Prerequisites

- **Node.js** (v16+)
- **npm** (v8+)
- **Python** (3.8+)
- **pip**
- **Git** (for cloning the repository)

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/GliaGrid.git
cd GliaGrid
```

### 2. Frontend Setup

```bash
# Install root-level dependencies
npm install

# Install frontend dependencies
cd frontend
npm install
cd ..
```

### 3. Backend Setup

```bash
# Navigate to backend directory
cd backend

# Create Python virtual environment
python -m venv venv

# Activate virtual environment
# On macOS/Linux:
source venv/bin/activate
# On Windows:
# venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt

# Return to root directory
cd ..
```

## Running the Application

You can start the entire application (frontend, backend, and Electron) with a single command:

```bash
npm start
```

This command concurrently runs:
- React frontend server (http://localhost:5173)
- Python backend server (http://localhost:8000)
- Electron desktop application wrapping both

## Development Details

### Project Structure

```
GliaGrid/
├── frontend/           # React frontend code
├── backend/            # Python FastAPI backend
├── main.js             # Electron main process
├── preload.js          # Electron preload script
└── package.json        # Node.js dependencies and scripts
```

### Scripts

- `npm start` - Start the complete application
- `npm run start-react` - Start only the frontend
- `npm run start-backend` - Start only the backend
- `npm run start-electron` - Start only the Electron app (requires frontend to be running)

## Troubleshooting

### Common Issues

1. **Backend server fails to start**
   - Ensure Python virtual environment is activated
   - Check that all dependencies are installed: `pip install -r requirements.txt`
   - Verify that port 8000 is not in use by another application

2. **Frontend development server fails**
   - Clear npm cache: `npm cache clean --force`
   - Delete node_modules and reinstall: `rm -rf node_modules && npm install`

3. **Large files handling issues**
   - Make sure you have the latest version with large files ignored in .gitignore
   
## Future Development

- Add downstream support for H5AD file handlers in rest of pipeline
- Add MIDcount across app for analysis and visualizations
- Implement comparisons with two projects