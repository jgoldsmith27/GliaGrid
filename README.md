## Installation and Setup

### Prerequisites
- Node.js (v16+)
- npm
- Python 3.8+
- pip

### Clone the Repository
```bash
git clone https://github.com/yourusername/GliaGrid.git
cd GliaGrid
```

### Setup
1. Install root dependencies:
```bash
npm install
```

2. Install frontend dependencies:
```bash
cd frontend
npm install
cd ..
```

3. Setup Python virtual environment:
```bash
cd backend
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
pip install -r requirements.txt
cd ..
```

### Running the Application
```bash
npm start
```
This will concurrently start:
- Frontend React server
- Backend Python server
- Electron desktop application

### Building for Production
```bash
npm run build
```

## Troubleshooting

### Large Files
If you encounter issues with large files when cloning or pushing:

1. Make sure you have the latest version with large files ignored in .gitignore
2. If you still have issues, you can download the release version directly from the releases page. 

### ToDos
- add layer boundary anotations to regular visualizations
- Complex receptor visualization support
- add downstream support for H5AD file handlers in rest of pipeline (they do not work currently)