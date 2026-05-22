# Running NGS Results Explorer

## API (port 8000)
```bash
cd ~/ngs-results-explorer && source .venv/bin/activate && uvicorn src.api:app --reload --port 8000
```
Open: http://localhost:8000/docs

## Dashboard (port 8052)
```bash
cd ~/ngs-results-explorer && source .venv/bin/activate && python3 -m src.dashboard
```
Open: http://localhost:8052
