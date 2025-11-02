# Primer‑BLAST3‑Demo 🧪
Containerized pipeline for automated primer design and BLAST validation.

## 🚀 Usage
Build and run locally:
```bash
docker build -t primer-blast3-demo .
docker run --rm -v ${PWD}:/data primer-blast3-demo python /app/src/primer_design.py
