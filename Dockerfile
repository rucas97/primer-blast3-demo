FROM python:3.12-slim

RUN apt-get update && apt-get install -y \
    ncbi-blast+ build-essential && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir biopython==1.85 primer3-py==2.2.0 matplotlib seaborn

CMD ["python", "primer_design.py"]
