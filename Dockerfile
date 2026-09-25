FROM python:3.12-slim

# minimap2 is the only external dependency; samtools ships inside pysam.
RUN apt-get update \
    && apt-get install --no-install-recommends -y minimap2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .

# Mount the working directory here; relative paths in -q/-r/-o resolve to it.
WORKDIR /data
ENTRYPOINT ["concatmap"]
