FROM debian:bookworm-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    ca-certificates \
    unzip \
    zlib1g-dev \
    perl \
    default-jre \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# install fastqc and symlink to path.
RUN curl -LsSf https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    -o fastqc.zip \
    && unzip fastqc.zip -d /usr/local/src/ \
    && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc \
    && rm fastqc.zip

# make our venv.
RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"
# Now install multiqc using uv :)
RUN uv pip install multiqc

# Verify installation
RUN for cmd in fastqc multiqc; do \
    command -v cmd &>/dev/null || exit 1; \
    done

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]
