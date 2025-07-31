FROM debian:bookworm-slim

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    make \
    wget \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev \
    pigz \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Download and install samtools
RUN mkdir -p /tmp/samtools \
    && wget -O /tmp/samtools/samtools.tar.bz2 \
    https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    && tar -xjf /tmp/samtools/samtools.tar.bz2 -C /tmp/samtools --strip-components=1 --no-same-owner \
    && cd /tmp/samtools \
    && ./configure --prefix=/usr/local \
    && make \
    && make test \
    && make install \
    && rm -rf /tmp/samtools

# Verify installation
RUN command -v samtools --help &>/dev/null || exit 1

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]