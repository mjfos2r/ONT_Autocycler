FROM ubuntu:22.04
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    zlib1g-dev \
    git \
    make \
    shellcheck \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

ARG USER_UID=1000
ARG USER_GID=1000
ARG USERNAME=devuser

RUN groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m ${USERNAME} -s /bin/bash \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

WORKDIR /workspace
RUN chown ${USER_UID}:${USER_GID} -R /workspace \
    && chown ${USER_UID}:${USER_GID} -R /opt

COPY . /opt/

RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN source /opt/venv/bin/activate \
    && uv pip install -r dev_deps.txt

USER ${USERNAME}

CMD [ "/bin/bash" ]







