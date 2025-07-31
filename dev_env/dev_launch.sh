#dontuse#!/usr/bin/env bash
# dev_launch.sh

# here's a handy helper script to auto launch the dev container.
export UID=$(id -u)
export GID=$(id -g)

docker-compose up -d
docker-compose exec dev bash