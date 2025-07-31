.PHONY: build start shell stop clean

build:
    UID=$(shell id -u) GID=$(shell id -g) docker-compose build

start:
    UID=$(shell id -u) GID=$(shell id -g) docker-compose up -d

shell:
    docker-compose exec dev bash

stop:
    docker-compose down

clean:
    docker-compose down --rmi all --volumes
