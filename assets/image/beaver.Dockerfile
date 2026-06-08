# syntax=docker/dockerfile:1.7
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG BEAVER_REPO=https://github.com/alejandrogzi/beaver.git
ARG BEAVER_REF=master

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        autoconf \
        automake \
        ca-certificates \
        g++ \
        gcc \
        git \
        libboost-all-dev \
        libhts-dev \
        make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone --depth 1 --branch "${BEAVER_REF}" "${BEAVER_REPO}" beaver-src

WORKDIR /opt/beaver-src
RUN aclocal \
    && autoconf \
    && autoheader \
    && automake -a \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && install -m 0755 src/beaver /usr/local/bin/beaver \
    && beaver \
    && cd /opt \
    && rm -rf /opt/beaver-src

WORKDIR /data
