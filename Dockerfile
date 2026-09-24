# PTESFinder v2 (PFv2).
#
# One image carrying PFv2 and the three tools it shells out to, so a run needs nothing from
# the host but the reads and a reference.
#
#   docker build -t conidiobolus/pfv2:2.2.0 .
#   docker run --rm conidiobolus/pfv2:2.2.0 --help
#   docker run --rm conidiobolus/pfv2:2.2.0 selftest
#
# The entrypoint is the `ptesfinder` CLI, also installed as `pfv2`.
#
# Builds for linux/amd64 and linux/arm64. The tools come from bioconda, which publishes
# aarch64 builds for all of them, so neither architecture needs anything compiled here.
# Aligner indexes are portable between the two: the STAR, Bowtie2, BWA and HISAT2 formats
# are all little-endian 64-bit with no architecture-specific packing, so an index built on
# x86_64 loads unchanged under the aarch64 binaries.

# ---------------------------------------------------------------- build
FROM eclipse-temurin:17-jdk AS build

WORKDIR /src
COPY manifest.mf setup.sh ./
COPY lib ./lib
COPY src ./src
COPY test/bio ./test/bio

# setup.sh runs the test suite and refuses to package the jar when it fails, so a broken
# build cannot become an image.
RUN bash setup.sh

# ---------------------------------------------------------------- runtime
FROM mambaorg/micromamba:1.5.8

ARG PFV2_VERSION=2.2.0

LABEL org.opencontainers.image.title="pfv2" \
      org.opencontainers.image.description="PTESFinder v2: annotation-free backsplice junction identification from RNA-seq" \
      org.opencontainers.image.version="${PFV2_VERSION}" \
      org.opencontainers.image.licenses="MIT" \
      pfv2.tool-version="${PFV2_VERSION}"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends tini procps \
 && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PFV2_HOME=/opt/pfv2 \
    PATH="/opt/conda/bin:${PATH}" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

USER root
COPY --from=build /src/PFv2.jar ${PFV2_HOME}/PFv2.jar
COPY --from=build /src/lib ${PFV2_HOME}/lib
COPY PFv2.sh ${PFV2_HOME}/PFv2.sh
COPY scripts ${PFV2_HOME}/scripts
COPY test/smoke ${PFV2_HOME}/test/smoke
COPY README.md CHANGELOG.md ${PFV2_HOME}/
COPY bin/ptesfinder /usr/local/bin/ptesfinder

RUN chmod +x /usr/local/bin/ptesfinder ${PFV2_HOME}/PFv2.sh \
 && ln -s ptesfinder /usr/local/bin/pfv2 \
 && ln -s ${PFV2_HOME}/PFv2.sh /usr/local/bin/PFv2.sh \
 && mkdir -p /work && chown $MAMBA_USER:$MAMBA_USER /work

# Non-root, so a bind-mounted output directory does not end up owned by root. Mount reads
# and references read-only and /work read-write.
USER $MAMBA_USER
WORKDIR /work

# Exercises stages 2, 3 and 5 against a synthetic dataset at build time, so a broken image
# fails the build rather than an analysis.
RUN ptesfinder selftest

ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/ptesfinder"]
CMD ["--help"]
