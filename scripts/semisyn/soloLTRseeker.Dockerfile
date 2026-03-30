# Dockerfile for soloLTRseeker
# Builds a self-contained image with all soloLTRseeker dependencies:
#   BLAST+, cd-hit, EMBOSS (water), BEDTools, and the soloLTRseeker scripts.
#
# Build:
#   docker build -t soloLTRseeker:latest -f scripts/semisyn/soloLTRseeker.Dockerfile .
#   OR (after copying this file into the soloLTRseeker repo):
#   docker build -t soloLTRseeker:latest tools/soloLTRseeker/
#
# Run:
#   docker run --rm \
#     -v /absolute/path/to/work:/work -w /work \
#     soloLTRseeker:latest \
#     /opt/soloLTRseeker/modules/soloLTRseeker ann_file.gff3 genome.fa

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    cd-hit \
    emboss \
    bedtools \
    gawk \
    bash \
    && rm -rf /var/lib/apt/lists/*

COPY modules/ /opt/soloLTRseeker/modules/
RUN chmod +x /opt/soloLTRseeker/modules/soloLTRseeker \
              /opt/soloLTRseeker/modules/*.sh 2>/dev/null || true

ENV PATH="/opt/soloLTRseeker/modules:${PATH}"

WORKDIR /work

ENTRYPOINT ["/opt/soloLTRseeker/modules/soloLTRseeker"]
