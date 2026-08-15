FROM mambaorg/micromamba:1.5-jammy

# Install all bioinformatics tools from both conda environments
# (align_map.yml + r_env.yml from the original Nextflow pipeline)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.9 \
    bioconda::bwa=0.7.17 \
    bioconda::samtools=1.15 \
    bioconda::deeptools=3.5.3 \
    bioconda::fastp=0.23.4 \
    conda-forge::r-base=4.1 \
    conda-forge::r-optparse \
    conda-forge::r-ape \
    conda-forge::r-stringr \
    conda-forge::r-ggplot2 \
    conda-forge::r-scales \
    conda-forge::r-rcolorbrewer \
    conda-forge::r-reshape2 \
    conda-forge::r-tibble \
    conda-forge::r-rsqlite \
    conda-forge::r-plyr \
    bioconda::bioconductor-rsubread=2.8.1 \
    bioconda::bioconductor-edger \
    bioconda::bioconductor-deseq2 \
    bioconda::bioconductor-enhancedvolcano \
    bioconda::r-ggbiplot \
    && micromamba clean --all --yes

USER root
RUN apt-get update && apt-get install -y curl wget && rm -rf /var/lib/apt/lists/*

ENV PYTHONUNBUFFERED=1
ENV PATH="/opt/conda/bin:${PATH}"