FROM quay.io/jupyter/r-notebook:latest

USER root
# Install required system libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libgdal-dev \
    libudunits2-dev \
    libsqlite3-dev \
    libpng-dev \
    libwebp-dev \
    libtiff-dev \
    libjpeg-dev \
    liblzma-dev \
    libz-dev \
    libcairo2-dev \
    libx11-dev

# Install some system utilities
RUN apt-get install -y --no-install-recommends \
    cmake
    
# handle some pacakges which doesn't look at system libraries
# and also umap-learn
RUN /opt/conda/bin/conda install -y -c conda-forge \
    zlib \
    cairo \
    libxml2-devel \
    xorg-xorgproto \
    freetype \
    libudunits2 \
    udunits2 \
    libgdal \
    armadillo \
    proj \
    libsqlite \
    umap-learn

# Install R packages for single-cell analysis
RUN /opt/conda/bin/conda install -y -c conda-forge \
    r-tidyverse \
    r-devtools \
    r-remotes \
    r-biocmanager \
    r-harmony \
    r-seurat \
    r-sf

# Install Bioconductor packages
RUN R -e 'options(warn=2);BiocManager::install(c("S4Vectors", "IRanges", "BiocGenerics"), Ncpus = parallel::detectCores())'
RUN R -e 'options(warn=2);BiocManager::install(c("XVector"), Ncpus = parallel::detectCores())'
RUN R -e 'options(warn=2);BiocManager::install(c("DelayedArray", "DelayedMatrixStats", "HDF5Array"), Ncpus = parallel::detectCores())'
RUN R -e 'options(warn=2);BiocManager::install(c("limma", "lme4", "batchelor", "ggrastr"), Ncpus = parallel::detectCores())'
RUN R -e 'options(warn=2);BiocManager::install(c("scran", "scater", "DropletUtils", "EnhancedVolcano", "DESeq2"), Ncpus = parallel::detectCores())'

# Install additional R packages from GitHub
RUN R -e 'options(warn=2);remotes::install_github("bnprks/BPCells/r")'
RUN R -e 'devtools::install_github("immunogenomics/presto")'
RUN R -e 'devtools::install_github("cole-trapnell-lab/monocle3")'
RUN R -e 'devtools::install_github("jokergoo/ComplexHeatmap")'
RUN R -e 'devtools::install_github("jinworks/CellChat")'
