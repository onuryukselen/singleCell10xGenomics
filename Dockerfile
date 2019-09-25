FROM nfcore/base
LABEL author="onur.yukselen@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/singlecell10xgenomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# Install dolphin-tools
RUN git clone https://github.com/dolphinnext/tools /usr/local/bin/dolphin-tools
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-singlecell-1.0/bin:/usr/local/bin/dolphin-tools/:$PATH