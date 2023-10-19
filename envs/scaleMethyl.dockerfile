FROM continuumio/miniconda3
# Path to yaml file containing conda dependencies
ARG CONDA_YML_PATH=.
# Install the conda environment
COPY $CONDA_YML_PATH/scaleMethyl.conda.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && conda clean -a --yes
# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH="/opt/conda/envs/scaleMethyl/bin:${PATH}"
# bc_parser, etc.
COPY $CONDA_YML_PATH/download-scale-tools.sh /
RUN /download-scale-tools.sh /tools
ENV PATH="/tools:${PATH}"