FROM continuumio/miniconda3
# Path to yaml file containing conda dependencies
ARG CONDA_YML_PATH=.
# Install the conda environment
COPY $CONDA_YML_PATH/scaleMethylPyQc.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && conda clean -a --yes
# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH="/opt/conda/envs/scaleMethylPyQc/bin:${PATH}"
