FROM nfcore/base:2.1
# Path to yaml file containing conda dependencies
ARG CONDA_YML_PATH=.
# Install the conda environment
COPY $CONDA_YML_PATH/scaleMethyl.conda.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && conda clean -a --yes

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/scaleMethylTssEnrichR/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name scaleMethylTssEnrichR > scaleMethylTssEnrichR.yml
