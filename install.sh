brew install git-lfs
# Install Anaconda Environment
eval "$(conda shell.bash hook)"
conda create -n g00xx python==3.10.10 -y
conda activate g00xx

# Install G00x
pip install poetry
poetry install --with dev

# Install R and R packages
conda install -c conda-forge -y r-base==4.3.3 r-remotes==2.5.0 r-git2r==0.33.0
conda install -c conda-forge -y r-tidyverse==2.0.0 r-coin==1.4_3 r-hmisc==5.1-3 r-survival=3.7-0
conda install -c conda-forge -y r-binom==1.1-1.1 r-exact==3.3 r-kableextra==1.4.0 r-sessioninfo==1.2.2 r-gee==4.13-29 r-latex2exp==0.9.6 r-confintr==1.0.2
conda install -c conda-forge -y r-testthat==3.2.1.1 r-cowplot==1.1.3
Rscript -e "remotes::install_github('FredHutch/VISCfunctions', ref = 'v1.2.4', dependencies = FALSE)"
