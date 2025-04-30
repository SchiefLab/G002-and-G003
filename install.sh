brew install git-lfs
# Install Anaconda Environment
eval "$(conda shell.bash hook)"
conda create -n g00x python==3.10.10 -y
conda activate g00x

# Install G00x
pip install poetry
poetry install --with dev
