brew install git-lfs
# Install Anaconda Environment
eval "$(conda shell.bash hook)"
conda create -n g00xx python==3.10.10 -y
conda activate g00xx

# Install G00x
pip install poetry
poetry install --with dev
