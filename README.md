# MidlineIdentifier -- Spatially quantify gene expression
*MidlineIdentifier* allows comparison of the spatial dynamics and polarization of gene profiles across budoids of varying sizes. A more detailed description of the method can be found in [our paper](https://doi.org/).

<!---[Schematics](./figs/Schematics.png)-->


## Installation

Download
```console
git clone https://github.com/AztekinLab/MidlineIdentifier.git
```

Create and activate a fresh conda environment:
```console
cd MidlineIdentifier
conda create -n MLI python=3.9
conda activate MLI
```

Finally, install `MidlineIdentifier`:
```console
pip install .
```


## Usage
Please see the [documentation](https://midlineidentifier.readthedocs.io/en/latest/) for full explanation.

### Quick start
```console
python example.py --dir ./Budoid_1A --sample Budoid_1A
```
