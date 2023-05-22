# AusSedThickRF

## Calculating the sedimentary thickness across Australia using passive seismic methods

This package contains notebooks required to estimate the sedimentary thickness across Australia using receiver functions.
It uses the [`hiperseis`](https://github.com/GeoscienceAustralia/hiperseis) package to download waveforms from select teleseismic events recorded at permanent and temporary seismic networks on the continent, and calculate the receiver functions.

Notebooks in `aussedthickrf` do all the processing.
In particular, `aussedthickrf/network_delays.ipynb` calculates the $T_{P_{sb}}$ delay time (see [Agrawal et al., 2022](https://academic.oup.com/gji/article/231/3/1850/6652500)).

All data needed or generated during processing is saved in `data`.

## Setup

`hiperseis` is still not fully pip installable because of some old dependencies so it is included as a submodule here.
This is in fact a particular fork of `hiperseis` that includes some minor modifications to allow it to work with more recent versions of python.
`hiperseis` itself also includes some submodules itself, which will all need to be installed.

The recommended python version for this package is python 3.9.
The package has not been used with other versions.
We recommend using `conda` to manage your environment, mainly because of the non-python third-party dependencies of `hiperseis`.

```bash
conda env create -f environment.yml
conda activate aussedthick
git clone --recurse-submodules https://github.com/auggiemarignier/AusSedThick.git
cd AusSedThick
cd hiperseis
pip install .
```