# EOEinBLG
[![DOI](https://zenodo.org/badge/501569913.svg)](https://zenodo.org/badge/latestdoi/501569913)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fernandopenaranda.github.io/EOEinBLG.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fernandopenaranda.github.io/EOEinBLG.jl/dev)

For the sake of reproducibility, in this repository we provide all the scripts and data for the simulations in: "Supercurrent mediated by helical edge modes in bilayer graphene"

The critical currents of the dense tight-binding model in Figs. 4(d,e) are computed in `simulations/computesupercurrents.jl`. The generated data files for panels d and e are stored in: `simulations/data/fig1vacuumedges/jcmat.csv` and `simulations/data/fig2vacuumedges/jcmat.csv`, respectively.

The script of the 4 corners model in Fig. 4(a,b) and Fig. S9 can be found in `Even odd.nb` in the `Mathematica` folder which also contains the plotting code used to build all theory plots in the main text and supplementary.
