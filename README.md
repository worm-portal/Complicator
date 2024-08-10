# The Complicator

Code Authors: Grayson Boyer, Everett Shock
Data Compilers: GEOPIG Lab, Apar Prasad

This ia a Python package for estimating standard state thermodynamic properties and Helgeson-Kirkham-Flowers (HKF) equation of state parameters for aqueous metal complexes with monovalent ligands using published methods from [Shock and Koretsky 1995](https://doi.org/10.1016/0016-7037(95)00058-8).

## Installation

This package can be installed through PyPI with:
`pip install Complicator`

## Usage

Supply a dataframe containing names of metals, ligands, stability constants for the first through fourth association at 25 °C, and optionally, standard state entropies of association for the first through fourth association complex. An example is assigned to the variable `df_input` in the example below.

```
import pandas as pd
from Complicator import complicate

df_input = pd.DataFrame({
        "Metal":["Ag+", "Al+3"],
        "Ligand":["OH-", "OH-"],
        "BETA_1":[2, 9.03], # stability constant for the first association complex
        "BETA_2":[3.97, 17.6], # second stability constant
        "BETA_3":[float('NaN'), 26.4], # third stability constant
        "BETA_4":[float('NaN'), 33.8], # fourth stability constant
        "Sass_1":[float('NaN'), 36.5], # entropy of association for the first complex (cal/mol/K)
        })

df_out, _ = complicate(df=df_input)

df_out
```

The output is a dataframe containing estimated thermodynamic properties, parameters for the revised Helgeson Kirkham Flowers (HKF) equation of state, and more. [The format of the output is explained in more detail here](https://worm-portal.asu.edu/docs/database/). This Water-Organic-Rock-Microbe (WORM) database format is designed to be compatible with the free and open source online geochemical modeling platform [WORM Portal](https://worm-portal.asu.edu/). Estimated thermodynamic properties of complexes from the Complicator can be used in conjunction with the rest of the WORM database to:
- calculate properties of reactions and create activity or predominance diagrams.
- expand the number of aqueous complexes available in geochemical speciation (equilibration) and mass transfer calculations.