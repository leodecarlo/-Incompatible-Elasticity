# Incompatible Elasticity — Firedrake implementation

This repository contains a Firedrake implementation of numerical experiments related to the model developed in:

> Samuel Amstutz and Nicolas van Goethem,  
> ***A second-order model of small-strain incompatible elasticity***, Mathematics and Mechanics of Solids, Volume 29, Issue 3
> <https://doi.org/10.1177/10812865231193427>.

The model describes small-strain continua in which the strain field may be incompatible. In this framework, the strain is treated as a primary geometric quantity, and the incompatibility of the strain plays a central role. This is motivated by materials with microscopic defects, such as dislocations, where incompatible deformations can arise at the macroscopic level.

The repository is related to the original FreeFEM implementation available at:

<https://github.com/samuel-amstutz/incompatibility>

The present code is not a literal line-by-line translation of the FreeFEM files. Instead, it gives a Firedrake implementation of the same family of numerical examples, using Firedrake finite element spaces, Firedrake mesh handling, and Firedrake output tools.

## Mathematical background

The paper develops a second-order small-strain model based on the strain field `E` and its incompatibility `inc(E)`. The incompatibility operator is the linearized curvature associated with the strain. In three dimensions, it is written abstractly as

```text
inc(E) = Curl Curl^T E.
```

In the two-dimensional setting used by the code, the scalar incompatibility has the form

```text
inc(E) = d_xx E_yy + d_yy E_xx - 2 d_xy E_xy.
```

The numerical formulation uses a decomposition of the total strain into compatible and incompatible parts:

```text
E = epsilon(u) + E_i
```

where:

- `u` is the displacement field,
- `epsilon(u)` is the compatible linearized strain,
- `E_i` is the incompatible part of the strain.

The energy contains both a classical elastic contribution and a second-order incompatibility contribution. In simplified notation, the model involves terms of the form

```text
A E · E + D inc(E) · inc(E),
```

where `A` is the usual elastic tensor and `D` is the incompatibility-related tensor.

The code also includes an internal scalar variable `theta`, called the compatibility modulus in the paper. This variable controls the effective tangent moduli and is updated during the alternating minimization procedure. The script outputs both mechanical fields and diagnostic quantities such as the von Mises stress squared and the yield-stress margin.

## Repository structure

```text
.
├── README.md
├── hct_traction.py
├── incomp.sh
├── .gitignore
└── meshes/
```

## Implementation overview

The main point of the Firedrake implementation is that the three numerical examples are handled by a single driver script:

```text
hct_traction.py
```

Although the file is named after the traction problem, it contains implementations of the three configurations:

```text
traction
necking
inclusion
```

These are selected with the command-line option:

```bash
--problem
```

The correspondence with the original FreeFEM repository is:

```text
Original FreeFEM file        Firedrake command
---------------------------------------------------------------
traction.edp                 python3 hct_traction.py --problem traction
necking.edp                  python3 hct_traction.py --problem necking
inhomogeneity.edp            python3 hct_traction.py --problem inclusion
```

Thus, despite the name `hct_traction.py`, the Firedrake code is not limited to the traction example. The same Python file also includes the necking and inclusion/inhomogeneity configurations.

## Main files

### `hct_traction.py`

Main Firedrake script.

This file contains the numerical implementation of the model. It defines:

- the Firedrake imports and finite element setup,
- the command-line interface,
- the mesh loading,
- the boundary labels,
- the strain and incompatibility operators,
- the elastic and incompatibility-dependent material coefficients,
- the mixed finite element formulation,
- the alternating minimization loop,
- the update of the internal variable `theta`,
- and the output of the computed fields.

The script supports the following problem choices:

```text
traction
necking
inclusion
```

The meaning of the choices is:

```text
Firedrake option        Mathematical / numerical example
--------------------------------------------------------
--problem traction      perforated plate under uniaxial traction
--problem necking       traction problem with necking geometry
--problem inclusion     plate with inclusion / inhomogeneity
```

Typical command-line usage is:

```bash
python3 hct_traction.py --problem traction
python3 hct_traction.py --problem necking
python3 hct_traction.py --problem inclusion
```

The script also accepts several options:

```text
--problem
--clscale
--n_iterates
--incomp_strain_elt
--incomp_strain_deg
--theta_deg
--theta_space
--method
--nonlinear_nmax
--scale_displacement
--print_all_args
```

For example:

```bash
python3 hct_traction.py \
    --problem traction \
    --clscale 0.04 \
    --n_iterates 20 \
    --scale_displacement 10.0 \
    -options_left 0
```

The option `-options_left 0` is passed through to PETSc/Firedrake; the script is written to tolerate this extra argument.

The incompatible strain can be discretized using either Hsieh-Clough-Tocher elements or Regge elements through:

```bash
--incomp_strain_elt HCT
--incomp_strain_elt Regge
```

The default is:

```text
--incomp_strain_elt HCT
```

The HCT formulation is the main implementation currently used by the example runs. The Regge option is present in the code, but some comments in the script indicate that parts of the Regge/DG formulation are still experimental or to be completed.

The update of `theta` can be performed with either:

```bash
--method Newton
--method bisection
```

The default is Newton.

### `incomp.sh`

Shell script collecting example runs.

The script assumes that the Firedrake virtual environment is already activated. It runs several representative cases:

- perforated plate under traction,
- necking example,
- refined necking example,
- inclusion example.

To run all examples in the script:

```bash
bash incomp.sh
```

or, after making the script executable,

```bash
chmod +x incomp.sh
./incomp.sh
```

The script currently contains runs of the form:

```bash
python3 hct_traction.py --problem "traction"  --clscale 0.04   --n_iterates 20 --scale_displacement 10.0 -options_left 0
python3 hct_traction.py --problem "necking"   --clscale 0.018  --n_iterates 50 --scale_displacement 5    -options_left 0
python3 hct_traction.py --problem "necking"   --clscale 0.0132 --n_iterates 50 --scale_displacement 5    -options_left 0
python3 hct_traction.py --problem "inclusion" --clscale 0.018  --n_iterates 10                         -options_left 0
python3 hct_traction.py --problem "inclusion" --clscale 0.1    --n_iterates 10                         -options_left 0
```

### `meshes/`

Directory containing the geometry and mesh files used by the simulations.

The code loads mesh files of the form:

```text
meshes/<meshname>_clscale<value>.msh
```

The three main mesh families are:

```text
perforated_plate_clscale*.msh
necking_plate_clscale*.msh
inclusion_plate_clscale*.msh
```

The parameter `clscale` controls the mesh size. For example, the traction example with

```bash
--problem traction --clscale 0.04
```

loads

```text
meshes/perforated_plate_clscale0.04.msh
```

The necking example with

```bash
--problem necking --clscale 0.018
```

loads

```text
meshes/necking_plate_clscale0.018.msh
```

The inclusion example with

```bash
--problem inclusion --clscale 0.1
```

loads

```text
meshes/inclusion_plate_clscale0.1.msh
```

### `.gitignore`

Git ignore file.

It excludes common temporary files, Python bytecode, editor backup files, LaTeX auxiliary files, and visualization outputs such as:

```text
*.pvd
*.vtu
```

These output files are generated by running the simulations and can be regenerated when needed.

## Relation with the original FreeFEM code

The original FreeFEM implementation is available at:

<https://github.com/samuel-amstutz/incompatibility>

The original repository contains the scripts:

```text
traction.edp
necking.edp
inhomogeneity.edp
```

In this Firedrake repository, the corresponding examples are handled by the single driver script `hct_traction.py` through the option `--problem`:

```text
FreeFEM file              Firedrake problem option
--------------------------------------------------
traction.edp              --problem traction
necking.edp               --problem necking
inhomogeneity.edp         --problem inclusion
```

Therefore, the repository can be described as a Firedrake implementation of the numerical examples from the FreeFEM repository, organized in one Python script rather than in three separate `.edp` files.

The Firedrake code uses Firedrake-specific finite element spaces, syntax, solver interfaces, mesh handling, and output routines. It should therefore be regarded as a Firedrake reimplementation of the same examples, not as a literal line-by-line translation.

## Requirements

The code requires a working Firedrake installation.

Python modules used by the main script include:

```text
firedrake
ufl
numpy
matplotlib
tabulate
argparse
warnings
time
```

The shell script `incomp.sh` assumes that the Firedrake environment has already been activated before running the examples.

A typical workflow is:

```bash
source /path/to/firedrake/bin/activate
cd /path/to/Incompatible-Elasticity
python3 hct_traction.py --problem traction --clscale 0.1 -options_left 0
```

Replace `/path/to/firedrake` by the location of your Firedrake installation.

## Output

The code writes visualization output using Firedrake's `VTKFile`. Output folders are named according to the finite element choice and the problem name, for example:

```text
HCT_traction/
HCT_necking/
HCT_inclusion/
Regge_traction/
```

The `.pvd` and `.vtu` files can be opened with ParaView.

The output fields include:

```text
incompatible strain Ei
displacement
pressure
theta
mu_A(theta)
full strain
inc of full strain
von Mises stress squared
yield stress squared minus von Mises stress squared
```

The script also saves an energy plot as:

```text
<element>_<problem>/energy.pdf
```

For example:

```text
HCT_traction/energy.pdf
```

## Example runs

### Perforated plate under uniaxial traction

```bash
python3 hct_traction.py \
    --problem traction \
    --clscale 0.04 \
    --n_iterates 20 \
    --scale_displacement 10.0 \
    -options_left 0
```

### Necking example

```bash
python3 hct_traction.py \
    --problem necking \
    --clscale 0.018 \
    --n_iterates 50 \
    --scale_displacement 5 \
    -options_left 0
```


```

### Inclusion / inhomogeneity example

```bash
python3 hct_traction.py \
    --problem inclusion \
    --clscale 0.1 \
    --n_iterates 10 \
    -options_left 0
```


## Reference

Samuel Amstutz and Nicolas van Goethem,  
**A second-order model of small-strain incompatible elasticity**.  
HAL: <https://hal.science/hal-03581050>

Original FreeFEM code:  
<https://github.com/samuel-amstutz/incompatibility>

Code attribution:  
The original forked implementation from which this code was adapted was written by  
**Francis Aznaran**, GitHub user **FAznaran**:  
<https://github.com/FAznaran>   



## Status

This is a research-code repository. The code is intended for experimentation, comparison with the FreeFEM implementation, and further development of Firedrake-based formulations for incompatible elasticity and incompatibility-driven plasticity.

The HCT-based implementation is the main current path. The Regge option is present as an alternative finite element choice, but should be treated as experimental.
