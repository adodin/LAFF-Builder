# Simple LAMMPS Force Field Builder
Written By: Amro Dodin (Geissler Group - UC Berkeley & LBL)

## Introduction

A python utility for constructing LAMMPS force field files.

Currently supports:

    - LJ Force fields
    - Easy to provide custom FF parameters via csv files
    - References for Force Field Parameters automatically bundled into file
    - Arithmetic & Geometric Mixing Rules (can be extended w/ custom rules)
    - Polarizable Drude Simulations including Thole Screening
    - Symmetrized Drude Force Fields (Ref. coming soon)

In the future, intended to be an all-in-one LAMMPS input file generator.
Similar to LEAP input generation in AMBER but in a python environment.

A small sample database (ff.csv) is provided with a few example ion force fields.

## Usage 

### Generating Force Field Files

In Python, it only takes three lines to create a lammps force field file.
For example, to create a non-polarizable NaCl force field from the default ff.csv parameters:

```python
from ForceField import ForceField
ff = ForceField(['spc', 'Na', 'Cl'])
ff.print('ff.spce-Na-Cl.lmp')
```

This creates a lammps input file called ```ff.spce-Na-Cl.lmp``` with all masses and pair coefficients.
The constructor for the ForceField object takes as a required argument a list  of Molecule and Atom names.
These are searched for by default in the provided ff.csv file. 
A list of custom csv files can also be provided with the ```db_fnames``` kwarg.

Or similarly, to create polarizable NaCO3 & NaCl (using the polarizable swm4-ndp water model) force fields we could run:

```python
from ForceField import ForceField
polNaClff = ForceField(['swm4-ndp', 'polNa', 'polCl'], polarizable=True)
polNaClff.print('ff.swm4-ndp-Na-Cl.lmp')
polNaCO3ff = ForceField(['swm4-ndp', 'polNa', 'CO3'], polarizable=True)
polNaCO3ff.print('ff.swm4-ndp-Na-CO3.lmp')
```

These files include Thole Screening as well as the ```fix drude``` lines required for Drude simulations in LAMMPS.
By default, all Drude Force Fields are symmetrized to prevent erroneous polarization due to non-Coulombic interactions.
This behavior can be modified by setting the kwarg ```symmetrized=False``` during consntruction.

### Including in LAMMPS Input File

Now that we have our LAMMPS force field file, we just need to include it where we normally specify pair coefficients.
For example, to include the polarizable NaCO3 force field we just add

```
include ff.swm4-ndp-Na-CO3.lmp
```

to our LAMMPS input script where we would normally specify masses and pair coefficients.

## Future Features

This software is a side project I'm building during my postdoc and so is unlikely to be regularly updated. 
New features will likely be added as they become useful to me.
Currently, I am planning to introduce the following features in future edits (in roughly this order of priority):

1. Molecule File Builder: Build LAMMPS .mol files from provided coordinates and topology information.
    - Ideally, this will allow convenient construction from SMILE strings or other data formats.
2. Initial Coordinate Builder: Build LAMMPS data files including initial coordinates.
    - Ideally, this will allow for conversion from xyz & pdb files to integrate with utilities like packmol
    - Alternatively, it would be nice to allow for simple random or on lattice constructions with a simple python input.
3. Input Script Generator: Build LAMMPS input files directly in python.
   - This is a bit of a stretch goal and the most involved to implement.
   - Ideally, this will allow simple lammps input files to be generated with only a few lines of python code.
   - Most likely, a few standard templates for e.g. minimization & npt/nvt simulations will be provided.
   - These can then be built on and modified with a few simple class methods.
   - Alternatively, input fixes and output settigns can be manually specified.
   - A simple syntax checker would be nice to have (e.g. checking that SHAKE fixes are specified in the right place)