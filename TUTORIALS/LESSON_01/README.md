# Basic Python and Atomic Simulation Environment

## :one: Purposes
- Know basic knowledge about python and atomic simulation environment (ASE)
- Make Ge bulk structure.

## :two: Hows
### Basic about Python 
- Tutorial can be found here: https://www.tutorialspoint.com/python/index.htm
- Start from: Python-Basic Syntax:  https://www.tutorialspoint.com/python/python_basic_syntax.htm
- Don't have to learn much. For now, just need to know the following:
```
1. How to write a python code
2. How to run a python code
3. How to do a simple math with python code
4. How to import a module/function
```
### Basic about Atomic Simulation Environment
- Document and tutorial can be found here: https://wiki.fysik.dtu.dk/ase/
- Don't have to learn much. For now, just need to learn:
```
1. How to make Ge bulk structure
2. How to visualize a structure
```
**Hint**:
* Start with section **`Modules/The Atoms object`**.
* Check out section **`Modules/Building things`**, `bulk` function on `ase` website.
* Ge has `face-centered diamond-cubic` structure.
* For visualization, Check out section **`Modules/visualize`**.

**Note**: Many versions of ASE has been installed in `challenger`. To use the default version:
```
module load ase
```

## :three: Final products
The following items must be delivered to `LESSON_01` folder before moving to `LESSON_02`
- A python code, name `bulkGe.py`, to make a primitive Ge bulk with experimental lattice parameter.

**Hint:** See **hints** in section :two:. We are looking for a 2 lines code: the first one is for loading `bulk` function from `ase.build`; the second one is for using `bulk` function to make `diamond` structure for `Ge`.

- A figure (*.png, *.jpg, or any other picture file) showing bulk Ge. 

**Hint:** Remeber to check with Ge bulk figure found online.
- A python code, name `supercellGe.py`, to make a 3x3x3 Ge bulk with lattice parameter that is 5% lager than experimental one.

**Hint:** Must do a simple math to calculate new lattice parameter. For making 3x3x3 of Ge bulk, use `atoms * ( 3, 3, 3 )`, where `atoms` is the `1x1x1` Ge bulk. This code should have 4 lines, two of which are the same as those in `bulkGe.py`, the other two are for calculating new lattice parameter and for making 3x3x3 supercell.
