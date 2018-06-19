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
4. How to import a module
```
### Basic about Atomic Simulation Environment
- Document and tutorial can be found here: https://wiki.fysik.dtu.dk/ase/
- Don't have to learn much. For now, just need to learn:
```
1. How to make Ge bulk structure
2. How to visualize a structure
```
**Hint**:
* Start with section **`The Atoms object`**.
* Check out section **`Building things`**, `bulk` module on `ase` website.
* Ge has `face-centered diamond-cubic` structure.
* For visualization, Check out section **`visualize`**.

**Note**: Many versions of ASE has been installed in `challenger`. To use the default version:
```
module load ase
```

## :three: Final products
The following items must be delivered to `LESSON_01` folder before moving to `LESSON_02`
- A python code, name `bulkGe.py`, to make a primitive Ge bulk with experimental lattice parameter.

**Hint:** See **hints** in section :two:
- A figure (*.png, *.jpg, or any other picture file) showing bulk Ge. 

**Hint:** Remeber to check with Ge bulk figure found online.
- A python code, name `supercellGe.py`, to make a 3x3x3 Ge bulk with lattice parameter that is 5% lager than experimental one.

**Hint:** Must do a simple math to calculate new lattice parameter. For making 3x3x3 of Ge bulk, use `Atoms * ( 3, 3, 3 )`, where `Atoms` is the `1x1x1` Ge bulk.
