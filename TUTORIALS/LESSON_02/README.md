# Basic Graph with matplotlib
Graph is very important. It is the most effective way to show idea and data. Knowing how to graph efficently will help imrpoving research performance. 

Why matplotlib? Because it will make graphing very easy. 

More infor: https://matplotlib.org/

## :one: Purposes
Using python for the following:
- Saving numerical data into a text file.
- Reading numerical data from a text file.
- Creating and doing simple math with an array.  
- Making line and scatter plots using matplotlib

## :two: Hows
### :large_blue_diamond: Saving numerical data into a text file 
- Google `write data to file python` to find good information/tutorials for learning.
- Below is a simple example:

To open a file, name `square.dat`, for writing:
```python
f = open( 'square.dat', 'w' )
```

To write two-column data, `x` (from 0 to 100) and `y` ( `y = x^2` ),  to that file:
```python
for x in range( 10 ):
    y = x**2
    print >>f, "%5i %7i" %( x, y )
```

To close the file (always close file after finshing):
```python
f.close( )
```
### :large_blue_diamond: Reading numerical data from a text file
- Learn about `numpy.loadtxt` at https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.loadtxt.html and other pages (use Google search)
- Below is a simple example of reading data from the `square.dat` to two array `u` and `v`:
```python
import numpy as np
data = np.loadtxt( 'square.dat' )
u = data[ :, 0 ]
v = data[ :, 1 ]
``` 
We can import only `loadtxt` module. However, for creating a good habit, always `import numpy as np`.

### :large_blue_diamond: Creating and doing simple math with an array
- For passing this `LESSON`, do focus on 1D `numpy.array`. Information about `numpy.array` can be found at: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.array.html
- This is an example of creating an 1D array `x`, the calculate `square` of all elements and save them to new array `y`.
```python
import numpy as np
x = np.array( [ 0, 1, 2, 3, 4, 5 ] )
y = x**2 
```
### :large_blue_diamond: Making line and scatter plots using matplotlib
- Tutorial for matplotlib: https://matplotlib.org/gallery/index.html. Focus on `Pyplot` and `Subplots, axes and figures` sections.
- Example: plot `v` vs `u` (see above for `u` and `v`) and save the plot as `plot.png`.
```python
import matplotlib.pyplot as plt
plt.plot( u, v )
plt.scatter( u, v )
plt.savefig( 'plot.png', dpi = 300, bbox_inches = 'tight' )
```

### :large_blue_diamond: More about `numpy`
Quote from http://www.numpy.org/
```
NumPy is the fundamental package for scientific computing with Python. 
It contains among other things:

- a powerful N-dimensional array object
- sophisticated (broadcasting) functions
- tools for integrating C/C++ and Fortran code
- useful linear algebra, Fourier transform, and random number capabilities

Besides its obvious scientific uses, NumPy can also be used as an efficient 
multi-dimensional container of generic data. Arbitrary data-types can be defined. 
This allows NumPy to seamlessly and speedily integrate with a wide variety of databases.
```
`numpy` tutorial: https://docs.scipy.org/doc/numpy/user/quickstart.html 

## :three: Final products
The following items must be delivered to `LESSON_02` folder before moving to `LESSON_03`

:heavy_check_mark: **A code, named `plot.py`, that does the following:**
- Create an `x1` array whose elements are intergers from -100 to 100. Hint: Using `np.array` and `range(m, n)`
- Calculate `y1 = sin(x1/10)`. Hint: use `np.sin()` 
- Write `x1` and `y1` to a two-colum datafile named `data.dat`. 
- Read data from `data.dat` to `x2` and `y2` arrays.
- Calculate `y3 = cos(x2/10) + 0.5`. Hint: use `np.cos()` 
- Plot `y1` vs `x1`, `y2 + 0.25` vs `x2`, and `y3` vs `x2` in the same plot.
- Save the plot as a figure named `plot.png`.

:heavy_check_mark: **The figure named `plot.png` obtained from the `plot.py` code.**
