# ipython 


## :one: Purposes
- Use `ipython` to perform simple data analysis and plots.

## :two: Hows

### Software needed
- ipython
- pandas
- numpy
- matplotlib

Before starting, load the modules:
```
module load pandas
```

### Exercise

`data.csv` is a dataset for a hyperthetical experiments on the resistivity (R) of a material as function of doping (N, P, O). There is also a column of data with label `T`. Experimentalist forgot what it is about. We will find the relationship between `R` and N, P, O doping concentrations and the unkown `T`.

Get the data set:
```
wget https://github.com/zoowe/dlePy/blob/master/TUTORIALS/IPYTHON/data.csv
```
If it does not work:
```
cp -a /shared/apps/dlePy/0.10/dlePy/TUTORIALS/IPYTHON/data.csv .
```
#### Step 1: Launch ipython

```
$ ipython
```

If you have not used `ipython`, do few simple maths to get familiar with it.

#### Step 2: Load nessesary modules for this exercise

```python
import pandas as pd                                # For data related task
from pandas.tools.plotting import scatter_matrix   # For plotting correlations
import matplotlib.pylot as plt                     # For graphing
%matplotlib                                        # Using matplotlib backend: TkAgg, iteractive mode for plotting
```

#### Step 3: Load data

Before reading data, use `less` or `vi` to inspect the `data.csv`. In `ipython`, use `!less` or `!vi`
```
!less data.csv
!vi data.csv
```

Read data
```python
data = pd.read_csv( 'data.csv', sep='\t' )
```

Run the following to see data
```python
data.columns
data.values.shape
data.head( 5 )
```

#### Step 4: Plot data

```python
scatter_matrix( data, alpha = 0.8 )
```

This plot shows the correlation, if any, between all variables. You will see that `T` is linearly depend on `N`. Thus, it is meaningless to have both. 

#### Step 5: Drop `T`

```python
data = data.drop( 'T', 1 )
```

Plot again to verify
```python
scatter_matrix( data, alpha = 0.8 )
```

There is no obvious correlation between `R` and doping concentrations. The next best guest is the cross terms, i.e. `N*P`, `N*O`, `O*P`. We do not know which one is good, thus we calculate all three and add to the dataset.

#### Step 6: Add cross terms to dataset

```python
data.insert( len( data.columns ), 'N*P', data[ 'N' ] * data[ 'P' ] )
data.insert( len( data.columns ), 'N*O', data[ 'N' ] * data[ 'O' ] )
data.insert( len( data.columns ), 'P*O', data[ 'P' ] * data[ 'O' ] )
```
If there are a lot more cross terms, we can use loops to do this task.

Plot them
```python
scatter_matrix( data, alpha = 0.8 )
```

The plot shows that `R` correlates with `N*P`.

#### Step 7: Fit `R` as a function of `N*P`

We first plot `R` and `N*P` in a nicer plot
```python
plt.figure(4, figsize = ( 3.5, 3.5 ) )    # Number 4 indicates the plot 4. 
                                          # The previous 3 should have been numbered as 1,2,3
plt.scatter( data[ 'N*P' ], data[ 'R' ], color = 'r', marker='o' , label  = 'data' )
```

The plot shows clearly a quadratic dependence. Now, do fitting.
```python
import numpy as np
# Fit to 2nd order polynomial
z2   = np.polyfit( data[ 'N*P' ], data[ 'R' ], 2 )
# Make it become a function
f2   = np.poly1d ( z2 )
# Calculate mean error
err2 = np.abs    ( f2( data[ 'N*P' ] ) - data[ 'R' ] ).mean( )

# Do the same for 3rd oder polynomial
z3   = np.polyfit( data[ 'N*P' ], data[ 'R' ], 2 )
f3   = np.poly1d ( z3 )
err3 = np.abs    ( f3( data[ 'N*P' ] ) - data[ 'R' ] ).mean( )
```

The two fits give similar errors. Thus, we should expect the their curves are very close. Plot them.
```
x = np,arange( np.min( data[ 'N*P' ] ), np.max( data[ 'N*P' ] ), 0.01 )
plt.plot( x, f2( x ), 'b-', label= '2nd order polyfit' )
plt.plot( x, f3( x ), 'g--', label= '3nd order polyfit' )
```

The plot shows identical curves.

#### Step 8: Find optimal value of doping concentrations.

```python
f2d  = np.polyder( f2, 1 )   # Take first derivative of f2
NP_2 = np.roots  ( f2d   )   # Optimal value of N*P from 2nd order fit

f3d  = np.polyder( f3, 1 )   # Take first derivative of f3
NP_3 = np.roots  ( f3d   )   # Optimal value of N*P from 3nd order fit.
                             # It should give 2 values, chose the appropriate one.
```

#### Step 9: Make figure nicer and save

```python
plt.figure ( 4 )               # working on figure 4
plt.legend ( loc = 1 )
plt.xlabel ( r'N$\times$P' )
plt.ylabel ( 'R' )
plt.savefig( 'fit.png', dpi = 600, bbox_inches = 'tight')

plt.figure ( 3 )               # working on figure 3
plt.savefig( 'data.png', dpi = 600, bbox_inches = 'tight')
```

#### Step 10: Save your working log

```python
%logstart -o -r -t working.log
```

