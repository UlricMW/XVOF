# XVOF
One dimensional hydro code for testing xfem method. 
This software is an experimental work trying to model spall fracture using enrichment and partition of unity 
method in a finite volume code [1, 2]. Hansbo & Hansbo enrichment is used [3].

## Code description
Code documentation is available in the local repository at xfv/doc/index.html
or at https://mariebichette.github.io/

#### Available physics

###### Bulk behavior
Material behavior includes simple models for hydrodynamics, elasticity, plasticity.
- For hydrodynamics, the equation of state is a MieGrüneisen type eos
- Linear elasticity and perfect plasticity is implemented
- **No** viscous plasticity is considered

For the equation of state, an external library C can be used in order to reduce computation time (faster operations + parallel computation). This library is computed using cmake and translated to a python module using SWIG. See below for how to use it.

###### Failure models
Two rupture models are available :

1) *imposed pressure* : once the cell reach rupture criterion, the pressure is imposed at the value specified in data file

2) *enrichment* : once the cell reach rupture criterion, kinematic and thermodynamic fields are enriched according to the Hansbo \& Hansbo idea.
Mass matrix can be lumped (Menouillard [4] or total lumping) to improve calculation time.
Penalty method is used to avoid the interpenetration of the discontinuity boundaries.

Damage can be modeled with a cohesive zone model with specific cohesive law. But no bulk damage is implemented

#### Non regression procedure
This code is versioned using GitHub. A non regression procedure is launched  by Travis at each *push* on the GitHub repository 
* A unittest collection
* Integration tests for hydrodynamics, elasticity, plasticity with and without enrichment (result should be compared to reference.hdf5)

Coverage is also tested within the non regression procedure.

Moreover, the documentation is automatically generated by Travis at each push on GitHub.

#### Future work
* Check CFL stability with respect to mass matrix.
* Implement other contact treatments.
* Implement simple damage models

#### Miscellaneous
This code is now ported in Python3.7, using type hints for local variables and method arguments.

## Installation
- Download the GitHub XVOF repository (`git clone` or download)
- Install pyenv with a python version > 3.7 (follow instructions in https://github.com/pyenv/pyenv#installation)
This will enable to modify the python of pyenv instead of the system python
- Install cmake and swig :
    ```
    apt-get install cmake
    apt-get install swig
    ```
Cmake will enable to build the C library for the eos computation and SWIG translates it as a python module
Cmake version should be > 3.12.
    
- To activate python of pyenv :
from XVOF repository type 
    ```
    pyenv local 3.7.7
    ``` 
    (python3.8 can currently create bug when installing the library C for equation of state, so it is recommended to use 3.7 version)
- To install the XVOF code and dependencies :
from XVOF repository type : 
    ```
    pip install --upgrade pip
    pip install -e .
    ```
This will also build the lib C for equation of state computation.

## Test case creation
Each case is composed of :
1) Input data, stored in the file XDATA.json.
2) A mesh file, which is a .txt file containing the coordinates of nodes.

To launch a test case :
- Create a data set "XDATA.json" and a meshfile "mesh.txt"
- From the test repository, type : 
```
    OMP_NUM_THREADS=2 python XtendedFiniteVolume <case-repository>
```
By default, the external lib C is used if it has been previously installed. 

To enforce the internal computation of the equation of state (with python module), type
```
    python XtendedFiniteVolume <case-repository> --use-internal-solver
```

## References
[1] Gorecki, M. (2019). Amélioration de la description physico-numérique de l’endommagement et de la rupture de la matière sous choc (Doctoral dissertation, École centrale de Nantes).

[2] Gorecki, M., Peillex, G., Pillon, L., & Moës, N. (2020). An enriched finite volume formulation for the simulation of ductile material failure under shock loading. Computational Mechanics, 1-22.

[3] Hansbo, A., & Hansbo, P. (2004). A finite element method for the simulation of strong and weak discontinuities in solid mechanics. Computer methods in applied mechanics and engineering, 193(33-35), 3523-3540.

[4] Menouillard, T., Rethore, J., Moes, N., Combescure, A., & Bung, H. (2008). Mass lumping strategies for X‐FEM explicit dynamics: application to crack propagation. International Journal for Numerical Methods in Engineering, 74(3), 447-474.