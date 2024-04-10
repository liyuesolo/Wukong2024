# WuKong Codebase

Depending on what you need you may choose to build the specific project by changing the CMakeLists.txt in the Project folder

Ideally basic simulation models such as,
-DiscreteShell, FEM2D/3D, EoLRods-
should have only the basic implementation such that they can be inherited whenever needed. 

### To-Dos
- [x] dynamics simulation (see Discrete Shell)
- [ ] Set up application base class
- [ ] Use curvature binormal for compute dihedral angles
- [ ] FEM2D
- [ ] IPC contact examples

### ThirdParty Library Dependencies
Polyscope and Libigl are from git fetch content
SuiteSparse has to be build locally

### Docker
docker build -t wukong_docker .

docker run -v local_foler/:/place_in_the_container -it --rm wukong_docker bash

Author: [Yue Li](https://liyuesolo.github.io/)

My dear CRL members please add your contribution!
