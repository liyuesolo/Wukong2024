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
Enable Docker to connect to host display to spawn GUI windows.
> $ xhost +

Navigate to repository home directory.
> $ cd Wukong2024

Build Docker image from Dockerfile in current working directory (.).
> $ docker build -t wukong_docker .

Open the Docker container. \
`-v ./:/{directory-name-in-container}` maps current working directory (.) to /{directory-name-in-container} in the Docker container. \
`--network=host -e DISPLAY=$DISPLAY --privileged` enables use of host display to spawn GUI windows.
> $ docker run -v ./:/{directory-name-in-container} -it --network=host -e DISPLAY=$DISPLAY --privileged --rm wukong_docker bash

Build the code.
> $ cd {directory-name-in-container} \
> $ ./build.sh

### More Info
Author: [Yue Li](https://liyuesolo.github.io/)

My dear CRL members please add your contribution! :heart:
