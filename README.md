# WuKong Codebase

Depending on what you need you may choose to build the specific project by changing the CMakeLists.txt in the Project folder

Ideally basic simulation models such as,
-DiscreteShell, FEM2D/3D, EoLRods-
should have only the basic implementation such that they can be inherited whenever needed. 

### Set up WuKong using Docker and VSCode

Enable Docker to connect to host display to spawn GUI windows.
> $ xhost +

Navigate to repository home directory.
> $ cd Wukong2024

Install Docker, Dev Containers extension in VSCode

Open VSCode type `control + p`, then type `>Reopen in Container`. This option will show up in the >< tab in the bottom left corner of vscode

The above command will open a dev container using the docker file we provided.

In this container, install C++ Extension Pack for dev container.

There you go!


### Set up WuKong using just Docker 

Build Docker image from Dockerfile in the .devcontainer directory (.).
> $ docker build -t wukong_docker .

Open the Docker container. \
`-v ./:/{directory-name-in-container}` maps current working directory (.) to /{directory-name-in-container} in the Docker container. \
`--network=host -e DISPLAY=$DISPLAY --privileged` enables use of host display to spawn GUI windows.
> $ docker run -v ./:/{directory-name-in-container} -it --network=host -e DISPLAY=$DISPLAY --privileged --rm wukong_docker bash

Build the code.
> $ cd {directory-name-in-container} \
> $ ./build.sh


### Docker Image Building Time
Building this docker image can take a while, for downloading MKL libraries and compiling SuiteSparse from the source code (just to remove a single print). 
In case you have a powerful workstation, considering changing all the `make -j8` to `make -j128`.

### Coding Convention

### Naming Convention

    NAMESPACE_EXAMPLE
    ClassExample
    functionExample
    variable_example
    TypenameExample

### To-Dos
- [x] dynamics simulation (see Discrete Shell)
- [ ] Set up application base class
- [ ] Use curvature binormal for compute dihedral angles
- [ ] FEM2D
- [ ] IPC contact examples

### More Info
Authors: [Yue Li](https://liyuesolo.github.io/) ...


My dear CRL members please add your contribution! :heart:
