# WuKong Codebase

Depending on what you need you may choose to build the specific project by changing the CMakeLists.txt in the Project folder

Ideally basic simulation models such as,
-DiscreteShell, FEM2D/3D, EoLRods-
should have only the basic implementation such that they can be inherited whenever needed. 

## Docker

Enable docker pull/push from command line by signing into our shared account "wukongsim".
> $ docker login --username=wukongsim

> $ Password: an_awesome_docker_account

Download the docker image. Change tag "linux" if needed.
> $ docker pull wukongsim/wukong:linux

If you wish to build the docker image from scratch from the Dockerfile (not recommended) or rebuild after modifying the dockerfile, run the following in the command line from the directory Wukong2024/.devcontainer.

LOGAN: Can we download cache layers using `docker pull` so that we can modify Dockerfile without building from scratch?
> $ docker build -t wukongsim/wukong:linux .

If finished modifying the Dockerfile, push to dockerhub.
> $ docker push wukongsim/wukong:linux

Enable Docker to connect to host display to spawn GUI windows.
> $ xhost +

### Run Docker in VSCode

Open the repository folder in VSCode. Install Docker and Dev Containers extensions in VSCode.

In VSCode, type `control + p`, then type `>Reopen in Container` (with the '>'). This option will show up in the >< tab in the bottom left corner of vscode.

The above command will open a dev container using the docker image we provided.

There you go! :100:

### Run Docker from Command Line

Navigate to repository home directory.
> $ cd Wukong2024

Open the Docker container. \
`-v ./:/{directory-name-in-container}` maps current working directory (.) to /{directory-name-in-container} in the Docker container. \
`--network=host -e DISPLAY=$DISPLAY --privileged` enables use of host display to spawn GUI windows.
> $ docker run -v ./:/{directory-name-in-container} -it --network=host -e DISPLAY=$DISPLAY --privileged --rm wukongsim/wukong:linux bash

Build the code.
> $ cd {directory-name-in-container}

> $ ./build.sh


### Docker Image Building Time
Building this docker image can take a while, for downloading MKL libraries and compiling SuiteSparse from the source code (just to remove a single print). 
In case you have a powerful workstation, considering changing all the `make -j8` to `make -j128`.

### Projects Tested Compiling
- Discrete Shell [x] Linux [] MacOs
- FEM3D  [x] Linux [] MacOs
- EoLRods  [x] Linux [] MacOs
- Isohedral Tiling  [x] Linux [] MacOs

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
Authors: [Yue Li](https://liyuesolo.github.io/), Logan Numerow, Peiyuan Xie ...


My dear CRL members please add your names if you contribute to the development of this codebase! :heart:
