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

Author: [Yue Li](https://liyuesolo.github.io/)

If WuKong contributes to an academic publication, you can cite it as:

```bib
@misc{wukong,
  title = {WuKong},
  author = {Yue Li},
  note = {https://github.com/liyuesolo/Wukong2024},
  year = {2024}
}
```