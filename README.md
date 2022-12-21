# phongdefo

### 1. About

This is a C++ Houdini SOP implementing [phong deformation](https://graphics.pixar.com/library/PhongDefo/paper.pdf).

### 2. Installation

Clone the repository, [source the Houdinin environment](https://www.sidefx.com/docs/hdk/_h_d_k__intro__compiling.html), and then run 'make install' from src/. The SOP depends on two header-only libraries, [eigen](https://gitlab.com/libeigen/eigen) for linear algebra operations and [nanoflann](https://github.com/jlblancoc/nanoflann) for its kdtree implementation. This has only been tested on a 64-bit linux machine, running the free version of houdini v19.5. Tinkering probably needed for compiling and installing on different systems.
