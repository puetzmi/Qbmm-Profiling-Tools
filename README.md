# Profiling Tools for the Basic Inversion Procedure of Quadrature-Based Moment Methods

The aim of this project is to provide tools to study the differences with respect to performance and accuracy using different implementations of the quadrature method of moments and derived methods.


#### Prerequisites

The tools/packages required for various purposes are listed below (the version numbers in parentheses indicate the recommended versions, which have been tested, though older versions may work as well).

 * **Building/running libraries and applications**
   * A C/C++ compiler. Compilers that have been tested are
     * the GNU C++ Compiler (11.2.0)
     * the Intel&reg; oneAPI DPC++/C++ Compiler (2022.2.0)
   * The Intel&reg; oneAPI Math Kernel Library (2022.2.0)
   * CMake (3.22.1)
 * **Building the documentation**
   * Doxygen (1.9.1)
 * **Data generation with Python scripts**
   * Python3 (3.10.6)
   * NumPy (1.21.5)
   * QuadMomPy (0.9.7)
 * **Testing**
   * CTest (3.22.1) for simple tests
   * Valgrind (3.18.1) for memory check

A Dockerfile to build a suitable image for data generation and building/running is provided in the project root directory.


#### Building The Project

For the following instructions, it is presumed that the binary files are not built in the source directory (which is strongly recommended) and that the build-directory name contains the compiler and build type such that the directory structure looks like this:
```
PROJECT_ROOT_DIR
│
├── applications
│   ├── ...
│
├── build
│   ├── gcc-release
│   ├── ...
│
├── CMakeLists.txt
│
├── doc
│   └── Doxyfile.in
│
├── ...
```
As indicated by the name of the build-directory, it is assumed that the default system compiler is the GNU C++ compiler and that the build type is a release build. 

For CMake to run properly the configuration files for the Intel&reg; oneAPI Math Kernel Library must be in a location where they can be found by CMake. In order to set the required variables source the appropriate script in the `scripts` directory. For example, in case of the GNU compiler
```bash
. scripts/load-gnu.sh
```

In order to build the project with the default parameters, change the current working directory to the build-directory and run CMake without any optional arguments.
```bash
> cd build/gcc-release
> cmake ../..
```
If the the Intel&reg; oneAPI Math Kernel Library is located in a non-default directory this step will fail unless one of the variables `INTEL_ONEAPI_ROOT` or `MKLROOT` is set to the right location.

If CMake finishes without errors the libraries and executables can be built by running
```bash
> make
```

For other configurations, suitable arguments must be passed to CMake. For example, to make a debug build of the project using the Intel&reg; oneAPI compiler as well as the Doxygen documentation, run the following command (assuming that the Intel&reg; compiler is installed as a non-default compiler and the commands `icx` and `icpx` are known):
```bash
> cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DBUILD_DOC ../..
> make
```
Alternatively, the CMake Curses GUI can be used to change the build configuration running 
```bash
> ccmake .
``` 
inside the build-directory. The relevant variables that can be passed to CMake or adjusted using the CMake Curses GUI are the following.
 * **CMAKE_BUILD_TYPE** : Basic type of the build, may be one of *Release*, *Debug*, *RelWithDebInfo*, *MinSizeRelease*.
 * **BUILD_DOC** : Build Doxygen documentation.
 * **WITH_TESTING** : Build tests and enable testing using CTest.
 * **USE_AVX2** : Use Advanced Vector Extensions (AVX) 2, disable if not supported by the processor.
 * **MALLOC_ALIGN_SIZE** : Specifies boundaries of aligned memory in bits.
 * **GDB_DEBUG_LEVEL** : The level of debug information for the GNU debugger. Levels other than 0 are only selectable if `CMAKE_BUILD_TYPE=Debug`.

When the project is built successfully the build-directory contains the subdirectories `lib` with shared libraries and `bin` with executable applications.


#### Building and using the documentation

If a supported version of Doxygen is installed and *BUILD_DOC* is enabled, an additional subdirectory `doc` will be created. The HTML documentation can then be opened by any suitable web browser, e.g. Mozilla Firefox:
```bash
> firefox doc/html/index.html
```

#### Testing

Provided that the required packages are installed, simple test applications can be run with the command
```bash
> make test
```
A complete check for memory leaks/errors is carried out by
```bash
> make test_memcheck
```

#### Data generation

Python scripts to generate input data for QMOM performance profiling are located in `scripts` in the source directory. For more information on the usage, see the comments in the Python files and the Doxygen documentation.


#### License

&copy; 2022 Michele Pütz

Open sourced under MIT license, the terms of which can be read here — [MIT License](http://opensource.org/licenses/MIT).