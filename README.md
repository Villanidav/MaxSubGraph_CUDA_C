# Maximum Common Substructure
Repository for the project on the Maximum Common Substructure problem.
The project is actually organized as a [Cmake](https://cmake.org/documentation/) project.

## Prerequisites
Other than GCC and Python, the project requires you to have installed:
- [RDKit](https://www.rdkit.org/docs/Install.html)

## Compiling
To compile it
```
$ cmake -S ./ -B ./build
$ cmake --build ./build -j
```


## Project Structure
The project is organized as follow:
- ***doc*** folder which contains the documentation of the previous project
- ***script*** containing the python file provided by the previous student
- ***src*** where you should place the source code of the C++ and CUDA porting of your application