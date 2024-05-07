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

## Important for adding new EXECUTABLE TESTS   ⭐
 Nel file CMakeLists.txt è presente una parte di add executable, ogni qual volta si voglia testare, 
 o aggiungere una nuova funzione al file di test, bisogna seguire questo procedimento:
 - inserisco l'header della funzione nella classe test.hpp o testFra.hpp
 - vado nel file cmake e per ogni funzione aggiunta scrivo "${source_path}/NOMEFUNZIONE.cpp" di fianco a "${source_path}/testFra.cpp"
 - salvo cmake file ctrl S
 - buildo nuovamente l'app ed eseguo ./test o ./testFra
NON AGGIUNGERE ALTRO NEL CMAKE FILE

## CMakeCache.txt Problems 
```
cd build
rm CMakeCache.txt
cd ..
$ cmake -S ./ -B ./build
```

