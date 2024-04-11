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

## Merge ⭐
```
$ git config pull.rebase false 
```

'''
git add .
git commit -m "aggiornamento"
git merge
git push 
'''

 - commit and synch

## Test ✅ / ❌
- gen_rotations : ✅
- Label : ✅
- select_vertex : ✅
- hood : ✅
- select_label : ✅
- gen_initial_label : ❌
- gen_ring_classes : ❌
- calc_bound : ✅
- gen_bond : ❌
## Rdkit
Data structure

    Basic data structure to describe a molecule: RDKit::RWMol or RDKit::ROMol

Parsing

    Parse a string that describe the protein (in PDB format): RDKit::PDBBlockToMol
    Parse a string that describe the ligand (in Mol2 format): RDKit::Mol2BlockToMol

Headers

    GraphMol/FileParsers/FileParsers.h <- to parse a file format GraphMol/GraphMol.h <- data structure definition

