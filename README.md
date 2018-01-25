Optimized Update/Prediction Assignment for Lifting Transforms on Graphs III: Multiresolution examples.
========
This is the source code related with the paper "Optimized Update/Prediction Assignment for Lifting Transforms on Graphs", 
Eduardo Martínez-Enríquez, Jesús Cid-Sueiro, Fenrando Díaz-de-María and Antonio Ortega.
========
MULTIRESOLUTION EXAMPLES

- To obtain the results with image examples in a multiresolution framework (Table I), run "Test_multiresolution_image.m". 
Test images used in the experiments are included. Change the "size_block" parameter depending on the size of 
the test image.


NOTE: The algorithms use parfor commands (execute loop iterations in parallel) from the Parallel Computing Toolbox. 
"Parfor" loops can be changed to "for" loops, but the execution time will be higher.
