# input format
a connected undirect graph, with each edge appears only once in the file

# run
* change the parameter in hi_treem.c, between "//PARAM" and "//PARAM END", and then make
* cat n10000.inp | ./hi_treem
* the output contains information for preprocessing and max-flow calculation

# how to generate input file
you can refer to igraph lib in python, for example, uniformly random graph (Erdos_Renyi), scale-free graph (Barabasi), unit disk graph (GRG).
