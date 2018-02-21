# TRMOD
RMOD-K Thermodynamic extension

This program will be used to implement material non-linearities on the RMOD-K mechanical.
The main ideia is to use the mechanical mesh as a coarse mesh that will be used to share 
nodal/element temperatures to the mechanical model.

** BUGS
- Realloc does not work. It seems to be a visual studio DEBUG HEAP problem(Currently allocation more memory with malloc)
