Simply run H2CO.m which is the main program.
Inside H2CO.m, the block

Vmax=15000               % cm-1

M = 15000;    Nsets=5;   % number of pts and number of sets to add into the collocation eq.

Ng = 15000

controls: 
Vmax - the maximum PES value to accept 
M*Nset - maximum number of collocation points
M - number of collocation points in one block, set M >= Ng
Nset - maximum number of blocks
Ng - number of basis functions 
