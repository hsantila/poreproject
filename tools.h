
/*------------------------------------------------------------

HEADER-FILE for tools.c

------------------------------------------------------------*/

#ifndef TOOLS
#define TOOLS

void gen_poreconf(double* polcoord,double rpol,double rpore, double* box, double* Qs,double* rs,double** xyz,char** types,char*** typeinfo, int Ntypes);
int overlap_pol(double* polcoord, double x, double y, double rpol, double r_ion, double* box);
int outside_pore(double* polcoord, double rpore, double x, double y);
int overlap_ion(int Np, double** xyz, int move_ndx, double* rs, double x, double y, double z, double* box);

#endif