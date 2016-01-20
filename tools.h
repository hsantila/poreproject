
/*------------------------------------------------------------

HEADER-FILE for tools.c

------------------------------------------------------------*/

#ifndef TOOLS
#define TOOLS

void gen_poreconf(double* polcoord,double rpol,double rpore, double* box, double* Qs,double* rs,double** xyz,char** types,char*** typeinfo, int Ntypes);
int overlap_pol(double* polcoord, double x, double y, double rpol, double r_ion, double* box);
int outside_pore(double* polcoord, double rpore, double x, double y);
int overlap_ion(int Np, double** xyz, int move_ndx, double* rs, double x, double y, double z, double* box);
int step(double** xyz,double* Qs, int Nmobile, double stepsize, double* box, double* polcoord, double* rs, double rpol,double rpore, double* energy);
int transl(double step, double** xyz, int* move_ndx, int Nmobile, double* box, double* x, double* y, double* z);
double calc_delta_e(double** xyz, double* Qs, double x, double y, double z, int move_ndx);
double calc_deltaE_switch(double** xyz, double* Qs, int move_ndx, int move_ndx2 );
double calc_energy(double** xyz, double* Qs);
void update(double** xyz, int move_ndx, double x, double y, double z);
int accept(double delta_u);

#endif
