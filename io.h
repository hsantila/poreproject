/*------------------------------------------------------------

HEADER-FILE for io.c

------------------------------------------------------------*/

#ifndef IO
#define IO

/* maximal k-vector */
#define Maxkmax 20

/*----------------------------------------------------------*/

int read_input(char* file_in,char* file_out,int* relax,int* nsteps, int* step_upd, int* gro_outp, int* var_outp, double* stepsize, double* alpha, double* rcut,int* kcut, char**** typeinfo, int* Ntypes,int* Ntot,double* rpore,double* rpol, double* tau, double* box);
int write_gro(char* file,int step, double** xyz, int Nall, double* polcoord, double* box, char** types, char*** typeinfo);
int write_energy(char* file,int step, double* energy);
int read_gro(char* file_gro,double** xyz,double* polcoord,double* Qs,double* rs, char*** typeinfo,int ntypes,char** types, int* gro_ntot,double* box);
#endif
