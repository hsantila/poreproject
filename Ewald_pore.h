/*------------------------------------------------------------

HEADER-FILE for ewald_pore.c

------------------------------------------------------------*/

#ifndef EWALD
#define EWALD

/* maximal k-vector */
#define Maxkmax 20

/*----------------------------------------------------------*/

extern void Ewald_init_line(double* length, int particlenumber, int polymern, double Ewaldalpha,
		double realspacecutoff, int* reciprocalspacecutoff, 
		double coulombprefactor, double epsilon, double tau, double* ccoord);
extern double Ewald_r_space_line(double**, double *,
			    double *, double *, double *);
extern double Ewald_k_space_line(double**, double *, 
			    double *, double *, double *);
extern double Ewald_r_space_part_line(double**, double *, double*,
			    double , int);
extern double Ewald_k_space_part_line(double** restrict, double* restrict, double* restrict,
			    double , int);
extern double Ewald_k_subset(double *, double *, double *, double *, int, int);
extern double Ewald_dipol(double**, double *,
			  double *, double *, double *);

double Ewald_r_space(double** xyz, double *Q, double** F);
//By hanne
void cylinder_to_cylinder_k(int pndx, double* fx, double* fy);
void cylinder_to_cylinder_r(double* fx, double* fy, int pndx);
void ion_to_cylinder_k(double** xyz, double* Q,int pndx, double* fx, double* fy);
void ion_to_cylinder_r(double** xyz,double* Q, int pndx, double* fx, double* fy);
double Ewald_r_space_line_parts(double**, double *,
			    double* , double* , double* );
double Ewald_k_space_line_parts(double**, double *, 
			    double* , double* , double*,double*);
double Ewald_k_space_doublemove(double** xyz, double *Q, double* MCcoord, double MCQ, double* MCcoord2, double MCQ2);
double Ewald_r_space_doublemove(double** xyz, double *Q, double *MCcoord, double MCQ, int MCion,double *MCcoord2, double MCQ2, int MCion2);

#endif