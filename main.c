#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "io.h"
int main(int argc, char **argv) {
    printf("Hello World!");

    //step parameters
int nstep=0;	//number of steps
double stepsize=0; 
int start_t=0;
int relax=0;
int acc=0;
int acc_sum=0;
//int lb;
int step_upd=0;
int gro_outp=0;
int var_outp=0;
int seed=0;
int temp;
//filenames

char* file_in;
char* file_out_base;
char* file_gro;
char file_out_gro[40];
char file_out_energ[40];
//indexes and indicators for input

int is_g=0;
int is_i=0;
int is_o=0;
int is_e=0;
int is_t=0;


int ndx_g=0;
int ndx_i=0;
int ndx_o=0;
int ndx_e=0;

//int ndx_t=0;



//Ewald parameters
double alpha=0;
//double coulomb=0;
double rcut=0;
int kcut[3];
//isotropic ewald convergence
int kmax_c=0;
int kmin_c=0;

//system spesification
//this is a uniformly charged cylinder in a pore. Pore can be made of beads or not

double polcoord[2];
//tells whether pore is of beads
int porebeads=0;
//total number of particles
int Ntot=0;
//total number of ions
int Nions=0;
int Ntypes=0;
//types indexed as in xyz
char**  types;
//contains amout, radius, charge. To pe parsed to separate matrixes
char*** typeinfo;
//charge per lenght of polymer
double tau;
double rpol;
double rpore;
//bead radiuses
double* rs;

double box[3];
double** xyz;
double* Qs; 

//analysis variables
double energy=0;
double etot=0;
int avg_count=0;
//------------------------------------------------------------------------

//                Reading input parameters

//------------------------------------------------------------------------



for (int i=1; i<argc; i++){
    if (strcmp(argv[i],"-g")==0){
      is_g = 1;
      ndx_g = i;
    }
    if (strcmp(argv[i],"-i")==0){
      is_i = 1;
      ndx_i = i;
    }
    if (strcmp(argv[i],"-o")==0){
      is_o = 1;
      ndx_o = i;
    }
    if (strcmp(argv[i],"-e")==0){
	  is_e = 1;
     	 ndx_e = i;
    }
    if (strcmp(argv[i],"-r")==0){
	seed =(int)strtod(argv[i+1], NULL); // optional seed for random number generator
    }		
    if (strcmp(argv[i],"-t")==0){
	is_t=1;
	start_t =(int)strtod(argv[i+1], NULL); // If gro file given
    }
    	
}

//check for input options
if (is_o==0 || is_i==0){
	printf("Invalid input\n");
	printf("Usage %s -i [datafile] -o [outfile] -s [MC steps] -g [grofile, optional] -t [starting time, if grofile specified] -e  [Rcut] [kmin] [kmax]\n", argv[0] );
	
	exit(0);
}

if((is_g==0 && is_t==1) || (is_g==1 && is_t==0))
{
	printf("Invalid input in optional parameters\n");
	printf("Usage %s -i [datafile] -o [outfile] -s [MC steps] -g [grofile, optional] -t [starting time, if grofile specified] -e [Rcut] [kmin] [kmax]\n", argv[0] );
	
	exit(0);
}


if(is_i==1)
{
	file_in=argv[ndx_i+1];
}
if(is_o==1)
{
	file_out_base=argv[ndx_o+1];
}
if(is_g==1)
{
	file_gro=argv[ndx_g+1];
}
strcpy(file_out_gro, file_out_base);
strcat(file_out_gro,".gro");
strcpy(file_out_energ, file_out_base);
strcat(file_out_energ,".enrg");

//------------------------------------------------------------------------------------------

//               Reading information from files, memory allocation for dynamic and system dependent matrixes

//------------------------------------------------------------------------------------------



	
	

temp=read_input(file_in,file_out_base,&relax,&nstep, &step_upd, &gro_outp, &var_outp, &stepsize, &alpha, &rcut, kcut, &typeinfo,&Ntypes,&Ntot,&rpore,&rpol,&tau, box);
 printf("%s %s %s %s %s\n", typeinfo[0][0],typeinfo[0][1], typeinfo[0][2],typeinfo[0][3],typeinfo[0][4]);
if(temp==0)
{
	printf("Cannot read input file\n");
	exit(0);
}


xyz=malloc(sizeof(double*)*Ntot);
	for (int i=0; i<Ntot;i++)
		xyz[i]=calloc(1,sizeof(double)*3);

Qs=calloc(1,sizeof(double)*Ntot);
rs=calloc(1,sizeof(double)*Ntot);

types=malloc(sizeof(char*)*Ntot);
	for (int i=0;i<Ntot;i++)
	    types[i]=malloc(sizeof(char)*5);

if (is_g==1)
{
//read in configuration, these are temporary variables for comparison between gro and input file

char*** gro_typeinfo;
int gro_ntot=0;
int gro_ntypes=0;
int chargesum=0;
char* dummyptr;

//read_gro(file_gro, xyz, gro_typeinfo, &gro_ntot,&gro_ntypes, box);

	if(gro_ntypes!=Ntypes)
	{
	printf("Missmatch in amount of particle types read");
	exit(0);
	}
	
printf("Read in configuration\n");
for(int i=0; i<gro_ntypes;i++)
{
  //hanne, fix the numer of stars in typeinfo
  printf("type: %s number: %s charge: %s radius: %s \n", gro_typeinfo[i][0],gro_typeinfo[i][1],gro_typeinfo[i][2], gro_typeinfo[i][3]);
  chargesum=chargesum+strtol(gro_typeinfo[i][1], &dummyptr,10)*strtol(gro_typeinfo[i][2],&dummyptr,10);
}

printf("Total charge of the system %d \n",chargesum);



}
else
{
  
//generate configuration
gen_poreconf(polcoord,rpol,rpore, box, Qs,rs,xyz,types, typeinfo, Ntypes);

}
write_gro(file_out_gro,0, xyz, Ntot, polcoord, box, types, typeinfo);
return 0;
}
