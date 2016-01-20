#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "io.h"
#include "Ewald_pore.h"
int main(int argc, char **argv) {
  

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
char* dummyptr;
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
//Number of Charged
int Ncharged=0;

//total number of ions (mobile particles)
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
    if (strcmp(argv[i],"-s")==0){
	nstep =(int)strtod(argv[i+1], NULL); // If parameter file contains steps, this is overwritten
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

if(temp==0)
{
	printf("Cannot read input file\n");
	exit(0);
}
printf("Ntot %d", Ntot);

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
//could use better failsafes

int gro_ntot=0;
int ret_gro=0;

ret_gro=read_gro(file_gro,xyz,polcoord,Qs,rs, typeinfo,Ntypes,types, &gro_ntot,box);

if(ret_gro==0)
{
  
  printf("Failure in gro read-in\n");
  exit(0);
  
}
	if(gro_ntot!=Ntot)
	{
	printf("Missmatch in amount of particles");
	exit(0);
	}
	
printf("Read in configuration\n");

for (int i=0;i<Nions;i++)
{
	printf("%s %d %lf %lf %lf %lf %lf\n",types[i], i, xyz[i][0],xyz[i][1],xyz[i][2], Qs[i], rs[i]);	

}


}
else
{
  
//generate configuration
gen_poreconf(polcoord,rpol,rpore, box, Qs,rs,xyz,types, typeinfo, Ntypes);

}

//Calculate the number of charged particles
for(int i=0;i<Ntot;i++)
{
  
  if (Qs[i]!=0)
    Ncharged++;
  
}
//Calculate number of ions
for(int i=0; i<Ntypes; i++)
{
  if(strcmp(typeinfo[i][4],"no")==0)
    Nions=Nions+strtol(typeinfo[i][1],&dummyptr,10);
  
}

write_gro(file_out_gro,0, xyz, Ntot, polcoord, box, types, typeinfo);

//------------------------------------------------------------------------

//               Ewald convergence test

//------------------------------------------------------------------------

if(is_e==1)
{
char file_abs[50];
char file_move[50];
char file_abs_parts[50];
char tmpstr[25];





/*FILE* move=fopen("move_energy.txt","w");
FILE* abs_ene=fopen("total_energy","w");
fprintf(move,"#0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9\n");
fprintf(abs_ene,"#0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9\n");*/
FILE* abs;
FILE* move;
FILE* abs_parts;


double e1=0;
double delta_e=0;
double delta_e1=0;
double delta_e2=0;
double delta_e3=0;




double ek_ii=0;
double ek_pi=0;
double ek_pp=0;
double ek_self=0;
double er_ii=0;
double er_pi=0;
double er_pp=0;
double ek_tot=0;
double er_tot=0;

int move_ndx;

double amin=0;
double amax=0;
double alpha=0;
char direct;
int k_const=0;
int karr[3];
int ndx=0;
move_ndx=Nions-1;



	if (ndx_e+8<argc)
	{
		printf("Not enough input parameters for ewald convergence test\n");
		printf("Usage %s -i [datafile] -o [outfile] -s [MC steps] -g [grofile, optional] -t [starting time, if grofile specified] -e [Rcut] [amin] [amax] [kmin] [kmax] [direction x, y or z] [kconst]\n", argv[0] );
		exit(0);
		//kconst=0 symmetry
	}	 	
	else
	{

		rcut = strtod(argv[ndx_e+1], NULL);
		amin=(double)strtod(argv[ndx_e+2], NULL);
		amax=(double)strtod(argv[ndx_e+3], NULL);
		kmin_c = (int)strtod(argv[ndx_e+4], NULL);
		kmax_c = (int)strtod(argv[ndx_e+5], NULL);
		direct=argv[ndx_e+6][0];

		if (direct=='x')
			ndx=0;
		if (direct=='y')
			ndx=1;
		if (direct=='z')
			ndx=2;

		k_const = (int)strtod(argv[ndx_e+7], NULL);
		printf("rcut %lf amin %lf amax %lf kmin %d kmax %d\n", rcut, amin, amax, kmin_c, kmax_c);	
		//perform ewald convergence test		
		for (int i=kmin_c; i<kmax_c+1; i++)
		{
			//fprintf(move,"%d", i);
			//fprintf(abs_ene,"%d", i);

			strcpy(file_abs,"tot_energy_conv");
			strcpy(file_move,"move_energy_conv");
			strcpy(file_abs_parts,"tot_energy_parts");

			sprintf(tmpstr,"_Rcut%1.2lfk%d%c%d",rcut,k_const,direct,i);
		
			strcat(file_move,tmpstr);
			strcat(file_move,".xvg");



			
			strcat(file_abs,tmpstr);
			strcat(file_abs,".xvg");
			strcat(file_abs_parts,tmpstr);
			strcat(file_abs_parts,".xvg");		

			printf("filuabs %s\n", file_abs);
			printf("filumove %s\n", file_move);


			abs=fopen(file_abs,"w");
			move=fopen(file_move,"w");
			abs_parts=fopen(file_abs_parts,"w");

			for (int j=0;j<amax*100;j++)
			{
				alpha = amin+(j+1)/100.0;
				if (k_const==0)
				{
				karr[0]=i;
				karr[1]=i;
				karr[2]=i;
				}
				else
				{

				karr[0]=i;
				karr[1]=i;
				karr[2]=k_const;
			
				}
				
	
				printf("kx %d ky %d kz %d\n", karr[0],karr[1], karr[2]);

				Ewald_init_line(box, Nions,1, alpha,rcut,karr , 1, 1, tau, polcoord);

				e1=calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]+0.1;
				xyz[move_ndx][1]=xyz[move_ndx][1]+0.5;
				xyz[move_ndx][2]=xyz[move_ndx][2]+0.5;
	
				delta_e=e1-calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]-0.1;
				xyz[move_ndx][1]=xyz[move_ndx][1]-0.5;
				xyz[move_ndx][2]=xyz[move_ndx][2]-0.5;


				xyz[move_ndx][0]=xyz[move_ndx][0]+0.5;
				xyz[move_ndx][1]=xyz[move_ndx][1]+0.5;
	
	
				delta_e1=e1-calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]-0.5;
				xyz[move_ndx][1]=xyz[move_ndx][1]-0.5;
	

	
				xyz[move_ndx][2]=xyz[move_ndx][2]+0.5;
	
				delta_e2=e1-calc_energy(xyz, Qs);
			
				xyz[move_ndx][2]=xyz[move_ndx][2]-0.5;

				xyz[move_ndx][0]=xyz[move_ndx][0]+7;
				xyz[move_ndx][1]=xyz[move_ndx][1]+7;
				xyz[move_ndx][2]=xyz[move_ndx][2]+0.1;
	
				delta_e3=e1-calc_energy(xyz, Qs);
				
				xyz[move_ndx][0]=xyz[move_ndx][0]-7;
				xyz[move_ndx][1]=xyz[move_ndx][1]-7;
				xyz[move_ndx][2]=xyz[move_ndx][2]-0.1;

				fprintf(move,"%.12lf %.12lf %.12lf %.12lf %.12lf\n", alpha, delta_e, delta_e1, delta_e2, delta_e3);	
				fprintf(abs,"%lf %lf\n",alpha, e1);

				//ek_tot=Ewald_r_space(xyz, Qs, 0);
				ek_tot=Ewald_k_space_line_parts(xyz, Qs, &ek_ii, &ek_pi,  &ek_pp, &ek_self);
				er_tot=Ewald_r_space_line_parts(xyz, Qs, &er_ii, &er_pi, &er_pp);
				fprintf(abs_parts,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", alpha, e1,ek_tot, er_tot,ek_pi, ek_ii, ek_pp, er_pi, er_ii, er_pp, ek_self);
			
	
										
			}			

				//fprintf(move,"\n");
				//fprintf(abs_ene,"\n");
			fclose(move);
			fclose(abs);
			fclose(abs_parts);		
		}
		exit(0);	
	}


}

//------------------------------------------------------------------------

//              Do monte carlo!

//------------------------------------------------------------------------

for (int i=0;i<Ntot;i++)
{
	printf("%s %d %lf %lf %lf %lf %lf\n",types[i], i, xyz[i][0],xyz[i][1],xyz[i][2], Qs[i], rs[i]);	

}
//initiate ewald
Ewald_init_line(box, Nions,1, alpha,rcut, kcut, 1, 1, tau, polcoord);

//initioate total energy
etot=calc_energy(xyz,Qs);

printf("askel %d energia %lf\n",-1, etot);

//set the random number generator to runtime, doesn't work with supercomputers
//srand ( (unsigned)time ( NULL ) );
//set the random number generator to input seed
srand(seed);

for (int i=start_t; i<=nstep;i++)
{
	//step function needs to be updated
	acc=step(xyz,Qs,Nions, stepsize,box, polcoord, rs, rpol,rpore, &energy);

	acc_sum=acc_sum+acc;
	etot=etot+energy;
		
	
	//if (acc==1)
//	{
		
		//printf("askel %d etot %lf deltaE %lf systemtot %lf \n",i, etot,energy,calc_energy(xyz,Qs));
			
//	}
	//else
	//	printf(" move rejected \n");
	//printf("step %d acc_sum %d\n", i, acc_sum);
	
	if(i<relax)
	{
		//do relaxation
		
		continue;

	}

	

	if((i+1)%step_upd == 0)
	{
		//update stepsize
		if((1.0*acc_sum/step_upd)>0.75)
		{
			stepsize=stepsize*2.0;
			if (stepsize>box[0] && stepsize>box[1] && stepsize>box[2])
			{
				printf("Attempt to increase step size over box size on timestep %d\n",i);
				stepsize=stepsize/10.0;
			}
			else 
				printf("Stepsize updated to %f on timestep %d\n", stepsize,i);			
							
		}
		if((1.0*acc_sum/step_upd)<0.25)
		{
						
			stepsize=stepsize/2.0;
			if(stepsize<0.1)
			{
				printf("Attempt to decrease steps size under limit 0.1 on timestep %d\n",i);
				stepsize=stepsize*10.0;
				
			}
			else
				printf("Stepsize updated to %f on timestep %d\n", stepsize,i);	

		}
		acc_sum=0.0;
	}
	
	if (i%gro_outp==0)
	{
	//basic printouts, energy, gro
		temp=write_gro(file_out_gro,i, xyz, Ntot, polcoord, box, types, typeinfo);
		
		if (temp==0)
			exit(0);
		printf("step %d acc_sum %d energy %f\n", i, acc_sum, etot);
		temp=write_energy(file_out_energ,i, &etot);
	}

		

}

//--------------------------------------------------------------------------------

//				Memory handling

//--------------------------------------------------------------------------------



for (int i=0; i<Ntot;i++)
	free(xyz[i]);
free(xyz);


free(Qs);
free(rs);


for (int i=0;i<Ntypes;i++)
{
	for(int j=0;j<5;j++)
	{
		free(typeinfo[i][j]);
	
	}
free(typeinfo[i]);	
}
free(typeinfo);


printf("Ntot %d", Ntot);



for (int i=0;i<Ntot;i++)
{
 
   free(types[i]);
}
free(types);


return 0;
}
