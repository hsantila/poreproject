#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <string.h>
#include "mymath.h"
#include "tools.h"
//#include "mymath.h"
#include "Ewald_pore.h"
//#include "io.h"
int step(double** xyz,double* Qs, int Nmobile, double stepsize, double* box, double* polcoord, double* rs, double rpol,double rpore, double* energy)
{
double delta_u=0;
int acc=0;
int move_ndx=0;
double x=0,y=0,z=0; //store trial move coordinates here
int ion_type=0;
int temp=0;

//choose move, assumes that there are several types of ions = switch is possible
double move=rand()/(double)RAND_MAX;
//printf("%lf\n",move);
if (move<0.6)
{
	transl(stepsize,xyz,&move_ndx, Nmobile,box,&x,&y,&z);
//printf("transl");

	//check for move overlapping with polymer
	
	
	if(overlap_ion(Nmobile, xyz, move_ndx, rs, x, y, z, box)==1)
	
	{
		*energy=0;
		return 0;
	}
	if (overlap_pol(polcoord, x, y, rpol, rs[move_ndx], box)==1)
	{
		*energy=0;
		return 0;
		
	}
	if(outside_pore(polcoord, rpore, x, y)==1)
	{
		*energy=0;
		return 0;

	}
	


	delta_u=calc_delta_e(xyz, Qs,x,y,z, move_ndx);
	acc=accept(delta_u);
	if (acc==1)
	{
		update(xyz,move_ndx,x,y,z);
		*energy=delta_u;
	
	}
	else
		*energy=0;
	return acc;
}
else
{
 //printf("switch"); 
move_ndx=(int)(Nmobile*(double)rand() / ((double)(RAND_MAX)+(double)(1) ));
int move_ndx2=(int)(Nmobile*(double)rand() / ((double)(RAND_MAX)+(double)(1) ));
delta_u=calc_deltaE_switch(xyz, Qs, move_ndx,move_ndx2);
if (abs(Qs[move_ndx])-abs(Qs[move_ndx2])<0.1)
{
*energy=0;
return 0;

}	
//printf(" ndx1 %d q1 %lf  ndx2 %d q2 %lf deltaU %lf \n", move_ndx, Qs[move_ndx],move_ndx2, Qs[move_ndx2], delta_u);
	
	acc=accept(delta_u);
	if (acc==1)
	{
		x=xyz[move_ndx2][0];
		y=xyz[move_ndx2][1];
		z=xyz[move_ndx2][2];
		double x2=xyz[move_ndx][0];
		double y2=xyz[move_ndx][1];
		double z2=xyz[move_ndx][2];
		update(xyz,move_ndx,x,y,z);
		update(xyz,move_ndx2,x2,y2,z2);
		*energy=delta_u;
	
	}
	else
		*energy=0;
	return acc;
}	


}

void gen_poreconf(double* polcoord,double rpol,double rpore, double* box, double* Qs,double* rs,double** xyz,char** types,char*** typeinfo, int Ntypes)
{
  
polcoord[0]=box[0]/2.0;
polcoord[1]=box[1]/2.0;
char* dummyptr;
int temp_N=0;
int indx=0;
int accept=0;
double r=0;
double x;
double y;
double z;
int res=100;
int Nnonpore=0;
int tmp=0;

for(int i=0; i<Ntypes;i++)
{
  //only static particles are generated first
  if(strcmp(typeinfo[i][4],"no")==0)
  {
     // printf("%s %d",typeinfo[i][0], strcmp(typeinfo[i][4],"yes"));
      temp_N=strtol(typeinfo[i][1],&dummyptr,10);
      for (int j=0;j<temp_N;j++)
      {
	strcpy(types[indx],typeinfo[i][0]);
	Qs[indx]=strtod(typeinfo[i][2], &dummyptr);
	rs[indx]=strtod(typeinfo[i][3], &dummyptr);
	
	while(accept==0)
	{
	 x=polcoord[0]+2.0*rand()/(double)(RAND_MAX)*rpore-rpore;
	  y=polcoord[1]+2.0*rand()/(double)(RAND_MAX)*rpore-rpore;
	  z=(rand()/(double)RAND_MAX*box[2]);
	  
	  if(overlap_pol(polcoord, x, y, rpol, rs[indx], box))
	    continue;
	  if(outside_pore(polcoord, rpore, x, y))
	    continue;
	  if(overlap_ion(indx, xyz,indx, rs, x, y, z,box))
	    continue;
	  accept=1;
	  
	}
	accept=0;
	xyz[indx][0]=x;
	xyz[indx][2]=z;
	xyz[indx][1]=y;
	indx=indx+1;
      }	
      
    
  } 
  
  
}

for(int i=0; i<Ntypes;i++)
{
  //polymer beads
  if(strcmp(typeinfo[i][0],"P")==0)
  {
      
      temp_N=strtol(typeinfo[i][1],&dummyptr,10);
      for (int j=0;j<temp_N;j++)
      {


	xyz[indx+j][0]=polcoord[0];
	xyz[indx+j][2]=box[2]/temp_N*(j);
	xyz[indx+j][1]=polcoord[1];
      	strcpy(types[indx+j],typeinfo[i][0]);
	Qs[indx+j]=0;
	rs[indx+j]=rpol;
	
      }	

    indx=indx+temp_N;
  } 
  
  
} 
Nnonpore=indx;
accept=0;

for(int i=0; i<Ntypes;i++)
{  
    //generates the pore by randomly distributing points based on artificial radius (scaled by 0.5) per point
 
  double theta=0;
  
    temp_N=strtol(typeinfo[i][1],&dummyptr,10);
    double Aperpoint=2*PI*rpore*box[2]/temp_N;
    double dummyr=sqrt(Aperpoint/(PI))*0.9;
    int reject=0;
   
    if(strcmp(typeinfo[i][0],"H")==0)
  {
      
     for (int j=0;j<temp_N;j++)
    {
      while (accept==0)
      {
       theta = 2.0*PI*rand()/(double)RAND_MAX;
       z = rand()*box[2]/(double)RAND_MAX;
        //generates the pore by randomly distributing points based on artificial radius (scaled by 0.25) per point x=rpore*cos(theta);
       y=polcoord[1]+rpore*sin(theta);
       x=polcoord[0]+rpore*cos(theta);
       for(int k=Nnonpore; k<Nnonpore+j;k++)
       { 
	 if(SQR(xyz[k][0]-x)+SQR(xyz[k][1]-y)+SQR(xyz[k][2]-z)<SQR(dummyr))
	 {
	   reject=reject+1;
	   tmp=1;
	    break;
	 }
       }
       if (tmp==0)
       {
	 xyz[indx+j][0]=x;
	 xyz[indx+j][1]=y;
	 xyz[indx+j][2]=z;
	 strcpy(types[indx+j],typeinfo[i][0]);
	 Qs[indx+j]=strtod(typeinfo[i][2], &dummyptr);
	rs[indx+j]=strtod(typeinfo[i][3], &dummyptr);
	 accept=1;
	 
       }
       tmp=0;
      }
      accept=0;
    }
      
    
  } 
  
}
} 

//---------------------------------------------------------------------------------------------------------------
int overlap_pol(double* polcoord, double x, double y, double rpol, double r_ion, double* box)
{
double dx=polcoord[0]-x;
double dy=polcoord[1]-y;

dx -= dround(dx/box[0])*box[0];
dy -= dround(dy/box[1])*box[1];

int temp=0;

if (sqrt((dy)*(dy)+(dx)*(dx))<(r_ion+rpol))
	temp=1;

return temp;
}

//-------------------------------------------------------------------------------------------------------------------
int overlap_ion(int Np, double** xyz, int move_ndx, double* rs, double x, double y, double z, double* box)
{


double dx;
double dy;
double dz;

/*printf("box %lf %lf %lf\n", box[0],box[1],box[2]);
printf("rions %lf %lf\n", rions[0],rions[1]);
printf("Np %d", Np);*/



for (int i=0;i<Np;i++)
{
		
	dx=xyz[i][0]-x;
	dy=xyz[i][1]-y;
	dz=xyz[i][2]-z;

	dx -= dround(dx/box[0])*box[0];
	dy -= dround(dy/box[1])*box[1];
	dz -= dround(dz/box[2])*box[2];

	if (move_ndx!=i && (SQR(dx)+SQR(dy)+SQR(dz)<SQR(rs[i]+rs[move_ndx])))
	{
		return 1;
	} 	
		
}
return 0;
} 
//--------------------------------------------------------------------------------------------------------------------------
int outside_pore(double* polcoord, double rpore, double x, double y)
{
double dx=polcoord[0]-x; 
double dy=polcoord[1]-y;

if((dx*dx+dy*dy)<rpore*rpore)
  return 0;
else
   return 1;
}
//--------------------------------------------------------------------------------------------------------------------------
int transl(double step, double** xyz, int* move_ndx, int Nmobile, double* box, double* x, double* y, double* z)
{
int temp=0;

//translate one aprticle

//choose particle

*move_ndx=(int)(Nmobile*((double)rand() / ((double)(RAND_MAX)+(double)(1)) ));
//printf("Nions %d attemp move on %d \n",Nions,*move_ndx);
//move
double r[3]={rand()/(double)RAND_MAX, rand()/(double)RAND_MAX, rand()/(double)RAND_MAX};
*x=xyz[*move_ndx][0]+(r[0]-0.5)*step*2;
*y=xyz[*move_ndx][1]+(r[1]-0.5)*step*2;
*z=xyz[*move_ndx][2]+(r[2]-0.5)*step*2;
//*d=5;

if (*z<0)
	*z=*z+box[2];	
else if (*z>box[2])
	*z=*z-box[2];
return temp;

}
//-----------------------------------------------------------------------------------------------------------------------
double calc_delta_e(double** xyz, double* Qs, double x, double y, double z, int move_ndx)
{
// calc ewald energy+other neccesary energies

double new_coord[3]={x,y,z};
double old_coord[3]={xyz[move_ndx][0],xyz[move_ndx][1],xyz[move_ndx][2]};
      double Er1=Ewald_r_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx); 
      double Ek1=Ewald_k_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx);
//      double Ed1=Ewald_dipol(xyz, Qs, 0,0,0);
//     double Old_e=calc_energy(xyz,Qs);
//double Old_er=Ewald_r_space_line(xyz, Qs,0,0,0);
//double Old_ek=Ewald_k_space_line(xyz, Qs,0,0,0);		
xyz[move_ndx][0]=new_coord[0];
xyz[move_ndx][1]=new_coord[1];
xyz[move_ndx][2]=new_coord[2];
      double Er2=Ewald_r_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);
      double Ek2=Ewald_k_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);      
//      double Ed2=Ewald_dipol(xyz, Qs, 0,0,0);
//    double New_e=calc_energy(xyz,Qs);
//double New_er=Ewald_r_space_line(xyz, Qs,0,0,0);
//double New_ek=Ewald_k_space_line(xyz, Qs,0,0,0);	
//printf("old_e %lf new_e %lf delta_e %lf partial_new %lf part_old %lf part_delta %lf\n", Old_e, New_e, New_e-Old_e,(Er2+Ek2),(Er1+Ek1), (Er2+Ek2)-(Er1+Ek1) );   
//printf("1 real %lf k %lf real_tot %lf k_tot %lf\n", Er1, Ek1, Old_er, Old_ek);
//printf("2 real %lf k %lf real_tot %lf k_tot %lf\n", Er2, Ek2, New_er, New_ek);
xyz[move_ndx][0]=old_coord[0];
xyz[move_ndx][1]=old_coord[1];
xyz[move_ndx][2]=old_coord[2];


//return (Er2+Ek2+Ed2)-(Er1+Ek1+Ed1);
//dipole off
return (Er2+Ek2)-(Er1+Ek1);
}
//-----------------------------------------------------
double calc_deltaE_switch(double** xyz, double* Qs, int move_ndx, int move_ndx2 )
{
// calc ewald energy+other neccesary energies
double old_coord[3]={xyz[move_ndx][0],xyz[move_ndx][1],xyz[move_ndx][2]};
double old_coord2[3]={xyz[move_ndx2][0],xyz[move_ndx2][1],xyz[move_ndx2][2]};
double new_coord[3]={xyz[move_ndx2][0],xyz[move_ndx2][1],xyz[move_ndx2][2]};
double new_coord2[3]={xyz[move_ndx][0],xyz[move_ndx][1],xyz[move_ndx][2]};
//energy in old conf for particle1
//      double Er11=Ewald_r_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx); 
//      double Ek11=Ewald_k_space_part_line(xyz, Qs, old_coord, Qs[move_ndx], move_ndx);
//energy in old conf for particle2
//      double Er12=Ewald_r_space_part_line(xyz, Qs, old_coord2, Qs[move_ndx2], move_ndx2); 
//      double Ek12=Ewald_k_space_part_line(xyz, Qs, old_coord2, Qs[move_ndx2], move_ndx2);
//      double Ed1=Ewald_dipol(xyz, Qs, 0,0,0);
//     double Old_e=calc_energy(xyz,Qs);	
double dm_k1=Ewald_k_space_doublemove(xyz, Qs, old_coord, Qs[move_ndx],old_coord2, Qs[move_ndx2]);
double dm_r1=Ewald_r_space_doublemove(xyz, Qs, old_coord, Qs[move_ndx], move_ndx,old_coord2,Qs[move_ndx2], move_ndx2);
	
xyz[move_ndx][0]=new_coord[0];
xyz[move_ndx][1]=new_coord[1];
xyz[move_ndx][2]=new_coord[2];
xyz[move_ndx2][0]=new_coord2[0];
xyz[move_ndx2][1]=new_coord2[1];
xyz[move_ndx2][2]=new_coord2[2];

//      double Er21=Ewald_r_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);
//      double Ek21=Ewald_k_space_part_line(xyz, Qs, new_coord, Qs[move_ndx], move_ndx);
//     double Er22=Ewald_r_space_part_line(xyz, Qs, new_coord2, Qs[move_ndx2], move_ndx2);
//      double Ek22=Ewald_k_space_part_line(xyz, Qs, new_coord2, Qs[move_ndx2], move_ndx2);        
//      double Ed2=Ewald_dipol(xyz, Qs, 0,0,0);
//    double New_e=calc_energy(xyz,Qs);
double dm_k2=Ewald_k_space_doublemove(xyz, Qs, new_coord, Qs[move_ndx],new_coord2, Qs[move_ndx2]);
double dm_r2=Ewald_r_space_doublemove(xyz, Qs, new_coord, Qs[move_ndx], move_ndx,new_coord2,Qs[move_ndx2], move_ndx2);	
//printf("old_e %lf new_e %lf delta_e %lf partial_new %lf part_old %lf part_delta %lf\n", Old_e, New_e, New_e-Old_e,(Er2+Ek2),(Er1+Ek1), (Er2+Ek2)-(Er1+Ek1) );   
//printf("real %lf k %lf dipole %lf\n", Er2-Er1, Ek2- Ek1, Ed2-Ed1);
xyz[move_ndx][0]=old_coord[0];
xyz[move_ndx][1]=old_coord[1];
xyz[move_ndx][2]=old_coord[2];
xyz[move_ndx2][0]=old_coord2[0];
xyz[move_ndx2][1]=old_coord2[1];
xyz[move_ndx2][2]=old_coord2[2];
//printf("%lf %lf %lf %lf\n", dm_r1, dm_r2,dm_k1, dm_k2);

//printf("optimized switch deltaE %lf switch delta_partE %lf delta_moveE_tots %lf\n",(dm_r2+dm_k2)-(dm_k1+dm_r1),(Er21+Ek21+Er22+Ek22)-(Er11+Ek11+Er12+Ek12),New_e-Old_e);
//printf("optimized switch deltaE %lf switch delta_partE %lf delta_moveE_tots %lf\n",dm_r2-dm_r1,(Er21+Er22)-(Er11+Er12),New_e-Old_e);
//printf("optimized switch deltaE %lf switch delta_partE %lf delta_moveE_tots %lf\n",dm_k2-dm_k1,(Ek21+Ek22)-(Ek11+Ek12),New_e-Old_e);
//return (Er2+Ek2+Ed2)-(Er1+Ek1+Ed1);
//dipole off
return (dm_r2+dm_k2)-(dm_k1+dm_r1);
}
//--------------------------------------------------------------

double calc_energy(double** xyz, double* Qs)
{
// calc ewald energy+other neccesary energies
double er=0;
double ek=0;
double ed=0;
er=Ewald_r_space_line(xyz, Qs,0,0,0);
ek=Ewald_k_space_line(xyz, Qs,0,0,0);
//ed=Ewald_dipol(xyz,Qs,0,0,0);
//dipole off;
ed=0;
//printf("total energy er %lf ek %lf ed %lf\n",er,ek,ed);
return er+ek+ed;
}

//--------------------------------------------------------------------

void update(double** xyz, int move_ndx, double x, double y, double z)
{

xyz[move_ndx][0]=x;
xyz[move_ndx][1]=y;
xyz[move_ndx][2]=z;


}
//---------------------------------------------------------------------
int accept(double delta_u)
{
int acc=0;
double tmp=rand()/(double)RAND_MAX;
//printf("delta u %lf, exp(-deltau) %lf, random %lf\n", delta_u,exp(-delta_u), tmp);
if (exp(-delta_u)>tmp)
	acc=1;
return acc;
}

