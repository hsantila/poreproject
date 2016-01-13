#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <string.h>
#include "mymath.h"
#include "tools.h"
//#include "mymath.h"
//#include "ewald_line.h"
//#include "io.h"


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
      printf("%s %d",typeinfo[i][0], strcmp(typeinfo[i][4],"yes"));
      temp_N=strtol(typeinfo[i][1],&dummyptr,10);
      for (int j=0;j<temp_N;j++)
      {
	types[indx]=typeinfo[i][0];
	Qs[indx]=strtol(typeinfo[i][2], &dummyptr,10);
	rs[indx]=strtol(typeinfo[i][3], &dummyptr,10);
	
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
      	types[indx+j]=typeinfo[i][0];
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
	 types[indx+j]=typeinfo[i][0];
	 Qs[indx+j]=strtol(typeinfo[i][2], &dummyptr,10);
	rs[indx+j]=strtol(typeinfo[i][3], &dummyptr,10);
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