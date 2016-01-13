#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include <errno.h>

int write_gro(char* file,int step, double** xyz, int Nall, double* polcoord, double* box, char** types, char*** typeinfo)
{
FILE *gro;

if ((gro = fopen(file,"a")) != NULL)
{


int res=100;
int atomnumber=1;
int resnumber=1;
char* str;
char atomn[5];
char atomtype[5];
double zkoord=0;


fprintf(gro,"Polymer and ions in a pore, t= %d\n", step);
fprintf(gro, "%d\n", Nall);

//write polymers
for(int i=0; i<Nall;i++)
{

      if(strcmp(types[i],"P")!=0 && strcmp(types[i],"H")!=0)
      {
	fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resnumber, types[i],types[i],atomnumber,xyz[i][0], xyz[i][1],xyz[i][2]);	
	resnumber=resnumber+1;
	atomnumber=atomnumber+1;
	continue;
      }
      else
      {
	fprintf(gro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resnumber, types[i],types[i],atomnumber,xyz[i][0], xyz[i][1],xyz[i][2]);
	atomnumber=atomnumber+1; 
	if(strcmp(types[i],types[i-1])==0)
	  continue;
	else
	  resnumber=resnumber+1;

      }

      
}
//write box

fprintf(gro,"%f %f %f\n", box[0], box[1], box[2]);
fclose(gro);
return 1;
}
else
{
	printf("Cannot open the gro file for writing, The reason *may* be %s\n", strerror(errno));
	return 0;
}
}

//------------------------------------------------------------------------------------------------------------
int read_input(char* file_in,char* file_out,int* relax,int* nsteps, int* step_upd, int* gro_outp, int* var_outp, double* stepsize, double* alpha, double* rcut,int* kcut, char**** typeinfo, int* Ntypes,int* Ntot,double* rpore,double* rpol, double* tau, double* box)
{

FILE *inp;
FILE *outp;
int ret=0;
char dataid[500];
double dataval=0.0; 
char retstring[60];
fpos_t types_begin;
char str1[10], str2[10], str3[10], str4[10], str5[10];
char* dummyptr;
char*** typeinfo_tmp;

if ((inp = fopen(file_in,"r")) != NULL)
{
	
	while(ret!=EOF)
	{
		ret=fscanf(inp,"%s", dataid);
	
		if(strstr(dataid,"r_pol")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*rpol=dataval;
		}
		if(strstr(dataid,"r_pore")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*rpore=dataval;
		}
		if(strstr(dataid,"taupol")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*tau=dataval;
		}

		if(strstr(dataid,"BOXX")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			box[0]=dataval;
		}
		if(strstr(dataid,"BOXY")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			box[1]=dataval;
		}
		if(strstr(dataid,"BOXZ")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			box[2]=dataval;
		}
		if(strstr(dataid,"stepsize")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*stepsize=dataval;
		}
		if(strstr(dataid,"stepupdateinterval")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*step_upd=(int)dataval;
		}
		if(strstr(dataid,"dataprintoutinterval")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*var_outp=(int)dataval;
		}
		if(strstr(dataid,"groprintoutinterval")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*gro_outp=(int)dataval;
		}
		if(strstr(dataid,"RELAXATIONSTEPS")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*relax=(int)dataval;
		}
		if(strstr(dataid,"MCsteps")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*nsteps=(int)dataval;
		}
		if(strstr(dataid,"Ealpha")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*alpha=dataval;
		}
		if(strstr(dataid,"Rcut")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			*rcut=dataval;
		}
		if(strstr(dataid,"EwaldKx")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			kcut[0]=(int)dataval;
		}
		if(strstr(dataid,"EwaldKy")!=NULL)
		{
			
			ret=fscanf(inp, " %lf\n", &dataval);
			kcut[1]=(int)dataval;
		}
		if(strstr(dataid,"EwaldKz")!=NULL)
		{
			ret=fscanf(inp, " %lf\n", &dataval);
			kcut[2]=(int)dataval;
		}
		if(strstr(dataid,"TYPES")!=NULL)
		{
		 
		  fgets(retstring, 60, inp);
		   fgetpos(inp,&types_begin); 
		
		
		  fgets(retstring, 60, inp);
		  
		  while (retstring!=NULL && retstring[0]!='\n')
		  {
		   
		
		    fgets(retstring, 60, inp);
		    
		     *Ntypes=*Ntypes+1;
		  }
		
		
		  typeinfo_tmp=malloc((*Ntypes)*sizeof(char**));
		    for(int i=0;i<(*Ntypes);i++)
		    {
		      typeinfo_tmp[i]=malloc(5*sizeof(char*));
		      for(int j=0;j<5;j++)
		      {
			typeinfo_tmp[i][j]=malloc(11*sizeof(char));
			
		      }
		    }
		  
		  fsetpos(inp, &types_begin);  
		  for (int i=0;i<(*Ntypes);i++)
		  {  
		  ret=fscanf(inp, "%s %s %s %s %s\n", str1,str2,str3,str4, str5);
		
		  strcpy(typeinfo_tmp[i][0],str1);
	
	
		  strcpy(typeinfo_tmp[i][1],str2);
		  strcpy(typeinfo_tmp[i][2],str3);
		  strcpy(typeinfo_tmp[i][3],str4);
		  strcpy(typeinfo_tmp[i][4],str5);
		  
		  *Ntot=*Ntot+strtol(typeinfo_tmp[i][1],&dummyptr,10);
		  
		  }
		  *typeinfo=typeinfo_tmp;
		  ret=fscanf(inp,"%s", dataid);
		  
		} 

		
	}




	fclose(inp);
	
outp=fopen(file_out,"w");
 fprintf(outp,"# BOX_X\nBOXX %lf\n\n",  box[0]);
 fprintf(outp,"# BOX_Y\nBOXY %lf\n\n",  box[1]);
 fprintf(outp,"# BOX_Z\nBOXZ %lf\n\n",  box[2]);
 fprintf(outp,"# line charge \ntaupol %lf\n\n",  *tau);
 fprintf(outp,"# MC steps\nMCsteps %d\n\n", *nsteps );
 fprintf(outp,"# GRO printout interval \ngroprintoutinterval %d\n\n", *gro_outp);
 fprintf(outp,"# Data printout interval \ndataprintoutinterval %d\n\n", *var_outp);
 fprintf(outp,"# Step update interval \nstepupdateinterval %d\n\n", *step_upd);
 fprintf(outp,"# Initial steps to be disregarded \nRELAXATIONSTEPS %d\n\n", *relax);
 fprintf(outp,"# stepsize\nstepsize %lf\n\n", *stepsize );
 fprintf(outp,"# K vectos in x\nEwaldKx %d\n\n",kcut[0]);
 fprintf(outp,"# K vectos in y\nEwaldKy %d\n\n",kcut[1]);
 fprintf(outp,"# K vectos in z\nEwaldKz %d\n\n",kcut[2]);
 fprintf(outp,"# Direct space cutoff\nRcut %lf\n\n",*rcut);
 fprintf(outp,"# Ewald summation alpha\nEalpha %lf\n\n",*alpha );
 fprintf(outp,"# TYPES amount charge radius static\n ");
 for (int i=0; i<*Ntypes;i++)
 {
   fprintf(outp,"%s %s %s %s %s\n", typeinfo_tmp[i][0],typeinfo_tmp[i][1], typeinfo_tmp[i][2],typeinfo_tmp[i][3],typeinfo_tmp[i][4]);

 }
 
// fprintf(outp,"# \n %d\n\n", );
fclose(outp);

	return 1;
}
else
{
	printf("Cannot read input parameters\n");
	return 0;
}

}

  
  


