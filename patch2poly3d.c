#include "interact.h"
#include "properties.h"

/* 

reads in a rectangular Okada style fault geometry from "geom.in" in
interact format and a "bc.in" boundary condition file and produces
output for Poly3D




*/


int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  FILE *in;
  int i,opcode,j,ecnt,dim,tri;
  COMP_PRECISION a[6],b[6],corner[MAX_NR_EL_VERTICES*3],lloc,wloc;
  char bcstring[200];
  if(argc > 2){
    fprintf(stderr,"%s tri[0] \n\t reads in interact patch format from %s\n",
	    argv[0],GEOMETRY_FILE);
    fprintf(stderr,"\t reads in interact boundary condition file from %s\n",
	    BC_FILE);
    fprintf(stderr,"\t and writes Poly3D output file to stdout\n");
    fprintf(stderr,"\t if tri is 1, will use two triangles, if zero will use quads\n");
    
    exit(-1);
  }
  if(argc==2)
    sscanf(argv[1],"%i",&tri);
  else
    tri=0;
  if(tri)
    fprintf(stderr,"%s: WARNING: using two triangles for each quad\n",argv[0]);
  else
    fprintf(stderr,"%s: using quads\n",argv[0]);
  /* test if we can deal with the boundary conditions */
  in=myopen(BC_FILE,"r");
  if(fscanf(in,"%i",&opcode) != 1){ /* op mode */
    fprintf(stderr,"%s: bc code read error\n",argv[0]);
    exit(-1);
  }
  switch(opcode){		/* test for implemented modes */
  case 1:			/* that's OK */
    break;
  default:
    fprintf(stderr,"%s: bc error: bc operational code %i not implemented\n",
	    argv[0],opcode);
    exit(-1);
    break;
  }
  fclose(in);
  /* read geometry, BCs and some other stuff */
  initialize(&medium, &fault,
	     READ_FAULT_PROPERTIES_DEF,MAX_NR_FLT_FILES_DEF,FALSE,TRUE,
	     COHESION_DEF,a,b, READ_STRESS_RELATION_FACTORS_DEF,
	     PRINT_SLIPLINE_DEF,TRUE,EPS_COMP_PREC,FALSE,EPS_COMP_PREC,
	     FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,SVD_SOLVER,
	     FALSE,PRESSURE_DEF,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,SVD_THRESHOLD,
	     FALSE,FALSE);
  /* 
     
     write Poly3D input file 

  */
  /* constants */
  printf("*********************************************************\n");
  printf("*               Section 1: CONSTANTS                    *\n");
  printf("*********************************************************\n\n");
  
  printf("title1 = \"automatically generated from interact input\"\n");
  printf("title2 = \"test\"\n\n");  
  
  printf("* elastic properties\n");
  printf("shear_mod = %g\n",SHEAR_MODULUS);
  printf("psn_ratio = %g\n",POISSON_NU);
  printf("youngs_mod = \nbulk_mod = \nlame_lambda = \n\n");

  printf("rem_bc_type = stress\n\n");

  /* prestress */
  printf("s11r = %g\ns22r = %g\ns33r = %g\ns12r = %g\ns13r = %g\ns23r = %g\n\n",
	 a[0]+medium->pressure,a[3]+medium->pressure,
	 a[5]+medium->pressure,
	 a[1],a[2],a[4]);
  
  printf("half_space = yes\n");
  printf("check_cond_num = no\n");
  printf("print_elt_geom = no\n");
  printf("elt_geom_csys  = global\n");
  printf("null_value     = -999.000000\n\n");

  printf("end *(CONSTANTS)\n\n");

  /* user coordinate system, which we won't use, i guess */
  printf("*********************************************************\n");
  printf("*               Section 2: USER COORDINATE SYSTEMS      *\n");
  printf("*********************************************************\n");
  printf("*\n");
  printf("* (1)        (2)     (3)    (4)    (5)     (6)     (7)     (8)       (9)\n");
  printf("*name      parent    x1o    x2o    x3o    rot1    rot2    rot3    rot order\n");
  printf("*--------- --------- ------ ------ ------ ------- ------- ------- ---------\n");
  printf("f1  global  0.000000  0.000000  0.000000  90.000000  0.000000  0.000000  123\n");
  printf("end *(USER COORDINATE SYSTEMS)\n\n");

  if(medium->print_bulk_fields){
    for(dim = i=0;i<3;i++)
      if(medium->n[i] > 1)
	dim ++;
    /* observational grid */
    printf("*******************************************************************************\n");
    printf("*                         Section 3: OBSERVATION GRIDS                        *\n");
    printf("*******************************************************************************\n");
    printf("* (1)    (2) (3)      (4)        (5)        (6)     (7)    (8)    (9)    (10)   (11)   (12) (13) (14) (15)\n");
    printf("*name    dim outp  endpt csys obspt csys outp csys x1beg  x2beg  x3beg  x1end  x2end  x3end  nx1 nx2 nx3\n");
    printf("*------- --- ----- ---------- ---------- --------- ------ ------ ------ ------ ------ ------ --- --- ---\n");
    printf("OGrid  %i  ds  global  global  global  %g %g %g %g %g %g %i %i %i\n\n",
	   dim,
	   medium->pxmin[INT_X], medium->pxmin[INT_Y],  medium->pxmin[INT_Z],
	   medium->pxmax[INT_X], medium->pxmax[INT_Y],  medium->pxmax[INT_Z],
	   medium->n[INT_X],medium->n[INT_Y],medium->n[INT_Z]);
	   
	   

    printf("end * (OBSERVATION GRIDS)\n");
  }

  /* vertices and elements */
  printf("************************************************************\n");
  printf("**          SECTION 4: OBJECTS/ELEMENTS/VERTICES          **\n");
  printf("**          (v = vertex, o = object, e = element)         **\n");
  printf("************************************************************\n");
#ifndef ALLOW_NON_3DQUAD_GEOM
  fprintf(stderr,"error, non quad geometry has to be allowed, set flags during compilation\n");
  exit(-1);
#endif
  /* check for special cases */
  for(i=0;i<medium->nrflt;i++)
    if(fabs(fault[i].dip - 90) < EPS_COMP_PREC)
      fprintf(stderr,"%s: WARNING: fault patch %i: dip of %g will freak out poly3d\n",
	      argv[0],i,fault[i].dip );

  /* loop through fault groups */
  for(ecnt=i=0;i<medium->nrgrp;i++){
    /* 
       use each fault group as an object: b (slip vector) or t
       (tractions) to be calculated
    */
    printf("* this is fault group %3i\n",i);
    printf("o  \"fg%03i\"  b  global\n\n",i);
    for(j=0;j<medium->nrflt;j++)
      if(fault[j].group == i){
#ifdef ALLOW_NON_3DQUAD_GEOM
	switch(fault[j].type){
	case OKADA_PATCH:
	  /* 
	     output of rectangular patch 
	  */
	  calculate_vertices(corner,(fault+j),&lloc,&wloc);
	  if(!tri){		/* quad */
	    /* vertices 
	       v name csys x1 x2 x3
	    */
	    ecnt++;
	    printf("* element %6i from patch %6i group %3i\n",
		   ecnt,j,i);
	    printf("v g%03if%06ia global %24.15e %24.15e %24.15e\n",
		   i,j,corner[0*3+INT_X],corner[0*3+INT_Y],corner[0*3+INT_Z]);
	    printf("v g%03if%06ib global %24.15e %24.15e %24.15e\n",
		   i,j,corner[1*3+INT_X],corner[1*3+INT_Y],corner[1*3+INT_Z]);
	    printf("v g%03if%06ic global %24.15e %24.15e %24.15e\n",
		   i,j,corner[2*3+INT_X],corner[2*3+INT_Y],corner[2*3+INT_Z]);
	    printf("v g%03if%06id global %24.15e %24.15e %24.15e\n",
		   i,j,corner[3*3+INT_X],corner[3*3+INT_Y],corner[3*3+INT_Z]);
	    
	    /* elements
	       e nrnd csys btype bval1 bval2 bval3 nnd1 nnd2 nnd3 ...
	       
	       btype can be b (for slip) or t (for traction)
	       
	       element local system in poly3d is 
	       down-dip, strike, normal
	    */
	    sprintf(bcstring,"%1s%1s%1s %15.8e %15.8e %15.8e",
		    (slip_type_bc(fault[i].mode[DIP]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[STRIKE]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[NORMAL]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[DIP]))?(-fault[i].u[DIP]):(-fault[i].s[DIP]),
		    (slip_type_bc(fault[i].mode[STRIKE]))?(fault[i].u[STRIKE]):(fault[i].s[STRIKE]),
		    (slip_type_bc(fault[i].mode[NORMAL]))?(fault[i].u[NORMAL]):(fault[i].s[NORMAL]));
	    printf("e 4 elocal %s g%03if%06ia g%03if%06ib g%03if%06ic g%03if%06id\n\n",
		   bcstring,i,j,i,j,i,j,i,j);
	  }else{		/* use two triangles */
	    ecnt++;
	    printf("* element %6i from patch %6i group %3i, triangle one\n",
		   ecnt,j,i);
	    printf("v g%03if%06ia global %24.15e %24.15e %24.15e\n",
		   i,j,corner[0*3+INT_X],corner[0*3+INT_Y],corner[0*3+INT_Z]);
	    printf("v g%03if%06ib global %24.15e %24.15e %24.15e\n",
		   i,j,corner[1*3+INT_X],corner[1*3+INT_Y],corner[1*3+INT_Z]);
	    printf("v g%03if%06ic global %24.15e %24.15e %24.15e\n",
		   i,j,corner[2*3+INT_X],corner[2*3+INT_Y],corner[2*3+INT_Z]);
	    sprintf(bcstring,"%1s%1s%1s %15.8e %15.8e %15.8e",
		    (slip_type_bc(fault[i].mode[DIP]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[STRIKE]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[NORMAL]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[DIP]))?(-fault[i].u[DIP]):(-fault[i].s[DIP]),
		    (slip_type_bc(fault[i].mode[STRIKE]))?(fault[i].u[STRIKE]):(fault[i].s[STRIKE]),
		    (slip_type_bc(fault[i].mode[NORMAL]))?(fault[i].u[NORMAL]):(fault[i].s[NORMAL]));
	    printf("e 3 elocal %s g%03if%06ia g%03if%06ib g%03if%06ic\n\n",
		   bcstring,i,j,i,j,i,j);

	    ecnt++;
	    printf("* element %6i from patch %6i group %3i, triangle two\n",
		   ecnt,j,i);
	    printf("v g%03if%06ia global %24.15e %24.15e %24.15e\n",
		   i,j,corner[0*3+INT_X],corner[0*3+INT_Y],corner[0*3+INT_Z]);
	    printf("v g%03if%06ic global %24.15e %24.15e %24.15e\n",
		   i,j,corner[2*3+INT_X],corner[2*3+INT_Y],corner[2*3+INT_Z]);
	    printf("v g%03if%06id global %24.15e %24.15e %24.15e\n",
		   i,j,corner[3*3+INT_X],corner[3*3+INT_Y],corner[3*3+INT_Z]);
	    sprintf(bcstring,"%1s%1s%1s %15.8e %15.8e %15.8e",
		    (slip_type_bc(fault[i].mode[DIP]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[STRIKE]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[NORMAL]))?("b"):("t"),
		    (slip_type_bc(fault[i].mode[DIP]))?(-fault[i].u[DIP]):(-fault[i].s[DIP]),
		    (slip_type_bc(fault[i].mode[STRIKE]))?(fault[i].u[STRIKE]):(fault[i].s[STRIKE]),
		    (slip_type_bc(fault[i].mode[NORMAL]))?(fault[i].u[NORMAL]):(fault[i].s[NORMAL]));
	    printf("e 3 elocal %s g%03if%06ia g%03if%06ic g%03if%06id\n\n",
		   bcstring,i,j,i,j,i,j);
	

	  }
	  break;
	default:
	  fprintf(stderr,"%s: fault type %i not implemented yet\n",
		  argv[0],fault[j].type);
	  exit(-1);
	}
#else

	calculate_vertices(corner,(fault+j),&lloc,&wloc);
	/* vertices 
	   v name csys x1 x2 x3
	*/
	ecnt++;
	printf("* element %6i from patch %6i group %3i\n",
	       ecnt,j,i);
	printf("v g%03if%06ia global %24.15e %24.15e %24.15e\n",
	       i,j,corner[0*3+INT_X],corner[0*3+INT_Y],corner[0*3+INT_Z]);
	printf("v g%03if%06ib global %24.15e %24.15e %24.15e\n",
	       i,j,corner[1*3+INT_X],corner[1*3+INT_Y],corner[1*3+INT_Z]);
	printf("v g%03if%06ic global %24.15e %24.15e %24.15e\n",
	       i,j,corner[2*3+INT_X],corner[2*3+INT_Y],corner[2*3+INT_Z]);
	printf("v g%03if%06id global %24.15e %24.15e %24.15e\n",
	       i,j,corner[3*3+INT_X],corner[3*3+INT_Y],corner[3*3+INT_Z]);
	
	/* elements
	   e nrnd csys btype bval1 bval2 bval3 nnd1 nnd2 nnd3 ...
	   
	   btype can be b (for slip) or t (for traction)
	   
	   element local system in poly3d is 
	   down-dip, strike, normal
	*/
	sprintf(bcstring,"%1s%1s%1s %15.8e %15.8e %15.8e",
		(slip_type_bc(fault[i].mode[DIP]))?("b"):("t"),
		(slip_type_bc(fault[i].mode[STRIKE]))?("b"):("t"),
		(slip_type_bc(fault[i].mode[NORMAL]))?("b"):("t"),
		(slip_type_bc(fault[i].mode[DIP]))?(-fault[i].u[DIP]):(-fault[i].s[DIP]),
		(slip_type_bc(fault[i].mode[STRIKE]))?(fault[i].u[STRIKE]):(fault[i].s[STRIKE]),
		(slip_type_bc(fault[i].mode[NORMAL]))?(fault[i].u[NORMAL]):(fault[i].s[NORMAL]));
	printf("e 4 elocal %s g%03if%06ia g%03if%06ib g%03if%06ic g%03if%06id\n\n",
	       bcstring,i,j,i,j,i,j,i,j);
#endif
      }
   
  }
  printf("end * (VERTICES,OBJECTS,ELEMENTS)\n\n");

  exit(0);
}

