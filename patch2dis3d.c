#include "interact.h"
#include "properties.h"

/* 

reads in a rectangular Okada style fault geometry from "geom.in" in
interact format and a "bc.in" boundary condition file and produces
output for DIS3D

based on dis3d from y. fialko as of aug. 2003

not fully tested

thorsten becker, thwbecker@post.harvard.edu

$Id: patch2dis3d.c,v 1.10 2011/01/09 02:02:43 becker Exp $


*/
void  get_dis3d_parameters(COMP_PRECISION , COMP_PRECISION ,COMP_PRECISION , COMP_PRECISION ,
			   COMP_PRECISION ,COMP_PRECISION ,COMP_PRECISION *, COMP_PRECISION *,
			   COMP_PRECISION *, COMP_PRECISION *);

int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  FILE *in;
  int i,opcode,j,nup,nsp,*nsbc,is_friction;
  COMP_PRECISION a[6],b[6],w1,w2,xd,yd;
  my_boolean plabel;
  if(argc != 1){
    fprintf(stderr,"%s\n\t reads in interact patch format from %s\n",
	    argv[0],GEOMETRY_FILE);
    fprintf(stderr,"\t reads in interact boundary condition file from %s\n",
	    BC_FILE);
    fprintf(stderr,"\t converts to triangular geometry and writes DIS3D output file to stdout\n");
    exit(-1);
  }
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
	     FALSE,PRESSURE_DEF,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,
	     SVD_THRESHOLD,FALSE,FALSE);
  /* check output options for conformity */
  if((medium->n[INT_X] <= 1) || medium->n[INT_Y] <= 1){
    fprintf(stderr,"%s: error nx and ny have to be > 1 \n",
	    argv[0]);
    exit(-1);
  }
  if(medium->n[INT_Z]!=1){
    fprintf(stderr,"%s: error, nz has to be unity, i.e. plane output\n",
	    argv[0]);
    exit(-1);
  }
  /* determine the number of patches with slip and stress boundary conditions */
  /* generate an array with entriues for each fault with number
     of stress bcs */
  if((nsbc=(int *)calloc(medium->nrflt,sizeof(int)))==NULL)
    MEMERROR("nsbc");
  for(nup=nsp=i=0;i < medium->nrflt;i++){
    /* 
       check geometry 
    */
    if((fabs(fault[i].dip) < EPS_COMP_PREC)||
       (fabs(fault[i].dip - 90) < EPS_COMP_PREC)){
      fprintf(stderr,"%s: error: dis3d can only deal with 0 < dip < 90\n",
	      argv[0]);
      exit(-1);
    }
#ifdef ALLOW_NON_3DQUAD_GEOM
    if(fault[i].type != RECTANGULAR_PATCH){
      fprintf(stderr,"%s: can only deal with rectangular patches\n",
	      argv[0]);
      exit(-1);
    }
#endif
    for(j=0;j<3;j++){		/* count the number of stress BCs */
      if(!slip_type_bc(fault[i].mode[j]))
	nsbc[i]++;
    }
    if((nsbc[i] != 0)&&(nsbc[i] != 3)){
      fprintf(stderr,"%s: BCs have to be all stress or all slip for a single patch for dis3d\n",
	      argv[0]);
      fprintf(stderr,"%s: if you want to specify a stress BC in bc.in, specify the other stress directions\n",argv[0]);
      fprintf(stderr,"%s: to be zero. For example, for a shear-stress in strike direction of unity on all patches\n",argv[0]);
      fprintf(stderr,"%s: specify:\n",argv[0]);
      fprintf(stderr,"%s: -1 -1 -1 10 1\n",argv[0]);
      fprintf(stderr,"%s: -1 -1 -1 20 0\n",argv[0]);
      fprintf(stderr,"%s: -1 -1 -1 30 0\n",argv[0]);
      exit(-1);
    }
    if(nsbc[i] != 3)	/* at least one slip BC */
      nup++;
    if(nsbc[i] > 0)	/* at least one stress BC */
      nsp++;
  }
  fprintf(stderr,"%s: determined %i slip and %i stress bcs for %i faults\n",
	  argv[0],nup,nsp,medium->nrflt);
  /* 
     
  write Dis3D input file 

  */
  printf("START\nautomatically produced by %s\n",argv[0]);
  /* major parameters */
  printf("CPARM\n%i %i %i %i %g %g %g\n",
	 medium->n[INT_X]*medium->n[INT_Y]*medium->n[INT_Z],
	 nup,nsp,nup+nsp,POISSON_NU, SHEAR_MODULUS,
	 STATIC_MU);
  printf("ITER\n1000 0.00001 1.0 1\nREMOTE\n1\n");
  /* stress boundary conditions, only one layer */
  printf("%g %g %g %g %g %g\n",
	 a[0],a[3],a[5],a[1],a[2],a[4]);
  printf("0. 0. 0. 0. 0. 0.\n0. 0. 0. 0. 0. 0.\n0. 0. 0. 0. 0. 0. 0.\n");
  /* why? */
  printf("PRINT\n");
  /* 
     displacement BCs first 
  */
  plabel = FALSE;
  for(i=0;i < medium->nrflt;i++){
    if(nsbc[i] == 0){
      if(!plabel){
	printf("FP1\n");
	plabel=TRUE;
      }
      get_dis3d_parameters(fault[i].x[INT_X],fault[i].x[INT_Y],fault[i].x[INT_Z], 
			   (COMP_PRECISION)fault[i].w,(COMP_PRECISION)fault[i].dip,
			   (COMP_PRECISION)fault[i].strike,&xd,&yd,&w1,&w2);

      /*
	K: Dummy variable; plane number (a carryover from FP1, dis3d)
	  FP(1,J)  = H = fault half-length                              
	  FP(2,J)  = DU = down-dip distance to upper fault edge         
	  FP(3,J)  = DL =         "            lower     "              
	  FP(4,J)  = THETA = fault dip                                  
	  FP(5,J)  = PHI = fault strike                                 
	  FP(6,J)  = X1C = global X1 coord of fault coord origin        
	  FP(7,J)  = X2C = global X2          "                         
	  FP(8,J)  = SS = strike-slip                                   
	  FP(9,J)  = DS = dip-slip                                      
	  FP(10,J) = OP = opening mode displacement across dislocation  
      */
      printf("1 %11g  %11g %11g  %11g %11g  %11g %11g  %11g %11g %11g\n",
	     fault[i].l,w1,w2,fault[i].dip,
	     270 - fault[i].strike,xd,yd,
	     -fault[i].u[STRIKE],fault[i].u[DIP], /* what about the dip??? */
	     fault[i].u[NORMAL]);
    }
  }
  /* 
     stress BCs next 
  */
  plabel = FALSE;
  for(i=0;i < medium->nrflt;i++){
    if(nsbc[i] == 3){
      if(!plabel){
	printf("FP3\n");
	plabel=TRUE;
      }
      get_dis3d_parameters(fault[i].x[INT_X],fault[i].x[INT_Y],fault[i].x[INT_Z], 
			   (COMP_PRECISION)fault[i].w,(COMP_PRECISION)fault[i].dip,
			   (COMP_PRECISION)fault[i].strike,&xd,&yd,&w1,&w2);
      /*
	comments from disdatex from y. fialko as of aug 2003:
	
	K: Dummy variable; plane number (a carryover from FP1, dis3d)
	H: Fault half-length (as in dis3d)
	DU: Down-dip distance to upper fault edge (as in FP1, dis3d)
	DL:           "          lower     "      (as in FP1, dis3d)
	THETA: Plane dip (as in FP1, dis3d)
	PHI: Plane strike (as in FP1, dis3d)
	X1C: Global X1 coord of fault coord origin (as in FP1, dis3d)
	X2C: 	"   X2	            "		   (as in FP1, dis3d)
	L: # of rows of subpatches that plane is divided into
	M: # of columns of subpatches that plane is divided into (fpsmod.f 
	    divides stress-specified planes into a uniform rectangular grid of
	    subpatches based on L and M)
	MCTYPE: Mohr-Coulomb flag.  0 if element is linear and b.c.'s are to be
		accepted literally; 1 if Mohr-Coulomb and must satisfy friction
		law (shear stress <= FMU times normal stress) and non-interpenetration.
		To handle effective stress properly, specified normal stress on Mohr-
		Coulomb elements should the pore pressure.  See Rubin 1992 or thesis 
		(Ch. 5, Appendix).
	OP,DOPDS,DOPDD: Normal stress b.c., and derivatives in strike and dip directions.
	SS,DSSDS,DSSDD: Strike-slip stress b.c., and derivatives in strike and dip 
			directions.
	DS,DDSDS,DDSDD: Dip-slip stress b.c., and derivatives in strike and dip 
			directions.
        My recollection is that OP, SS, and DS refer to the stress at the plane
	    corner with the lowest natural coordinates; that is, it is one of the
	    two top corners, and is the corner in the negative PHI direction.  This
	    is of course irrelevant if there are no gradients in plane tractions.
	    This can be checked using the IFSTAT=1 or 2 option.
	Sign conventions for OP, SS, DS should be checked (By looking at surface 
	    displacements, for example).  My recollection is that negative stresses give
	    rise to plane displacements that would be reckoned negative by DIS3D; 
	    positive stresses to positive displacements.
      */
      for(j=0;j<3;j++){
	if((fault[i].mode[j] !=  STRIKE_SLIP)&&
	   (fault[i].mode[j] !=  DIP_SLIP)&&
	   (fault[i].mode[j] !=  NORMAL_SLIP)){
	  fprintf(stderr,"%s: error patch %i: dir: %i bc %i not implemented yet\n",
		  argv[0],i,j,fault[i].mode[j]);
	  exit(-1);
	}
      }
      is_friction = 0;
      printf("1 %11g  %11g %11g  %11g %11g  %11g %11g  %i %i %i\n",
	     fault[i].l,w1,w2,fault[i].dip,
	     270 - fault[i].strike,xd,yd,1,1,is_friction);
      printf("%11g %11g %11g  %11g %11g %11g  %11g %11g %11g\n",
	     fault[i].s[NORMAL],0.,0.,
	     -fault[i].s[STRIKE],0.,0.,
	     fault[i].s[DIP],0.,0.);
    }
  }
  printf("PRINT\nCO3\nLL\n");
  printf("%g %g %i %g %g %i %g\n",
	 medium->pxmin[INT_X],(medium->pxmax[INT_X]-medium->pxmin[INT_X])/
	 (COMP_PRECISION)(medium->n[INT_X]-1),
	 medium->n[INT_X],
	 medium->pxmin[INT_Y],(medium->pxmax[INT_Y]-medium->pxmin[INT_Y])/
	 (COMP_PRECISION)(medium->n[INT_Y]-1),
	 medium->n[INT_Y],
	 medium->pxmin[INT_Z]);
  printf("DISP\nPRINT\nSTRESS\nPRINT\nSTOP\n");

  exit(0);
}
/* 

calculate some of the geometrical parameters for dis3d

input: 
w: interact patch half width
x,y,z: x, y, and z (< 0) coordinates of interact patch center 
dip: dip in degrees in interact convention
strike: in deg in interact convention

output:

xd,yd: dis3d style shifted patch coordinates
w1,w2: distance to upper and lower edge of fault


*/
void get_dis3d_parameters(COMP_PRECISION x, COMP_PRECISION y,COMP_PRECISION z, COMP_PRECISION w,
			  COMP_PRECISION dip,COMP_PRECISION strike,
			  COMP_PRECISION *xd, COMP_PRECISION *yd,
			  COMP_PRECISION *w1, COMP_PRECISION *w2)
{
  COMP_PRECISION sind,xp,betarad,tmp;

  sind = sin(DEG2RADF(dip));
  if(fabs(sind)<EPS_COMP_PREC){
    fprintf(stderr,"get_dis3d_parameters: error, dip too small for w1/w2 calculation\n");
    exit(-1);
  }
  *w1 = (-z - w * sind)/sind; /* upper fault continuation edge*/
  *w2 = *w1 + 2.0 * w; /* lower fault continuation edge */

  /* surface projection of distance of fault patch center along dip to surface */
  tmp = w + *w1;
  xp = sqrt(SQUARE(tmp) - SQUARE(z));
  /* strike */
  betarad = DEG2RADF(strike - 90.0);
  /* offset of dip projected something */
  // shift coordinate from patch center to surface projection
  *xd = x + xp * sin(betarad);
  *yd = y + xp * cos(betarad);
}




