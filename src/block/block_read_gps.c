
/*
  part of blockinvert

  
  read in GPS velocities, assign blocks, block codes and 

  set up a general projection for all observationa sites
  
  also assigns mod->block[].fixed array (boolean for blocks fixed with respect
  to rigid motions, and the mod->block[].center array which holds the
  average coordinates of each block (block center of gravity) in Cartesian

  also returns the RMS of the velocity values

  constrain_euler will be set to 2 (lock last block) by some remove_net_trans settings

  $Id: block_read_gps.c,v 1.2 2003/12/23 04:05:02 becker Exp $
    
*/
#include "interact.h"
#include "blockinvert.h"

#ifndef DEBUG
#define DEBUG
#endif


void read_gps_velocities(struct bmd *mod,struct prj *projection,
			 COMP_PRECISION *velrms,char **argv,
			 my_boolean verbose, 
			 int *remove_net_trans, int *constrain_euler)
{
  int i,j,loc_bcode,npdim;
  size_t isize;
  COMP_PRECISION xcmean[3],maxlat,minlat,latrange;
  my_boolean in_fixed_block, remove_global_net_rotation;
  FILE *in,*in2;
#ifdef BLOCK_SPHERICAL
  int n,np3d;
  COMP_PRECISION *aloc,*net_block,*loc_vc;
#endif
  mod->fblock_nsites = 0;
  static COMP_PRECISION gfac = BLOCK_GFAC;
  /* if remove_net_trans is set to -2 on entry,
     will determine the global net rotation

     if remove_global_net_rotation is set to TRUE, 
     will remove this rotation from the input velocities

  */
  if(*remove_net_trans == -2 )
    remove_global_net_rotation = TRUE;
  else
    remove_global_net_rotation = FALSE;
  
  //
  // initialize 
  mod->nrgp = mod->nrb = 0;
  //
  xcmean[INT_X]=xcmean[INT_Y]=xcmean[INT_Z]=0.0;
  minlat = FLT_MAX;maxlat = FLT_MIN;
  //
  // input file
  in = fopen(GPS_VEL_FILE,"r");
  if(!in){
    fprintf(stderr,"read_gps_velocities: WARNING: no GPS velocity file %s found\n",
	    GPS_VEL_FILE);
    return;
  }else{
    fprintf(stderr,"read_gps_velocities: expecting velocities in mm/yr from %s\n",
	    GPS_VEL_FILE);
  }
  //
  // node dependent arrays
  //
  my_ivecrealloc(&mod->bcode,BLOCK_DIM,"block_read_gps: 1a");
  my_vecrealloc(&mod->gx,BLOCK_DIM,"block_read_gps: 1b");
  my_vecrealloc(&mod->gcx,3,"block_read_gps: 1c");
  my_vecrealloc(&mod->v,BLOCK_DIM,"block_read_gps: 1d");
  my_vecrealloc(&mod->sigv,BLOCK_DIM,"block_read_gps: 1e");
  my_vecrealloc(&mod->rho,1,"block_read_gps: 1f");
  /* 
     block dependent 
  */
  init_blocks(mod,(int)MAX_NR_BLOCK);
  /*
    
    expecting 
    
    lon lat ve[mm/yr] vn[mm/yr] sigve[mm/yr] sigvn[mm/yr] rho_corr bcode
    
    format
    
    on inbput, bcode has to increase from 1 if bcode is negative, will
    fix the rigid velocities of this plate (internally, bcode will go
    from 0 ... nrb-1)

    rho_corr just gets passed through
    
  */
  npdim=0;
  /* 
     start velocity input loop
  */
  while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %i",
	       (mod->gx+npdim+INT_X),(mod->gx+npdim+INT_Y), /* lon lat in degrees */
	       (mod->v+npdim+INT_X),(mod->v+npdim+INT_Y), /* ve vn in mm/yr */
	       (mod->sigv+npdim+INT_X),(mod->sigv+npdim+INT_Y),	/* uncertainties */
	       (mod->rho+mod->nrgp),&loc_bcode) == 8){ /* rho (correlation) and 
							  block code */
    /*

      geographical coordinates. convert to cartesian to find
      mean, and for possible use in spherical calculation
      
    */
    lonlat2xyz_deg(mod->gx[npdim+INT_X],mod->gx[npdim+INT_Y],(mod->gcx+3*mod->nrgp)); 
    add_b_to_a_vector_3d(xcmean,(mod->gcx+3*mod->nrgp)); /* add to mean */
    if(mod->gx[npdim+INT_Y] > maxlat)
      maxlat = mod->gx[npdim+INT_Y];
    if(mod->gx[npdim+INT_Y] < minlat)
      minlat = mod->gx[npdim+INT_Y];
    /*
      
      sigma uncertainties
      
    */
    if((fabs(mod->sigv[npdim+INT_X]) < EPS_COMP_PREC)||
       (fabs(mod->sigv[npdim+INT_Y]) < EPS_COMP_PREC)){
      fprintf(stderr,"%s: read_gps_velocities: error: at least one sigma value is zero\n",
	      argv[0]);
      exit(-1);
    }
    /*
      
      deal with block codes. we use two of them for x and y
      components, this is redundant but simpler for the loops below
      
    */
    if(loc_bcode == 0){
      fprintf(stderr,"%s: read_gps_velocities: error, codes have to be != 0\n",
	      argv[0]);
      exit(-1);
    }else if(loc_bcode < 0){	/* the block around this point is
				   fixed for rigid block velocities */
      in_fixed_block = TRUE;
      loc_bcode = -loc_bcode;
    }else{
      in_fixed_block = FALSE;
    }
    if(loc_bcode > mod->nrb){
      if(loc_bcode >= MAX_NR_BLOCK){
	fprintf(stderr,"%s: read_gps_velocities: too many blocks (%i), increase MAX_NR_BLOCK (%i)\n",
		argv[0],loc_bcode,MAX_NR_BLOCK);
	exit(-1);
      }
      mod->nrb = loc_bcode;
    }
    /*

      internally, we have codes go from 0 .. nrcode-1 for assignments
      
    */
    loc_bcode--;
    /* 
       assign the block code to all bcode entries for this node
    */
    for(i=0;i < BLOCK_DIM;i++)
      mod->bcode[npdim + i] = loc_bcode; 
    mod->block[loc_bcode].pcnt++;
    /* this block is fixed */
    if(in_fixed_block){
      if(!mod->block[loc_bcode].fixed){
	mod->block[loc_bcode].fixed = TRUE;
	fprintf(stderr,"%s: WARNING: node %i lead to locked block %i\n",argv[0],
		mod->nrgp+1,loc_bcode+1);
      }
    }
    /*

      increment and allocate more space

    */
    mod->nrgp++;
    npdim += BLOCK_DIM;
    isize = BLOCK_DIM*(mod->nrgp+1);
    my_vecrealloc(&mod->v, isize,"v");my_vecrealloc(&mod->sigv,isize,"sigv");
    my_vecrealloc(&mod->gx,isize,"gx");my_vecrealloc(&mod->rho,mod->nrgp+1,"rho");
    my_vecrealloc(&mod->gcx,(mod->nrgp+1)*3,"gcx");
    my_ivecrealloc(&mod->bcode,(mod->nrgp+1)*BLOCK_DIM,"bcode");
    /* 
       end velocity data loop, now we know the number of blocks and number
       of GPS velocities
    */
  }
  fclose(in);	
  if(npdim != (mod->nrgp * BLOCK_DIM)){
    fprintf(stderr,"read_gps_velocities: assignment error, npdim: %i np * BLOCK_DIM: %i\n",
	    mod->nrgp,BLOCK_DIM);
    exit(-1);
  }
  /* determine RMS velocity on input */
  *velrms = rms(mod->v,npdim);
  /* 
     reallocate block structure 
  */
  init_blocks(mod,mod->nrb);
  /* 
     report to stdout
  */
  fprintf(stderr,"%s: read_gps_velocities: read %i original vel. (absmean/rms: %g(%g)/%g(%g) and %i blocks)\n",
	  argv[0],mod->nrgp,
	  mean_abs(mod->v,1,npdim),mean_abs(mod->sigv,1,npdim),
	  *velrms,rms(mod->sigv,npdim),mod->nrb);
  //
  // check if more than one velocity observation was assigned to each
  // block
  //
  for(i=j=0;i < mod->nrb;i++){
    if(mod->block[i].pcnt <= 1){
      fprintf(stderr,"%s: WARNING: block code %i was assigned to only %i velocities\n",
	      argv[0],i+1,mod->block[i].pcnt);
    }else
      if(verbose)
	fprintf(stderr,"nrp(%i): %i ",i+1,mod->block[i].pcnt);
    if(mod->block[i].fixed)
      j++;			/* count the fixed blocks */
  }
  if(verbose)
    fprintf(stderr,"\n");
  if(j)
    fprintf(stderr,"%s: read_gps_velocities: WARNING: %i blocks are fixed for rigid motions based on block codes\n",
	    argv[0],j);
  if(mod->nrgp){
    my_vecrealloc(&mod->gpx,BLOCK_DIM * mod->nrgp,"rgps");
    /*
      
    determine geographical center of data and use as projection center
    also determine standard parallels to go from geographic to plane
    approximation
    
    */
    normalize_3d(xcmean);
    xyz2lonlat_deg(xcmean,&projection->clon,&projection->clat);
    latrange = maxlat - minlat;
    projection->lat1 = minlat + (1./3.) * latrange;
    projection->lat2 = minlat + (2./3.) * latrange;
    if(verbose)
      fprintf(stderr,"%s: read_gps_velocities: using projection type %i: plon: %g plat: %g azi: %g lat1: %g lat2: %g\n",
	      argv[0],projection->type,projection->clon,
	      projection->clat,projection->azi,
	      projection->lat1,projection->lat2);
    //
    // convert observational site location to projected coordinates,
    // sums up contributions for the estimate of the block centroids,
    // and creates the average
    //
    project_gps_coordinates(mod->gx,mod->gpx,mod->nrgp,(int)BLOCK_DIM,
			    mod->bcode,mod->block,mod->nrb,projection,
			    argv);
    
  } /* end loop for nrgp != 0 */
  /* 
     deal with changes in reference frame
  */
  /* 
     adjust to remove net translation of last block 
  */
  switch(*remove_net_trans){
  case -3:			/* remove the net-rotation of the last block based on all sites
				   in this block */
    *remove_net_trans = mod->nrb-1;
    mod->select_fblock_sites = FALSE;
    //*constrain_euler = 2;	/* lock last block */
    break;
  case -4:
    *remove_net_trans = mod->nrb-1;
    mod->select_fblock_sites = TRUE;
    //*constrain_euler = 2;	/* lock last block */
    break;
  default:
    mod->select_fblock_sites = FALSE;
    break;
  }
  if(mod->select_fblock_sites){
    /* select the sites */
    in2=fopen(RIGIDBLOCKSITES_FILE,"r");
    if(!in2){			/* default sites */
      fprintf(stderr,"read_gps_velocities: constraining default sites for ref frame\n");
      mod->fblock_nsites = 3;
      my_vecrealloc(&mod->fblock_sites,BLOCK_DIM*mod->fblock_nsites,"read_gps");
      mod->fblock_sites[0*BLOCK_DIM+0] = 245.396; mod->fblock_sites[0*BLOCK_DIM+1] =  34.8066; /* NEED */
      mod->fblock_sites[1*BLOCK_DIM+0] = 244.673; mod->fblock_sites[1*BLOCK_DIM+1] =  34.8063; /* 0809 */
      mod->fblock_sites[2*BLOCK_DIM+0] = 244.577; mod->fblock_sites[2*BLOCK_DIM+1] = 35.5410;
    }else{			/* read form file */
      fprintf(stderr,"read_gps_velocities: reading individual sites for ref frame from %s\n",
	      RIGIDBLOCKSITES_FILE);
      mod->fblock_nsites = 0;
      my_vecrealloc(&mod->fblock_sites,BLOCK_DIM*(mod->fblock_nsites+1),"read_gps");
      while(fscanf(in2,TWO_CP_FORMAT,
		   &mod->fblock_sites[mod->fblock_nsites*BLOCK_DIM+0],
		   &mod->fblock_sites[mod->fblock_nsites*BLOCK_DIM+1])==2){
	mod->fblock_nsites++;
	my_vecrealloc(&mod->fblock_sites,BLOCK_DIM*(mod->fblock_nsites+1),"read_gps");
      }
      fclose(in2);
    }
    for(i=0;i<mod->fblock_nsites;i++)
      fprintf(stderr,"read_gps_velocities: sites %2i for ref frame: %11g, %11g\n",
	      i+1,mod->fblock_sites[i*BLOCK_DIM+0],mod->fblock_sites[i*BLOCK_DIM+1]);
  }
  mod->changed_reference_frame = FALSE;	/* this might get reset below */
  mod->omega_corr[INT_X] = mod->omega_corr[INT_Y] = mod->omega_corr[INT_Z] = 0.0;
  if(*remove_net_trans >= 0){		
    /* 
       
    FIND NET ROTATION FOR A SINGLE BLOCK AND REMOVE IT FROM THE
    VELOCITINT_Y FIELD
    
    */
    if(*remove_net_trans >= mod->nrb){
      fprintf(stderr,"read_gps_velocities: ERROR: rem. net trans. block %i out of nrb: %i\n",
	      *remove_net_trans,mod->nrb);
      exit(-1);
    }
#ifndef BLOCK_SPHERICAL
    fprintf(stderr,"%s: read_gps_velocities: ERROR, non-spherical not working anymore\n",
	    argv[0]);
    exit(-1);
#else
    /* 
       find best fit block motion for a single block from input
       velocities this will initialize the cartesian velocities, if
       not init already 
    */
    find_spherical_rotation(mod,*remove_net_trans,
			    mod->omega_corr,FALSE,TRUE,FALSE, /* orig vel, verbose, non-global */
			    TRUE,mod->select_fblock_sites, /* compute misfit */
			    mod->fblock_sites,mod->fblock_nsites);
    fprintf(stderr,"%s: read_gps_velocities: removing best-fit omega: (%g, %g, %g) (deg/Myr) from block %c\n",
	    argv[0],mod->omega_corr[INT_X]/gfac,
	    mod->omega_corr[INT_Y]/gfac,
	    mod->omega_corr[INT_Z]/gfac,
	    bname(*remove_net_trans));
    /* 

    compute global velocities from best fit block motion 

    */
    aloc = NULL;
    /* get the modified A matrix with all block codes set to
       remove_net_trans as a block */
    assemble_block_a(&aloc,mod->gpx,mod->bcode,mod->nrgp,
		     mod->nrb,mod->block,*remove_net_trans,mod);
    /*

    multiply with rotation vector for all blocks

    */
    n = mod->nrb * BLOCK_NBASE;
    np3d = mod->nrgp * 3;
    /* create global rotation vector, only entries are for 
       remove_net_trans block  */
    my_vecalloc(&net_block,n,"net_block");
    for(i=0;i<n;i++)		/* init with zeroes */
      net_block[i] = 0.0;
    /* 
       assign to global rotation vector 
    */
    a_equals_b_vector((net_block+(*remove_net_trans)*BLOCK_NBASE),
		      mod->omega_corr,3);
    /* 
       create global rotation velocities 
    */
    my_vecalloc(&loc_vc,np3d,"net_block");
    /* 
       calculate velocities in Cartesian 
    */
    calc_Ax_ftn(aloc,np3d,n,net_block,loc_vc);
    /* 
       remove the cartesian velocity field loc_vc from the GPS data
       velocity field, both in cartesian and polar systems and write
       the corrected and correction velocity fields to files. this
       overwrites the input velocities
    */
    remove_cvel_from_vel(mod,loc_vc,np3d);
    /* 
       recompute RMS
    */
    *velrms = rms(mod->v,mod->nrgp * BLOCK_DIM);
    fprintf(stderr,"%s: read_gps_velocities: modified rms of corrected field: %g\n",
	    argv[0],*velrms);
    free(aloc);free(loc_vc);free(net_block);
#endif
  }else{
    /* 

    DETERMINE THE GLOBAL NET ROTATION

    and remove it from velocity data, if remove_global_net_rotation is
    true
    
    */
#ifdef BLOCK_SPHERICAL
    fprintf(stderr,"%s: read_gps_velocities: global, best-fit rotation is\n",
	    argv[0]);
    find_spherical_rotation(mod,-1,mod->omega_corr,FALSE,TRUE, /* global, orig vel, verbose */
			    remove_global_net_rotation,
			    FALSE,FALSE,	/* no misfit, no fblock sites */
			    mod->fblock_sites,mod->fblock_nsites);
    if(!remove_global_net_rotation){
      /* reset omega_corr */
      mod->omega_corr[INT_X] = mod->omega_corr[INT_Y] =
	mod->omega_corr[INT_Z] = 0.0;
    }
#else
    fprintf(stderr,"%s: read_gps_velocities: error, not implemented for cartesian\n",
	    argv[0]);
    exit(-1);
#endif
  }
  /* 

     find the spherical Euler pole of all input blocks individually
     after all corrections are done.

  */
  fprintf(stderr,"%s: read_gps_velocities: individual block net rotations:\n",
	  argv[0]);
  /* 
     find the best fit rigid motion at input 
  */
  for(i=0;i< mod->nrb;i++)
    find_spherical_rotation(mod,i,mod->block[i].xrigid,
			    FALSE,TRUE,FALSE, /* orig vel, verbose, non remove global */
			    FALSE,FALSE, /* no misfits, no fblock sites */
			    mod->fblock_sites,mod->fblock_nsites);
}
/*
  
  project the gps velocity coordinates to mercator projection psace
  and average to block center array

  input: x[n*BLOCK_DIM] px[n*BLOCK_DIM]


*/
void project_gps_coordinates(COMP_PRECISION *x, COMP_PRECISION *px,
			     int n, int dim, int *bcode,
			     struct bck *block,int nrb,
			     struct prj *projection,char **argv)
{
  COMP_PRECISION gin[3],gout[3],dummy=0;
  int i,j;
  if(projection->azi != 90){
    fprintf(stderr,"%s: read_gps_velocities: azimuth != 90 not implemented yet (%g)\n",
	    argv[0],projection->azi);
    exit(-1);
  }
  if((projection->type != PROJECT_AZI)&&
     (projection->type != OMERC_AZI)&&
     (projection->type != LCONFORM)){
    fprintf(stderr,"%s: read_gps_velocities: projection type %i not implemented\n",
	    argv[0],projection->type);
    exit(-1);
  }
  fprintf(stderr,"%s: read_gps_velocities: projection: %i %g %g %g %g %g\n",
	  argv[0],projection->type,projection->clon,
	  projection->clat,projection->azi,
	  projection->lat1,projection->lat2);  
  gin[INT_Z] = 0.0;
  for(i=j=0;i < n;i++,j+=dim){
    /* 
       project lon lat to projected coordinates
       
       (routine assumes a Z - value)
       
       gx: geographic coordinates
       gpx: local plane projection
       
       also assign mean projected coordinates for the blocks
       
    */
    a_equals_b_vector(gin,(x+j),dim); /* have to do this  since gin[3] */
    /* 
       this is projected into the local, oblique mercator system  
    */
    geoproject(gin,gout,projection->type,projection->clon,
	       projection->clat,projection->azi,dummy,dummy,
	       projection->lat1,projection->lat2,(int)FALSE);
    a_equals_b_vector((px+j),gout,dim);	/* same here, gout[3] */
    /* 
       
    average the block centroid coordinate in cartesian

    */				/* gout gets reused */
    lonlat2xyz(gin[INT_X],gin[INT_Y],gout);
    add_b_to_a_vector_3d(block[bcode[j]].center,gout);
  }
  for(i=0;i < nrb;i++)
    normalize_3d(block[i].center);
}
/*
  
read in block centroids from a file with simple
lon lat format

*/
void block_read_centroids(COMP_PRECISION **cx, int *nrb)
{
  FILE *in;
  int onrb,i;
  onrb = *nrb;
  in=myopen(BLOCK_CENTROID_FILE,"r");
  i=0;
  *cx=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*BLOCK_DIM);
  while(fscanf(in,"%lf %lf",(*cx+BLOCK_DIM*i+INT_X),(*cx+BLOCK_DIM*i+INT_Y))==2){
    i++;
    *cx=(COMP_PRECISION *)
      realloc(*cx,sizeof(COMP_PRECISION)*BLOCK_DIM*(i+1));
    if(!*cx)
      MEMERROR("");
  }
  fclose(in);
  *nrb = i;    
  if(*nrb != onrb){
    fprintf(stderr,"block_read_centroids: mismatch in block number from GPS and centroids\n");
    exit(-1);
  }else{
    fprintf(stderr,"block_read_centroids: read %i centroid locations\n",
	    *nrb);
  }
}
/*  
    
find the spherical rotation vector (Euler vector) for a block with
code "code". if use_model_vel is set to TRUE, will use model
predictions, else the input GPS velocity data

if code is set to -1, will determine the best fit rotation for all
blocks. in this case, if remove_global is set, the velocities will be
corrected for this

*/
void find_spherical_rotation(struct bmd *mod,int code, 
			     COMP_PRECISION *omega,
			     my_boolean use_model_vel,
			     my_boolean verbose,
			     my_boolean remove_global,
			     my_boolean compute_misfit,
			     my_boolean select_fblock_sites,
			     COMP_PRECISION *fblock_sites,
			     int fblock_nsites)
{
  int i,j,k,m,nconstrained;
  static int n = 3;	
  COMP_PRECISION *a,*b,*aorig,*borig,*bbestfit,
    aloc[3][3],*p1=NULL,wmax,*a2,tmp,lon,lat,
    wcutoff=1e-8,mag,*loc_vc,vp[3];
  char cstring[2];
#ifdef DEBUG
  FILE *out;
#endif 
  static COMP_PRECISION gfac = BLOCK_GFAC;
  if(!mod->pbase_init)		/* initialize basis vectors */
    init_gps_pbase(mod);
  if(!use_model_vel){
    if(!mod->cvel_init)		/* initialize cartesian 
				   GPS velocities */
      init_cart_gps_velsig(mod,FALSE);
  }else{
    /* assign cartesian model velocities */
    init_cart_gps_velsig(mod,TRUE);
  }
  if(remove_global){
    if(code >= 0){
      fprintf(stderr,"find_spherical_rotation: logic error, can only remove global if code < 0\n");
      exit(-1);
    }
    if(use_model_vel){
      fprintf(stderr,"find_spherical_rotation: logic error, remove global only makes sense for data\n");
      exit(-1);
    }
    /* make sure that we are verbose in this case */
    verbose = TRUE;
  }
  
  a=b=aorig=borig=NULL;
  /* 

  this will produce a matrix with m by 3 entries
  where m is the number of points found for block 
  "code"

  */
  nconstrained=0;
  for(m=i=0;i < mod->nrgp;i++){		/* loop through all points */
    if((code < 0)||		/* either select all nodes */
       ((!select_fblock_sites)&&
       (mod->bcode[i*BLOCK_DIM] == code))|| /* or those in a block */
       (select_fblock_sites && 
	node_in_geo_list(fblock_sites,fblock_nsites,
			 (mod->gx+BLOCK_DIM*i)))){	
      /* 
	 select this velocity is point is in block or if all points
	 should be selected
	 |    0  r_z  -r_y |
	 | -r_z    0   r_x |
	 |  r_y -r_x     0 |

      */
      nconstrained++;
      aloc[0][0]=              0.0;aloc[0][1]=  mod->gcx[i*3+INT_Z];aloc[0][2]= -mod->gcx[i*3+INT_Y];
      aloc[1][0]= -mod->gcx[i*3+INT_Z];aloc[1][1]=              0.0;aloc[1][2]=  mod->gcx[i*3+INT_X];
      aloc[2][0]=  mod->gcx[i*3+INT_Y];aloc[2][1]= -mod->gcx[i*3+INT_X];aloc[2][2]=              0.0;
      //print_matrix_C(&aloc[0][0],3,3,stderr,FALSE);fprintf(stderr,"\n");
      my_vecrealloc(&a,(m+3)*n,"a loc");my_vecrealloc(&b,(m+3),"b loc");
      my_vecrealloc(&aorig,(m+3)*n,"aorig loc");my_vecrealloc(&borig,(m+3),"borig loc");
      for(j=0;j < 3;j++){	/* assign to A and b vectors */
	tmp = 1.0 / mod->sigvc[i*3+j];
	if(use_model_vel){	/* switch here */
	  borig[m+j] = mod->vmodc[i*3+j]; /* original */
	  b[m+j] = borig[m+j] * tmp; /* scaled */
	}else{
	  borig[m+j] = mod->vc[i*3+j];
	  b[m+j] = borig[m+j] * tmp;
	}
	for(k=0;k < n;k++){
	  aorig[(m+j)*n+k] = aloc[j][k]; /* original */
	  a[(m+j)*n+k] = aorig[(m+j)*n+k] * tmp; /* scaled */
	}
      }
      m += 3;
    }
  } /* end nrgp loop */
  if((select_fblock_sites) && 
     (nconstrained != fblock_nsites)){
    fprintf(stderr,"find_spherical_rotation: error: not all %i sites for ref frame found, only %i\n",
	    fblock_nsites,nconstrained);
    exit(-1);
  }
  if(m < 6){
    fprintf(stderr,"find_spherical_rotation: error, less than two datapoints found\n");
    exit(-1);
  }else{
    /* flip the A matrix */
    my_vecalloc(&a2,m*n,"a2");
    for(i=0;i<m;i++)
      for(j=0;j<n;j++)
	a2[j*m+i] = a[i*n+j];
    /* this solves for omega BUT OVERWRITES a2 */
    svd_driver_lapack(a2,omega,b,m,n,
		      &wcutoff,0,&wmax,&p1,0,TRUE,FALSE);
    if(compute_misfit){
      /* 
	 compute and report misfit between rigid motion and velocity at
	 sites which were used to define frame
      */
      my_vecalloc(&bbestfit,m,"fsr");
      calc_Ax(aorig,m,n,omega,bbestfit);
      //print_matrix_C(aorig,m,n,stderr,FALSE);fprintf(stderr,"\n");
      //print_vector(omega,n,stderr);fprintf(stderr,"\n");
      //print_vector(bbestfit,m,stderr);fprintf(stderr,"\n");
      fprintf(stderr,"fsr: RMS misfit (elast def.) for %i sites for ref frame: %g\n",
	      nconstrained,sqrt(distance_squared(bbestfit,borig,m)/(COMP_PRECISION)m));
      free(bbestfit);
    }
    if(verbose){
      /* 
	 convert omega vector from solution 
	 to Euler pole in different units
      */
      calc_geo_euler_pole(omega,&lon,&lat,&mag);
      (code >= 0)?(sprintf(cstring,"%c",bname(code))):(sprintf(cstring,"%2i",code));
      if(use_model_vel)
	fprintf(stderr,"fsr: model w: %6.3f, %6.3f, %6.3f (deg/Myr) lon lat r: %8.3f %8.3f %6.3f (%3i pts, bn: %s) %s %s\n",
		reformat_small(omega[INT_X]/gfac),reformat_small(omega[INT_Y]/gfac),
		reformat_small(omega[INT_Z]/gfac),lon,lat,reformat_small(mag),m/3,cstring,
		(code < 0)?("global"):(""),
		(mod->changed_reference_frame)?("mod. RF"):("orig. RF"));
      else
	fprintf(stderr,"fsr: input w: %6.3f, %6.3f, %6.3f (deg/Myr) lon lat r: %8.3f %8.3f %6.3f (%3i pts, bs: %s) %s %s\n",
		reformat_small(omega[INT_X]/gfac),reformat_small(omega[INT_Y]/gfac),
		reformat_small(omega[INT_Z]/gfac),lon,lat,reformat_small(mag),m/3,cstring,
		(code < 0)?("global"):(""),
		(mod->changed_reference_frame)?("mod. RF"):("orig. RF"));
    }
    if(code < 0){
      /* 
	 calculate the best fit velocity field
	 
      */
      if((n != 3) || (m != mod->nrgp * 3)){
	fprintf(stderr,"find_spherical_rotation: error: for remove_global, m (%i) should be nrgp (%i) * 3\n",
		m,mod->nrgp);
	exit(-1);
      }
      my_vecalloc(&loc_vc,m,"loc_vc");
      /* 

	 the following is a bit clumsy, but leave for now 

      */
      /* reassign A without uncertainties 

	 |    0  r_z  -r_y |
	 | -r_z    0   r_x |
	 |  r_y -r_x     0 |
      
      */
      for(i=j=0;i<m;i+=3,j++){
	a[    i*3+0]=  0.0;            a[    i*3+1]=  mod->gcx[j*3+INT_Z];a[    i*3+2]= -mod->gcx[j*3+INT_Y];
	a[(i+1)*3+0]= -mod->gcx[j*3+INT_Z];a[(i+1)*3+1]=  0.0;            a[(i+1)*3+2]=  mod->gcx[j*3+INT_X];
	a[(i+2)*3+0]=  mod->gcx[j*3+INT_Y];a[(i+2)*3+1]= -mod->gcx[j*3+INT_X];a[(i+2)*3+2]=  0.0;
      }
      /* flip */
      for(i=0;i<m;i++)
	for(j=0;j<n;j++)
	  a2[j*m+i] = a[i*n+j];
      /* 
	 calculate velocities corresponding to net rotation 
      */
      calc_Ax_ftn(a2,m,n,omega,loc_vc);
      fprintf(stderr,"find_spherical_rotation: misfit between GPS and best-fit net rotation: %g\n",
	      sqrt(distance_squared(loc_vc,mod->vc,m)/(COMP_PRECISION)m));

      if(remove_global){
	/* 
	   remove the best-fit net rotation velocity field 
	*/
	fprintf(stderr,"find_spherical_rotation: WARNING: removing rotation component of (%g, %g, %g)\n",
		omega[INT_X]/gfac,omega[INT_Y]/gfac,omega[INT_Z]/gfac);
	remove_cvel_from_vel(mod, loc_vc, m);
#ifdef DEBUG
	/* 
	   
	output of velocity field which would be used to correct,
	in polar system
	
	*/
	/* 
	   convert to polar system 
	*/
	fprintf(stderr,"find_spherical_rotation: writing correction net rotation  to %s\n",
		VEL_GR_FILE);
	out = myopen(VEL_GR_FILE,"w");
	for(i=j=0;i < mod->nrgp;i++,j+=BLOCK_DIM){
	  cv2pv((loc_vc+i*3),vp, (mod->pbase+i*9));
	  fprintf(out,"%12g %12g %12g %12g %12g %12g %12g\n",
		  mod->gx[j+INT_X],mod->gx[j+INT_Y],
		  vp[INT_PHI],-vp[INT_THETA],mod->sigv[j+INT_X], 
		  mod->sigv[j+INT_Y],mod->rho[i]);
	}
	fclose(out);
#endif	
      }else{
	fprintf(stderr,"find_spherical_rotation: net rotation will be left untouched in GPS data.\n");
      }
      free(loc_vc);
    }
    free(a2);
  }
  free(a);free(b);free(aorig);free(borig);
}


/* 
   initialize the polar basis vectors and cartesian location 
   vectors
*/
void init_gps_pbase(struct bmd *mod)
{
  int i;
  my_vecrealloc(&mod->pbase,9*mod->nrgp,"init_gps_pbase 1");
  for(i=0;i<mod->nrgp;i++){
    // spherical basis vectors
    calculate_polar_base(mod->gx[i*BLOCK_DIM+INT_X],mod->gx[i*BLOCK_DIM+INT_Y],
			 (mod->pbase+i*9));
  }
  mod->pbase_init = TRUE;
}
/* 
   convert the spherical system velocities and uncertainties to
   cartesian. if use_model is TRUE, will use model predictions, else
   input GPS
*/
void init_cart_gps_velsig(struct bmd *mod,my_boolean use_model)
{
  int i;
  COMP_PRECISION xp[3];
  if(BLOCK_DIM != 2){
    fprintf(stderr,"init_cart_gps_velsig: this routine assumes BLOCK_DIM == 2\n");
    exit(-1);
  }
  if(!mod->pbase_init)		/* initialize basis vectors */
    init_gps_pbase(mod);
  if(use_model){
    /* 
       predicted velocities
    */
    my_vecrealloc(&mod->vmodc,3*mod->nrgp,"init_cart_gps_velsig 1");
    for(i=0;i<mod->nrgp;i++){
      /* convert the velocities and uncertainties to cartesian 
	 vectors */
      // velocities
      xp[INT_R] = 0.0;xp[INT_THETA] = -mod->vmod[i*BLOCK_DIM+INT_Y];
      xp[INT_PHI]=mod->vmod[i*BLOCK_DIM+INT_X];
      pv2cv(xp,(mod->vmodc+i*3),(mod->pbase+i*9));
    }
    mod->cvelmod_init = TRUE;
  }else{
    /* 
       
    observed velocities

    */
    my_vecrealloc(&mod->sigvc,3*mod->nrgp,"init_cart_gps_velsig 2");
    my_vecrealloc(&mod->vc,3*mod->nrgp,"init_cart_gps_velsig 1");
    for(i=0;i<mod->nrgp;i++){
      xp[INT_R] = 0.0;
      xp[INT_THETA] = -mod->v[i*BLOCK_DIM+INT_Y]; /* south */
      xp[INT_PHI]=mod->v[i*BLOCK_DIM+INT_X]; /* west */
      pv2cv(xp,(mod->vc+i*3),(mod->pbase+i*9));
      convert_sig_p2c((mod->sigv+i*BLOCK_DIM),(mod->sigvc+i*3),
		      (mod->pbase+i*9));
    }
    mod->cvel_init = TRUE;
  }
}
/* 
   
   using a stored cartesian velocity field, recalculate the polar
   component, which is defined as ve, vn
   
*/
void calculate_c2p_gps(COMP_PRECISION *vc, COMP_PRECISION *vp,
		       struct bmd *mod)
{
  COMP_PRECISION xp[3];
  int i;
  if(BLOCK_DIM != 2){
    fprintf(stderr,"calculate_c2p_gps: this routine assumes BLOCK_DIM == 2\n");
    exit(-1);
  }
  if(!mod->pbase_init)		/* initialize basis vectors */
    init_gps_pbase(mod);
  /* 
     make room for polar velocities 
  */
  my_vecrealloc(&mod->v,BLOCK_DIM * mod->nrgp,"calc_c2p");
  for(i=0;i<mod->nrgp;i++){
    cv2pv((vc+i*3),xp, (mod->pbase+i*9));
    vp[i*BLOCK_DIM+INT_X] =  xp[INT_PHI];  /* v_West */
    vp[i*BLOCK_DIM+INT_Y] = -xp[INT_THETA];	/* v_North */
  }
}
/*

convert polar (lon lat) UNCERTAINTIES  to cartesian ones

not really clear how that should properly be done

*/
void convert_sig_p2c(COMP_PRECISION *vsig,COMP_PRECISION *vsigx,
		     COMP_PRECISION *pbase)
{
  COMP_PRECISION xp[3];
  int i;
  xp[INT_THETA] = -vsig[INT_Y];
  xp[INT_PHI]   =  vsig[INT_X];
  /* this turns out to be a well behaved solution */
  xp[INT_R]     = (vsig[INT_X]+vsig[INT_Y])/2.0;
  //xp[R] = hypot(vsig[INT_X],vsig[INT_Y]); /* that's OK, too */
  pv2cv(xp,vsigx,pbase);
  for(i=0;i<3;i++){
    if(vsigx[i] < 0)
      vsigx[i] = -vsigx[i];
  }
}
/* 

remove the cartesian velocity field loc_vc from the input GPS data
velocity field, both in cartesian and polar systems. write the
corrected and correction velocity fields to files, if debuggin output
is requested. 

also assigns the vcorp array, the spherical representation of the
correction


given a cartesian velocity correction vector loc_vc (nrgp * 3),
remove this vector from the GPS data 

np3d = nrgp * 3

*/
void remove_cvel_from_vel(struct bmd *mod, 
			  COMP_PRECISION *loc_vc,
			  int np3d)
{
  /* remove global rotation velocities from cartesian */
  sub_b_from_a_vector(mod->vc,loc_vc,np3d);
  /* recalculate the polar component GPS velocities */
  calculate_c2p_gps(mod->vc,mod->v,mod);
  /* 
     assign correction velocities in polar system 
  */
  mod->changed_reference_frame = TRUE;
  my_vecrealloc(&mod->vcorp,mod->nrgp*BLOCK_DIM,"vcorp");
  /* save the polar components of the correction velocities  */
  calculate_c2p_gps(loc_vc,mod->vcorp,mod);
#ifdef DEBUG
  /* 
     output of modified velocities in polar system 
  */
  fprintf(stderr,"remove_cvel_from_vel: writing corrected velocities to %s\n",
	  VELCOR_FILE);
  print_simple_vel(mod->gx,mod->v,mod->sigv,mod->rho,mod->nrgp,
		   VELCOR_FILE);
#endif
}

char bname(int i)
{
  char bname_string[26]="ABCDEFGHIJKLMNOPQRSTUVXYZ";
  if((i<0)||(i > 25)){
    fprintf(stderr,"bname: code %i out of range\n",i);
    exit(-1);
  }
  return bname_string[i];
}

/*  
    
test if a lon/lat coordinate x is in a list

*/
my_boolean node_in_geo_list(COMP_PRECISION *fblock_sites,
			    int fblock_nsites,COMP_PRECISION *x)
{
  int i;
  if(BLOCK_DIM != 2){
    fprintf(stderr,"node_in_geo_list: error, expected BLOCK_DIM to be 2, but is %i\n",
	    BLOCK_DIM);
    exit(-1);
  }
  for(i=0;i < fblock_nsites;i++){
    if(dist_on_sphere(x[0],x[1],fblock_sites[i*BLOCK_DIM+0],
		      fblock_sites[i*BLOCK_DIM+1]) < 1e-4)
      return TRUE;
  }
  return FALSE;
}
