/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: output.c,v 2.60 2011/01/07 07:19:58 becker Exp $

  
  various output formats, see also input.c and matrixio.c


*/
#include "interact.h"
#include "properties.h"

/* 
   
   slip line files, ie. slip on a fault at a certain time

*/
void print_slip_line(struct med *medium,struct flt *fault)
{
  if(medium->slip_line_init)
    if(medium->time - medium->slip_line_time >= medium->slip_line_dt-EPS_COMP_PREC)
      flush_slipline(medium,fault);
}
void flush_slipline(struct med *medium,
		    struct flt *fault)
{
  int i,j;
  for(i=0;i<medium->nrgrp;i++){
    fprintf(medium->slip_line_out[i],"# time: %g group: %i\n",medium->time,i);
    for(j=0;j<medium->nrflt;j++)
      if(fault[j].group == i)
	fprintf(medium->slip_line_out[i],"%g %g %g %g %g %i\n",
		fault[j].pos[STRIKE],fault[j].pos[DIP],fault[j].u[STRIKE],
		fault[j].u[DIP],fault[j].u[NORMAL],j);
    fprintf(medium->slip_line_out[i],"\n\n");
    fflush(medium->slip_line_out[i]);
  }
  medium->slip_line_time=medium->time;
}
/* 
   print stress and displacements to individual files
   for each fault group of several, averaged patches 

*/
#define NR_VAL 8
#define SEARCH_VAL_EXTREMA

void print_fault_stress_and_slip(struct med *medium,struct flt *fault,
				 my_boolean print_now)
{
  int i,j,n[MAX_NR_FLT_FILES_DEF];
  static COMP_PRECISION print_time=0.0;
  COMP_PRECISION val[MAX_NR_FLT_FILES_DEF][NR_VAL],
    ctmp,abs_tau;
#ifdef SEARCH_VAL_EXTREMA
  COMP_PRECISION val_min[MAX_NR_FLT_FILES_DEF][NR_VAL],
    val_max[MAX_NR_FLT_FILES_DEF][NR_VAL];
#endif
  if((print_now) || (medium->time - print_time >= 
		     medium->print_interval-EPS_COMP_PREC)){
    if(medium->flt_stress_init){
      /* average, init */
      for(i=0;i<MAX_NR_FLT_FILES_DEF;i++){
	for(j=0;j<NR_VAL;j++){
	  val[i][j]=0.0;
#ifdef SEARCH_VAL_EXTREMA
	  val_min[i][j]= 1e20; val_max[i][j]=-1e20;
#endif
	}
	n[i]=0;
      }
      /* 
	 loop through faults 
	 and average over patches
      */
      for(i=0;i<medium->nrflt;i++){
	// obtain the appropriate absolute shear stress
	calc_absolute_shear_stress(&abs_tau,i,fault);
	ctmp = coulomb_stress(abs_tau,(COMP_PRECISION)fault[i].mu_s,
			      fault[i].s[NORMAL],medium->cohesion);
	// coulomb stress based on s_s and s_d
	val[fault[i].group][0] += ctmp;
#ifdef SEARCH_VAL_EXTREMA
	if(val_min[fault[i].group][0] > ctmp)
	  val_min[fault[i].group][0] = ctmp;
	if(val_max[fault[i].group][0] < ctmp)
	  val_max[fault[i].group][0] = ctmp;
#endif
	// stress drop
	val[fault[i].group][1] += 
	  ((COMP_PRECISION)(fault[i].mu_s-fault[i].mu_d)*
	   fault[i].s[NORMAL]);
	//
	// slip components
	//
	val[fault[i].group][2] += fault[i].u[STRIKE];
	val[fault[i].group][3] += fault[i].u[NORMAL];
	val[fault[i].group][4] += fault[i].u[DIP];
	//
	// stress components
	//
	val[fault[i].group][5] += fault[i].s[STRIKE];
	val[fault[i].group][6] += fault[i].s[NORMAL];
	val[fault[i].group][7] += fault[i].s[DIP];
	n[fault[i].group]++;
      }
      // output here
      for(i=0;i<medium->nrgrp;i++){
	// time
	fprintf(medium->flt_stress_out[i],"%g ",medium->time);
	
	for(j=0;j < NR_VAL;j++){
	  /* get average */
	  if(n[i])
	    val[i][j] /= (COMP_PRECISION)n[i]; /* find mean */
	  fprintf(medium->flt_stress_out[i],"%g ",val[i][j]);
	}
#ifdef SEARCH_VAL_EXTREMA
	fprintf(medium->flt_stress_out[i],"%g %g",val_min[i][0],
		val_max[i][0]);
#endif
	fprintf(medium->flt_stress_out[i],"\n");
      }
    }
    if(medium->events_init)
      fflush(medium->events_out);
    print_time=medium->time;
  }
}
/* 

   print stress and displacements on patches, sorted by fault groups

*/
void print_fault_data(char *filename,struct med *medium,struct flt *fault)
{
  int i,j,grp;
  FILE *out;
  COMP_PRECISION val[NR_VAL],abs_tau;
#ifdef ALLOW_NON_3DQUAD_GEOM
  COMP_PRECISION global_dip_rad,sin_global_dip_rad,cos_global_dip_rad;
  COMP_PRECISION gstrike[3],gdip[3],gnormal[3],slip[3],trac[3],old_rdip = FLT_MAX,old_rstrike = FLT_MAX;
#endif
  out=myopen(filename,"w");
  fprintf(stderr,
	  "print_fault_data: writing to \"%s\"\n",filename);
  fprintf(stderr,
	  "print_fault_data: format: pos_1 pos_2 area s_c m_d*s_n u_s u_n u_d s_s s_n s_d patch_nr group_nr\n");
  fprintf(out,"# fault state: note that stresses are only those due to slip and do not include background stress\n");
  fprintf(out,"%10s %10s %10s %10s %10s %13s %13s %13s %13s %13s %13s %5s %5s\n",
	  "#    pos_1/W","pos_2/W","area","s_c", "m_d*s_n", "u_s", "u_n", "u_d", 
	  "s_s", "s_n", "s_d","pach#","grp#");

  /* 
     loop through faults and sort by group
  */
  for(grp=0;grp<medium->nrgrp;grp++)
    for(i=0;i<medium->nrflt;i++){
      if(fault[i].group == grp){
	calc_absolute_shear_stress(&abs_tau,i,fault);
	val[0] = // Coulomb stress
	  coulomb_stress(abs_tau,(COMP_PRECISION)fault[i].mu_s,
			 fault[i].s[NORMAL],medium->cohesion);
	val[1] = ((COMP_PRECISION)fault[i].mu_d)*fault[i].s[NORMAL];
#ifdef ALLOW_NON_3DQUAD_GEOM
	if(fault[i].type==TRIANGULAR){
	  /* triangular element will have local stress and slip
	     rotated into global system */
	  if((fault[i].strike != old_rstrike)||(fault[i].dip != old_rdip)){
	    /* need to compute new global basis vectors */
	    fprintf(stderr,"patch %03i triangular, rotating to global strike %g dip %g (not repeating for same angles)\n",
		    i,fault[i].strike,fault[i].dip);
	    /* compute appropriate projection vectors */
	    global_dip_rad   = DEG2RADF((COMP_PRECISION)fault[i].dip);
	    my_sincos(&sin_global_dip_rad,&cos_global_dip_rad,global_dip_rad);
	    calc_quad_base_vecs(gstrike, gnormal, gdip,
				fault[i].sin_alpha, fault[i].cos_alpha,
				sin_global_dip_rad,   cos_global_dip_rad);
	    
	    old_rdip = fault[i].dip;
	    old_rstrike = fault[i].strike;
	  }
	 
	  for(j=0;j<3;j++){	/* shear components only */
	    slip[j] = fault[i].t_strike[j]*fault[i].u[STRIKE]+fault[i].t_dip[j]*fault[i].u[DIP];
	    trac[j] = fault[i].t_strike[j]*fault[i].s[STRIKE]+fault[i].t_dip[j]*fault[i].s[DIP];
	  }
	  val[2] = project_vector(slip,gstrike);
	  val[3] = fault[i].u[NORMAL];//project_vector(slip,gnormal);
	  val[4] = project_vector(slip,gdip);

	  val[5] = project_vector(trac,gstrike);
	  val[6] = fault[i].s[NORMAL];	  //val[6] = project_vector(trac,gnormal);
	  val[7] = project_vector(trac,gdip);
	  
	}else{
	  val[2] = fault[i].u[STRIKE];// slip values
	  val[3] = fault[i].u[NORMAL];
	  val[4] = fault[i].u[DIP];
	  val[5] = fault[i].s[STRIKE];// resolved stress values
	  val[6] = fault[i].s[NORMAL];
	  val[7] = fault[i].s[DIP];
	}
#else
	val[2] = fault[i].u[STRIKE];// slip values
	val[3] = fault[i].u[NORMAL];
	val[4] = fault[i].u[DIP];
	val[5] = fault[i].s[STRIKE];// resolved stress values
	val[6] = fault[i].s[NORMAL];
	val[7] = fault[i].s[DIP];
#endif
	// local position of patch in fault group
	fprintf(out,"%10.3e %10.3e %10.3e ",fault[i].pos[0],fault[i].pos[1],
		fault[i].area);
	fprintf(out,"%10.3e %10.3e ",val[0],val[1]);
	//
	// remove tiny numbers from output?
	//
	for(j=2;j<NR_VAL;j++){
	  // nah, leave those in 
	  //val[j] = reformat_small(val[j]);
	  fprintf(out,"%13.6e ",val[j]);
	}
	fprintf(out,"%5i %5i",i,grp);
	fprintf(out,"\n");
      }
    }
  fclose(out);
}

#undef NR_VAL

/* 

   print mean, standard dev and correlation lengths for stresses on a fault group
   every medium->spcorr_timeinterval
*/
void print_fault_stress_stat(FILE *out,int grp,struct med *medium,
			     struct flt *fault)
{
  int i,j,n;
  COMP_PRECISION stat[3][2]={{0,0},{0,0},{0,0}};
  static my_boolean init=FALSE,pcorr;
  static int nrbin=30;
  static COMP_PRECISION last_spcorr_time=-1.0;
  static COMP_PRECISION xcr[3],range_frac = 0.3,cr=0.1;/* crit_val, range and 
							  cutoff for spatial corr */
  if(!init){// init phase
    if((medium->nrgrp == 1) && (grp == 0))
      pcorr=TRUE;
    else
      pcorr=FALSE;
    init=TRUE;
  }
  // get simple stress field stats
  n=0;
  for(i=0;i<medium->nrflt;i++){
    if(fault[i].group == grp){
      n++;
      for(j=0;j<3;j++){
	stat[j][0] += fault[i].s[j];
	stat[j][1] += SQUARE(fault[i].s[j]);
      }
    }
  }
  if(n==0){
    fprintf(stderr,"print_fault_stress_stat: no patch found for group %i\n",
	    grp);
  }else if(n==1){
    fprintf(out,"%g %5i %g %g %g %g %g %g 0 0 0\n",
	    medium->time,grp,stat[STRIKE][0],0.0,stat[DIP][0],0.0,
	    stat[NORMAL][0],0.0);
  }else{
    fprintf(out,"%12.5e %5i ",medium->time,grp);
    for(i=0;i<3;i++){
      stat[i][1] = sqrt((((COMP_PRECISION)n) * stat[i][1] - SQUARE(stat[i][0])) /
			((((COMP_PRECISION)n)*(((COMP_PRECISION)n)-1.0))));
      stat[i][0] /= (COMP_PRECISION)n;
      fprintf(out,"%12.5e %12.5e ",stat[i][0],stat[i][1]);
    }
    // spatial correlation of stresses 
    if(pcorr){
      //
      // if more than spcorr_interval time elapsed, 
      // recalculate the spatial correlation of the stress field
      // 
      if(medium->time - last_spcorr_time > medium->spcorr_interval){
	for(i=0;i<3;i++)
	  calc_spatial_correlation(fault,medium->nrflt,(int)2,i,nrbin,
				   &medium->cr_r,&medium->cr_xr,&medium->cr_nr,
				   range_frac,cr,(xcr+i),&medium->cr_w1,
				   &medium->cr_w2,&medium->cr_w3);
	last_spcorr_time = medium->time;
      }
      fprintf(out,"%12.5e %12.5e %12.5e\n",xcr[STRIKE],xcr[DIP],xcr[NORMAL]);
    }else{
      fprintf(out,"0 0 0\n");
    }
  }
}
#undef NR_VAL
/* 

   print patch quantities for a certain fault group to geomview file

*/
void print_group_data_geom(char *filename,struct med *medium,
			   struct flt *fault,
			   int grp,int mode,
			   COMP_PRECISION fixed_range)
{
  int i,j,k;
  FILE *out;
  COMP_PRECISION val,corner[4][3],col[3],lloc,wloc;
  static COMP_PRECISION stats[4][5],range[4];
  static my_boolean init[4]={FALSE,FALSE,FALSE,FALSE};
  if(mode>3){
    fprintf(stderr,"print_group_data_geom: error: mode: %i has to be smaller than 4\n",
	    mode);
    exit(-1);
  }
  out=myopen(filename,"w");
  fprintf(stderr,"print_group_data_geom: writing type %i data in COFF format to \"%s\"\n",
	  mode,filename);
  // color QUAD file
  fprintf(out,"CQUAD\n");
  if(!init[mode]){
    //
    // check min,max, and mean values based on all faults
    //
    stats[mode][0] =  stats[mode][4]=0.0;
    stats[mode][1] =  FLT_MAX;
    stats[mode][2] = -FLT_MAX;
    for(i=0;i<medium->nrflt;i++){
      stats[mode][4] += 1.0;
      val = select_val_for_print((fault+i),mode);
      stats[mode][0] += val;
      if(val < stats[mode][1])
	stats[mode][1]=val;
      if(val > stats[mode][2])
	stats[mode][2]=val;
    }
    stats[mode][0] /= stats[mode][4];
    range[mode] = stats[mode][2] - stats[mode][1];
    fprintf(stderr,"print_group_data_geom: min/mean/max: %g/%g/%g, %g\n",
	    stats[mode][1],stats[mode][0],
	    stats[mode][2],range[mode]);
    if(range[mode] == 0.0)
      range[mode]=1.0;
    init[mode] = TRUE;
  }
  for(i=0;i<medium->nrflt;i++){
    if(fault[i].group == grp){
      // go to 0...1 range
      if(fixed_range==0.0)
	val= (select_val_for_print((fault+i),mode) -
	      stats[mode][1])/range[mode];
      else
	val = (select_val_for_print((fault+i),mode) -
	       stats[mode][1])/fixed_range;
      if(val < 0)
	val=0.0;
      if(val > 1.0)
	val=1.0;
      // select color
      col[0]=val;col[1]=0.;col[2]=1.0-val;
      // output
      calculate_corners(corner,(fault+i),&lloc,&wloc);
      for(j=0;j<4;j++){
	for(k=0;k<3;k++)
	  fprintf(out,"%g ",corner[j][k]/CHAR_FAULT_DIM);
	fprintf(out," %g %g %g 1.0\n",
		col[0],col[1],col[2]);
      }
    }
  }
  fclose(out);
}
//
// choose a certain fault property for output
//
COMP_PRECISION select_val_for_print(struct flt *fault,int mode)
{
  switch(mode){// length of slip vector
  case 0:// slip vector length
    return norm_3d(fault->u);
    break;
  case 1:
    //strike component
    return fault->u[STRIKE];
    break;
  case 2:// dip
    return fault->u[DIP];
    break;
  case 3:// normal
    return fault->u[NORMAL];
    break;
  default:
    fprintf(stderr,"select_val_for_print: error: mode %i undefined\n",
	    mode);
    exit(-1);
    break;
  }
  return NAN;
}
/* 

   print stresses as stress tensor components 
   1 to 6 to file, either for grid or for the locations as specified in the 
   output location file

*/
void print_stress(struct med *medium,struct flt *fault)
{
  COMP_PRECISION x[3],dx[3];
  int i,j,k,l,m,nz,nxy;
  FILE *out;
  my_boolean use_fault_plane;
  if(medium->print_bulk_fields){
    //
    // output on grid
    //
    fiddle_with_limits_for_plot(medium,&nz,&use_fault_plane,dx,TRUE);
    nxy=medium->n[INT_X] * medium->n[INT_Y];// for array access in POSS
    out=myopen(STRESS_OUT_FILE,"w");
    for(i=0,x[INT_X]=medium->pxmin[INT_X];i<medium->n[INT_X];x[INT_X]+=dx[INT_X],i++)
      for(j=0,x[INT_Y]=medium->pxmin[INT_Y];j<medium->n[INT_Y];x[INT_Y]+=dx[INT_Y],j++)
	for(k=0,x[INT_Z]=medium->pxmin[INT_Z];k<nz;x[INT_Z]+=dx[INT_Z],k++){
	  // check if all stresses are finite
	  if(medium->suppress_nan_output){
	    for(l=0;l<6;l++)
	      if(!FINITE_TEST(medium->s[POSS(i, j, k, l)]))
		break;
	  }else{
	    l=6;
	  }
	  if(l == 6){
	    if(!use_fault_plane || medium->ok[i*medium->n[INT_Y]+j]){
	      for(m=0;m<3;m++){
		//x[m]=reformat_small(x[m]);
		fprintf(out,"%14.7e ",x[m]);
	      }
	      /* if stresses are in a 6 component vector, then
		 
		   s[0] = sm[INT_X][INT_X];
		   s[1] = sm[INT_X][INT_Y];
		   s[2] = sm[INT_X][INT_Z];
		   s[3] = sm[INT_Y][INT_Y];
		   s[4] = sm[INT_Y][INT_Z];
		   s[5] = sm[INT_Z][INT_Z];
	      */
	      for(l=0;l<6;l++)
		fprintf(out,"%14.7e ",medium->s[POSS(i, j, k, l)]);
	      fprintf(out,"\n");
	    }
	  }
	}
    fclose(out);
    fprintf(stderr,"print_stress: stress output in \"%s\"\n",STRESS_OUT_FILE);
    fprintf(stderr,"print_stress: format: x y z s_xx s_xy s_xz s_yy s_yz s_zz\n");

    out=myopen(STRESS_HDR_FILE,"w");
    fprintf(out,"%g %g %i %g %g %i %g %g %i\n",
	    medium->pxmin[INT_X],medium->pxmax[INT_X],medium->n[INT_X],
	    medium->pxmin[INT_Y],medium->pxmax[INT_Y],medium->n[INT_Y],
	    medium->pxmin[INT_Z],medium->pxmax[INT_Z],medium->n[INT_Z]);
    fclose(out);
    fprintf(stderr,"print_stress: dimension information in \"%s\"\n",
	    STRESS_HDR_FILE);
  }
  if(medium->read_oloc_from_file){
    /*

      output on xyz tripels

    */
    out=myopen(STRESS_OUT_FILE,"w");
    for(k=j=i=0;i<medium->olocnr;i++,j+=6,k+=3){
      if(medium->suppress_nan_output){
	for(l=0;l<6;l++)
	  if(!FINITE_TEST(medium->s[j+l]))
	    break;
      }else{
	l=6;
      }
      if(l==6){
	for(m=0;m<3;m++){
	  //medium->xoloc[k+m] = reformat_small(medium->xoloc[k+m]);
	  fprintf(out,"%14.7e ",medium->xoloc[k+m]);
	}
	for(l=0;l<6;l++)
	  fprintf(out,"%14.7e ",medium->s[j+l]);
	fprintf(out,"\n");
      }
    }
    fclose(out);
    fprintf(stderr,"print_stress: stress output in \"%s\"\n",STRESS_OUT_FILE);
    fprintf(stderr,"print_stress: format: x y z s_xx s_xy s_xz s_yy s_yz s_zz\n");
  }
}
/* 
   print stress fields resolved on fault 
   orientation to file 
*/
void print_stress_on_fault(struct med *medium,struct flt *fault,int nrf)
{
  COMP_PRECISION x[3],dx[3],sm[3][3],st,sn,sd;
  int i,j,k,nz,nxy;
  FILE *out;
  my_boolean use_fault_plane;
  fiddle_with_limits_for_plot(medium,&nz,&use_fault_plane,dx,TRUE);
  nxy = medium->n[INT_X] * medium->n[INT_Y];
  out=myopen(FAULT_STRESS_OUT_FILE,"w");
  for(i=0,x[INT_X]=medium->pxmin[INT_X];i<medium->n[INT_X];x[INT_X]+=dx[INT_X],i++)
    for(j=0,x[INT_Y]=medium->pxmin[INT_Y];j<medium->n[INT_Y];x[INT_Y]+=dx[INT_Y],j++)
      for(k=0,x[INT_Z]=medium->pxmin[INT_Z];k<nz;x[INT_Z]+=dx[INT_Z],k++){
	if(!use_fault_plane || medium->ok[i*medium->n[INT_Y]+j]){
	  fprintf(out,"%14.7e %14.7e %14.7e ",x[INT_X],x[INT_Y],x[INT_Z]);
	  sm[INT_X][INT_X]=medium->s[POSS(i, j, k, 0)];
	  sm[INT_X][INT_Y]=sm[INT_Y][INT_X]=medium->s[POSS(i, j, k, 1)];
	  sm[INT_X][INT_Z]=sm[INT_Z][INT_X]=medium->s[POSS(i, j, k, 2)];
	  sm[INT_Y][INT_Y]=medium->s[POSS(i, j, k, 3)];
	  sm[INT_Y][INT_Z]=sm[INT_Z][INT_Y]=medium->s[POSS(i, j, k, 4)];
	  sm[INT_Z][INT_Z]=medium->s[POSS(i, j, k, 5)];
	  calc_three_stress_components(sm,fault[nrf].normal,fault[nrf].t_strike,
				       fault[nrf].normal,fault[nrf].t_dip,
				       &st,&sn,&sd);
	  //st=reformat_small(st);;
	  //sd=reformat_small(sd);
	  //sn=reformat_small(sn);
	  fprintf(out,"%14.7e %14.7e %14.7e",st,sd,sn);
	  fprintf(out,"\n");
	}
      }
  fclose(out);
  fprintf(stderr,"print_stress_on_fault: stress output in %s\n",
	  FAULT_STRESS_OUT_FILE);
  fprintf(stderr,"print_stress_on_fault: format: x y z s_strike s_dip s_normal\n");

	
}
/* 

   print the displacement field, either from grid or from
   spotted observations

*/
void print_displacement(struct med *medium,struct flt *fault)
{
  COMP_PRECISION x[3],dx[3];
  int i,j,k,l,m,nz,nxy;
  FILE *out;
  my_boolean use_fault_plane;
  if(medium->print_bulk_fields){
    // on grid
    fiddle_with_limits_for_plot(medium,&nz,&use_fault_plane,dx,TRUE);
    nxy = medium->n[INT_X] * medium->n[INT_Y];// for POSU
    if(nxy * nz == 0){
      fprintf(stderr,"print_displacement: error: nx: %i ny: %i nz: %i\n",
	     medium->n[INT_X],medium->n[INT_Y],nz);
      exit(-1);
    }
    out=myopen(DISP_OUT_FILE,"w");
    for(i=0,x[INT_X]=medium->pxmin[INT_X];i<medium->n[INT_X];x[INT_X]+=dx[INT_X],i++)
      for(j=0,x[INT_Y]=medium->pxmin[INT_Y];j<medium->n[INT_Y];x[INT_Y]+=dx[INT_Y],j++)
	for(k=0,x[INT_Z]=medium->pxmin[INT_Z];k<nz;x[INT_Z]+=dx[INT_Z],k++){
	  if(medium->suppress_nan_output){
	    for(l=0;l<3;l++)
	      if(!FINITE_TEST(medium->u[POSU(i, j, k, l)]))
		break;
	  }else{
	    l=3;
	  }
	  if(l==3){
	    if(!use_fault_plane || medium->ok[i*medium->n[INT_Y]+j]){
	      for(m=0;m<3;m++){
		//x[m]=reformat_small(x[m]);
		fprintf(out,"%14.7e ",x[m]);
	      }
	      for(l=0;l<3;l++)
		fprintf(out,"%14.7e ",medium->u[POSU(i, j, k, l)]);
	      fprintf(out,"\n");
	    }
	  }
	}
    fclose(out);
    fprintf(stderr,"print_displacement: displacement output in %s\n",
	    DISP_OUT_FILE);
    fprintf(stderr,"print_displacement: format: x y z u_x u_y u_z\n");
    out=myopen(DISP_HDR_FILE,"w");
    fprintf(out,"%g %g %i %g %g %i %g %g %i\n",
	    medium->pxmin[INT_X],medium->pxmax[INT_X],medium->n[INT_X],
	    medium->pxmin[INT_Y],medium->pxmax[INT_Y],medium->n[INT_Y],
	    medium->pxmin[INT_Z],medium->pxmax[INT_Z],medium->n[INT_Z]);
    fclose(out);
    fprintf(stderr,"print_displacement: dimension information in \"%s\"\n",
	    DISP_HDR_FILE);
  }
  if(medium->read_oloc_from_file){
    // on xyz tripels
    out=myopen(DISP_OUT_FILE,"w");
    for(i=j=0;i<medium->olocnr;i++,j+=3){
      if(medium->suppress_nan_output){
	for(l=0;l<3;l++)
	  if(!FINITE_TEST(medium->u[j+l]))
	    break;
      }else{
	l=3;
      }
      if(l == 3){
	for(m=0;m<3;m++){
	  //medium->xoloc[m+j] = reformat_small(medium->xoloc[m+j]);
	  fprintf(out,"%14.7e ",medium->xoloc[m+j]);
	}
	for(l=0;l<3;l++)
	  fprintf(out,"%14.7e ",medium->u[j+l]);
	fprintf(out,"\n");
      }
    }
    fclose(out);
    fprintf(stderr,"print_displacement: displacement output in %s\n",DISP_OUT_FILE);
    fprintf(stderr,"print_displacement: format: x y z u_x u_y u_z\n");
  }
}
/*
  print equations as set up for stress
*/
void print_equations(int naflt,my_boolean *sma,
		     int *nameaf,A_MATRIX_PREC *b,
		     int nreq, char *title,struct flt *fault)
{
  int i,j,eqc1;
  for(eqc1=i=0;i<naflt;i++)
    for(j=0;j<3;j++)
      if(sma[i*3+j]){
	fprintf(stderr,
		"%s: eq %3i/%3i flt %5i/g:%5i m: %1i o_s: %12.4e sd: %12.4e target_s: %12.4e o_slip: %12.4e\n",
		title, eqc1, nreq, nameaf[i], fault[nameaf[i]].group, j,
		fault[nameaf[i]].s[j], b[eqc1], b[eqc1]+fault[nameaf[i]].s[j],
		fault[nameaf[i]].u[j]);
	eqc1++;
      }
}
//
// print slip solutions
//
void print_solutions(int naflt, int *nameaf,struct flt *fault,struct med *medium,
		     char *title)
{
  int i;
  COMP_PRECISION abs_tau;
  for(i=0;i<naflt;i++){
    calc_absolute_shear_stress(&abs_tau,nameaf[i],fault);
    fprintf(stderr,
	    "after t:%11g flt:%3i ss/d/n/c/as/sfn: %10.3e/%10.3e/%10.3e/%10.3e/%10.3e/%10.3e us/d/f: %10.3e/%10.3e/%5.2f, %s\n",  
	    medium->time,nameaf[i],fault[nameaf[i]].s[STRIKE],
	    fault[nameaf[i]].s[DIP],fault[nameaf[i]].s[NORMAL],
	    coulomb_stress(abs_tau,(COMP_PRECISION)fault[nameaf[i]].mu_s,
			   fault[nameaf[i]].s[NORMAL],medium->cohesion),
	    abs_tau,fabs(fault[nameaf[i]].s[NORMAL])*fault[nameaf[i]].mu_d,
	    fault[nameaf[i]].u[STRIKE],fault[nameaf[i]].u[DIP],
	    fault[nameaf[i]].u[STRIKE]/
	    ((fabs(fault[nameaf[i]].u[DIP]) > 0)?(fault[nameaf[i]].u[DIP]):(1.0)),
					title);  
  }
}
/*

  writes the content of the moment stack to a file when time has
  progressed after all activations for this event

  output is in the cevents.dat file

*/
void flush_moment_stack(struct med *medium)
{
  int i,j;
  static unsigned short pline=0;
  pline ++;
  for(i=0;i<medium->nrgrp;i++)
    if(medium->fault_group[i].moment != 0.0){
      fprintf(medium->cevents_out,
	      "%20.15e %5i %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\n",
	      medium->old_moment_time,i,// time of old group activation
	      medium->fault_group[i].moment,// temporaty group moment release
	      medium->total_moment,// total moment at old time
	      // weighted slip according to current slip area
	      // i.e. the slip was added according to patches area
	      // that's why we renormalize here
	      medium->fault_group[i].slip[STRIKE]/medium->fault_group[i].sarea,
	      medium->fault_group[i].slip[DIP]/medium->fault_group[i].sarea,
	      medium->fault_group[i].slip[NORMAL]/medium->fault_group[i].sarea,
	      medium->fault_group[i].slip[3+INT_X]/medium->fault_group[i].sarea,
	      medium->fault_group[i].slip[3+INT_Y]/medium->fault_group[i].sarea,
	      medium->fault_group[i].slip[3+INT_Z]/medium->fault_group[i].sarea,
	      medium->fault_group[i].cent[INT_X]/medium->fault_group[i].sarea,
	      medium->fault_group[i].cent[INT_Y]/medium->fault_group[i].sarea,
	      medium->fault_group[i].cent[INT_Z]/medium->fault_group[i].sarea);
      if(medium->debug)
	fprintf(stderr,
		"event: %20.15e %5i %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
		medium->old_moment_time,i,// time of old group activation
		medium->fault_group[i].moment,// temporary group moment release
		// weighted slip according to current slip area
		// i.e. the slip was added according to patches area
		// that's why we renormalize here
		medium->fault_group[i].slip[STRIKE]/medium->fault_group[i].sarea,
		medium->fault_group[i].slip[DIP]/medium->fault_group[i].sarea,
		medium->fault_group[i].slip[NORMAL]/medium->fault_group[i].sarea,
		medium->fault_group[i].slip[3+INT_X]/medium->fault_group[i].sarea,
		medium->fault_group[i].slip[3+INT_Y]/medium->fault_group[i].sarea,
		medium->fault_group[i].slip[3+INT_Z]/medium->fault_group[i].sarea,
		medium->fault_group[i].cent[INT_X]/medium->fault_group[i].sarea,
		medium->fault_group[i].cent[INT_Y]/medium->fault_group[i].sarea,
		medium->fault_group[i].cent[INT_Z]/medium->fault_group[i].sarea);
      
    }
  // flush every couple of lines
  if(pline > 10){
    fflush(medium->cevents_out);
    pline=0;
  }
  // reset group oriented output
  medium->old_moment_time = medium->time;
  for(i=0;i<medium->nrgrp;i++){
    medium->fault_group[i].moment = 
      medium->fault_group[i].sarea = 0.0;
    for(j=0;j < 6;j++)		/*  */
      medium->fault_group[i].slip[j]=0.0;
    for(j=0;j < 3;j++)		/*  */
      medium->fault_group[i].cent[j]=0.0;
    
  }
#ifdef USE_PGPLOT
  if(medium->op_mode == SIMULATE_LOADING_AND_PLOT)
    plot_moment_array(medium,medium->momrel,2);
#endif
}

/*

  adjust volume displacement and stress plot limits to allow for
  output of displacements and stresses on a plane aligned with the
  average fault plane. x and y then go either in strike/dip or in
  strike/normal directions

  returns new bounds and the use_fault_plane flag

*/
void fiddle_with_limits_for_plot(struct med *medium,int *nz,
				 my_boolean *use_fault_plane,COMP_PRECISION *dx,
				 my_boolean check_for_stress_calc)
{
  if(check_for_stress_calc)
    if(!medium->bulk_field_init){
      fprintf(stderr,"fiddle_with_limits_for_plot: bulk fields were not initialized, failed");
      exit(-1);
    }
  if(medium->n[INT_Z]<0){
    *nz=1;// uses fault plane set up
    *use_fault_plane =TRUE;
  }else{
    *nz=medium->n[INT_Z];
    *use_fault_plane = FALSE;
  }
  if(medium->n[INT_X]<2)
    dx[INT_X]=0.0;
  else
    dx[INT_X]=(medium->pxmax[INT_X]-medium->pxmin[INT_X])/
      ((COMP_PRECISION)(medium->n[INT_X]-1));
  if(medium->n[INT_Y]<2)
    dx[INT_Y]=0.0;
  else
    dx[INT_Y]=(medium->pxmax[INT_Y]-medium->pxmin[INT_Y])/
      ((COMP_PRECISION)(medium->n[INT_Y]-1));
  if(*nz<2)
    dx[INT_Z]=0.0;
  else
    dx[INT_Z]=(medium->pxmax[INT_Z]-medium->pxmin[INT_Z])/
      ((COMP_PRECISION)(*nz-1));
}

void time_report(char *sub,char *out_string,struct med *medium)
{
  struct timespec current_time;
  double tisec,tsec;
  clock_gettime(CLOCK_REALTIME,&current_time);

  tisec = (double)medium->init_time.tv_sec  + (double)medium->init_time.tv_nsec/1e9;
  tsec =  (double)current_time.tv_sec  +      (double)current_time.tv_nsec/1e9;
  
  fprintf(stderr,"%s: %s, %.3f s since init\n",sub,out_string,tsec-tisec);

}
