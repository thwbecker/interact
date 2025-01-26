/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: print_patch_geometry.c,v 1.23 2003/01/07 02:52:29 tbecker Exp $

  output of patch gemeotry and other stuff to *out stream
  in various formats

*/
#include "interact.h"
#include "properties.h"

/* pass a single fault patch here, hence */
/* flt_offset should normally be zero */
void print_patch_geometry_and_bc(int flt_offset,struct flt *fault,
				 int opmode,
				 COMP_PRECISION stop_time,
				 my_boolean calculate_base_vecs,
				 FILE *out,
				 my_boolean shrink_patches,
				 COMP_PRECISION *scalar)
{
  int k,l;
  static my_boolean bc_init=FALSE,g_init=FALSE;
  static int nrf,bc_code;
  COMP_PRECISION corner[4][3],alpha,sin_dip,cos_dip,x[3],lfac,leeway;
#ifdef ALLOW_NON_3DQUAD_GEOM
  COMP_PRECISION tmpdbl;
#endif
  //
  // shrink patches for easier viewing
  //
  if(shrink_patches)
    leeway = 0.9;
  else
    leeway = 1.0;

  if(calculate_base_vecs){
#ifdef ALLOW_NON_3DQUAD_GEOM
    if(fault[flt_offset].type != TRIANGULAR){
      // normal (rectangular, 2-D, or point source)
      alpha=90.0-(COMP_PRECISION)fault[flt_offset].strike;
      my_sincos_deg(&fault[flt_offset].sin_alpha,&fault[flt_offset].cos_alpha,(COMP_PRECISION)alpha);
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[flt_offset].dip);
      calc_base_vecs(fault[flt_offset].t_strike,fault[flt_offset].normal,fault[flt_offset].t_dip,
		     fault[flt_offset].sin_alpha,fault[flt_offset].cos_alpha,sin_dip,cos_dip);
    }else{// triangular element
      get_alpha_dip_tri_gh((fault+flt_offset)->xt,&(fault+flt_offset)->sin_alpha,
			   &(fault+flt_offset)->cos_alpha,&tmpdbl,&(fault+flt_offset)->w);
      (fault+flt_offset)->dip=(float)tmpdbl;
      (fault+flt_offset)->area = (fault+flt_offset)->w;
      alpha=RAD2DEGF(asin((fault+flt_offset)->sin_alpha));
      (fault+flt_offset)->strike= 90.0 - alpha;
    }
#else
    //
    // rectangular patch
    //
    alpha=90.0-(COMP_PRECISION)fault[flt_offset].strike;
    my_sincos_deg(&fault[flt_offset].sin_alpha,&fault[flt_offset].cos_alpha,(COMP_PRECISION)alpha);
    my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[flt_offset].dip);
    calc_base_vecs(fault[flt_offset].t_strike,fault[flt_offset].normal,fault[flt_offset].t_dip,
		   fault[flt_offset].sin_alpha,fault[flt_offset].cos_alpha,sin_dip,cos_dip);
#endif
  }
  switch(opmode){
  case PATCH_OUT_MODE:{
    /*
      
      output in the PATCH format as read by interact

    */
#ifdef ALLOW_NON_3DQUAD_GEOM
    switch(fault[flt_offset].type){
    case POINT_SOURCE:{
      // L' = -aspect_ratio = - L/W , W' = area = 4 L W
      fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i\n",
	      fault[flt_offset].x[X], fault[flt_offset].x[Y], 
	      fault[flt_offset].x[Z],fault[flt_offset].strike,
	      fault[flt_offset].dip,-fault[flt_offset].l/fault[flt_offset].w,
	      4.0*fault[flt_offset].w*fault[flt_offset].l,fault[flt_offset].group);
      break;
    }
    case TWO_DIM_SEGMENT_PLANE_STRAIN:
    case TWO_DIM_SEGMENT_PLANE_STRESS:
    case RECTANGULAR_PATCH:{// normal rectangular fault, this is the default
      fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i\n",
	      fault[flt_offset].x[X], fault[flt_offset].x[Y], 
	      fault[flt_offset].x[Z],fault[flt_offset].strike,
	      fault[flt_offset].dip,
	      fault[flt_offset].l,fault[flt_offset].w,fault[flt_offset].group);
      break;
    }
    case TRIANGULAR:{// xt has to be assigned and allocated before
      fprintf(out,"%15.8e %15.8e %15.8e %10.6f %10.6f %15.8e %15.8e %6i ",
	      999.,999.,999.,999.,999.,-1.,-1.,fault[flt_offset].group);
      fprintf(out,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	      fault[flt_offset].xt[X],fault[flt_offset].xt[Y],fault[flt_offset].xt[Z],
	      fault[flt_offset].xt[3+X],fault[flt_offset].xt[3+Y],fault[flt_offset].xt[3+Z],
      	      fault[flt_offset].xt[6+X],fault[flt_offset].xt[6+Y],fault[flt_offset].xt[6+Z]);
      break;
    }
    default:
      fprintf(stderr,"print_fault_patch: patch type %i undefined\n",
	      fault[flt_offset].type);
      exit(-1);
      break;
    }
#else
    // only rectangular
    fprintf(out,"%15.8e %15.8e %15.8e %10.6f %10.6f %15.8e %15.8e %6i\n",
	    fault[flt_offset].x[X], fault[flt_offset].x[Y], fault[flt_offset].x[Z],fault[flt_offset].strike,
	    fault[flt_offset].dip,fault[flt_offset].l,fault[flt_offset].w,fault[flt_offset].group);
#endif
    break;
  }
  case GEOMVIEW_MODE:{
    /*

      output in GEOMVIEW format for 3-D viewing

     */
    if(!g_init){
      fprintf(out,"QUAD\n");
      g_init=TRUE;
    }
#ifdef ALLOW_NON_3DQUAD_GEOM
    switch(fault[flt_offset].type){
    case POINT_SOURCE:
    case TWO_DIM_SEGMENT_PLANE_STRAIN:
    case TWO_DIM_SEGMENT_PLANE_STRESS:
    case RECTANGULAR_PATCH:{
      calculate_bloated_corners(corner,(fault+flt_offset),leeway);
      for(k=0;k<4;k++){
	for(l=0;l<3;l++)
	  if(fabs(corner[k][l]/CHAR_FAULT_DIM)>EPS_COMP_PREC)
	    fprintf(out,"%g ",corner[k][l]/CHAR_FAULT_DIM);
	  else
	    fprintf(out,"0.0 ");
	fprintf(out,"\n");
      }
      break;
    }
    case TRIANGULAR:{
      for(k=0;k<3;k++)
	for(l=0;l<3;l++)
	  fprintf(out,"%g ",fault[flt_offset].xt[k*3+l]);
      for(k=l=0;l<3;l++)
	fprintf(out,"%g ",fault[flt_offset].xt[k*3+l]);
      fprintf(out,"\n");
      break;
    }}
#else
    // only rectangular faults
    calculate_bloated_corners(corner,(fault+flt_offset),leeway);
    for(k=0;k<4;k++){
      for(l=0;l<3;l++)
	if(fabs(corner[k][l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%g ",corner[k][l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"0.0 ");
      fprintf(out,"\n");
    }
#endif
    break;
  }
  /*

    output in XYZ format for plotting with pxxyz from
    GMT

  */
  case PSXYZ_SCALAR_MODE:
  case PSXYZ_MODE:{
    if(opmode == PSXYZ_SCALAR_MODE)
      fprintf(out,"> -Z%g\n",scalar[flt_offset]);
#ifdef ALLOW_NON_3DQUAD_GEOM
    switch(fault[flt_offset].type){
    case TWO_DIM_SEGMENT_PLANE_STRAIN:
    case TWO_DIM_SEGMENT_PLANE_STRESS:{
      calculate_bloated_corners(corner,(fault+flt_offset),leeway);
      // draw segment with endbars
      lfac = fault[flt_offset].l * 0.2;
      for(l=0;l<3;l++)
	x[l] = corner[0][l] + fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%g %g %g\n",x[X]/CHAR_FAULT_DIM,x[Y]/CHAR_FAULT_DIM,x[Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = corner[0][l] - fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%g %g %g\n",x[X]/CHAR_FAULT_DIM,x[Y]/CHAR_FAULT_DIM,x[Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = corner[0][l];
      fprintf(out,"%g %g %g\n",x[X]/CHAR_FAULT_DIM,x[Y]/CHAR_FAULT_DIM,x[Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = corner[1][l];
      fprintf(out,"%g %g %g\n",x[X]/CHAR_FAULT_DIM,x[Y]/CHAR_FAULT_DIM,x[Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = corner[1][l] + fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%g %g %g\n",x[X]/CHAR_FAULT_DIM,x[Y]/CHAR_FAULT_DIM,x[Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = corner[1][l] - fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%g %g %g\n",x[X]/CHAR_FAULT_DIM,x[Y]/CHAR_FAULT_DIM,x[Z]/CHAR_FAULT_DIM);
      break;
    }
    case POINT_SOURCE:
    case RECTANGULAR_PATCH:{
      calculate_bloated_corners(corner,(fault+flt_offset),leeway);
      for(k=0;k<4;k++){
	for(l=0;l<3;l++){
	  if(fabs(corner[k][l]/CHAR_FAULT_DIM)>
	     EPS_COMP_PREC)
	    fprintf(out,"%g ",corner[k][l]/CHAR_FAULT_DIM);
	  else
	    fprintf(out,"0.0 ");
	}
	fprintf(out,"\n");
      }
      break;
    }
    case TRIANGULAR:{
      for(k=0;k<3;k++){
	for(l=0;l<3;l++)
	  fprintf(out,"%g ",fault[flt_offset].xt[k*3+l]);
	fprintf(out,"\n");
      }
      for(k=l=0;l<3;l++)
	fprintf(out,"%g ",fault[flt_offset].xt[k*3+l]);
      fprintf(out,"\n");
      break;
    }
    }
#else
    calculate_bloated_corners(corner,(fault+flt_offset),leeway);
    for(k=0;k<4;k++){
      for(l=0;l<3;l++){
	if(fabs(corner[k][l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%g ",corner[k][l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"0.0 ");
      }
      fprintf(out,"\n");
    }
#endif
    if(opmode != PSXYZ_SCALAR_MODE)
      fprintf(out,">\n");
    break;
  }
  case CORNEROUT_MODE:{
#ifdef ALLOW_POINT_SOUCE
    if((fault[flt_offset].type != POINT_SOURCE)&&
       (fault[flt_offset].type != RECTANGULAR_PATCH)&&
       (!patch_is_2d(fault[flt_offset].type))){
      fprintf(stderr,"print_patch_geometry: corner out defined for rectangular, point source and 2-Dpatches only\n");
      exit(-1);
    }
#endif
    calculate_bloated_corners(corner,(fault+flt_offset),leeway);
    for(k=0;k<4;k++){
      for(l=0;l<3;l++){
	if(fabs(corner[k][l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%g ",corner[k][l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"0 ");
      }
    }
    fprintf(out,"\n");
    break;
  }
  case BC_OUT_MODE:{
    if(!bc_init){
      fprintf(out,"%i\n",SIMULATE_LOADING_AND_PLOT);
      fprintf(out,"%g %g %g\n",0.01,0.01,stop_time);
      fprintf(out,"%g %g\n",0.01,0.1); 
      bc_init=TRUE;
      nrf=0;
    }
    if(fault[flt_offset].strike > 45 && 
       fault[flt_offset].strike <= 135)
      bc_code=STRIKE_SLIP_RIGHTLATERAL;
    else
      bc_code=STRIKE_SLIP_LEFTLATERAL;
    fprintf(out,"-1 %i %i %i\n",nrf,nrf,bc_code);
    nrf++;
    break;
  }
  case VEC_OUT_MODE:{
    fprintf(out,"%12.5e %12.5e %12.5e\t %12.5e %12.5e %12.5e\t %12.5e %12.5e %12.5e\n",
	   fault[flt_offset].t_strike[X],fault[flt_offset].t_strike[Y],
	   fault[flt_offset].t_strike[Z],
	   fault[flt_offset].t_dip[X],fault[flt_offset].t_dip[Y],
	   fault[flt_offset].t_dip[Z],
	   fault[flt_offset].normal[X],fault[flt_offset].normal[Y],
	   fault[flt_offset].normal[Z]);
    break;
  }
  case XYZ_AND_VEC_MODE:{
#ifdef ALLOW_POINT_SOUCE
    if((fault[flt_offset].type != POINT_SOURCE)&&
       (fault[flt_offset].type != RECTANGULAR_PATCH)&&
       (!patch_is_2d(fault[flt_offset].type))){
      fprintf(stderr,"print_patch_geometry: corner out defined for rectangular, point source and 2-Dpatches only\n");
      exit(-1);
    }
#endif
    calculate_bloated_corners(corner,(fault+flt_offset),leeway);
    for(k=0;k<4;k++){
      for(l=0;l<3;l++){
	if(fabs(corner[k][l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%12.5e ",corner[k][l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"%12.5e ",0.0);
      }
    }
    fprintf(out,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	    fault[flt_offset].t_strike[X],fault[flt_offset].t_strike[Y],
	    fault[flt_offset].t_strike[Z],
	    fault[flt_offset].t_dip[X],fault[flt_offset].t_dip[Y],
	    fault[flt_offset].t_dip[Z],
	    fault[flt_offset].normal[X],fault[flt_offset].normal[Y],
	    fault[flt_offset].normal[Z]);
    break;
  }
  default:{
    fprintf(stderr,"print_patch_geometry: mode %i undefined\n",
	    opmode);
    exit(-1);
  }}
}



