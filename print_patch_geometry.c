/*
  interact: model fault interactions using dislocations in a halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu

  output of patch gemeotry and other stuff to *out stream
  in various formats

*/
#include "interact.h"
#include "properties.h"

/* pass a single fault patch here, hence */
/* flt_offset should normally be zero */
int print_patch_geometry_and_bc(int flt_offset,struct flt *fault,
				int opmode,
				COMP_PRECISION stop_time,
				MODE_TYPE calculate_base_vecs,
				FILE *out,
				my_boolean shrink_patches,
				COMP_PRECISION *scalar)
{
  int k,l,i,ncon,ielmul,nvert;
  static my_boolean bc_init=FALSE;
  static int nrf,bc_code;
  COMP_PRECISION vertex[MAX_NR_EL_VERTICES*3],sin_dip,cos_dip,leeway;
  double alpha;
#ifdef ALLOW_NON_3DQUAD_GEOM
  COMP_PRECISION lfac,x[3];
#endif
  //
  // shrink patches for easier viewing
  //
  if(shrink_patches)
    leeway = 0.9;
  else
    leeway = 1.0;

  if(calculate_base_vecs){
    /* 
       
       recompute the base vectors?

     */
#ifdef ALLOW_NON_3DQUAD_GEOM
    if((!is_triangular(fault[flt_offset].type))&&(fault[flt_offset].type != IQUAD)){
      // normal (rectangular, 2-D, or point source)
      alpha=90.0-(double)fault[flt_offset].strike;
      my_sincos_degd(&fault[flt_offset].sin_alpha,&fault[flt_offset].cos_alpha,alpha);
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[flt_offset].dip);
      calc_quad_base_vecs(fault[flt_offset].t_strike,fault[flt_offset].normal,fault[flt_offset].t_dip,
			  fault[flt_offset].sin_alpha,fault[flt_offset].cos_alpha,sin_dip,cos_dip);
    }else{// triangular element
      get_tri_prop_based_on_gh((fault+flt_offset));
    }
#else
    //
    // rectangular patch
    //
    alpha=90.0-(COMP_PRECISION)fault[flt_offset].strike;
    my_sincos_degd(&fault[flt_offset].sin_alpha,&fault[flt_offset].cos_alpha,(double)alpha);
    my_sincos_deg(&sin_dip,&cos_dip,fault[flt_offset].dip);
    calc_quad_base_vecs(fault[flt_offset].t_strike,fault[flt_offset].normal,fault[flt_offset].t_dip,
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
	      fault[flt_offset].x[INT_X], fault[flt_offset].x[INT_Y], 
	      fault[flt_offset].x[INT_Z],fault[flt_offset].strike,
	      fault[flt_offset].dip,-fault[flt_offset].l/fault[flt_offset].w,
	      4.0*fault[flt_offset].w*fault[flt_offset].l,fault[flt_offset].group);
      break;
    }
    case TWO_DIM_SEGMENT_PLANE_STRAIN:
    case TWO_DIM_SEGMENT_PLANE_STRESS:{// normal rectangular fault, this is the default
      fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i\n",
	      fault[flt_offset].x[INT_X], fault[flt_offset].x[INT_Y], 
	      fault[flt_offset].x[INT_Z],fault[flt_offset].strike,
	      fault[flt_offset].dip,
	      fault[flt_offset].l,fault[flt_offset].w,fault[flt_offset].group);
      break;
    }
    case TRIANGULAR_M244:
    case TRIANGULAR_M236:
    case TRIANGULAR_HYBR:
    case TRIANGULAR:{// xt has to be assigned and allocated before
      fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i ",
	      fault[flt_offset].x[INT_X],	      fault[flt_offset].x[INT_Y],	      fault[flt_offset].x[INT_Z],
	      fault[flt_offset].strike,fault[flt_offset].dip,
	      -fault[flt_offset].l,-fault[flt_offset].w,
	      fault[flt_offset].group);
      fprintf(out,"%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n",
	      fault[flt_offset].xn[INT_X  ],fault[flt_offset].xn[INT_Y  ],fault[flt_offset].xn[INT_Z],
	      fault[flt_offset].xn[3+INT_X],fault[flt_offset].xn[3+INT_Y],fault[flt_offset].xn[3+INT_Z],
      	      fault[flt_offset].xn[6+INT_X],fault[flt_offset].xn[6+INT_Y],fault[flt_offset].xn[6+INT_Z]);
      break;
    }
    case IQUAD:{// xt has to be assigned and allocated before
      fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i ",
	      fault[flt_offset].x[INT_X],	      fault[flt_offset].x[INT_Y],	      fault[flt_offset].x[INT_Z],
	      fault[flt_offset].strike,fault[flt_offset].dip,
	      fault[flt_offset].l,-fault[flt_offset].w,
	      fault[flt_offset].group);
      fprintf(out,"%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n",
	      fault[flt_offset].xn[9+INT_X],fault[flt_offset].xn[9+INT_Y],fault[flt_offset].xn[9+INT_Z],
	      fault[flt_offset].xn[12+INT_X],fault[flt_offset].xn[12+INT_Y],fault[flt_offset].xn[12+INT_Z],
      	      fault[flt_offset].xn[3+INT_X],fault[flt_offset].xn[3+INT_Y],fault[flt_offset].xn[3+INT_Z],
	      fault[flt_offset].xn[6+INT_X],fault[flt_offset].xn[6+INT_Y],fault[flt_offset].xn[6+INT_Z]);
      break;
    }
    case OKADA_PATCH:{
      fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i\n",
	      fault[flt_offset].x[INT_X], fault[flt_offset].x[INT_Y], fault[flt_offset].x[INT_Z],fault[flt_offset].strike,
	      fault[flt_offset].dip,fault[flt_offset].l,fault[flt_offset].w,fault[flt_offset].group);
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
    fprintf(out,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i\n",
	    fault[flt_offset].x[INT_X], fault[flt_offset].x[INT_Y], fault[flt_offset].x[INT_Z],fault[flt_offset].strike,
	    fault[flt_offset].dip,fault[flt_offset].l,fault[flt_offset].w,fault[flt_offset].group);
#endif
    break;
  }
  /*

    output in XYZ format for plotting with pxxyz from
    GMT

  */
  case PSXYZ_SCALAR_MODE:
  case PSXYZ_STRIKE_DISP_OUT_MODE:
  case PSXYZ_MODE:{
    calculate_bloated_vertices(vertex,(fault+flt_offset),leeway);
#ifdef ALLOW_NON_3DQUAD_GEOM
    switch(fault[flt_offset].type){
    case TWO_DIM_SEGMENT_PLANE_STRAIN:
    case TWO_DIM_SEGMENT_PLANE_STRESS:{
      if(opmode == PSXYZ_SCALAR_MODE)
	fprintf(out,"> -Z%e\n",scalar[flt_offset]);
      else if(opmode == PSXYZ_STRIKE_DISP_OUT_MODE)
	fprintf(out,"> -Z%e\n",fault[flt_offset].u[STRIKE]);
      // draw segment with endbars
      lfac = fault[flt_offset].l * 0.2;
      for(l=0;l<3;l++)
	x[l] = vertex[0*3+l] + fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%20.15e %20.15e %20.15e\n",x[INT_X]/CHAR_FAULT_DIM,x[INT_Y]/CHAR_FAULT_DIM,x[INT_Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = vertex[0*3+l] - fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%20.15e %20.15e %20.15e\n",x[INT_X]/CHAR_FAULT_DIM,x[INT_Y]/CHAR_FAULT_DIM,x[INT_Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = vertex[0*3+l];
      fprintf(out,"%20.15e %20.15e %20.15e\n",x[INT_X]/CHAR_FAULT_DIM,x[INT_Y]/CHAR_FAULT_DIM,x[INT_Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = vertex[1*3+l];
      fprintf(out,"%20.15e %20.15e %20.15e\n",x[INT_X]/CHAR_FAULT_DIM,x[INT_Y]/CHAR_FAULT_DIM,x[INT_Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = vertex[1*3+l] + fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%20.15e %20.15e %20.15e\n",x[INT_X]/CHAR_FAULT_DIM,x[INT_Y]/CHAR_FAULT_DIM,x[INT_Z]/CHAR_FAULT_DIM);
      for(l=0;l<3;l++)
	x[l] = vertex[1*3+l] - fault[flt_offset].normal[l] * lfac;
      fprintf(out,"%20.15e %20.15e %20.15e\n",x[INT_X]/CHAR_FAULT_DIM,x[INT_Y]/CHAR_FAULT_DIM,x[INT_Z]/CHAR_FAULT_DIM);
      break;
    }
    case IQUAD:
      ielmul = number_of_subpatches((fault+flt_offset));
      for(i=0;i < ielmul;i++){
	if(opmode == PSXYZ_SCALAR_MODE)
	  fprintf(out,"> -Z%e\n",scalar[flt_offset]);
	else if(opmode == PSXYZ_STRIKE_DISP_OUT_MODE)
	  fprintf(out,"> -Z%e\n",fault[flt_offset].u[STRIKE]);
	ncon = ncon_of_subpatch((fault+flt_offset),i);
	for(k=0;k < ncon;k++){
	  for(l=0;l<3;l++){
	    if(fabs(vertex[node_number_of_subelement((fault+flt_offset),k, i)*3+l]/CHAR_FAULT_DIM) > EPS_COMP_PREC)
	      fprintf(out,"%20.15e ",vertex[node_number_of_subelement((fault+flt_offset),k, i)*3+l]/CHAR_FAULT_DIM);
	    else
	      fprintf(out,"0.0 ");
	  }
	  fprintf(out,"\n");
	}
      }
      break;
    case POINT_SOURCE:
    case TRIANGULAR_M244:
    case TRIANGULAR_M236:
    case TRIANGULAR_HYBR:
    case TRIANGULAR:
    case OKADA_PATCH:{
      if(opmode == PSXYZ_SCALAR_MODE)
	fprintf(out,"> -Z%e\n",scalar[flt_offset]);
      else if(opmode == PSXYZ_STRIKE_DISP_OUT_MODE)
	fprintf(out,"> -Z%e\n",fault[flt_offset].u[STRIKE]);
      ncon = ncon_of_subpatch((fault+flt_offset),0);
      for(k=0;k < ncon;k++){
	//fprintf(stderr,"%i/%i %20.15e %20.15e %20.15e\n",k,ncon,vertex[k*3+0],vertex[k*3+1],vertex[k*3+2]);
	for(l=0;l<3;l++){
	  if(fabs(vertex[k*3+l]/CHAR_FAULT_DIM) > EPS_COMP_PREC)
	    fprintf(out,"%20.15e ",vertex[k*3+l]/CHAR_FAULT_DIM);
	  else
	    fprintf(out,"0.0 ");
	}
	fprintf(out,"\n");
      }
      break;
    }}
#else
    if(opmode == PSXYZ_SCALAR_MODE)
      fprintf(out,"> -Z%e\n",scalar[flt_offset]);
    else if(opmode == PSXYZ_STRIKE_DISP_OUT_MODE)
      fprintf(out,"> -Z%e\n",fault[flt_offset].u[STRIKE]);
    for(k=0;k<4;k++){
      for(l=0;l<3;l++){
	if(fabs(vertex[k*3+l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%20.15e ",vertex[k*3+l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"0.0 ");
      }
      fprintf(out,"\n");
    }
#endif
    if((opmode != PSXYZ_SCALAR_MODE)&&(opmode !=  PSXYZ_STRIKE_DISP_OUT_MODE))
      fprintf(out,">\n");
    break;
  }
  case VERTEXOUT_MODE:{
    calculate_bloated_vertices(vertex,(fault+flt_offset),leeway);
    nvert = nvert_of_patch((fault+flt_offset));
    for(k=0;k < nvert;k++){
      for(l=0;l<3;l++){
	if(fabs(vertex[k*3+l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%20.15e ",vertex[k*3+l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"0 ");
      }
      fprintf(out,"\n");
    }
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
	   fault[flt_offset].t_strike[INT_X],fault[flt_offset].t_strike[INT_Y],
	   fault[flt_offset].t_strike[INT_Z],
	   fault[flt_offset].t_dip[INT_X],fault[flt_offset].t_dip[INT_Y],
	   fault[flt_offset].t_dip[INT_Z],
	   fault[flt_offset].normal[INT_X],fault[flt_offset].normal[INT_Y],
	   fault[flt_offset].normal[INT_Z]);
    break;
  }
  case XYZ_AND_VEC_MODE:{
#ifdef ALLOW_NON_3DQUAD_GEOM
    if((fault[flt_offset].type != POINT_SOURCE)&&
       (fault[flt_offset].type != OKADA_PATCH)&&
       (!patch_is_2d(fault[flt_offset].type))){
      fprintf(stderr,"print_patch_geometry: vertex out defined for rectangular, point source and 2-Dpatches only\n");
      exit(-1);
    }
#endif
    calculate_bloated_vertices(vertex,(fault+flt_offset),leeway);
    for(k=0;k<4;k++){
      for(l=0;l<3;l++){
	if(fabs(vertex[k*3+l]/CHAR_FAULT_DIM)>
	   EPS_COMP_PREC)
	  fprintf(out,"%12.5e ",vertex[k*3+l]/CHAR_FAULT_DIM);
	else
	  fprintf(out,"%12.5e ",0.0);
      }
    }
    fprintf(out,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	    fault[flt_offset].t_strike[INT_X],fault[flt_offset].t_strike[INT_Y],
	    fault[flt_offset].t_strike[INT_Z],
	    fault[flt_offset].t_dip[INT_X],fault[flt_offset].t_dip[INT_Y],
	    fault[flt_offset].t_dip[INT_Z],
	    fault[flt_offset].normal[INT_X],fault[flt_offset].normal[INT_Y],
	    fault[flt_offset].normal[INT_Z]);
    break;
  }
  default:{
    fprintf(stderr,"print_patch_geometry: mode %i undefined\n",
	    opmode);
    exit(-1);
  }}
  return nvert_of_patch((fault+flt_offset));
}



