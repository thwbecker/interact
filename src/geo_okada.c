/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: geo_okada.c,v 1.6 2003/03/19 20:06:13 becker Exp $


  driver routine for eval_geo_green


  evaluate Okada's rectangular slip in half-space formula
  for geographical coordinates

*/
#include "interact.h"
#include "blockinvert.h"

int main(int argc, char **argv)
{
  COMP_PRECISION fx[2][2],pole[2]={0,90},x[3],px[3],u[3],pu[3],ps[3][3],
    s[3][3],center[3],l,w,tdepth,depth,dummy=0,dip,
    disp[3]={0,0,0},azi,smazi,cmazi;
  FILE *out1,*out2;
  int iret;

  my_boolean dxfout = TRUE;
  if(argc<6){
    fprintf(stderr,"%s lon1 lat1 lon2 lat2 total_depth_extent\n",
	    argv[0]);
    exit(-1);
  }
  /*
    read in lon,lat of endpoints 
  */
  sscanf(argv[1],"%lf",&fx[0][0]);
  sscanf(argv[2],"%lf",&fx[0][1]);
  sscanf(argv[3],"%lf",&fx[1][0]);
  sscanf(argv[4],"%lf",&fx[1][1]);
  sscanf(argv[5],"%lf",&tdepth);
  //
  // get the projection poles and fault length
  get_projected_fault_parameters(fx,tdepth,center,&azi,&dip,&l,&w,
				 &depth);
  //
  // sin and cos of azimuth for back-rotation from projected frame
  // have to change to alpha frame
  my_sincos(&smazi,&cmazi,-(PIHALF - DEG2RADF(azi)));
  dip = 90.0;disp[STRIKE]=1.0;

  fprintf(stdout,"# f1: %g %g  f2: %g %g a: %g C: %g %g P: %g %g l: %g w: %g z: %g\n",
	  fx[0][0],fx[0][1],fx[1][0],fx[1][1],azi,
	  center[0],center[1],pole[0],pole[1],l,w,-depth);
  if(dxfout){
    // output files for testing
    out1 = myopen("tmp.dx.p","w");
    out2 = myopen("tmp.dx.c","w");
  }
  
  //
  // loop through points to be projected 
  //
  // input format is: lon lat z
  //
  while(fscanf(stdin,"%lf %lf %lf",(x+INT_X),(x+INT_Y),(x+INT_Z))==3){
    geoproject(x, px, FLT_ROT_PROJECTION, center[0], center[1], azi,
	       dummy, dummy, dummy, dummy, (int)FALSE);
    //
    // evaluate the displacement and stresses in the projected 
    // frame
    //
    eval_okada_basic(px,l,w,dip,depth,disp,pu,ps,&iret);
    /*

      since the oblique mercator projection is conformal, we can just
      rotate the stresses and displacements back in the east-north-up
      frame
      
    */
    rotate_vec(pu, u,(double)  cmazi, (double)smazi);
    rotate_mat_z(ps,s, (double) cmazi, (double) smazi);
    if(dxfout){
      //
      // geographic output of displacements in lon lat azi length format
      fprintf(out1,"%g %g %g %g\n",
	      x[INT_X],x[INT_Y],vec_to_strike(u),hypot(u[INT_X],u[INT_Y]));
      // cartesian output of displacements in projected frame
      fprintf(out2,"%g %g %g %g\n",
	      px[INT_X],px[INT_Y],vec_to_strike(pu),hypot(pu[INT_X],pu[INT_Y]));
    }
    // output of lon lat ve vn to stdout
    fprintf(stdout,"%g %g %g %g\n",x[INT_X],x[INT_Y],u[INT_X],u[INT_Y]);
  }
  if(dxfout){
    fclose(out1);fclose(out2);
  }
  GMT_end (1, argv);
  return 0;
}
