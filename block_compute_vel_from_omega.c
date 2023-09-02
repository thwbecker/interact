/*


reads in faults in the blockinvert format
reads in an omega file in 
w_x sig_wx ... 
format,
and write velocities to stdout


$Id: block_compute_vel_from_omega.c,v 1.2 2004/03/25 23:48:47 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"

int main(int argc, char **argv)
{
  COMP_PRECISION dummy=0,*cx,velrms,x[3],px[3],dir[2],lfac;
  FILE *out;
  int i,j,minus_one=-1;
  struct prj projection;
  struct bmd *mod;
  init_block_mods(&mod);
  projection.type = OMERC_AZI;projection.azi=90;
  // read in velocities to get an appropriate projection
  read_gps_velocities(&mod,&projection,&velrms,argv,FALSE,
		      &minus_one);
  //
  // read in centroids
  block_read_centroids(&cx,&mod->nrb);
  //
  // simply read in faults
  read_bflt(&mod,15.,projection,FALSE,FLT_MAX,
	    FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,
	    FALSE,FALSE);
  /*
    check the fault assignment by comparing the 
    centroid locations
  */
  x[INT_Z] = 0.0;
  fprintf(stderr,"%s: read centroids ok, testing %i faults\n",
	  argv[0],mod->nflt);
  for(i=0;i < mod->nflt;i++){
    for(j=0;j < 2;j++){
      /* 

	project left centroid in fault-local system

      */
      a_equals_b_vector(x,(cx+(mod->fault+i)->block[j]*2),2); 
      geoproject(x, px, FLT_ROT_PROJECTION,(mod->fault+i)->x[INT_X],
		 (mod->fault+i)->x[INT_Y],(mod->fault+i)->azi,dummy, 
		 dummy, dummy, dummy, (int)FALSE);
      dir[j] = -px[INT_Y] / norm(px,2);
    }
    if(0)
      fprintf(stderr,"flt %3i: (%12g, %12g) a: %12g lc: %12g, %12g dir: %12g rc: %12g, %12g dir: %12g\n",
	      i+1,mod->fault[i].x[INT_X],mod->fault[i].x[INT_Y],mod->fault[i].azi,
	      cx[mod->fault[i].block[0]*2+INT_X],
	      cx[mod->fault[i].block[0]*2+INT_Y],dir[0],
	      cx[mod->fault[i].block[1]*2+INT_X],
	      cx[mod->fault[i].block[1]*2+INT_Y],dir[1]);
    /*

      check if block are on the correct side

      the left (block[0]) block should have a positive direction
      while the right (block[0]) should have a negative direction
      (dir[]) in the definition from above ( dir[j] = -px[Y] / norm(px,2);)
      
    */
    if(fabs(dir[0]) > fabs(dir[1])){
      if(dir[0] < 0)
	flip_block_code((mod->fault+i),dir);
    }else{
      if(dir[1] > 0)
	flip_block_code((mod->fault+i),dir);
    }
    if((dir[0] < -0.2) || (dir[1] > 0.2)){
      fprintf(stderr,"%s: WARNING: flt %i: blocks: %i %i dirs: %g %g \n",
	      argv[0],i,mod->fault[i].block[0]+1,mod->fault[i].block[1]+1,dir[0],dir[1]);
    }
  }
  /*
    
    output of resorted faults

  */
  out = myopen("fltcodes.gmt","w");
  for(i=0;i < mod->nflt;i++){
    fprintf(stdout,"%15.10f %15.10f %15.10f %15.10f %12g %5i %5i %12g %12g\n",
	    mod->fault[i].ex[0][INT_X],mod->fault[i].ex[0][INT_Y],
	    mod->fault[i].ex[1][INT_X],mod->fault[i].ex[1][INT_Y],
	    mod->fault[i].dip,mod->fault[i].block[0]+1,mod->fault[i].block[1]+1,
	    mod->fault[i].lfac,mod->fault[i].ld);
    //
    lfac = (mod->fault[i].l/1.2) / 111.1949;
    for(j=0;j < 2;j++){
      // output of left/right block codes for plotting
      // mod->fault number, l/r, center location offset-location, block code
      if(j==1)lfac = - lfac;
      fprintf(out,"%i %i %g %g %g %g %i\n",
	      i+1,j+1,mod->fault[i].x[INT_X],mod->fault[i].x[INT_Y],
	      mod->fault[i].x[INT_X] + lfac * mod->fault[i].evec[NORMAL*3+INT_X],
	      mod->fault[i].x[INT_Y] + lfac * mod->fault[i].evec[NORMAL*3+INT_Y],
	      mod->fault[i].block[j]+1);
    }
  }
  fclose(out);
  fprintf(stderr,"%s: written fault block codes to fltcodes.gmt\n",
	  argv[0]);
  return 0;
}
