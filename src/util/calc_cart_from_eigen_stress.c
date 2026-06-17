/*
  interact: model fault interactions using dislocations in a 
            halfspace

  (c) Thorsten Becker, thwbecker@post.harvard.edu


*/
#include "interact.h"
/*

  given stresses in the principal component reference frame, calculate
  cartesian stresses

*/
int main(int argc, char **argv)
{
  COMP_PRECISION sr[3][3],x,y,z,eval[3],evec[9],
    s[3][3],r[3][3];
  int i,j,i3,imode=0;
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      ccfes_help(argv);
    }else if(strcmp(argv[i],"-a")==0){// read in strike and dip
      imode=1;
    }else if(strcmp(argv[i],"-ah")==0){// read in horizontal system
      imode=2;
    }else{
      fprintf(stderr,"%s: option \"%s\" not understood\n\n",
	      argv[0],argv[i]);
      ccfes_help(argv);
    }
  }
  while(read_vecs(imode,&x,&y,&z,evec,eval)){
    for(i=i3=0;i<3;i++,i3+=3){
      // normalize eigenvectors
      normalize((evec+i3),3);
      // check orthogonality
      for(j=i+1;j<3;j++)
	if(fabs(dotp_3d((evec+i3),(evec+j*3)))>EPS_COMP_PREC){
	  fprintf(stderr,"%s: error, evectors %i and %i not orthogonal, dotp:  %g\n",
		  argv[0],i+1,j+1,dotp_3d((evec+i3),(evec+j*3)));
	  fprintf(stderr,"%s: evec %i: %g %g %g\n",
		  argv[0],i+1,evec[i3+INT_X],evec[i3+INT_Y],evec[i3+INT_Z]);
	  fprintf(stderr,"%s: evec %i: %g %g %g\n",
		  argv[0],j+1,evec[j*3+INT_X],evec[j*3+INT_Y],evec[j*3+INT_Z]);
	  exit(-1);
	}
      // assemble rotion matrix, prepare stress matrix
      for(j=0;j<3;j++){
	r[j][i] = evec[i3+j];
	s[i][j] = 0.0;
      }
      s[i][i] = eval[i];
    }
    // get cartesian stresses by rotation
    rotate_mat(s,sr,r);
    fprintf(stdout,"%g %g %g \t%g %g %g %g %g %g\n",
	    x,y,z,sr[INT_X][INT_X],sr[INT_X][INT_Y],sr[INT_X][INT_Z],sr[INT_Y][INT_Y],sr[INT_Y][INT_Z],sr[INT_Z][INT_Z]);
    
  }
  exit(0);
}

// read in eigensystem
my_boolean read_vecs(int imode, COMP_PRECISION *x,COMP_PRECISION *y,
		  COMP_PRECISION *z,COMP_PRECISION *evec,COMP_PRECISION *eval)
{
  int rcnt,i;
  COMP_PRECISION strike[3],dip[3],s1h,s2h,azi;
  if(imode == 0){// read in cartesian format
    rcnt = fscanf(stdin,FIFTEEN_CP_FORMAT,x,y,z,
		  (eval),  (evec+INT_X),(evec+INT_Y),(evec+INT_Z),
		  (eval+1),(evec+3+INT_X),(evec+3+INT_Y),(evec+3+INT_Z),
		  (eval+2),(evec+6+INT_X),(evec+6+INT_Y),(evec+6+INT_Z));
    if(rcnt == 15)return TRUE; else return FALSE;
  }else if(imode == 1){// read in strike dip format
    rcnt = fscanf(stdin,TWELVE_CP_FORMAT,x,y,z,
		  (eval),(strike),(dip),
		  (eval+1),(strike+1),(dip+1),
		  (eval+2),(strike+2),(dip+2));
    if(rcnt == 12){
      // convert from strike/dip to vector
      for(i=0;i<3;i++)
	angles_to_vec(dip[i],strike[i],(evec+i*3));
      return TRUE;
    }else
      return FALSE;
  }else if(imode == 2){// read in horizontal stress state format
    rcnt = fscanf(stdin,THREE_CP_FORMAT,x,y,z) +
      fscanf(stdin,THREE_CP_FORMAT,&s1h,&s2h,&azi);
    if(rcnt == 6){
      fprintf(stderr,"read_vecs: converting horizontal (plane stress) state, mean stress: %g max shear: %g\n",
	    (s1h+s2h)/2,(s1h-s2h)/2);
      stress_vec_from_hstate(s1h,s2h,0.0,azi,0,evec,eval);
      return TRUE;
    }else
      return FALSE;
  }else{
    fprintf(stderr,"read_vecs: imode %i undefined\n",imode);
    exit(-1);
  }
}

void ccfes_help(char **argv)
{
  fprintf(stderr,"%s\n",argv[0]);
  fprintf(stderr,"\treads stresses in eigensystem format\n");
  fprintf(stderr,"\tx y z s1 s1x s1y s1z s2 s2x s2y s2z s3 s3x s3y s3z\n");
  fprintf(stderr,"\tand converts stress state into cartesian reference frame, output format\n");
  fprintf(stderr,"\t x y z s_xx s_xy s_xz s_yy s_yz s_zz\n\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"\toptions:\n");
  fprintf(stderr,"\t-a\t reads in x y z s1 strike1 dip1 s2 strike2 dip2 s3 strike3 dip3 format instead of cartesian components\n");
  fprintf(stderr,"\t\there, strike and dip are defined as in interact\n");
  fprintf(stderr,"\t-ah\tstress state is assumed to be horizontal. in this case, the input format is\n");
  fprintf(stderr,"\t\tx y z s1h s2h azi\n");
  fprintf(stderr,"\twhere s1h and s2h are the major and minor stresses in the horizontal plane, and azi is the azimuth of s1h\n");
  fprintf(stderr,"\n");
  exit(-1);
}
