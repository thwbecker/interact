#include <stdio.h>
#include <string.h>
#include "converter.h"

void openBC( char * path, int *x1, int *x2, int *xn, int *y1, int *y2, 
             int *yn, int *z1, int *z2, int *zn, int *fp, int *dis )
{
  FILE * inFile;

  char fileStream[50];
  char * opMode;
  char * printMode;
  char * minx;
  char * maxx;
  char * xnum;
  char * miny;
  char * maxy;
  char * ynum;
  char * minz;
  char * maxz;
  char * znum;
  char * inc;
  char * start;
  char * stop;
  char * code;
  char * value;

  char * space = " ";
  char * enter = "\n";

  inFile = fopen( path, "r" );

  if( inFile == NULL )
  {
    fprintf( stderr, "Error opening inFile %s!  Program Terminating!!!\n", 
	     path );
    exit(1);
  }

  fgets( fileStream, 50, inFile );
  opMode = strtok( fileStream, space );
  
  fgets( fileStream, 50, inFile );
  printMode = strtok( fileStream, space );

  fgets( fileStream, 50, inFile );
  minx = strtok( fileStream, space );
  maxx = strtok( NULL, space );
  xnum = strtok( NULL, space );
  miny = strtok( NULL, space );
  maxy = strtok( NULL, space );
  ynum = strtok( NULL, space );
  minz = strtok( NULL, space );
  maxz = strtok( NULL, space );
  znum = strtok( NULL, space );

  *x1 = atoi( minx );
  *x2 = atoi( maxx );
  *xn = atoi( xnum );
  *y1 = atoi( miny );
  *y2 = atoi( maxy );
  *yn = atoi( ynum );
  *z1 = atoi( minz );
  *z2 = atoi( maxz );
  *zn = atoi( znum );

  *z1 = 0 - *z1;
  *z2 = 0 - *z2;

  fgets( fileStream, 50, inFile );
  inc = strtok( fileStream, space );
  start = strtok( NULL, space );
  stop = strtok( NULL, space );
  code = strtok( NULL, space );
  value = strtok( NULL, space);

  // printf( "%s %s %s %s %s\n", inc, start, stop, code, value );

  *fp = atoi( code );
  *dis = atoi( value );

  fclose( inFile );
}
