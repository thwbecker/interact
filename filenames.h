/* 
   FILENAMES AND OUTPUT OPTIONS 
*/
#define GEOMETRY_FILE "geom.in" // patch geometry input
#define BC_FILE "bc.in" // boundary conditions
#define FAULT_PROP_FILE "fp.in" // fault properties
#define FAULT_STRESS_INIT_FILE "fsi.in" // stress initialization for faults
#define STRESS_OUT_FILE "stress.out" // stress field output
#define STRESS_HDR_FILE "stress.hdr" // header file
#define FAULT_STRESS_OUT_FILE "flt_stress.out" // stress on fault output
#define DISP_OUT_FILE "disp.out" // displacement field output
#define DISP_HDR_FILE "disp.hdr" // header 
#define FAULT_DATA_PREFIX "flt"
#define ONE_STEP_FAULT_DATA_FILE "flt.dat" // stress and displacement on fault data
#define EVENT_FILE_BINARY  "events.bin"    // binary event log for loading simulation
#define EVENT_FILE_ASCII  "events.dat"     // ASCII version
#define RESTART_EVENT_FILE_BINARY  "restart.events.bin" // restart file, where code looks for old EVENT_FILE_BINARY
#define RESTART_EVENT_FILE_ASCII  "restart.events.dat"  // same for ASCII
#define CEVENT_FILE "cevents.dat" // ASCII log of fault events (as opposed to patch events which are logged in EVENT_FILE_BINARY 
#define STRESS_STAT_FILE "sstat.dat" // stress heterogeneity statistics
#define OLOC_FILE "oloc.dat" // 
#define RES_STRESS_FILE "rstress.dat" // resolved stress on fault segments for one-step loading 
#define PLANE_COORD_FILE "plane.xyz" // file with real coordinates of fault plane local output

// program will atuomatically append the PID and ".dat" for the matrix
// and PID and ".hdr" for the header file
#define INTERACTION_MATRIX_FILE "/tmp/i"

#define A_MATRIX_FILE "a.dat"
#define B_VECTOR_FILE "b.dat"
#define A_CON_MATRIX_FILE "ac.dat"
#define DEBUG_A_CON_MATRIX_ASCII_OUT "ac.asc"
#define DEBUG_B_CON_VECTOR_ASCII_OUT "bc.asc"
#define DEBUG_X_CON_VECTOR_ASCII_OUT "xc.asc"
#define DEBUG_A_MIX_MATRIX_ASCII_OUT "am.asc"
#define DEBUG_B_MIX_VECTOR_ASCII_OUT "bm.asc"
#define DEBUG_X_MIX_VECTOR_ASCII_OUT "xm.asc"
#define SLIP_LINE_PREFIX "slipline"
#define STRESS_RELATION_FILE "smlin.in"

#define SV_FILE "sv.dat"
