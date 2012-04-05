#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>

//
// auxiliary functions
//
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
double read_timer( )
{
  static bool initialized = false;
  static struct timeval start, end;
  if( !initialized )
    {
      gettimeofday( &start, NULL );
      initialized = true;
    }
  gettimeofday( &end, NULL );
  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

int find_option( int argc, char **argv, const char *option )
{
  for( int i = 1; i < argc; i++ )
    if( strcmp( argv[i], option ) == 0 )
      return i;
  return -1;
}


char *read_string( int argc, char **argv, const char *option, char *default_value )
{
  int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
      return argv[iplace+1];
    return default_value;
}



//
//  solvers
//
int build_table( int nitems, int cap, int *T, int *w, int *v )
{
  int wj = w[0], vj = v[0];

  for( int i = 0;  i <  wj;  i++ ) T[i] = 0;
  for( int i = wj; i <= cap; i++ ) T[i] = vj;
    
  for( int j = 1; j < nitems; j++ ) 
    {
      int wj = w[j], vj = v[j];
      for( int i = 0;  i <  wj;  i++ ) T[i+cap+1] = T[i];
      for( int i = wj; i <= cap; i++ ) T[i+cap+1] = max( T[i], T[i-wj] + vj );
        
      T += cap+1;
    }
    
  return T[cap]; // Max value attained
}

void backtrack( int nitems, int cap, int *T, int *w, int *u )
{
  int i = nitems*(cap+1) - 1;
  for( int j = nitems-1; j > 0; j-- )
    {
      u[j] = T[i] != T[i-cap-1]; // Is the jth item used?
      i -= cap+1 + (u[j] ? w[j] : 0 );
    }
  u[0] = T[i] != 0;
}

//
//  benchmarking program
//
int main (int argc, char** argv)
{

  char *savename = read_string( argc, argv, "-o", NULL );
  FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
  
  //srand48( (unsigned int)time(NULL) ); // Generates a sequence of 48-bit ints
  srand48( 1000 ); // Modified srand 

  //these constants have little effect on runtime
  int max_value  = 1000;
  int max_weight = 1000;
    
  //these set the problem size
  int capacity   = 2999; //9; //999; // Max weight that bag can hold 
  int nitems     = 5000; // 100; //5000;
    
  //allocate arrays
  int *weight = (int*)malloc( nitems * sizeof(int) );
  int *value  = (int*)malloc( nitems * sizeof(int) );
  int *used   = (int*)malloc( nitems * sizeof(int) );
  int *total  = (int*)malloc( nitems * (capacity+1) * sizeof(int) );
  if( !weight || !value || !used || !total )
    {
      fprintf( stderr, "Failed to allocate memory" );
      return -1;
    }
    
  // init
  max_weight = min( max_weight, capacity );//don't generate items that don't fit
  for( int i = 0; i < nitems; i++ )
    {
      weight[i] = 1 + (lrand48()%max_weight);
      value[i]  = 1 + (lrand48()%max_value);
    }
    
  // time the solution
  double seconds = read_timer( );

  int best_value = build_table( nitems, capacity, total, weight, value );
  backtrack( nitems, capacity, total, weight, used ); // total,weight,used arrays may be changed

  seconds = read_timer( ) - seconds;
        
  printf( "%d items, capacity: %d, time: %g\n", nitems, capacity, seconds );

  //sanity check
  int total_weight = 0, nused = 0, total_value = 0;
  for( int i = 0; i < nitems; i++ )
    if( used[i] )
      {
	nused++;
	total_weight += weight[i];
	total_value += value[i];
      }
  if( best_value != total_value || total_weight > capacity )
    printf( "INVALID SOLUTION\n" );
    
  printf( "%d items used, value %d, weight %d\n", nused, total_value, total_weight );

  if( fsave ) {
    fprintf(fsave, "%d items used, value %d, weight %d\n", nused, total_value, total_weight );

    for (int j=0; j<nitems; j++) {
      fprintf( fsave, "Index %d: %d %d\n", j, weight[j], value[j]);
    }
//    for (int j=0; j<(nitems * (capacity+1)); j++) {
//      fprintf( fsave, "Index %d: %d\n", total[j]); // Print total array
//    }

    fclose( fsave );

  }


  //release resources
  free( weight );
  free( value );
  free( used );
  free( total );

  return 0;
}
