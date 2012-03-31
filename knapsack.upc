#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>
#include <string.h>

//
// auxiliary functions
//
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
double read_timer( )
{
    static int initialized = 0;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = 1;
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
int build_table( int nitems, int cap, shared int *T, shared int *w, shared int *v )
{
    int wj, vj;
    
    wj = w[0];
    vj = v[0];
    upc_forall( int i = 0;  i <  wj;  i++; &T[i] ) T[i] = 0;
    upc_forall( int i = wj; i <= cap; i++; &T[i] ) T[i] = vj;
    upc_barrier;
    
    for( int j = 1; j < nitems; j++ )
    {
        wj = w[j];
        vj = v[j];
        upc_forall( int i = 0;  i <  wj;  i++; &T[i] ) T[i+cap+1] = T[i];
        upc_forall( int i = wj; i <= cap; i++; &T[i] ) T[i+cap+1] = max( T[i], T[i-wj]+vj );
        upc_barrier;
        
        T += cap+1;
    }
    
    return T[cap];
}

void backtrack( int nitems, int cap, shared int *T, shared int *w, shared int *u )
{
    int i, j;
    
    if( MYTHREAD != 0 )
        return;
    
    i = nitems*(cap+1) - 1;
    for( j = nitems-1; j > 0; j-- )
    {
        u[j] = T[i] != T[i-cap-1];
        i -= cap+1 + (u[j] ? w[j] : 0 );
    }
    u[0] = T[i] != 0;
}

//
//  serial solver to check correctness
//
int solve_serial( int nitems, int cap, shared int *w, shared int *v )
{
    int i, j, best, *allocated, *T, wj, vj;
    
    //alloc local resources
    T = allocated = malloc( nitems*(cap+1)*sizeof(int) );
    if( !allocated )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }
    
    //build_table locally
    wj = w[0];
    vj = v[0];
    for( i = 0;  i <  wj;  i++ ) T[i] = 0;
    for( i = wj; i <= cap; i++ ) T[i] = vj;
    for( j = 1; j < nitems; j++ ) 
    {
        wj = w[j];
        vj = v[j];
        for( i = 0;  i <  wj;  i++ ) T[i+cap+1] = T[i];
        for( i = wj; i <= cap; i++ ) T[i+cap+1] = max( T[i], T[i-wj]+vj );
        T += cap+1;
    }
    best = T[cap];
    
    //free resources
    free( allocated );
    
    return best;
}


//shared int* sharedMem = upc_all_alloc(100, sizeof(int));


//
//  benchmarking program
//
shared [3] int* foo;

int main( int argc, char** argv )
{
   char *savename = read_string( argc, argv, "-o", NULL );
   FILE *fsave = savename ? fopen( savename, "w" ) : NULL;

    printf("1: No. of threads is %d \n", THREADS);
    printf("2: Current thread is %d \n", MYTHREAD);
 
    int i, best_value, best_value_serial, total_weight, nused, total_value;
    double seconds;
    shared int *weight;
    shared int *value;
    shared int *used;
    shared int *total;
    
    //these constants have little effect on runtime
    int max_value  = 1000;
    int max_weight = 1000;
    
    //these set the problem size
    int capacity   = 9; //999;
    int nitems     = 100; //5000;
    
    //srand48( (unsigned int)time(NULL) + MYTHREAD );
    srand48( 1000 );    

    //allocate distributed arrays, use cyclic distribution
    weight = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    value  = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    used   = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    total  = (shared int *) upc_all_alloc( nitems * (capacity+1), sizeof(int) );
    if( !weight || !value || !total || !used )
    {
fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }



    //FIX: 
    upc_barrier;
    
    // init
    max_weight = min( max_weight, capacity );//don't generate items that don't fit into bag
    upc_forall( i = 0; i < nitems; i++; i )
    {

        printf("2: No. of threads is %d \n", THREADS);

        weight[i] = 1 + (lrand48()%max_weight);
        value[i]  = 1 + (lrand48()%max_value);
    }



    
    // time the solution
    seconds = read_timer( );
    
    best_value = build_table( nitems, capacity, total, weight, value );
    backtrack( nitems, capacity, total, weight, used );
    
    seconds = read_timer( ) - seconds;





    
    // check the result
    if( MYTHREAD == 0 )
    {
        printf( "%d items, capacity: %d, time: %g\n", nitems, capacity, seconds );
        
        best_value_serial = solve_serial( nitems, capacity, weight, value );
        
        total_weight = nused = total_value = 0;
        for( i = 0; i < nitems; i++ )
            if( used[i] )
            {
                nused++;
                total_weight += weight[i];
                total_value += value[i];
            }
        
        printf( "%d items used, value %d, weight %d\n", nused, total_value, total_weight );
        
        if( best_value != best_value_serial || best_value != total_value || total_weight > capacity )
            printf( "WRONG SOLUTION\n" );


    if( fsave ) {
      fprintf(fsave, "%d items used, value %d, weight %d\n", nused, total_value, total_weight );
      //for (int j=0; j<nitems; j++) {
      //  fprintf( fsave, "Index %d: %d %d %d\n", j, used[j], weight[j], value[j]);
      //}
      //for (int j=0; j<(nitems * (capacity+1)); j++) {
      //  fprintf( fsave, "Index %d: %d\n", total[j]); // Print total array
      //}
      fclose( fsave );
    }




        
        //release resources
        upc_free( weight );
        upc_free( value );
        upc_free( total );
        upc_free( used );
    }
    
    return 0;
}
