#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <upc.h>
#include <upc_collective.h>
#include <string.h>

#define BLK_SIZE 42 

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


// Once row 1 is done, row 2 can start.
// Split based on capacity
// Each thread will work on cap/THREADS no. of elements
int build_table_local( int nitems, int cap, int padCap, shared [BLK_SIZE] int *T, int *Tlocal, int *w, int *v )
{
    int wj, vj;
    
    wj = w[0];
    vj = v[0];

    // Initialization stage    
    //int interval = (padCap/THREADS)+1; // Replaced with BLK_SIZE 
    int startIdx = BLK_SIZE*MYTHREAD;

    memset ((void*) (Tlocal+startIdx), 0, max(min(BLK_SIZE,wj-startIdx), 0)*sizeof(int));


    if (wj < startIdx) {
        memset ((void*) (Tlocal+startIdx), vj, BLK_SIZE*sizeof(int));
    }
    else if (wj < startIdx + BLK_SIZE) {
         memset ((void*) (Tlocal+wj), vj, (startIdx+BLK_SIZE-wj)*sizeof(int));
    }

    upc_memput((shared void *) (T+startIdx), (void *) (Tlocal+startIdx), BLK_SIZE*sizeof(int));

/*
    for( int j = 1; j < nitems-1; j+=2 )
    {
        int w1 = w[j];
        int v1 = v[j];

	int w2 = w[j+1];
	int v2 = v[j+1];
	
	upc_memget( (void*) (Tlocal+max(startIdx-w1-w2, 0)), (shared void *) (T+max(startIdx-w1-w2,0)), min(w1+w2,startIdx)*sizeof(int)); 

	int cpylength = max(min(BLK_SIZE+w2,w1-startIdx+w2), 0);
	memcpy ((void*) (Tlocal+max(startIdx-w2,0)+padCap+1), (void *) (Tlocal+max(startIdx-w2,0)), cpylength*sizeof(int));

        for( int i = max(startIdx-w2, w1); i < startIdx+BLK_SIZE; i++ ) {
	    Tlocal[i+padCap+1] = max( Tlocal[i], Tlocal[i-w1]+v1 );
	}

	upc_memput((shared void *) (T+startIdx+padCap+1), (void *) (Tlocal+startIdx+padCap+1), BLK_SIZE*sizeof(int));

        memcpy ((void*) (Tlocal+startIdx+2*padCap+2), (void *) (Tlocal+startIdx+padCap+1), max(min(BLK_SIZE,w2-startIdx), 0)*sizeof(int));

        for( int i = max(startIdx, w2); i < startIdx+BLK_SIZE; i++ ) {
            Tlocal[i+2*padCap+2] = max( Tlocal[i+padCap+1], Tlocal[i-w2+padCap+1]+v2 );
        }



	upc_memput((shared void *) (T+startIdx+2*padCap+2), (void *) (Tlocal+startIdx+2*padCap+2), BLK_SIZE*sizeof(int));

        upc_barrier;
        
        T += 2*padCap+2;
        Tlocal += 2*padCap+2;
    }
*/

    for( int j = 1; j < nitems; j++ )
    {
        wj = w[j];
        vj = v[j];

        // 1st UPC for loop
        memcpy ((void*) (Tlocal+startIdx+padCap+1), (void *) (Tlocal+startIdx), max(min(BLK_SIZE,wj-startIdx), 0)*sizeof(int));

        //upc_barrier;

        // 2nd UPC for loop
        upc_memget( (void*) (Tlocal+max(startIdx-wj, 0)), (shared void *) (T+max(startIdx-wj,0)), min(wj,startIdx)*sizeof(int));

        for( int i = max(startIdx, wj); i < startIdx+BLK_SIZE; i++ ) {
            Tlocal[i+padCap+1] = max( Tlocal[i], Tlocal[i-wj]+vj );
        }

        upc_memput((shared void *) (T+startIdx+padCap+1), (void *) (Tlocal+startIdx+padCap+1), BLK_SIZE*sizeof(int));

        upc_barrier;

        T += padCap+1;
        Tlocal += padCap+1;
    }


    return T[cap];
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
        upc_forall( int i = 0;  i <  wj;  i++; &T[i] ) { 
	  T[i+cap+1] = T[i];
	  //printf("Thread no. is %d \n", upc_threadof(&T[i]));
	}
        upc_forall( int i = wj; i <= cap; i++; &T[i] ) {
	  T[i+cap+1] = max( T[i], T[i-wj]+vj );
	}
        upc_barrier;
        
        T += cap+1;
    }
    
    return T[cap];
}

void backtrack( int nitems, int cap, shared [BLK_SIZE] int *T, shared int *w, shared int *u )
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


//
//  benchmarking program
//

int main( int argc, char** argv )
{

    shared int *weight;
    shared int *value;
    shared int *used;
    shared [BLK_SIZE] int *total=NULL;

    int* local;
    shared [BLK_SIZE] int *global=NULL;

   char *savename = read_string( argc, argv, "-o", NULL );
   FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
 
    int i, best_value, best_value_serial, total_weight, nused, total_value;
    double seconds;
    
    //these constants have little effect on runtime
    int max_value  = 1000;
    int max_weight = 1000;
    
    //these set the problem size
    int capacity   = 999; //9; //999;
    int nitems     = 5000; //100; //5000;
    double mult = ceil(((double)capacity+1)/(double)THREADS);
    int padCapacity = (int)mult*THREADS - 1; // Padded capacity should be a multiple of thread no.

    if (MYTHREAD == 0) {
       printf ("mult is % d\n", (int)mult);
       printf ("No. of threads is % d\n", THREADS);
       printf ("Padded capacity is %d \n", padCapacity);   
       printf ("Block size should be %d \n", (padCapacity+1)/THREADS);
    }

    //srand48( (unsigned int)time(NULL) + MYTHREAD );
    srand48( 1000 );    

    // Allocate distributed arrays, use cyclic distribution
    weight = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    value  = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    used   = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    total  = (shared [BLK_SIZE] int *) upc_all_alloc( nitems * (padCapacity+1), sizeof(int) );
    if( !weight || !value || !total || !used )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    upc_barrier;
    
    // Allocate local arrays, use blocking 
    int *weightLoc = (int*) malloc( nitems * sizeof(int) );
    int *valueLoc  = (int*) malloc( nitems * sizeof(int) );
    int *usedLoc   = (int*) malloc( nitems * sizeof(int) );
    int *totalLoc  = (int*) malloc( nitems * (padCapacity+1) * sizeof(int) );
    if( !weightLoc || !valueLoc || !totalLoc || !usedLoc )
    {
        fprintf( stderr, "Failed to allocate local memory" );
        upc_global_exit( -1 );
    }

    upc_barrier;

 
    // 
    // Init. Prepare arrays in thread 0
    //
    max_weight = min( max_weight, padCapacity );  //do not generate items that don't fit into bag
    if (MYTHREAD == 0) {
       for( int i = 0; i < nitems; i++ ) {
         weight[i] = 1 + (lrand48()%max_weight);
         value[i]  = 1 + (lrand48()%max_value);
       }
    }

    upc_barrier;

    for (int j=0; j<nitems; j++) {
      weightLoc[j] = weight[j];
      valueLoc[j] = value[j];
    }

    upc_barrier;

    
    // time the solution
    // best_value = build_table( nitems, padCapacity, total, weight, value );
    seconds = read_timer( );

    best_value = build_table_local(nitems, capacity, padCapacity, total, totalLoc, weightLoc, valueLoc );
    backtrack( nitems, padCapacity, total, weight, used );
    
    seconds = read_timer( ) - seconds;

    upc_barrier;
 
    // check the result
    if( MYTHREAD == 0 )
    {
        printf( "%d items, capacity: %d, time: %g\n", nitems, padCapacity, seconds );
        
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
        
        if( best_value != best_value_serial )
            printf( "WRONG SOLUTION 1: Best val frm parallel not equal to best val from serial... best_value %d, best_serial %d \n", best_value, best_value_serial );
	
        if( best_value != total_value )
            printf( "WRONG SOLUTION 2: Best val not equal to total val (Reasonable for padded capacity case).. best_value %d, total_value %d\n", best_value, total_value );

        if( total_weight > padCapacity )
            printf( "WRONG SOLUTION 3: Total weight not larger than padded capacity.. total_weight %d, padCapacity %d \n", total_weight, padCapacity );
    if( fsave ) {
      fprintf(fsave, "%d items used, value %d, weight %d\n", nused, total_value, total_weight );

      if (MYTHREAD == 0) {
        for (int j=0; j<nitems; j++) {
          fprintf( fsave, "Index %d: %d %d\n", j, weightLoc[j], valueLoc[j]);
        }
      }


      //for (int j=0; j<nitems; j++) {
      //  fprintf( fsave, "Index %d: %d %d\n", j, weight[j], value[j]);
      //}
      //for (int j=0; j<(nitems * (padCapacity+1)); j++) {
      //  fprintf( fsave, "Index %d: %d\n", total[j]); // Print total array
      //}
      fclose( fsave );
    }




        
        //release resources
        upc_free( weight );
        upc_free( value );
        upc_free( total );
        upc_free( used );

	free( weightLoc );
	free( valueLoc );
	free( totalLoc );
	free( usedLoc );
    }
    
    return 0;
}
