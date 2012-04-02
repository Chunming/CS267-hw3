#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>
#include <upc_collective.h>
#include <string.h>

#define COUNT_PER_PE 4

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
int build_table_local( int nitems, int cap, shared [250] int *T, int *Tlocal, int *w, int *v )
{
    int wj, vj;
    
    wj = w[0];
    vj = v[0];

    // Initialization stage    
    int interval = (cap/THREADS)+1; 
    int startIdx = interval*MYTHREAD;
    int count = 0;

    memset ((void*) (Tlocal+startIdx), 0, min(startIdx+interval,wj)*sizeof(int));
    //for (int i=startIdx; i < min(wj,startIdx+interval); i++) Tlocal[i] = 0;


    if (wj < startIdx + interval) {    
	memset ((void*) (Tlocal+startIdx+wj), vj, (startIdx+interval)*sizeof(int));
    }

    //for (int i=wj; i <= min(cap,startIdx+interval); i++) Tlocal[i] = vj;

    upc_memput((shared void *) (T+startIdx), (void *) (Tlocal+startIdx), interval*sizeof(int));
    //for (int i=startIdx; i<(startIdx+interval); i++) T[i] = Tlocal[i];
    upc_barrier;

/*
    for( int j = 1; j < nitems; j++ )
    {
        wj = w[j];
        vj = v[j];
        upc_forall( int i = 0;  i <  wj;  i++; &T[i] ) T[i+cap+1] = T[i];
        upc_forall( int i = wj; i <= cap; i++; &T[i] ) T[i+cap+1] = max( T[i], T[i-wj]+vj );
        upc_barrier;

        T += cap+1;
    }
*/




    for( int j = 1; j < nitems; j++ )
    {
        wj = w[j];
        vj = v[j];
	
	// 1st UPC for loop

        //for( int i = startIdx; i <  min(startIdx+interval,wj); i++ ) {
	//  Tlocal[i+cap+1] = Tlocal[i];
	//}

	memcpy ((void*) (Tlocal+startIdx+cap+1), (void *) (Tlocal+startIdx), min(startIdx+interval,wj)*sizeof(int));

    	upc_barrier;

       // for( int i = startIdx; i <  min(startIdx+interval,wj); i++ ) {
//	  Tlocal[i+cap+1] = Tlocal[i];
//	}

	
	// 2nd UPC for loop
	upc_memget( (void*) (Tlocal+max(startIdx-wj, 0)), (shared void *) (T+max(startIdx-wj,0)), min(wj,startIdx)*sizeof(int)); 

        for( int i = startIdx; i <= min(startIdx+interval, cap); i++ ) {
	    if (wj <= i) {
	       Tlocal[i+cap+1] = max( Tlocal[i], Tlocal[i-wj]+vj );
	    }
	}


/*
	for( int i = wj; i <= min(startIdx+interval, cap); i++ ) {
          if((i-wj)/interval != MYTHREAD ) {
            Tlocal[i-wj] = T[i-wj];
            Tlocal[i+cap+1] = max( Tlocal[i], Tlocal[i-wj]+vj );
          }
          else {
            Tlocal[i+cap+1] = max( Tlocal[i], Tlocal[i-wj]+vj );
          }
        }

*/
	upc_memput((shared void *) (T+startIdx+cap+1), (void *) (Tlocal+startIdx+cap+1), interval*sizeof(int));
	
        //for (int i=startIdx; i<(startIdx+interval); i++) T[i+cap+1] = Tlocal[i+cap+1];
        //upc_memput( (shared void*) (global+MYTHREAD*COUNT_PER_PE), (void*) local, COUNT_PER_PE*sizeof(int) );
        //upc_memput( (shared void*) (&T[cap+1]+MYTHREAD*interval), (void*) Tlocal, interval*sizeof(int) );

        //for (int i=startIdx; i<(startIdx+interval); i++) T[i+cap+1] = Tlocal[i+cap+1];
        upc_memput( (shared void*) (T+startIdx+cap+1), (void*) (Tlocal+startIdx+cap+1), 250*sizeof(int) );

        upc_barrier;
        
        T += cap+1;
        Tlocal += cap+1;
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

void backtrack( int nitems, int cap, shared [250] int *T, shared int *w, shared int *u )
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
    shared [250] int *total=NULL;

    int* local;
    shared [250] int *global=NULL;




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
    
    //srand48( (unsigned int)time(NULL) + MYTHREAD );
    srand48( 1000 );    

    // Allocate distributed arrays, use cyclic distribution
    weight = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    value  = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    used   = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    total  = (shared [250] int *) upc_all_alloc( nitems * (capacity+1), sizeof(int) );
    if( !weight || !value || !total || !used )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    upc_barrier;

    //
    // Allocate local arrays, use blocking
    //
    int *weightLoc = (int*) malloc( nitems * sizeof(int) );
    int *valueLoc  = (int*) malloc( nitems * sizeof(int) );
    int *usedLoc   = (int*) malloc( nitems * sizeof(int) );
    int *totalLoc  = (int*) malloc( nitems * (capacity+1) * sizeof(int) );

    upc_barrier;

    // 
    // Test segment
    //

    local = (int *)malloc(sizeof(int)*5000*1000);
    for (int i=0;i<5000*1000;i++) { 
      local[i] = MYTHREAD*10000000+i;
    }
    upc_barrier;

    size_t nBytes = sizeof(int);
    global  = (shared [250] int *) upc_all_alloc( 5000*1000, nBytes ); // numBlocks, numBytes
    upc_barrier;

    //
    // Copy data from local to global
    //
    //upc_memput( (shared void*) (global+MYTHREAD*250), (void*) local, 250*sizeof(int) );
    //upc_barrier;

    for (int i=MYTHREAD*250; i<(MYTHREAD*250+250); i++) global[i] = 0;
    for( int j = 1; j < 10; j++ ) // 10 instead of nitems
    {
        for (int i=MYTHREAD*250; i<(MYTHREAD*250+250); i++) global[i+999+1] = local[i+999+1];
        //upc_memput( (shared void*) (global+MYTHREAD*COUNT_PER_PE), (void*) local, COUNT_PER_PE*sizeof(int) );
        //upc_memput( (shared void*) (&global[cap+1]+MYTHREAD*250), (void*) &local[cap+1], 250*sizeof(int) );

        upc_barrier;
        
        global += 999+1;
        local += 999+1;
    }

/*
    if (MYTHREAD == 0) {
       for( int i = 0; i < 10*1000; i++ ) { //10 instead of nitems
         //fprintf(fsave, "global at %d is %d \n", i, global[i]);
       }
    }

    upc_barrier; //FIX: 
  */
 
    // 
    // Init. Prepare arrays in thread 0
    //
    max_weight = min( max_weight, capacity );  //do not generate items that don't fit into bag
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
    seconds = read_timer( );
    
    best_value = build_table_local(nitems, capacity, total, totalLoc, weightLoc, valueLoc );
    //best_value = build_table( nitems, capacity, total, weight, value );
    backtrack( nitems, capacity, total, weight, used );
    
    seconds = read_timer( ) - seconds;

    upc_barrier;
    
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

      if (MYTHREAD == 0) {
        for (int j=0; j<nitems; j++) {
          fprintf( fsave, "Index %d: %d %d\n", j, weightLoc[j], valueLoc[j]);
        }
      }


      //for (int j=0; j<nitems; j++) {
      //  fprintf( fsave, "Index %d: %d %d\n", j, weight[j], value[j]);
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

	free( weightLoc );
	free( valueLoc );
	free( totalLoc );
	free( usedLoc );
    }
    
    return 0;
}
