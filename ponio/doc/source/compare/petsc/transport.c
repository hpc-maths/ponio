#include <petsc.h>

PetscErrorCode
transport( TS ts, double t, Vec vec_y, Vec vec_dy, void* params )
{
    double a  = ( (double*)params )[0];
    double dx = ( (double*)params )[1];
    int n_x   = ( (double*)params )[2];

    double const* y;
    double* dy;

    VecGetArrayRead( vec_y, &y );
    VecGetArray( vec_dy, &dy );

    dy[0] = -( fmax( a, 0. ) * ( y[0] - y[n_x - 1] ) + fmin( a, 0. ) * ( y[1] - y[0] ) ) / dx;

    for ( int i = 1; i < n_x - 1; ++i )
    {
        dy[i] = -( fmax( a, 0. ) * ( y[i] - y[i - 1] ) + fmin( a, 0. ) * ( y[i + 1] - y[i] ) ) / dx;
    }

    dy[n_x - 1] = -( fmax( a, 0. ) * ( y[n_x - 1] - y[n_x - 2] ) + fmin( a, 0. ) * ( y[0] - y[n_x - 1] ) ) / dx;

    VecRestoreArray( vec_dy, &dy );
    VecRestoreArrayRead( vec_y, &y );

    return PETSC_SUCCESS;
}

int
main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, NULL, NULL );

    // space parameter
    int n_x   = 500;
    double* x = (double*)malloc( n_x * sizeof( double ) );
    for ( int i = 0; i < n_x; ++i )
    {
        x[i] = (double)i / n_x;
    }
    double dx = x[1] - x[0];

    // velocity
    double a = 1.0;

    // time parameter
    double t0 = 0.;
    double tf = 0.3;
    double dt = dx / a;

    // initial condition
    double* y0 = (double*)malloc( n_x * sizeof( double ) );
    for ( int i = 0; i < n_x; ++i )
    {
        y0[i] = 0;
        if ( 0.25 <= x[i] && x[i] < 0.5 )
        {
            y0[i] = x[i] - 0.25;
        }
        if ( 0.5 <= x[i] && x[i] < 0.75 )
        {
            y0[i] = -x[i] + 0.75;
        }
    }

    double params[3] = { a, dx, (double)n_x };

    int n_sol    = ( tf - t0 ) / dt;
    double** sol = (double**)malloc( n_sol * sizeof( double* ) );
    int i        = 0;
    for ( i = 0; i < n_sol; ++i )
    {
        sol[i] = (double*)malloc( n_x * sizeof( double ) );
    }

    TS ts;
    TSCreate( PETSC_COMM_SELF, &ts );
    TSSetType( ts, TSRK );
    TSRKSetType( ts, TSRK1FE );

    TSSetTime( ts, t0 );
    TSSetMaxTime( ts, tf );
    TSSetTimeStep( ts, dt );

    TSSetRHSFunction( ts, NULL, transport, params );

    Vec vec_y0;
    VecCreateSeqWithArray( PETSC_COMM_SELF, 1, n_x, y0, &vec_y0 );
    TSSetSolution( ts, vec_y0 );

    int n_step;
    double tn;
    TSGetTime( ts, &tn );
    TSGetStepNumber( ts, &n_step );

    Vec vec_yn;
    double const* yn;

    while ( tn < tf )
    {
        TSGetStepNumber( ts, &n_step );
        TSGetSolution( ts, &vec_yn );
        VecGetArrayRead( vec_yn, &yn );

        memcpy( sol[n_step], yn, n_x * sizeof( double ) );

        VecRestoreArrayRead( vec_yn, &yn );

        TSStep( ts );

        TSGetTime( ts, &tn );
    }

    // save solution
    char const* filename = "transport.txt";
    FILE* output_file    = fopen( filename, "w" );
    if ( !output_file )
    {
        perror( "fopen" );
        exit( EXIT_FAILURE );
    }

    for ( int i = 0; i < n_x; ++i )
    {
        fprintf( output_file, "%f", x[i] );
        for ( int n = 0; n < n_sol; ++n )
        {
            fprintf( output_file, " %f", sol[n][i] );
        }
        fprintf( output_file, "\n" );
    }
    fclose( output_file );

    for ( int i = 0; i < n_sol; ++i )
    {
        free( sol[i] );
    }
    free( sol );

    VecDestroy( &vec_yn );
    VecDestroy( &vec_y0 );

    PetscFinalize();

    return 0;
}
