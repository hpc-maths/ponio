#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

static int
stepper_set_driver_null( void* vstate, gsl_odeiv2_driver const* d )
{
    return GSL_SUCCESS;
}

typedef struct
{
    double* k;
    double* y0;
    double* y_onestep;
} euler_state_t;

static void*
euler_alloc( size_t dim )
{
    euler_state_t* state = (euler_state_t*)malloc( sizeof( euler_state_t ) );

    if ( state == 0 )
    {
        GSL_ERROR_NULL( "failed to allocate space for euler_state", GSL_ENOMEM );
    }

    state->k = (double*)malloc( dim * sizeof( double ) );

    if ( state->k == 0 )
    {
        free( state );
        GSL_ERROR_NULL( "failed to allocate space for k", GSL_ENOMEM );
    }

    state->y0 = (double*)malloc( dim * sizeof( double ) );

    if ( state->y0 == 0 )
    {
        free( state->k );
        free( state );
        GSL_ERROR_NULL( "failed to allocate space for y0", GSL_ENOMEM );
    }

    state->y_onestep = (double*)malloc( dim * sizeof( double ) );

    if ( state->y_onestep == 0 )
    {
        free( state->y0 );
        free( state->k );
        free( state );
        GSL_ERROR_NULL( "failed to allocate space for ytmp", GSL_ENOMEM );
    }

    return state;
}

static int
euler_apply( void* vstate,
    size_t dim,
    double t,
    double h,
    double y[],
    double yerr[],
    double const dydt_in[],
    double dydt_out[],
    gsl_odeiv2_system const* sys )
{
    euler_state_t* state = (euler_state_t*)vstate;

    size_t i;

    double* const k = state->k;

    if ( dydt_in != NULL )
    {
        // DBL_MEMCPY (k, dydt_in, dim);
        for ( i = 0; i < dim; i++ )
        {
            k[i] = dydt_in[i];
        }
    }
    else
    {
        int s = GSL_ODEIV_FN_EVAL( sys, t, y, k );

        if ( s != GSL_SUCCESS )
        {
            return s;
        }
    }

    for ( i = 0; i < dim; i++ )
    {
        y[i] = y[i] + h * k[i];
    }

    /* Derivatives at output */

    if ( dydt_out != NULL )
    {
        int s = GSL_ODEIV_FN_EVAL( sys, t + h, y, dydt_out );

        if ( s != GSL_SUCCESS )
        {
            return s;
        }
    }

    /* Error estimation */

    for ( i = 0; i < dim; i++ )
    {
        yerr[i] = 0.;
    }

    return GSL_SUCCESS;
}

static int
euler_reset( void* vstate, size_t dim )
{
    euler_state_t* state = (euler_state_t*)vstate;

    size_t i;
    for ( i = 0; i < dim; ++i )
    {
        state->k[i]         = 0.;
        state->y0[i]        = 0.;
        state->y_onestep[i] = 0.;
    }

    return GSL_SUCCESS;
}

static unsigned int
euler_order( void* vstate )
{
    euler_state_t* state = (euler_state_t*)vstate;
    state                = 0; /* prevent warnings about unused parameters */
    return 1;
}

static void
euler_free( void* vstate )
{
    euler_state_t* state = (euler_state_t*)vstate;
    free( state->k );
    free( state->y0 );
    free( state->y_onestep );
    free( state );
}

static gsl_odeiv2_step_type const euler_type = { "euler", // name
    1,                                                    // can use dydt_in
    1,                                                    // gives exact dydt_out
    &euler_alloc,
    &euler_apply,
    &stepper_set_driver_null,
    &euler_reset,
    &euler_order,
    &euler_free };

gsl_odeiv2_step_type const* gsl_odeiv2_step_euler = &euler_type;

int
upwind( double t, double const y[], double dy[], void* params )
{
    double a  = ( (double*)params )[0];
    double dx = ( (double*)params )[1];
    int n_x   = ( (double*)params )[2];

    dy[0] = -( fmax( a, 0. ) * ( y[0] - y[n_x - 1] ) + fmin( a, 0. ) * ( y[1] - y[0] ) ) / dx;

    for ( int i = 1; i < n_x - 1; ++i )
    {
        dy[i] = -( fmax( a, 0. ) * ( y[i] - y[i - 1] ) + fmin( a, 0. ) * ( y[i + 1] - y[i] ) ) / dx;
    }

    dy[n_x - 1] = -( fmax( a, 0. ) * ( y[n_x - 1] - y[n_x - 2] ) + fmin( a, 0. ) * ( y[0] - y[n_x - 1] ) ) / dx;

    return GSL_SUCCESS;
}

int
main()
{
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

    gsl_odeiv2_system transport_pb = { upwind, NULL, n_x, &params };
    gsl_odeiv2_driver* meth        = gsl_odeiv2_driver_alloc_y_new( &transport_pb, gsl_odeiv2_step_euler, dt, 1e-6, 0.0 );

    int n_sol    = ( tf - t0 ) / dt;
    double** sol = (double**)malloc( n_sol * sizeof( double* ) );
    int i        = 0;
    for ( i = 0; i < n_sol; ++i )
    {
        sol[i] = (double*)malloc( n_x * sizeof( double ) );
    }

    double* yn = y0;

    int n     = 0;
    double tn = t0;
    while ( tn < tf )
    {
        memcpy( sol[n], yn, n_x * sizeof( double ) );

        double tnp1 = tn + dt;
        int status  = gsl_odeiv2_driver_apply( meth, &tn, tnp1, yn );
        if ( status != GSL_SUCCESS )
        {
            printf( "error, return value=%d\n", status );
            break;
        }

        ++n;
    }

    // save solution
    char const* filename = "transport.txt";
    FILE* output_file    = fopen( filename, "w" );
    if ( !output_file )
    {
        perror( "fopen" );
        exit( EXIT_FAILURE );
    }

    for ( i = 0; i < n_x; ++i )
    {
        fprintf( output_file, "%f", x[i] );
        for ( n = 0; n < n_sol; ++n )
        {
            fprintf( output_file, " %f", sol[n][i] );
        }
        fprintf( output_file, "\n" );
    }
    fclose( output_file );

    for ( i = 0; i < n_sol; ++i )
    {
        free( sol[i] );
    }
    free( sol );

    return 0;
}
