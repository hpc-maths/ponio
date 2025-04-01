#include <stdio.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int upwind(double t, const double y[], double dy[], void *params)
{
  double a = ((double *)params)[0];
  double dx = ((double *)params)[1];
  int n_x = ((double *)params)[2];

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
  int n_x = 500;
  double * x = (double *)malloc(n_x * sizeof(double));
  for ( int i = 0 ; i<n_x; ++i)
  {
    x[i] = (double)i/n_x;
  }
  double dx = x[1] - x[0];

  // velocity
  double a = 1.0;

  // time parameter
  double t0 = 0.;
  double tf = 0.3;
  double dt = dx / a;

  // initial condition
  double * y0 = (double *)malloc(n_x * sizeof(double));
  for ( int i = 0 ; i<n_x; ++i)
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

  gsl_odeiv2_system transport_pb = {upwind, NULL, n_x, &params};
  gsl_odeiv2_driver * meth = gsl_odeiv2_driver_alloc_y_new(&transport_pb, gsl_odeiv2_step_rk2, dt, 1e-6, 0.0);

  int n_sol = (tf-t0)/dt;
  double ** sol = (double **)malloc(n_sol*sizeof(double *));
  int i=0;
  for (i=0; i<n_sol; ++i)
  {
    sol[i] = (double *)malloc(n_x * sizeof(double));
  }

  double * yn = y0;

  int n = 0;
  double tn = t0;
  while( tn < tf )
  {
    memcpy(sol[n], yn, n_x*sizeof(double));

    double tnp1 = tn + dt;
    int status = gsl_odeiv2_driver_apply(meth, &tn, tnp1, yn);
    if (status != GSL_SUCCESS)
    {
      printf ("error, return value=%d\n", status);
      break;
    }

    ++n;
  }

  // save solution
  const char* filename = "transport.txt";
  FILE* output_file = fopen(filename, "w");
  if (!output_file) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }

  for (i=0; i<n_x; ++i)
  {
    fprintf(output_file, "%f", x[i]);
    for (n=0; n<n_sol; ++n)
    {
      fprintf(output_file, " %f", sol[n][i]);
    }
    fprintf(output_file, "\n");
  }
  fclose(output_file);

  for (i=0; i<n_sol; ++i)
  {
    free(sol[i]);
  }
  free(sol);

  return 0;
}
