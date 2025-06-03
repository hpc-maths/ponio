#include <stdio.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int arenstorf(double t, const double y[], double dy[], void *params)
{
  double mu = *(double *)params;

  double const y1 = y[0];
  double const y2 = y[1];
  double const y3 = y[2];
  double const y4 = y[3];

  double const r1 = sqrt( ( y1 + mu ) * ( y1 + mu ) + y2 * y2 );
  double const r2 = sqrt( ( y1 - 1. + mu ) * ( y1 - 1. + mu ) + y2 * y2 );

  dy[0] = y3;
  dy[1] = y4;
  dy[2] = y1 + 2. * y4 - ( 1. - mu ) * ( y1 + mu ) / ( r1 * r1 * r1 ) - mu * ( y1 - 1. + mu ) / ( r2 * r2 * r2 );
  dy[3] = y2 - 2. * y3 - ( 1. - mu ) * y2 / ( r1 * r1 * r1 ) - mu * y2 / ( r2 * r2 * r2 );

  return GSL_SUCCESS;
}

int
main()
{
  double mu = 0.012277471;
  double y0[4] = { 0.994, 0., 0., -2.00158510637908252240537862224 };
  double t0 = 0.0;
  double tf = 17.0652165601579625588917206249;
  double dt = 1e-5;

  gsl_odeiv2_system arenstorf_pb = {arenstorf, NULL, 4, &mu};
  const gsl_odeiv2_step_type * meth = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(meth, 4);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-5, 1e-5);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(4);

  int max_n_sol = 100;
  double ** sol = (double **)malloc(max_n_sol*sizeof(double[6]));
  int i=0;
  for (i=0; i<max_n_sol; ++i)
  {
    sol[i] = (double*)malloc(sizeof(double[6]));
  }

  double yn[4];
  memcpy(yn, y0, 4*sizeof(double));

  int n = 0;
  double tn = t0;
  while( tn < tf )
  {
    sol[n][0] = tn;
    sol[n][1] = yn[0];
    sol[n][2] = yn[1];
    sol[n][3] = yn[2];
    sol[n][4] = yn[3];
    sol[n][5] = dt;

    int status = gsl_odeiv2_evolve_apply(e, c, s, &arenstorf_pb, &tn, tf, &dt, yn);

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value=%d\n", status);
      break;
    }

    ++n;
  }

  const char* filename = "arenstorf.txt";
  FILE* output_file = fopen(filename, "w");
  if (!output_file) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  for (i=0; i<n; ++i)
  {
    fprintf(output_file, "%f %f %f %f %f %f\n", sol[i][0], sol[i][1], sol[i][2], sol[i][3], sol[i][4], sol[i][5] );
  }
  fclose(output_file);

  for (i=0; i<max_n_sol; ++i)
  {
    free(sol[i]);
  }
  free(sol);

  return 0;
}
