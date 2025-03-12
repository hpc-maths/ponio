#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int lorenz(double t, const double y[], double dy[], void *params)
{
  double sigma = ((double *)params)[0];
  double rho = ((double *)params)[1];
  double beta = ((double *)params)[2];

  dy[0] = sigma * (y[1] - y[0]);
  dy[1] = y[0] * (rho - y[2]) - y[1];
  dy[2] = y[0] * y[1] - beta * y[2];

  return GSL_SUCCESS;
}

int
main()
{
  double params[3] = { 10.0, 28.0, 8.0 / 3.0 };
  double y0[3] = {1.0, 1.0, 1.0};
  double t0 = 0.0;
  double tf = 20.0;
  double dt = 0.01;

  gsl_odeiv2_system lorenz_pb = {lorenz, NULL, 3, &params};
  gsl_odeiv2_driver * meth = gsl_odeiv2_driver_alloc_y_new (&lorenz_pb, gsl_odeiv2_step_rk4, dt, 1e-6, 0.0);

  int n_sol = (tf-t0)/dt;
  double ** sol = (double **)malloc(n_sol*sizeof(double[4]));
  int i=0;
  for (i=0; i<n_sol; ++i)
  {
    sol[i] = (double*)malloc(sizeof(double[4]));
  }

  double yn[3];
  memcpy(yn, y0, 3*sizeof(double));

  int n = 0;
  double tn = t0;
  while( tn < tf )
  {
    sol[n][0] = tn;
    sol[n][1] = yn[0];
    sol[n][2] = yn[1];
    sol[n][3] = yn[2];

    double tnp1 = tn + dt;
    int status = gsl_odeiv2_driver_apply(meth, &tn, tnp1, yn);
    if (status != GSL_SUCCESS)
    {
      printf ("error, return value=%d\n", status);
      break;
    }

    ++n;
  }

  const char* filename = "lorenz.txt";
  FILE* output_file = fopen(filename, "w");
  if (!output_file) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  for (i=0; i<n_sol; ++i)
  {
    fprintf(output_file, "%f %f %f %f\n", sol[i][0], sol[i][1], sol[i][2], sol[i][3] );
  }
  fclose(output_file);

  for (i=0; i<n_sol; ++i)
  {
    free(sol[i]);
  }
  free(sol);

  return 0;
}
