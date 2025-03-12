#include <petsc.h>

PetscErrorCode lorenz(TS ts, double t, Vec vec_y, Vec vec_ydot, void* params)
{
    double sigma = ((double *)params)[0];
    double rho = ((double *)params)[1];
    double beta = ((double *)params)[2];

    const double* y;
    double* ydot;

    VecGetArrayRead(vec_y, &y);
    VecGetArray(vec_ydot, &ydot);

    ydot[0] = sigma*( y[1] - y[0] );
    ydot[1] = rho*y[0] - y[1] - y[0]*y[2];
    ydot[2] = y[0]*y[1] - beta*y[2];

    VecRestoreArray(vec_ydot, &ydot);
    VecRestoreArrayRead(vec_y, &y);

    return PETSC_SUCCESS;
}

int
main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, NULL, NULL);

  double params[3] = { 10.0, 28.0, 8.0 / 3.0 };
  double y0[3] = {1.0, 1.0, 1.0};
  double t0 = 0.0;
  double tf = 20.0;
  double dt = 0.01;

  int n_sol = (tf-t0)/dt;
  double ** sol = (double **)malloc(n_sol*sizeof(double[4]));
  for (int i=0; i<n_sol; ++i)
  {
    sol[i] = (double*)malloc(sizeof(double[4]));
  }

  TS ts;
  TSCreate(PETSC_COMM_SELF, &ts);
  TSSetType(ts, TSRK);
  TSRKSetType(ts, TSRK4);

  TSSetTime(ts, 0.);
  TSSetMaxTime(ts, tf);
  TSSetTimeStep(ts, dt);

  TSSetRHSFunction(ts, NULL, lorenz, params);

  Vec vec_y0;
  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 3, y0, &vec_y0);
  TSSetSolution(ts, vec_y0);

  int n_step;
  double tn;
  TSGetTime(ts, &tn);
  TSGetStepNumber(ts, &n_step);

  Vec vec_yn;
  double * yn;

  while( tn < tf )
  {
    TSGetStepNumber(ts, &n_step);
    TSGetSolution(ts, &vec_yn);
    VecGetArray(vec_yn, &yn);

    sol[n_step][0] = tn;
    sol[n_step][1] = yn[0];
    sol[n_step][2] = yn[1];
    sol[n_step][3] = yn[2];

    VecRestoreArray(vec_yn, &yn);

    TSStep(ts);

    TSGetTime(ts, &tn);
  }

  const char* filename = "lorenz.txt";
  FILE* output_file = fopen(filename, "w");
  if (!output_file) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<n_sol; ++i)
  {
    fprintf(output_file, "%f %f %f %f\n", sol[i][0], sol[i][1], sol[i][2], sol[i][3] );
  }
  fclose(output_file);

  for (int i=0; i<n_sol; ++i)
  {
    free(sol[i]);
  }
  free(sol);

  VecDestroy(&vec_yn);
  VecDestroy(&vec_y0);

  PetscFinalize();

  return 0;
}
