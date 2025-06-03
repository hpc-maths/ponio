#include <petsc.h>

PetscErrorCode arenstorf(TS ts, double t, Vec vec_y, Vec vec_ydot, void* params)
{
    double mu = *(double *)params;

    const double* y;
    double* ydot;

    VecGetArrayRead(vec_y, &y);
    VecGetArray(vec_ydot, &ydot);

    double const y1 = y[0];
    double const y2 = y[1];
    double const y3 = y[2];
    double const y4 = y[3];

    double const r1 = sqrt( ( y1 + mu ) * ( y1 + mu ) + y2 * y2 );
    double const r2 = sqrt( ( y1 - 1. + mu ) * ( y1 - 1. + mu ) + y2 * y2 );

    ydot[0] = y3;
    ydot[1] = y4;
    ydot[2] = y1 + 2. * y4 - ( 1. - mu ) * ( y1 + mu ) / ( r1 * r1 * r1 ) - mu * ( y1 - 1. + mu ) / ( r2 * r2 * r2 );
    ydot[3] = y2 - 2. * y3 - ( 1. - mu ) * y2 / ( r1 * r1 * r1 ) - mu * y2 / ( r2 * r2 * r2 );

    VecRestoreArray(vec_ydot, &ydot);
    VecRestoreArrayRead(vec_y, &y);

    return PETSC_SUCCESS;
}

int
main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, NULL, NULL);

  double mu = 0.012277471;
  double y0[4] = { 0.994, 0., 0., -2.00158510637908252240537862224 };
  double t0 = 0.0;
  double tf = 17.0652165601579625588917206249;
  double dt = 1e-5;

  int max_n_sol = 100;
  double ** sol = (double **)malloc(max_n_sol*sizeof(double[6]));
  for (int i=0; i<max_n_sol; ++i)
  {
    sol[i] = (double*)malloc(sizeof(double[6]));
  }

  TS ts;
  TSCreate(PETSC_COMM_SELF, &ts);
  TSSetType(ts, TSRK);
  TSRKSetType(ts, TSRK5DP);
  TSSetTolerances(ts, 1e-5, NULL, 1e-5, NULL);

  TSSetTime(ts, 0.);
  TSSetMaxTime(ts, tf);
  TSSetTimeStep(ts, dt);

  TSSetRHSFunction(ts, NULL, arenstorf, &mu);

  Vec vec_y0;
  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 4, y0, &vec_y0);
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
    TSGetTimeStep(ts, &dt);
    TSGetSolution(ts, &vec_yn);
    VecGetArray(vec_yn, &yn);

    sol[n_step][0] = tn;
    sol[n_step][1] = yn[0];
    sol[n_step][2] = yn[1];
    sol[n_step][3] = yn[2];
    sol[n_step][4] = yn[3];
    sol[n_step][5] = dt;

    VecRestoreArray(vec_yn, &yn);

    TSStep(ts);

    TSGetTime(ts, &tn);
  }

  const char* filename = "arenstorf.txt";
  FILE* output_file = fopen(filename, "w");
  if (!output_file) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<n_step; ++i)
  {
    fprintf(output_file, "%f %f %f %f\n", sol[i][0], sol[i][1], sol[i][2], sol[i][3] );
  }
  fclose(output_file);

  for (int i=0; i<max_n_sol; ++i)
  {
    free(sol[i]);
  }
  free(sol);

  VecDestroy(&vec_yn);
  VecDestroy(&vec_y0);

  PetscFinalize();

  return 0;
}
