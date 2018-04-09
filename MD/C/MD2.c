/*
 *  Simple molecular dynamics code.
 *  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
 *
 * This program implements:
 *     long range inverse square forces between particles. F = G * m1*m2 / r**2
 *     viscosity term     F = -u V
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * Coordinates are relative to a large central mass and the entire system is moving relative to the
 * viscous media.
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * This program was developed as part of a code optimisation course
 * and is therefore deliberately inefficient.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"

void vis_force(int N, double *f, double *vis, double *vel);
void add_norm(int N, double *r, double *delta);
double force(double W, double delta, double r);
void wind_force(int N, double *f, double *vis, double vel);

void evolve(int count, double dt)
{
  int step;
  int i, j, k, l;
  int collided;
  double Size;
  /*
 * Loop over timesteps.
 */
  for (step = 1; step <= count; step++)
  {
    printf("timestep %d\n", step);
    printf("collisions %d\n", collisions);

    /* Optimization 2 */
    for (j = 0; j < Ndim; j++)
    {
      for (i = 0; i < Nbody; i++)
      {
        f[j][i] = -vis[i] * (velo[j][i] + wind[j]); 
        if (j == 0) r[i] = 0.0;
        r[i] += (pos[j][i] * pos[j][i]);
        if (j == Ndim-1) r[i] = sqrt(r[i]);
        
      }
    }

    for (l = 0; l < Ndim; l++)
    {
      for (i = 0; i < Nbody; i++)
      {

        f[l][i] = f[l][i] -
                  G * mass[i] * M_central * pos[l][i] / pow(r[i],3.0);
      }
    }

    /* calculate pairwise separation of particles */
    for (l = 0; l < Ndim; l++)
    {
      k = 0;
      for (i = 0; i < Nbody; i++)
      {
        for (j = i + 1; j < Nbody; j++)
        {
          delta_pos[l][k] = pos[l][i] - pos[l][j];
          k = k + 1;
        }
      }
    }

    for (i = 0; i < Ndim; i++)
    {
      for ( k = 0; k < Npair; k++)
      {
        if (i == 0)
          delta_r[k] = 0.0;
        delta_r[k] += (delta_pos[i][k] * delta_pos[i][k]);
        if (i == Ndim - 1)
          delta_r[k] = sqrt(delta_r[k]);
      }
    }

    /*
 * add pairwise forces.
 */
    int collis[Nbody][Nbody] = {0};
    for (l = 0; l < Ndim; l++)
    {
      k = 0;
      for (i = 0; i < Nbody; i++)
      {
        for (j = i + 1; j < Nbody; j++)
        {
          Size = radius[i] + radius[j];
          collided = 0;

          /*  flip force if close in */
          if (delta_r[k] >= Size)
          {
            f[l][i] = f[l][i] -
                      G * mass[i] * mass[j] * delta_pos[l][k] / pow(delta_r[k], 3.0);
            f[l][j] = f[l][j] +
                      G * mass[i] * mass[j] * delta_pos[l][k] / pow(delta_r[k], 3.0);
          }
          else
          {
            f[l][i] = f[l][i] +
                      G * mass[i] * mass[j] * delta_pos[l][k] / pow(delta_r[k], 3.0);
            f[l][j] = f[l][j] -
                      G * mass[i] * mass[j] * delta_pos[l][k] / pow(delta_r[k], 3.0);
            collis[i][j] = 1;
          }
          if ((l == Ndim - 1) && collis[i][j] == 1)
          {
            collisions++;
          }
          k = k + 1;
        }
      }
    }

    /* update positions */
    for (j = 0; j < Ndim; j++)
    {
      for (i = 0; i < Nbody; i++)
      {
        pos[j][i] = pos[j][i] + dt * velo[j][i];
        velo[j][i] = velo[j][i] + dt * (f[j][i] / mass[i]);
      }
    }
  }
}
