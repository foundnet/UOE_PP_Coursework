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

    /* set the viscosity term in the force calculation */
    for (j = 0; j < Ndim; j++)
    {
          for(int i=0;i<Nbody; i++)
          {
            f[j][i] = -vis[i] * velo[j][i];
            f[j][i] = f[j][i] -vis[i] * wind[j];
            r[i] = 0.0;
            r[i] += (pos[j][i] * pos[j][i]);
            if (j == Ndim-1) r[i] = sqrt(r[i]);
            f[j][i] = f[j][i] - force(G * mass[i] * M_central, pos[j][i], r[i]); 
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

    /* calculate norm of seperation vector */
    for (k = 0; k < Npair; k++)
    {
      delta_r[k] = 0.0;
    }
    for (i = 0; i < Ndim; i++)
    {
      add_norm(Npair, delta_r, delta_pos[i]);
    }
    for (k = 0; k < Npair; k++)
    {
      delta_r[k] = sqrt(delta_r[k]);
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
          /*  flip force if close in */
          if (delta_r[k] >= Size)
          {
            f[l][i] = f[l][i] -
                      force(G * mass[i] * mass[j], delta_pos[l][k], delta_r[k]);
            f[l][j] = f[l][j] +
                      force(G * mass[i] * mass[j], delta_pos[l][k], delta_r[k]);
          }
          else
          {
            f[l][i] = f[l][i] +
                      force(G * mass[i] * mass[j], delta_pos[l][k], delta_r[k]);
            f[l][j] = f[l][j] -
                      force(G * mass[i] * mass[j], delta_pos[l][k], delta_r[k]);
            collis[i][j] = 1;
          }

          if ((l == Ndim - 1) && (collis[i][j] = 1))
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