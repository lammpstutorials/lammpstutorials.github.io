/*
 * WHAM -- histogram manipulation routines
 * $Revision$
 * $Author: alan $
 * $Date: 2003/12/16 19:11:20 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "wham-2d.h"

//  Allocate memory for a histogram
struct histogram *hist_alloc(int min_nonzero_x, int max_nonzero_x,
                             int min_nonzero_y, int max_nonzero_y,
                             int num_points, int num_mc_samples)
{
struct histogram *hp;
int size, i;

hp = (struct histogram *) malloc(sizeof(struct histogram));
if (!hp)
    {
    printf("malloc failed allocating histogram\n");
    exit(-1);
    }

hp->first_x = min_nonzero_x;
hp->last_x = max_nonzero_x;
hp->first_y = min_nonzero_y;
hp->last_y = max_nonzero_y;
hp->num_points = num_points;
hp->num_mc_samples = num_mc_samples;

size = (max_nonzero_x - min_nonzero_x + 1);
hp->data = (double **) malloc(sizeof(double*)*size);
if (!hp->data)
    {
    printf("failure allocating histogram data\n");
    exit(-1);
    }
for (i=0; i<size; i++)
    {
    hp->data[i] = (double *)malloc(sizeof(double)*
                                    (max_nonzero_y - min_nonzero_y + 1));
    if (!(hp->data[i]))
        {
        printf("failure allocating histogram data[%d]\n", i);
        exit(-1);
        }
    }

#if 0
// Copy the data into the local histogram
for (i=min_nonzero_x; i<=max_nonzero_x; i++)
    {
    for (j=min_nonzero_y; j<=max_nonzero_y; j++)
        {
        hp->data[i-min_nonzero_x][j-min_nonzero_y] = HISTOGRAM[i][j];
        }
    }
#endif

// allocate space for the cumulative
size *= max_nonzero_y - min_nonzero_y + 1;
hp->cum = (double *)malloc(sizeof(double)*(size+1));
if (!(hp->cum))
    {
    printf("failure allocating cumulative distribution\n");
    exit(-1);
    }

return hp;
}

/* Get a value from a histogram structure
 * Given i and j, the indices into the global histogram, return 
 * the correct histogram value.  Since we only store range of histogram
 * indices containing the nonzero values, we have to check the index value 
 * against the range in the structure.
 */
double get_histval(struct histogram *hist, int i, int j)
{
if ( (i < hist->first_x) || (i > hist->last_x) 
   ||(j < hist->first_y) || (j > hist->last_y))
    {
    return 0.0;
    }
else
    {
    return hist->data[i - hist->first_x][j - hist->first_y];
    }
// can't happen
return 0.0;
}


// Allocate memory for a hist_group structure
struct hist_group *hist_group_alloc(void)
{
struct hist_group *hp;
hp =  (struct hist_group *) malloc(sizeof(struct hist_group));
if (!hp)
    {
    printf("failure allocating memory for hist_group\n");
    exit(-1);
    }
return hp;
}

// Allocate memory for all the pieces of a hist_group structure
struct hist_group *make_hist_group(int num_windows)
{
int i;
struct hist_group *h;

h = hist_group_alloc();
h->num_windows = num_windows;
// allocate space for the bias locations
h->bias_locations = (double **) malloc(sizeof(double *)*num_windows);
for(i=0;i<num_windows;i++)
    {
    h->bias_locations[i] = (double *) malloc(sizeof(double) * 2);
    }

if (!h->bias_locations)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->springX = (double *) malloc(sizeof(double)*num_windows);
if (!h->springX)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->springY = (double *) malloc(sizeof(double)*num_windows);
if (!h->springY)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->F = (double *) malloc(sizeof(double)*num_windows);
if (!h->F)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->F_old = (double *) malloc(sizeof(double)*num_windows);
if (!h->F_old)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->kT = (double *) malloc(sizeof(double)*num_windows);
if (!h->kT)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->partition = (double *) malloc(sizeof(double)*num_windows);
if (!h->partition)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->hists = (struct histogram *) malloc(sizeof(struct histogram)*num_windows);
if (!h->hists)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }
return h;
}


void clear_HISTOGRAM(void)
{
int i,j;
for (i=0; i<NUM_BINSx; i++)
    {
    for (j=0; j<NUM_BINSy; j++)
        {
        HISTOGRAM[i][j] = 0.0e0;
        }
    }
}

// Store the current F values
void save_free(struct hist_group *h)
{
int i;
for (i=0; i<h->num_windows; i++)
    {
    h->F_old[i] = h->F[i];
    h->F[i] = 0.0;
    }
}

// Check if the F's have stopped changing
int is_converged(struct hist_group *h)
{
double error;
int i;

for (i=0; i<h->num_windows; i++)
    {
    error = fabs(h->F[i] - h->F_old[i]);
    if (error > TOL) return 0;
    }
return i;
}

// Calculate the average degree of convergence of the F's
double average_diff(struct hist_group *h)
{
double error;
int i;

error = 0.0;
for (i=0; i<h->num_windows; i++)
    {
    error += fabs(h->F[i] - h->F_old[i]);
    }
error /= (double) h->num_windows;
return error;
}


// Calculate the free energy, setting the minimum value to 0
void calc_free(double **free, double **prob, double kT, 
               int use_mask, int **mask)
{
int i,j;
double offset;
double min = 1e50;

for (i=0; i<NUM_BINSx; i++)
    {
    for (j=0;j<NUM_BINSy; j++)
        {
        if (use_mask && !(mask[i][j]))
            {
            prob[i][j] = 0.0;
            free[i][j] = MASKED;
            }
        else
            {
            free[i][j] = -kT * log(prob[i][j]);
            if (free[i][j] < min)
                {
                min = free[i][j];
                }
            }
        }
    }

offset = min;
for (i=0; i<NUM_BINSx; i++)
    {
    for (j=0;j<NUM_BINSy; j++)
        {
        if (!use_mask || mask[i][j]) // leave masked values where they are
            {
            free[i][j] -= offset;
            }
        }
    }

}

double calc_bias(struct hist_group *h, int index, double *coor)
{
double springx, springy;
double dx,dy;
double half_px=0.0;
double half_py=0.0;

if (PERIODICx)
    {
    half_px = PERIODx/2.0;
    }
if (PERIODICy)
    {
    half_py = PERIODy/2.0;
    }

springx = h->springX[index];
springy = h->springY[index];
dx = coor[0] - h->bias_locations[index][0];
dy = coor[1] - h->bias_locations[index][1];
if (PERIODICx)
    {
    dx = fabs(dx);
    if (dx > half_px)
        {
        dx -= PERIODx;
        }
    }
if (PERIODICy)
    {
    dy = fabs(dy);
    if (dy > half_py)
        {
        dy -= PERIODy;
        }
    }
return 0.5*(springx*(dx*dx) + springy*(dy*dy));
}

void calc_coor(int i, int j, double *coor)
{
coor[0] = HIST_MINx + BIN_WIDTHx*((double)i+0.5); 
coor[1] = HIST_MINy + BIN_WIDTHy*((double)j+0.5); 
}


