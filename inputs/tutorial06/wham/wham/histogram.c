/*
 * WHAM -- histogram manipulation routines
 * $Revision$
 * $Author: alan $
 * $Date: 2003/10/06 13:38:42 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "wham.h"


//  Allocate memory for a histogram
struct histogram *hist_alloc(int first, int last, int num_points, 
                             int num_mc_samples)
{
struct histogram *hp;

hp = (struct histogram *) malloc(sizeof(struct histogram));
if (!hp)
    {
    printf("malloc failed allocating histogram\n");
    exit(-1);
    }

hp->first = first;
hp->last = last;
hp->num_points = num_points;
hp->num_mc_samples = num_mc_samples;
hp->data = (double *) malloc(sizeof(double)*(last - first + 1) );
hp->cum = (double *) malloc(sizeof(double)*(last - first + 1) );

if (!(hp->data) || !(hp->cum))
    {
    printf("Failure allocating data or cum for histogram\n");
    exit(-1);
    }

return hp;
}

/* Get a value from a histogram structure
 * Given index, the index into the global histogram, return 
 * the correct histogram value.  Since we only store range of histogram
 * indices containing the nonzero values, we have to check the index value 
 * against the range in the structure.
 */
double get_histval(struct histogram *hist, int index)
{
if ( (index < hist->first) || (index > hist->last) )
    {
    return 0.0;
    }
else
    {
    return hist->data[index - hist->first];
    }
// can't happen
return 0.0;
}

double get_cumval(struct histogram *hist, int index)
{
if (index < hist->first)
    {
    return 0.0;
    }
else if (index > hist->last)
    {
    return 1.0;
    }
else
    {
    return hist->cum[index - hist->first];
    }
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
struct hist_group *h;

h = hist_group_alloc();
h->num_windows = num_windows;
h->bias_locations = (double *) malloc(sizeof(double)*num_windows);
if (!h->bias_locations)
    {
    printf("allocation error in make_hist_group: %s\n", strerror(errno));
    exit(errno);
    }

h->spring_constants = (double *) malloc(sizeof(double)*num_windows);
if (!h->spring_constants)
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
int i;
for (i=0; i<NUM_BINS; i++)
    {
    HISTOGRAM[i] = 0.0e0;
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


// Calculate the free energy, and returning the bin index at that point.
int calc_free(double *free, double *prob, double kT)
{
int i,bin_min;
double offset;
double min = 1e50;

bin_min = 0;
for (i=0; i<NUM_BINS; i++)
    {
    free[i] = -kT * log(prob[i]);
    if (free[i] < min) 
        {
        min = free[i];
        bin_min = i;
        }
    }
offset = min;
for (i=0; i<NUM_BINS; i++)
    {
    free[i] -= offset;
    }

return bin_min;
}

double calc_bias(struct hist_group *h, int index, double coor)
{
double spring, loc;
double dx;

spring = h->spring_constants[index];
loc = h->bias_locations[index];
dx = coor - loc;
// minimum image if periodic -- assumes angles in degrees
if (PERIODIC)
    {
    dx = fabs(dx);
    if (dx > PERIOD/2.0)
        {
        dx -= PERIOD;
        }
    }
return 0.5*dx*dx*spring;
}


double calc_coor(int i)
{
return HIST_MIN + BIN_WIDTH*((double)i+0.5); 
}
