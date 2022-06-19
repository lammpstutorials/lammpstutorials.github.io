/* Fumction prototypes and data types for WHAM
 *
 * $Revision$
 * $Author: alan $
 * $Date: 2003/10/06 13:39:49 $
 */
#include <stdio.h>

// A few global variables
extern double HIST_MIN;
extern double HIST_MAX;
extern double BIN_WIDTH;
extern double TOL;
extern double kT;
extern int    NUM_BINS;
extern int    PERIODIC;
extern double PERIOD;

// Some predefined periodic units
#define DEGREES   360.0
#define RADIANS   6.28318530717959

#define k_B 0.001982923700 // Boltzmann's constant in kcal/mol K
//#define k_B  0.0083144621 // Boltzmann's constant kJ/mol-K
//#define k_B 1.0  // Boltzmann's constant in reduced units


// global (untrimmed) histogram, global to prevent reallocation
extern double *HISTOGRAM;


/*
 * histogram datatype
 *
 * Here's the logic -- instead of storing the whole histogram, we just
 *                     store the range containing the nonzero values.
 *                     So, after we've read the data from a given window
 *                     into a fullsized array, we'll locate the index of
 *                     the first and last nonzero values in the histogram,
 *                     and just keep what's in between.  This will lead to
 *                     a large memory savings, especially if we have a lot
 *                     of windows (and when we change over to 2 or 3
 *                     dimensions)
 */

struct histogram
{
double *data; // the actual histogram array
double *cum; // the cumulative histogram array
int first; // index of the 0th element of data in the global histogram
int last;  // index of the last element of data in the global histogram
int num_points;
int num_mc_samples;
};



// This will be a data structure containing everything we need, so
// the routine which reads the metadata can return it all as 1 pointer
struct hist_group
{
int num_windows; // number of separate biased trajectory windows
double *bias_locations; // array of locations of the bias for each window
double *spring_constants; // array of spring constants for the biases for
                          // each window
double *F; // array of free energy perturbations due to restraint
double *F_old; // array of free energy perturbations due to restraint
               // from previous iteration
double *kT; // array of sampling temperatures, in kcal/mol
double *partition; // configurational integral for each window, needed
                   // for normalization purposes if sampling at different
                   // temperatures
struct histogram *hists; // array of histograms for each window
};


// Function prototypes
// file_read.c
int get_numwindows(FILE *file);
int is_metadata(char *line);
int read_metadata(FILE *file, struct hist_group *hist_group);
int read_data(char *filename, int have_energy);



// histogram.c
double get_histval(struct histogram *hist, int index);
struct histogram *hist_alloc(int first, int last, int num_points,
                             int num_mc_samples);
struct hist_group *hist_group_alloc(void);
struct hist_group *make_hist_group(int num_windows);
void clear_HISTOGRAM(void);
void save_free(struct hist_group *h);
int is_converged(struct hist_group *h);
double average_diff(struct hist_group *h);
int calc_free(double *free, double *prob, double kT);
double calc_bias(struct hist_group *h, int index, double coor);
double calc_coor(int i);

// wham.c
void wham_iteration(struct hist_group* hist_group, double *prob,
                   int have_energy);

// bootstrap.c
int get_rand_bin(double *cum, int num_bins, long *idum);
void mk_new_hist(double *cum, double *dist, int num_bins, int num_points,
                 long *idum);

// Numerical Recipes routines
double ran2(long *idum);
double locate(double xx[], int n, double x, int *j);
