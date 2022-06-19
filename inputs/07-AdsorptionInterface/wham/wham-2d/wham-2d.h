/* Fumction prototypes and data types for WHAM
 *
 * $Revision$
 * $Author: alan $
 * $Date: 2003/12/16 19:11:51 $
 */
#include <stdio.h>

// A few global variables
extern double HIST_MINx;
extern double HIST_MAXx;
extern double BIN_WIDTHx;
extern double HIST_MINy;
extern double HIST_MAXy;
extern double BIN_WIDTHy;
extern double TOL;
extern double kT;
extern int    NUM_BINSx;
extern int    NUM_BINSy;

extern int    PERIODICx,PERIODICy;  // flags to turn on periodicity
extern double PERIODx, PERIODy;     // flags to control periodic interval

// A couple of predefined periodic units
#define DEGREES   360.0
#define RADIANS   6.28318530717959

#define k_B 0.001982923700 // Boltzmann's constant in kcal/mol K
//#define k_B  0.0083144621 // Boltzmann's constant kJ/mol-K
//#define k_B 1.0  // Boltzmann's constant in reduced units



// Value inserted for the free energy of masked values
#define MASKED 9999999.

// global (untrimmed) histogram, global to prevent reallocation
extern double **HISTOGRAM;


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
 *                     dimensions) */

struct histogram
{
double **data; // the actual histogram array
double *cum; // cumulative histogram array
int first_x; // index of the 0th element of data in the global histogram
int last_x;  // index of the last element of data in the global histogram
int first_y;
int last_y;
int num_points;
int num_mc_samples;
int shift; // location of first nonzero element of cum
};

// This will be a data structure containing everything we need, so
// the routine which reads the metadata can return it all as 1 pointer
struct hist_group
{
int num_windows; // number of separate biased trajectory windows
double **bias_locations; // array of locations of the bias for each window
                         // 2d for x and y
double *springX; // array of spring constants for the biases for
                          // each window
double *springY; // array of spring constants for the biases for
                          // each window
double *F; // array of free energy perturbations due to restraint
double *F_old; // array of free energy perturbations due to restraint
               // from previous iteration
double *kT; // array of sampling temperatures, in kcal/mol
double *partition; // configurational integral for each window, needed for
                   // normalization purposes if sampling at different
                   // temperatures
struct histogram *hists; // array of histograms for each window
};

// The mask datatype specifies a rectangular region of the reaction coordinate
// surface which will be ignored when computing the free energy.  If using this
// functionality, we'll store an array of mask rectangles.
struct mask
{
double xmin;
double xmax;
double ymin;
double ymax;
};

// Function prototypes
// file_read.c
int get_numwindows(FILE *file);
int is_metadata(char *line);
int read_metadata(FILE *file, struct hist_group *hist_group,
                  int use_mask, int **mask);
//int read_maskfile(FILE *file, struct mask *mask_array);
int read_data(char *filename, int have_energy, int use_mask, int **mask);
//int build_mask(int num_masks, struct mask *mask_array, int **mask);
void build_mask(int **mask);




// histogram.c
double get_histval(struct histogram *hist, int i, int j);
struct histogram *hist_alloc(int min_nonzero_x, int max_nonzero_x,
                             int min_nonzero_y, int max_nonzero_y,
                             int num_points, int num_mc_samples);
struct hist_group *hist_group_alloc(void);
struct hist_group *make_hist_group(int num_windows);
void clear_HISTOGRAM(void);
void save_free(struct hist_group *h);
int is_converged(struct hist_group *h);
double average_diff(struct hist_group *h);
void calc_free(double **free, double **prob, double kT,
               int use_mask, int **mask);
double calc_bias(struct hist_group *h, int index, double *coor);

void calc_coor(int i,int j, double *coor);


// wham-2d.c

int parse_periodic(char *c, double *period);
void wham_iteration(struct hist_group* hist_group, double **prob,
                            int have_energy, int use_mask, int **mask);

// bootstrap.c
int get_rand_bin(double *cum, int num_bins, long *idum);
void mk_new_hist(double *cum, double *dist, int num_bins, int num_points,
                 long *idum);

// Numerical Recipes routines
double ran2(long *idum);
double locate(double xx[], int n, double x, int *j);
