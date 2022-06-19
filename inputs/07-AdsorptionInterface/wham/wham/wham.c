/* 
 * WHAM -- Weighted Histogram Analysis Method
 *
 * Reference 1: Computer Physics Communications, 91(1995) 275-282, Benoit Roux
 * Body of code references equation numbers from this paper.
 *
 * $Revision$
 * $Author: agrossf $
 * $Date: 2005/06/24 15:05:03 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "wham.h"


#define COMMAND_LINE "Command line: wham [P|Ppi|Pval] hist_min hist_max num_bins tol temperature numpad metadatafile freefile [num_MC_trials randSeed]\n"

double HIST_MAX,HIST_MIN,BIN_WIDTH,TOL;
double *HISTOGRAM;
double kT;
int  NUM_BINS;
int PERIODIC;
double PERIOD;

int main(int argc, char *argv[])
{
/* moved to global
double kT; // temperature
*/
int i,j;
int len;
int first;
int bin_min;
int have_energy;
char *freefile;
FILE *METAFILE, *FREEFILE; 
struct hist_group *hist_group;
struct histogram  *hp;
double coor;
//double num, denom;
//double bias, bf;
double error;
double *free_ene;

double *prob,*final_prob;
double *ave_p;
double *ave_p2;
double *ave_pdf;
double *ave_pdf2;
double *ave_F;
double *ave_F2;
double sum;
double *final_f;
int iteration;
int max_iteration = 100000;
int numpad;
int num_mc_runs;
int num_used;
char *c;
long idum;
double pdf;

if (argc < 2)
    {
    printf( COMMAND_LINE );
    exit(-1);
    }

// Print the command line out into the output file
printf("#");
for (i=0; i<argc; i++)
    {
    printf(" %s", argv[i]);
    }
printf("\n");

if (toupper(argv[1][0]) == 'P')
    {
    PERIODIC = 1;
    len = strlen(argv[1]);
    if (len == 1)
        {
        PERIOD = DEGREES;  // 360
        }
    else
        {
        c= &(argv[1][1]);
        for (i=0; i<len-1;i++)
            {
            c[i] = toupper(c[i]);
            }
        if (strncmp(c,"PI",2) == 0)
            {
            PERIOD = RADIANS;  // 2 pi
            }
        else
            {
            PERIOD = atof(c);
            }
        }
    printf("#Turning on periodicity with period = %f\n", PERIOD);

    // now shift down the other command line arguments
    argc--;
    argv++;
    }
else
    {
    PERIODIC = 0;
    }

// Parse command line arguments
if (argc != 9 && argc !=11)
    {
    printf( COMMAND_LINE );
    exit(-1);
    }
    
HIST_MIN = atof(argv[1]);
HIST_MAX = atof(argv[2]);
NUM_BINS = atoi(argv[3]);
BIN_WIDTH = (HIST_MAX - HIST_MIN) / (double) NUM_BINS;
TOL = atof(argv[4]);
kT = atof(argv[5]) * k_B;

numpad = atoi(argv[6]);

METAFILE = fopen(argv[7], "r");
if (METAFILE == (FILE *)NULL)
    {
    printf("couldn't open metadatafile %s: %s\n", argv[7], strerror(errno));
    exit(errno);
    }

i = strlen(argv[8]);
freefile = (char *) malloc(i * sizeof(char));
freefile = argv[8];
if (!freefile)
    {
    printf("couldn't allocate space for freefile name: %s\n", strerror(errno));
    exit(errno);
    }

if (argc == 11)
    {
    num_mc_runs = atoi(argv[9]);
    idum = atol(argv[10]);
    if (idum > 0)
        {
        idum = -idum;
        }
    // initialize the random number generator
    ran2(&idum);
    }
else
    {
    num_mc_runs = 0;
    }

HISTOGRAM = (double *) malloc(sizeof(double) * NUM_BINS);
if (!HISTOGRAM)
    {
    printf("couldn't allocate space for HISTOGRAM: %s\n", strerror(errno));
    exit(errno);
    }

prob = (double *) malloc(sizeof(double) * NUM_BINS);
if (!prob)
    {
    printf("couldn't allocate space for prob: %s\n", strerror(errno));
    exit(errno);
    }

final_prob = (double *) malloc(sizeof(double) * NUM_BINS);
if (!final_prob)
    {
    printf("couldn't allocate space for final_prob: %s\n", strerror(errno));
    exit(errno);
    }

ave_p = (double *) malloc(sizeof(double) * NUM_BINS);
if (!ave_p)
    {
    printf("couldn't allocate space for ave_p: %s\n", strerror(errno));
    exit(errno);
    }

ave_p2 = (double *) malloc(sizeof(double) * NUM_BINS);
if (!ave_p2)
    {
    printf("couldn't allocate space for ave_p2: %s\n", strerror(errno));
    exit(errno);
    }

ave_pdf = (double *) malloc(sizeof(double) * NUM_BINS);
if (!ave_pdf)
    {
    printf("couldn't allocate space for ave_pdf: %s\n", strerror(errno));
    exit(errno);
    }

ave_pdf2 = (double *) malloc(sizeof(double) * NUM_BINS);
if (!ave_pdf2)
    {
    printf("couldn't allocate space for ave_pdf2: %s\n", strerror(errno));
    exit(errno);
    }
 
free_ene = (double *) malloc(sizeof(double) * NUM_BINS);
if (!free_ene)
    {
    printf("couldn't allocate space for free_ene: %s\n", strerror(errno));
    exit(errno);
    }  

i = get_numwindows(METAFILE);
printf("#Number of windows = %d\n", i);

ave_F = (double *) malloc(sizeof(double) * i);
if (!ave_pdf)
    {
    printf("couldn't allocate space for ave_F: %s\n", strerror(errno));
    exit(errno);
    }

ave_F2 = (double *) malloc(sizeof(double) * i);
if (!ave_F2)
    {
    printf("couldn't allocate space for ave_F2: %s\n", strerror(errno));
    exit(errno);
    }


hist_group = make_hist_group(i);
//printf("From hist_group: %d\n", hist_group->num_windows);

i = read_metadata(METAFILE, hist_group);
assert(i == hist_group->num_windows);

// Figure out if we have trajectories at different temperatures.
// Missing temperatures are set to -1 in read_metadata, and
// since we require that either all trajectories specify a temperature
// or all trajectories are assumed to be at the wham temperature, we only 
// have to check one of them
if (hist_group->kT[0] > 0)
    {
    have_energy = 1;
    }
else
    {
    have_energy = 0;
    for (i=0; i< hist_group->num_windows; i++) 
        {
        hist_group->kT[i] = kT;
        }
    }

// allocate memory to store the final F values (for when we do MC bootstrap)
final_f = (double *) malloc(sizeof(double) * hist_group->num_windows);
if (!final_f)
    {
    printf("couldn't allocate space for final_f: %s\n", strerror(errno));
    exit(errno);
    } 

free(HISTOGRAM);


// for each window, zero out the estimated perturbation due to the restraints
for (i=0; i< hist_group->num_windows; i++)
    {
    hist_group->F[i]=0.0; // nonzero ensures not converged first time
    hist_group->F_old[i]=0.0;
    }

// Do the actual WHAM stuff, iterate to self consistency
iteration = 0;
first = 1;
while (! is_converged(hist_group) || first )
    {
    first = 0;
    save_free(hist_group);
    wham_iteration(hist_group, prob,have_energy);

    // Dump out some info
    iteration++;
    if (iteration % 10 == 0)
        {
        error = average_diff(hist_group);
        printf("#Iteration %d:  %f\n", iteration, error);
        }

    // Dump out the histogram and free energy
    if (iteration % 100 == 0)
        {
        calc_free(free_ene,prob,kT);
        for (i=0; i< NUM_BINS; i++)
            {
            coor = calc_coor(i);
            printf("%f\t%f\t%f\n", coor, free_ene[i], prob[i]);
            }
        printf("\n");

        // Write the bias values to stdout
        printf("# Dumping simulation biases, in the metadata file order \n");
        printf("# Window  F (free energy units)\n");
        for (j=0; j<hist_group->num_windows;j++)
            {
            printf("# %d\t%f\n", j, hist_group->F[j]);
            final_f[j] = hist_group->F[j];
            }
        }
    // Cheesy bailout if we're going on too long
    if (iteration >= max_iteration) 
        {
        printf("Too many iterations: %d\n", iteration);
        break;
        }
    }


// We're done, write out the free energy, histogram, and bias values

// Write the bias values to stdout
printf("# Dumping simulation biases, in the metadata file order \n");
printf("# Window  F (free energy units)\n");
for (j=0; j<hist_group->num_windows;j++)
    {
    printf("# %d\t%f\n", j, hist_group->F[j]);
    final_f[j] = hist_group->F[j];
    }


// Normalize the probability, store it
sum = 0.0;
for (i=0; i < NUM_BINS; i++)
    {
    sum += prob[i];
    }
for (i=0; i < NUM_BINS; i++)
    {
    prob[i] /= sum;;
    final_prob[i] = prob[i];
    }

// Compute the free energy from the normalized probability
bin_min = calc_free(free_ene, prob,kT);

// Do the requested number of bootstrap monte carlo error analysis runs.
if (num_mc_runs <= 0)
    {
    printf("# No MC error analysis requested\n");
    }

// initialize averaging arrays
for (i=0; i< NUM_BINS; i++)
    {
    ave_p[i] = 0.0;
    ave_p2[i] = 0.0;
    ave_pdf[i] = 0.0;
    ave_pdf2[i] = 0.0;
    }
for (i=0; i< hist_group->num_windows; i++)
    {
    ave_F[i] = 0.0;
    ave_F2[i] = 0.0;
    }

for (i=0; i< num_mc_runs; i++)
    {
      // pick a set of fake data sets
      for (j=0; j<hist_group->num_windows;j++)
        {
        hp = &hist_group->hists[j];
        //printf("Faking %d: %d  %d\n", i,j,hp->num_points);
        num_used = hp->last - hp->first + 1;
        mk_new_hist(hp->cum, hp->data, num_used, hp->num_mc_samples, &idum);
      
        hist_group->F_old[j] = 0.0;
        hist_group->F[j] = 0.0;
        }
      
      // perform WHAM iterations on the fake data sets
      iteration = 0;
      first = 1;
      while (! is_converged(hist_group) || first )
        {
        first = 0;
        save_free(hist_group);
        wham_iteration(hist_group, prob,have_energy);
        iteration++;
        // Cheesy bailout if we're going on too long
        if (iteration >= max_iteration)
          {
          printf("Too many iterations: %d\n", iteration);
          break;
          }
        }  
      printf("#MC trial %d: %d iterations\n", i, iteration);
      printf("#PMF values\n");
      // accumulate the average and stdev of the resulting probabilities
      sum = 0.0;
      for (j=0; j < NUM_BINS; j++)
        {
        sum += prob[j];
        }
      for (j=0; j < NUM_BINS; j++)
        {
        prob[j] /= sum;
        }
      for (j=0; j < NUM_BINS; j++)
          {
          pdf = -kT*log(prob[j]);
      
          ave_p[j] += prob[j];
          ave_pdf[j] += pdf;
          ave_p2[j] += prob[j] * prob[j];
          ave_pdf2[j] += pdf*pdf;
          }
      for (j=0; j<hist_group->num_windows;j++) 
          {
          ave_F[j] += hist_group->F[j] - hist_group->F[0];
          ave_F2[j] += hist_group->F[j]*hist_group->F[j] ; 
          }
    } 
 for (i=0; i < NUM_BINS; i++)
   {
   ave_p[i] /= (double)num_mc_runs;
   ave_p2[i] /= (double)num_mc_runs;
   ave_p2[i] = sqrt(ave_p2[i] - ave_p[i]*ave_p[i]);
   ave_pdf[i] /= (double)num_mc_runs;
   ave_pdf2[i] /= (double)num_mc_runs;
   ave_pdf2[i] = sqrt(ave_pdf2[i] - ave_pdf[i]*ave_pdf[i]);
   }

 for (i=0; i < hist_group->num_windows; i++)
   {
   ave_F[i] /= (double)num_mc_runs;
   ave_F2[i] /= (double)num_mc_runs;
   ave_F2[i] = sqrt(ave_F2[i] - ave_F[i]*ave_F[i]);
   }


FREEFILE = fopen(freefile, "w");
if (!FREEFILE)
    {
    printf("couldn't open %s: %s\n", freefile, strerror(errno));
    printf("dumping free energy,probability, and window free energies to stdout\n");
    for (i=0; i< NUM_BINS; i++)
        {
        coor = calc_coor(i);
        printf("%f\t%f\t%f\t%f\t%f\n", coor, free_ene[i], ave_pdf2[i], 
                    prob[i], ave_p2[i]);
        }
    for (i=0; i<hist_group->num_windows; i++)
        {
        fprintf(FREEFILE,"%d\t%f\t%f\n", i, final_f[i],ave_F2[i]);  
        }

    exit(errno);
    }
else
    {
    // write out header
    fprintf(FREEFILE, "#Coor\t\tFree\t+/-\t\tProb\t\t+/-\n");
    // write out the leading padded values
    for (i=-numpad; i<0; i++)
        {
        coor = calc_coor(i);
        fprintf(FREEFILE,"%f\t%f\t%f\t%f\t%f\n", coor, free_ene[NUM_BINS+i], 
        ave_pdf2[NUM_BINS+i], 
        final_prob[NUM_BINS+i],
        ave_p2[NUM_BINS+i]);
        }
    // write out the center values
    for (i=0; i<NUM_BINS; i++)
        {
        coor = calc_coor(i);
        fprintf(FREEFILE,"%f\t%f\t%f\t%f\t%f\n", coor, free_ene[i],
        ave_pdf2[i],final_prob[i], 
        ave_p2[i]);  
        }

    // write out the trailing padded values
    for (i=0; i<numpad; i++)
        {
        coor = calc_coor(NUM_BINS+i);
        fprintf(FREEFILE,"%f\t%f\t%f\t%f\t%f\n", coor, free_ene[i], 
        ave_pdf2[i],final_prob[i], 
        ave_p2[i]); 
        }

    fprintf(FREEFILE, "#Window\t\tFree\t+/-\t\n");
    for (i=0; i<hist_group->num_windows; i++)
        {
        fprintf(FREEFILE,"#%d\t%f\t%f\n", i, final_f[i],ave_F2[i]);  
        }
    }


exit(0);
}


/*******************************************************************/

/*
 *  Perform a single WHAM iteration
 */
void wham_iteration(struct hist_group* hist_group, double *prob, 
                   int have_energy)
{
int i,j;
double num, denom, bias, bf, coor;
// loop over bins of global histogram
for (i=0; i<NUM_BINS; i++)
    {
    coor = calc_coor(i);
    num = 0.0;
    denom = 0.0;
    /*
     *   use previous biases to estimate probability
     *   Equation 8 in Reference 1
     */
    for (j=0; j<hist_group->num_windows;j++)
        {
        num += (double) get_histval( &(hist_group->hists[j]),i);
        bias = calc_bias(hist_group,j,coor);
        bf = exp((hist_group->F_old[j] - bias) / hist_group->kT[j]);
        /* If we have energies, then the histograms contain boltzmann 
         * factors.  If not, they contain integer counts.  Accordingly, 
         * we either use a partition function to normalize, or simply the 
         * number of points.
         */
        if (have_energy)
            {
            denom += (double) hist_group->partition[j] * bf;
            }
        else
            {
            denom += (double) hist_group->hists[j].num_points * bf;
            }
        }
    prob[i] = num / denom;
    /*
     *   use new probability to update bias estimate
     *   Equation 9 from Reference 1
     */
    for (j=0; j<hist_group->num_windows;j++)
        {
        bias = calc_bias(hist_group,j,coor);
        bf = exp(-bias/hist_group->kT[j]) * prob[i];
        hist_group->F[j] += bf;
        }
    }
/*
 *   take natural log of Equation 9 from Reference 1
 *   because we store F rather than exp(-F/kT)
 */
for (j=0; j<hist_group->num_windows;j++)
    {
    hist_group->F[j] = -hist_group->kT[j] * log(hist_group->F[j]);
    }
for (j=hist_group->num_windows - 1;j>=0;j--)
    {
    hist_group->F[j] -= hist_group->F[0];
    }
}
