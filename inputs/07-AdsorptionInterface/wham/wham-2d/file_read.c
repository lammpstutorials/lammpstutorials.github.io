/*
 * WHAM -- Weighted Histogram Analysis Method
 * Utility Code -- read metadata and data files
 *
 *
 * $Revision$
 * $Author: alan $
 * $Date: 2003/12/11 18:16:00 $
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include "wham-2d.h"
#define LINESIZE 1024

int get_numwindows(FILE *file)
{
char *line;
int  num_windows;

num_windows = 0;
line = malloc(sizeof(char)*LINESIZE);
if (!line)
    {
    printf("couldn't allocate space for line\n");
    exit(-1);
    }

// make sure we're at the beginning of the file
rewind(file);
line = fgets(line,LINESIZE,file);
while (line != NULL)
    {
    if (is_metadata(line)) num_windows++;
    line = fgets(line,LINESIZE,file);
    }
free(line);
return num_windows;
}



// check if this is a legitimate line, and not a comment or blank line
int is_metadata(char *line)
{
int i,length; 
int not_space;

// check if it's a comment
if (line[0] == '#') return 0;
// verify it's not a blank line
length = strlen(line);
not_space = 0;
i = 0;
while (!not_space && i < length)
    {
    if (!isspace(line[i])) not_space = 1;
    i++;
    }
return not_space;
}


// given a metadata file point and the hist_group structure,
// read the data files from the metadata file, read the data, build the
// histograms.
// return the number of windows allocated

int read_metadata(FILE *file, struct hist_group *hist_group, 
                  int use_mask, int **mask)
{
char *line;
char filename[LINESIZE];
int vals;
double locx, locy, springx, springy, temp;
double sum;
int num_points=0; 
int current_window=0;
int min_nonzero_x = 0;
int max_nonzero_x = NUM_BINSx -1;
int min_nonzero_y = 0;
int max_nonzero_y = NUM_BINSy -1;
int mc_samples;
double correl_time;
int i,j;
int bin, last_bin;
int found_nonzero;
int have_temp = 0; // flag metadata lines which specify a temperature
int have_notemp = 0; // flag metadata lines which don't specify a temperature
struct histogram *h;

line = (char *) malloc(sizeof(char) * LINESIZE);
if (!line)
    {
    printf("couldn't allocate space for line\n");
    exit(-1);
    }

rewind(file);
line = fgets(line,LINESIZE,file);
while (line != NULL)
    {
    if (is_metadata(line))
        {
        vals = sscanf(line, "%s %lf %lf %lf %lf %lf %lf", filename, &locx, 
                                              &locy, &springx, &springy, 
                                              &correl_time, &temp);
        if (vals >= 5)
            {
            //printf("%s,  %f,   %f\n", filename, loc, spring);
            hist_group->bias_locations[current_window][0] = locx;
            hist_group->bias_locations[current_window][1] = locy;
            hist_group->springX[current_window] = springx;
            hist_group->springY[current_window] = springy;
            if (vals == 5)
                {
                correl_time = 1.0;
                }
            if (vals > 6)
                {
                hist_group->kT[current_window] = k_B * temp;
                have_temp = 1;
                }
            else
                {
                hist_group->kT[current_window] = -1.0;
                have_notemp = 1;
                }
            // Either all metadata lines must specify a temperature, or
            // none may have one specified
            if ( (have_temp && have_notemp) || (!have_temp && !have_notemp) )
                {
                printf(
                    "Some but not all metadata lines specify a temperature\n");
                exit(-1);
                }
            // Read data into the global histogram
            num_points = read_data(filename, have_temp, use_mask, mask);
            // make sure we didn't screw up trying to read
            if (num_points < 0)
                {
                printf("Error trying to read %s: %s\n", filename,
                                                    strerror(errno));
                exit(errno);
                }

            mc_samples = (int)(num_points/correl_time);
            if (mc_samples < 1)
                {
                printf("# Correl time is too big for %s:\n", filename);
                printf("# You have %d points, correl time %f\n",
                                            num_points, correl_time);
                printf("# Bootstrap error analysis will crash\n");
                }
            
#if 0
            if (mc_samples < 0)
                {
                printf("# Number of points per mc sample unset\n");
                printf("# for %s, using total number of point (%d)\n",
                        filename, num_points);
                }
            else if (num_points < mc_samples)
                {
                printf("# Number of points %s (%d) is less than specified\n",
                        filename, num_points);
                printf("# number of points per MC sample (%d)\n", mc_samples);
                printf("# Resetting the MC sample size to %d\n", num_points);
                mc_samples = num_points;
                }
#endif

            found_nonzero = 0;
            for (i=0; i<NUM_BINSx;i++)
                {
                for(j=0; j<NUM_BINSy;j++)
                    {
                    if (HISTOGRAM[i][j] > 0)
                        {
                        min_nonzero_x = i;
                        found_nonzero = 1;
                        break;
                        }
                    }
                if (found_nonzero) break;
                }

            found_nonzero = 0;
            for (i=NUM_BINSx-1; i>=0;i--)
                {
                for(j=0; j<NUM_BINSy;j++)
                    {
                    if (HISTOGRAM[i][j] > 0)
                        {
                        max_nonzero_x = i;
                        found_nonzero = 1;
                        break;
                        }
                    }
                if (found_nonzero) break;
                }

            min_nonzero_y=NUM_BINSy-1;
            max_nonzero_y=0;
            for(i=min_nonzero_x;i<=max_nonzero_x;i++)
                {
                for (j=0;j<NUM_BINSy;j++)
                    {
                    if (HISTOGRAM[i][j] > 0)
                        {
                        if (j <= min_nonzero_y)
                            {
                            min_nonzero_y = j;
                            }
                        if (j >= max_nonzero_y)
                            {
                            max_nonzero_y = j;
                            }
                        }
                    }
                }

            if ((min_nonzero_x > max_nonzero_x) || 
                (min_nonzero_y > max_nonzero_y) )
                {
                // Change the error message to match the one I just
                // added to 1d wham
                printf("Error reading time series file %s\n", filename);
                printf("No data points within the histogram bounds:\n");
                printf("X\t[%f, %f]\n", HIST_MINx, HIST_MAXx);
                printf("Y\t[%f, %f]\n", HIST_MINy, HIST_MAXy);
                printf("You need to either change the bounds of your ");
                printf("histogram or remove this file from the metadata");
                printf(" file\n");
                exit(-2);

                }
            // Unnecessary, since we now exit if either of these conditions are
            // met.  This should save me 5-10 emails a year.  Ha...
            //assert(min_nonzero_x <= max_nonzero_x);
            //assert(min_nonzero_y <= max_nonzero_y);

            // allocate a trimmed histogram
            h = hist_alloc(min_nonzero_x, max_nonzero_x,
                           min_nonzero_y, max_nonzero_y,
                           num_points, mc_samples);
            hist_group->hists[current_window] = *h;


            // copy the data into the histogram and compute cum
            bin = 0;
            hist_group->hists[current_window].cum[0] = 0.0;
            for (i=min_nonzero_x; i<=max_nonzero_x; i++)
                {
                for (j=min_nonzero_y; j<=max_nonzero_y; j++)
                    {
                    bin++;
                    hist_group->hists[current_window].cum[bin] = 
                            hist_group->hists[current_window].cum[bin-1] +
                            HISTOGRAM[i][j];
   hist_group->hists[current_window].data[i-min_nonzero_x][j-min_nonzero_y] =
                                 HISTOGRAM[i][j];
                    }
                }

            sum = hist_group->hists[current_window].cum[bin];

            last_bin = bin;
            // normalize the cumulative distribution
            for (bin=0; bin<=last_bin; bin++)
                {
                hist_group->hists[current_window].cum[bin] /= sum;
                }

            // compute the normalization (partition function) if we used 
            // energies
            // Compute even if we didn't use energies -- useful for MC
            hist_group->partition[current_window] = sum;
            current_window++;
            }
        else
            {
            printf("Error parsing metafile: %s\n", filename);
            }
        }
    line = fgets(line,LINESIZE,file);
    }
free(line);
return(current_window);
}

/*
// Read a mask file and build the array of masking rectangles
int read_maskfile(FILE *file, struct mask *mask_array)
{
char *line;
int numvals;
int current_maskline = 0;
double xmin, xmax, ymin, ymax;


line = (char *) malloc(sizeof(char) * LINESIZE);
if (!line)
    {
    printf("couldn't allocate space for line\n");
    exit(-1);
    } 

rewind(file);
line = fgets(line,LINESIZE,file);
while (line != NULL)
    {
    if (is_metadata(line))
        {
        numvals = sscanf(line, "%lf %lf %lf %lf", &xmin, &xmax, &ymin, &ymax);
        if (numvals < 4)
            {
            printf("Error processing mask file line:\n");
            printf("%s", line);
            printf("Only found %d values, need 4\n", numvals);
            exit(-1);
            }
        mask_array[current_maskline].xmin = xmin;
        mask_array[current_maskline].xmax = xmax;
        mask_array[current_maskline].ymin = ymin;
        mask_array[current_maskline].ymax = ymax;
        current_maskline++;
        //printf("%d\n", current_maskline);
        }
    line = fgets(line,LINESIZE,file);
    }
//printf("%d\n", current_maskline);
free(line);
return(current_maskline);
}
*/

// Read a datafile, dump it into the global histogram
// If the second argument, have_energy, is nonzero, this tells us we should 
// be looking for fourth column in the data file, containing the energy of the
// molecule for each timepoint, in addition to its location.
// Return negative value on error, number of data points included
// in the histogram otherwise
//
// clears HISTOGRAM as a side effect
int read_data(char *filename, int have_energy, int use_mask, int **mask)
{
FILE *file;
char *line;
double time, value_x, value_y, energy;
int vals;
int index_x, index_y;
int num_points;

clear_HISTOGRAM();
num_points = 0;

line = (char *) malloc(sizeof(char) * LINESIZE);
if (!line)
    {
    printf("couldn't allocate space for line\n");
    exit(-1);
    }

file = fopen(filename, "r");
if (!file) 
    {
    free(line);
    return -1;
    }

line = fgets(line,LINESIZE,file);
while (line != NULL)
    {
    if (line[0] != '#')
        {
        if (have_energy)
            {
            vals = sscanf(line,"%lf %lf %lf %lf", &time, &value_x, 
                                                  &value_y, &energy);
            if (vals != 4)
                {
                printf("Failure reding %s: missing energy value\n",
                        filename);
                exit(-1);
                }
            }
        else
            {
            vals = sscanf(line,"%lf %lf %lf", &time, &value_x, &value_y);
            }
        if ( (value_x < HIST_MAXx) && (value_x >HIST_MINx) 
        && (value_y < HIST_MAXy) && (value_y >HIST_MINy) ) 
            {
            index_x = (int) ((value_x - HIST_MINx) / BIN_WIDTHx);
            index_y = (int) ((value_y - HIST_MINy) / BIN_WIDTHy);
            if (have_energy)
                {
                HISTOGRAM[index_x][index_y] += exp(-energy/kT);
                }
            else
                {
                HISTOGRAM[index_x][index_y] += 1.0e0;
                }
            num_points++;
            if (use_mask)
                {
                mask[index_x][index_y] = 1;
                }
            }
        }
    line = fgets(line,LINESIZE,file);
    }

fclose(file);
free(line);
return num_points;
}























