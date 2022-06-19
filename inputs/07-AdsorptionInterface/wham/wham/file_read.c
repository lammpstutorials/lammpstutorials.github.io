/*
 * WHAM -- Weighted Histogram Analysis Method
 * Utility Code -- read metadata and data files
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include "wham.h"
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

int read_metadata(FILE *file, struct hist_group *hist_group)
{
char *line;
char filename[LINESIZE];
int vals, mc_samples;
double loc, spring, temp;
double sum;
double correl_time;
int num_points=0; 
int current_window=0;
int min_nonzero = 0;
int max_nonzero = NUM_BINS -1;
int num_used;
int still_zero;
int i;
int have_temp = 0;  // flag metadata lines which specify a temperature
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
        vals = sscanf(line, "%s %lf %lf %lf %lf", filename, &loc, &spring, 
                                              &correl_time, &temp);
        if (vals >= 3)
            {
            //printf("%s,  %f,   %f\n", filename, loc, spring);
            hist_group->bias_locations[current_window] = loc;
            hist_group->spring_constants[current_window] = spring;

            if (vals == 3)
                {
                correl_time = 1.0;
                }

            if (vals >4)
                {
                hist_group->kT[current_window] = k_B * temp;
                have_temp = 1;
                }
            else
                {
                hist_group->kT[current_window] = -1.0;
                have_notemp = 1;
                }
            // Either all lines must have a temperature specified or
            // none may have one specified
            if ( (have_temp && have_notemp) || (!have_temp && !have_notemp) )
                {
                printf(
               "Some but not all metadata lines specify a temperature\n");
                exit(-1);
                }
            // Read data into the global histogram
            num_points = read_data(filename, have_temp);
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

            // find first non_zero value
            still_zero = 1;
            i = 0;
            while (still_zero && (i < NUM_BINS))
                {
                if (HISTOGRAM[i] > 0)
                    {
                    still_zero = 0;
                    }
                else
                    {
                    i++;
                    }
                }
            min_nonzero = i;

            // find last non_zero value
            still_zero = 1;
            i = NUM_BINS - 1;
            while (still_zero && (i >= min_nonzero))
                {
                if (HISTOGRAM[i] > 0)
                    {
                    still_zero = 0;
                    }
                else
                    {
                    i--;
                    }
                }
            max_nonzero = i;

            // This assert confuses the hell out people, so I'm replacing
            // it with a readable error message.
            //assert(min_nonzero <= max_nonzero);
            if (min_nonzero > max_nonzero)
                {
                printf("Error reading time series file %s\n",  filename);
                printf("No data points within histogram bounds [%f, %f]\n",
                        HIST_MIN, HIST_MAX);
                printf("You need to either change the bounds of your ");
                printf("histogram or remove this file from the metadata ");
                printf("file.\n");
                exit(-2);
                }

            //printf("Min = %d\tMax = %d\n", min_nonzero, max_nonzero);

            // allocate a trimmed histogram
            h = hist_alloc(min_nonzero, max_nonzero, num_points, 
                           mc_samples);
            hist_group->hists[current_window] = *h;

            // Copy data into histogram
            for (i=min_nonzero; i<=max_nonzero; i++)
                {
                hist_group->hists[current_window].data[i-min_nonzero] =
                                                                  HISTOGRAM[i];
                // store the cumulative distribution
                if (i == min_nonzero)
                    {
                    hist_group->hists[current_window].cum[0] = 0.0;
                    }
                else
                    {
                    hist_group->hists[current_window].cum[i-min_nonzero] =
                       hist_group->hists[current_window].cum[i-min_nonzero-1] 
                       + HISTOGRAM[i-1];
                    }
                }
           
            num_used = max_nonzero - min_nonzero;

            sum = hist_group->hists[current_window].cum[num_used] +
                        HISTOGRAM[max_nonzero];

            for (i=0; i<=num_used; i++)
                {
                hist_group->hists[current_window].cum[i] /= sum;
                /*
                printf("%d\t%f\t%f\n", i, 
                        hist_group->hists[current_window].cum[i],
                        hist_group->hists[current_window].cum[num_used]);
                 */
                }
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

// Read a datafile, dump it into the global histogram
// If the second argument, have_energy, is nonzero, this tells us we should be 
// looking for a third column in the datafile, containing the energy of the
// molecule for each timepoint, in addition to its location.  
// Return negative value on error, number of data points included
// in the histogram otherwise
//
// clears HISTOGRAM as a side effect
int read_data(char *filename, int have_energy)
{
FILE *file;
char *line;
double time, value,energy;
int vals;
int index;
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
    //fclose(file);
    return -1;
    }

line = fgets(line,LINESIZE,file);
while (line != NULL)
    {
    if (line[0] != '#')
        {
        if (have_energy)
            {
            vals = sscanf(line,"%lf %lf %lf", &time, &value, &energy);
            if (vals != 3)
                {
                printf("failure reading %s: missing energy value\n",
                        filename);
                exit(-1);
                }
            }
        else
            {
            vals = sscanf(line,"%lf %lf", &time, &value);
            if (vals != 2)
                {
                printf("failure reading %s: missing position value\n",
                        filename);
                exit(-1);
                }
            }
        if ( (value < HIST_MAX) && (value > HIST_MIN) )
            {
            index = (int) ((value - HIST_MIN) / BIN_WIDTH);
            if (have_energy)
                {
                HISTOGRAM[index] += exp(-energy/kT);
                }
            else
                {
                HISTOGRAM[index]+= 1.0e0;
                }
            num_points++;
            }
        }
    line = fgets(line,LINESIZE,file);
    }

fclose(file);
free(line);
return num_points;
}























