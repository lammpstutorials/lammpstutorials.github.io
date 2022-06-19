void locate(double xx[], int n, double x, int *j)
{
	int ju,jm,jl;
	int ascnd;

    n -= 1;  // n is the index of the last entry

	jl=0;
	ju=n+1;
	//ascnd=(xx[n] > xx[1]);  // original version
	ascnd=(xx[n] > xx[0]);  // I think this makes it zero based
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
    // Added by AMG: If there are several equal bins, the standard
    //               version puts everything in the highest valued bin.
    //               We need it to always put things in the lowest valued
    //               of the equal bins.  This is critical in 2d wham, where
    //               there are lots of zero-density bins.
    // Linear scan is ugly, but works well enough here because we're always
    // scanning a small number of values
    while ( (xx[jl] <= xx[jl-1]) && (jl > 0) && (jl != n))
        {
        //printf(":::%d\t%f\t%f\n", jl, x, xx[jl]);
        jl--;
        } 
	*j=jl;
}
/* (C) Copr. 1986-92 Numerical Recipes Software "k=($5.&12$o1j. */
