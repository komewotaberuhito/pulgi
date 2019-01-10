#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//#define WINDOWS (1)  // Comment out this when using on Unix-like system.

#ifdef WINDOWS
#include <Windows.h>
#endif

#define NBINSX (16384)  // <- Modify this if necessary

int data[NBINSX];
double time[NBINSX];

/* read data from file and store them to data[NBINSX] */
int read_data(const char * fname);

double time_to_us(const double t) {
	/* time (in micro sec) per x value */
	return (t*1.e+06);
}

/* main routine */
int main(int argc, char ** argv){
  char fname[1024];
  FILE *fp, *fp2;
  int i, j, k, l;

  int bin = 32;

  /* n: histogram with reduced bins: data[NBINSX]->n[64] */
  double n[bin];
  /* y: natural logarithm of n[64]: y[i]=log(n[i]) */
  double y[bin];
  /* ey: error of y */
  double ey[bin];
  /* time */
  double t[bin];

  double E[2][2];
  /* error matrix */
  double sigma[2][2];
  /* solution of the fit: f(t) = alpha[0] + alpha[1] t */
  double alpha[2];
  double beta[2];
  double det;
  double tmin_bin, tmax_bin, dt_bin;
  int nbins_merge;

  int merge;

  int bin_min, bin_max;
  int threshold = 10;
  double chi2;

  /* read file name from command line */
  if (argc > 1){
    strncpy(fname, argv[1], 1024);
  }else{
    fprintf(stderr, "no input file!\n");
    exit(1);
  }

  /* read data from file and store data to a array */
  read_data(fname);

  merge = (NBINSX+bin-1) / bin;
  if (NBINSX >= bin) {
	  nbins_merge = bin;
  }
  else {
	  nbins_merge = NBINSX;
  }

  printf("nbins_merge = %d\n", nbins_merge);
  
  /* merge merge=(NBINSX+63)/64 bins to one bin to inclease statistics */
  /* number of events for merged bins: n[j] */
  /* \Sigma_{i=j*merge}^{(j+1)*merge} data[i] = n[j] */
  for (j = 0; j < bin; j++) n[j] = 0;
  for (i = 0; i < NBINSX; i++) n[i/merge] += data[i];
  /* centeral value of the bin */
  /* time for merged bins: t[j] */
  if (nbins_merge > 1) {
	  dt_bin = time[1] - time[0];
	  for (j = 0; j < nbins_merge; j++) {
	    tmin_bin = time_to_us(time[merge*j] - 0.5*dt_bin); //第二項中にいれた
	    tmax_bin = time_to_us(time[merge*(j + 1) - 1] + 0.5*dt_bin);
	    t[j] = 0.5*(tmin_bin + tmax_bin); //binの中間
	  }
  }
  else {
	  t[0] = time_to_us(time[0]);
  }

  /*同時の信号を除くため一つ目を除去*/
  bin_min = 1;

  printf("n[j] > %d\n", threshold);
  /*イベント数がthreshold以下のものがあればそれ以降は取らない*/
  for (j=0; j < nbins_merge; j++) {
    if (n[j] <= threshold) {
      bin_max = j - 1;
      break;
    }
  }

  printf("number of merged data = %d\n", bin_max - bin_min + 1);
  printf("tmin = %lf, tmax = %lf\n",
	 time_to_us(time[merge*bin_min] - 0.5*dt_bin),
	 time_to_us(time[merge*(bin_max + 1) -1] + 0.5*dt_bin));
  
  /* log(n[j]) -> y[j] +- ey[j] */
  fp=fopen("life-merged-bins.dat", "w");
  for (j = bin_min; j <= bin_max; j++){
    if (n[j]>0){ /* number of events > 0 */
      /* calculate logarithm of the number of events */
      y[j] = log(n[j]);
      /* caluculate its error */
      /* ey[j] = ... */
      ey[j] = 1./sqrt(n[j]);
      /* write values to the file (life-merged-bins.txt) for gnuplot */
      fprintf(fp,"%f %f %f\n", t[j], y[j], ey[j]);
    }else{ /* number of events = 0 */
      y[j] = -1;
      ey[j] = -1;
    }
  }
  fclose(fp);

  /* linear minimum chi square fit                                   */
  /* f(t) = \Sigma_k alpha[k] t^k (k=0,1 in this fit)                */
  /* ( f(t) = alpha[0] + alpha[1] t )                                */
  /* d chi^2/d a[k] = \Sigma_i (f(t[j])-y[j])/(ey[j])^2 * t[j]^k = 0 */
  /* expression in matrix: beta = E alpha                            */
  /* where                                                           */
  /* alpha[k]-> solution of the fit                                  */
  /* beta[k] = \Sigma_j t[j]^k * y[j] / ey[j]^2                      */
  /* E[k][l] = \Sigma_j t[j]^(k+l) / ey[j]^2                         */

  /* initialize beta[k] and E[m][n] */
  for (k = 0; k < 2; k++){
    beta[k] = 0.;
    for (l = 0; l < 2; l++) E[k][l] = 0.;
  }

  /* calculate beta[k] and E[m][n] */
  for (j = bin_min; j <= bin_max; j++){
    if (ey[j] < 0) continue; /* ignore empty bin */
    const double w=1./ey[j]/ey[j];
    beta[0] +=      y[j]*w;
    beta[1] += t[j]*y[j]*w;
    E[0][0] +=           w;
    E[1][0] +=      t[j]*w;
    E[0][1] +=      t[j]*w;
    E[1][1] += t[j]*t[j]*w;
  }

  /* invert covariant matrix E */
  /* alpha = E^{-1} beta = sigma beta */
  /* error matrix: sigma = E^{-1} */
  det = E[0][0]*E[1][1]-E[0][1]*E[1][0];
  sigma[0][0] =  1./det*E[1][1];
  sigma[0][1] = -1./det*E[0][1];
  sigma[1][0] = -1./det*E[1][0];
  sigma[1][1] =  1./det*E[0][0];

  /* calculate alpha: alpha = sigma beta */
  /* error of alpha: diagonal part of the error matrix, sigma */
  alpha[0] = sigma[0][0]*beta[0]+sigma[1][0]*beta[1];
  alpha[1] = sigma[1][0]*beta[0]+sigma[1][1]*beta[1];
  printf("f(t) = alpha[0] + alpha[1] t\n");
  printf("alpha[0] = %f +- %f\n", alpha[0], sqrt(sigma[0][0]));
  printf("alpha[1] = %f +- %f\n", alpha[1], sqrt(sigma[1][1]));
  
  /* print tau and its error here */
  printf("tau = %lf +- %lf\n", -1./alpha[1], sqrt(sigma[1][1])/pow(alpha[1],2));

  /* calculate chi square here */
  chi2 = 0.;
  for(j = bin_min; j <= bin_max; ++j) {
    chi2 += pow((y[j] - (alpha[0]+alpha[1]*t[j]))/ey[j], 2);
  }

  printf("chi2 = %lf, chi2/n = %lf\n",
	 chi2, chi2/(1.0*(bin_max-bin_min+1 - 2)));

  /* for gnuplot */
  fp2 = fopen("tmp.scr", "w");
  fprintf(fp2,"set terminal postscript enhanced color lw 3.0\n");
  fprintf(fp2,"set output \"life.ps\" \n");
  fprintf(fp2,"set xlabel \"TIME ({/Symbol m}s)\" \n");
  fprintf(fp2,"set ylabel \"ln(COUNTS/bin)\" \n");
  fprintf(fp2,"plot \"life-merged-bins.dat\" with yerrorbar, %f*x+%f\n", 
	  alpha[1], alpha[0]);
  fclose(fp2);

#ifdef WINDOWS
  Sleep(1000);
#else
  sleep(1);
#endif

  system("gnuplot tmp.scr");
  return 0;
}

int read_data(const char * fname){
  FILE * fp;
  int i;
  char buff[1024];
  fp=fopen(fname, "r");
  if (!fp){
    fprintf(stderr, "wrong input file: %s!\n", fname);
    exit(1);
  }
  i = 0;
  while (1) {
	  if (i < NBINSX) {
		  if (fgets(buff, sizeof(buff), fp) != NULL) {
			  if (2 == sscanf(buff, "%lf %d\n", time + i, data + i)) {
				  ++i;
			  }
		  }
		  else {
			  fprintf(stderr, "Error. NBINSX (%d) is larger than the input data. Aborted.", NBINSX);
			  exit(1);
		  }
	  }
	  else {
		  if (fgets(buff, sizeof(buff), fp) == NULL) {
			  break;
		  }
		  else {
			  fprintf(stderr, "Error. NBINSX (%d) is smaller than the input data. Aborted.", NBINSX);
			  exit(1);
		  }
	  }
  }
  return i;
}
