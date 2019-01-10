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

double t_min = 1.024;
double t_max = 10.24;

/* read data from file and store them to data[NBINSX] */
int read_data(const char * fname);

double time_to_us(const double t) {
	/* time (in micro sec) per x value */
	return (t*1.e+06);
}

double lnL(double tau) {
  double sum = 0.;
  double t[NBINSX];
  for (int i = 0; i < NBINSX; i++){
    t[i] = time_to_us(time[i]);
    if ((t[i] > t_min) && (t[i] < t_max)) {
      sum += data[i]*(-log(tau) - t[i]/tau
		      - log(exp(-t_min/tau) - exp(-t_max/tau)));
    }
  }
  return sum;
}

double golden_section(double (*f)(double));

double dfdx(double (*f)(double), double x) {
  double h = 1.0e-8;
  return (f(x+0.5*h)-f(x-0.5*h))/h;
}
  
/* main routine */
int main(int argc, char ** argv){
  FILE *fp, *fp2;
  char fname[1024];
  int i, j, k, l, sum;
  double mean;
  double tau;
  double d_tau;
  double t[NBINSX];

  double lnL_max;
  double tau_p, tau_m;

  double epsilon;

  /* read file name from command line */
  if (argc > 1){
    strncpy(fname, argv[1], 1024);
  }else{
    fprintf(stderr, "no input file!\n");
    exit(1);
  }

  /* read data from file and store data to a array */
  read_data(fname);

  /* caluculate mean value, most probable value and so on */

  sum = 0;
  mean = 0.;
  for (i = 0; i < NBINSX; i++){
    t[i] = time_to_us(time[i]);
    if ((t[i] > t_min) && (t[i] < t_max)) { 
      sum += data[i];
      mean += data[i]*t[i];
    }
  }
  mean /= (double)sum;
  /* print calculation result */
  printf("tmin = %lf, tmax = %lf\n", t_min, t_max);
  printf("total count = %d, mean value = %f\n", sum, mean);
  /* unbinned maximum likelihood method */
  /* likelihood L: */
  /*  L = Pi_i 1/tau e^{-t[i]/tau} (i: loop over events) */
  /*  lnL = Sigma_i -ln(tau) - t[i]/tau; */
  /*  d lnL / d tau = Sigma_i {-1/tau + t[i]/tau^2 } = 0 */
  /*  -> tau = Sigma_i t[i] / Sigma_i 1  */
  /*         = mean                      */
  tau = golden_section(lnL);
  lnL_max = lnL(tau);

  epsilon = 1.0;
  tau_p = tau + 0.1;
  for (j = 0; fabs(epsilon) > 1.0e-12; ++j) {
    epsilon = (lnL(tau_p)-lnL_max+0.5)/dfdx(lnL, tau_p);
    tau_p = fmax(tau_p - epsilon, tau + 0.001);
    if (j > 1000) {
      printf("tau_p error\n");
      return -1;
    }
  }
  
  epsilon = 1.0;
  tau_m = tau - 0.1;
  for (j = 0; fabs(epsilon) > 1.0e-12; ++j) {
    epsilon = (lnL(tau_m)-lnL_max+0.5)/dfdx(lnL, tau_m);
    tau_m = fmin(fmax(tau_m - epsilon, 0.01), tau - 0.001);
    if (j > 1000) {
      printf("tau_m error\n");
      return -1;
    }
  }

  printf("lnL_max = %lf\n", lnL_max);
  printf("tau = %lf + %lf - %lf\n", tau, tau_p - tau, tau - tau_m);
      
  /* calculate the error of tau here */
  /* Delta_tau = 1/sqrt(-{d^2 lnL / (d tau')^2}|_{tau'=tau}) */

  /* scan likelihood around the mean and write to file for gnuplot */
  /* estimate positive and negartive error separately */
  /* 1 sigma error: the value(time) for lnL_max-1 */
  fp = fopen("life-likelihood.dat", "w");
  /* adjust range and step size in below */
  for (j = 0; j < 400; j++){
    const double delta = j*0.0005;
    if (delta < 1.1 * fmax(tau_p - tau, tau - tau_m)) {
      fprintf(fp, "%f %f\n", tau+delta, lnL(tau+delta));
      fprintf(fp, "%f %f\n", tau-delta, lnL(tau-delta));
    }
  }
  fclose(fp);

  fp2 = fopen("tmp2.scr", "w");
  fprintf(fp2,"set terminal postscript enhanced color lw 3.0\n");
  fprintf(fp2,"set output \"likelihood.ps\" \n");
  fprintf(fp2,"set xlabel \"{/Symbol t} ({/Symbol m}s)\" \n");
  fprintf(fp2,"set ylabel \"LOG LIKELIHOOD (ln(L))\" \n");
  fprintf(fp2,"plot \"life-likelihood.dat\", ");
  fprintf(fp2,"%.10lf title \"lnL_{max} - 0.5\"\n", lnL_max - 0.5);
  fclose(fp2);

#ifdef WINDOWS
  Sleep(1000);
#else
  sleep(1);
#endif

  system("gnuplot tmp2.scr");
  return 0;
}

int read_data(const char * fname) {
	FILE * fp;
	int i;
	char buff[1024];
	fp = fopen(fname, "r");
	if (!fp) {
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

double golden_section(double (*f)(double)) {
  double x = 2.0;
  double dx = 0.01;
  double a, b, c;
  double w = (3.0 - sqrt(5)) / 2.0; /* golden ratio */

  /* 初期囲い込み開始 */
  if (f(x-dx) < f(x) && f(x+dx) < f(x)) {
    a = x - dx;
    b = x;
    c = x + dx;
  } else if (f(x-dx) > f(x)) {
    /* 左へ進む */
    c = x;
    b = x - dx;
    dx = 2 * dx;
    while (f(b-dx) > f(b)) {
      c = b;
      b = b - dx;
      dx = 2 * dx;
    }
    a = b - dx;
  } else {
    /* 右へ進む */
    a = x;
    b = x + dx;
    dx = 2 * dx;
    while (f(b+dx) > f(b)) {
      a = b;
      b = b + dx;
      dx = 2 * dx;
    }
    c = b + dx;
  }
  /* 初期囲い込み完了 */

  /* 黄金分割法 */
  while ((c-a) > 1.0e-6 || (f(b)-f(a) > 1.0e-12) || (f(b)-f(c) > 1.0e-12)) {
    if ((c-b) > (b-a)) {
      x = b + w * (c-b);
      if (f(x) > f(b)) {
        a = b;
        b = x;
      } else {
        c = x;
      }
    } else {
      x = b + w * (a-b);
      if (f(x) > f(b)) {
        c = b;
        b = x;
      } else {
        a = x;
      }
    }
  }
  return x;
}
