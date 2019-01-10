#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NBINSX (256)  // <- Modify this if necessary

int data[NBINSX];
double adc[NBINSX];

/* read data from file and store them to data[NBINSX] */
int read_data(const char * fname);

double adc_to_energy(const double x){
  /* energy (in MeV) per x value (MODIFY THIS) */
  return x*0.0096 - 0.048;
}
  
/* main routine */
int main(int argc, char ** argv){
  FILE* fp;
  /* input file name */
  char fname[1024];
  int i;
  /* sum of the number of events */
  int sum;
  /* maximum number of events */
  int max;
  /* most probable value */
  double most_probable;
  /* <- add lines below for the calculation of the mean value */
  double mean; 

  /* read file name from command line */
  if (argc > 1){
    strncpy(fname, argv[1], 1024);
  }else{
    fprintf(stderr, "no input file!\n");
    exit(1);
  }

  /* read data from file and store data to a array: data[NBINSX] */
  /* data[NBINSX] contain number of events for corresponding ADC value */
  read_data(fname);

  /* caluculate mean value, peak value and so on */
  /* caluculate mean value here */
  sum = 0;
  max = 0;
  mean = 0.;
  most_probable = -1.;
  for (i = 0; i < NBINSX; i++){
    /* convert: ADC count i -> energy */
    const double energy = adc_to_energy(adc[i]);
    /* i-> ADC count    data[i] -> number of events for "i" */
    sum += data[i];
    mean += energy * data[i];
    if (data[i] > max){
      max = data[i];
      most_probable = energy;
    }
  }

  mean /= sum;
  
  /* print calculation result */
  printf("number of events = %d\n", sum);
  printf("most probable energy loss = %e MeV, mean = %e MeV\n", 
	 most_probable, mean);
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
			  if (2 == sscanf(buff, "%lf %d\n", adc + i, data + i)) {
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
