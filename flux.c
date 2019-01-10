#include <stdio.h>
#include <math.h>

int main(){
  /* modify the following three constants */
  const double x_max = 35.0;
  const double y_max = 35.0;
  const double step = 1.;
  /* number of events (MODIFY THIS) */
  const long int n_events = 4787146;
  /* distance between the two scintillator */
  const double d = 8.0;
  /* total time */
  const double time = 498721.8;
  int isconst = 1; /*if flux is const*/
  double x1, y1, x2, y2;
  /* integral of  dOmega dS' */
  for (isconst = 0; isconst < 2; isconst++) {
  double integ = 0.;
  for (x1 = 0.; x1 < x_max; x1+= step){
    for (y1 = 0.; y1 < y_max; y1+= step){
      for (x2 = 0.; x2 < x_max; x2+= step){
	for (y2 = 0.; y2 < y_max; y2+= step){
	  /* write function here to calculate integral */
	  
	  integ += pow(step, 4.)* pow(d, 4 - 2*isconst) / pow((pow((x1-x2),2.)+pow((y1-y2),2.)+d*d),3-isconst);
	}
      }
    }
  }
  if (isconst) printf("\nwhen flux is const.\n");
  else printf("\nwhen flux depends on cos^2 theta.\n");
  printf ("total time = %f sec \n", time);
  printf ("number of events = %d\n", n_events);
  printf ("integral dOmega dS' = %f sr*cm^2 \n", integ);
  printf ("flux = %f particles/sr/cm^2/s \n", (double)n_events/time/integ);
  }
}
