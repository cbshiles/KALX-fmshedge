/*************************** ex-stoc.cpp **********************************
* Author:        Agner Fog
* Date created:  2001-11-11
* Last modified: 2010-03-09 by P. Frazier for ORIE 5582
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Example showing how to use non-uniform random variate generator library.   *
* Necessary files are in stocc.zip.                                          *
*                                                                            *
* Instructions:
* Compile for console mode and run.
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <stdio.h>                     // define printf()
#include <time.h>                      // define time()
#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"                     // define random library classes

int main() {
   int seed = (int)time(0);            // random seed
   StochasticLib1 sto(seed);           // make instance of random library
   int i;                              // loop counter
   double r;                           // random number
   int ir;                             // random integer number

   // make random numbers with uniform distribution
   printf("Random numbers with uniform distribution:\n");
   for (i=0; i<16; i++) {
      ir = sto.IRandom(0, 20);
      printf("%8i  ", ir);
   }

   // make random numbers with normal distribution
   printf("\n\nRandom numbers with normal distribution:\n");
   for (i=0; i<16; i++) {
      r = sto.Normal(10, 4);
      printf("%8.5f  ", r);
   }

   // make random numbers with poisson distribution
   printf("\n\nRandom numbers with poisson distribution:\n");
   for (i=0; i<16; i++) {
      ir = sto.Poisson(10);
      printf("%8i  ", ir);
   }

   // make random numbers with binomial distribution
   printf("\n\nRandom numbers with binomial distribution:\n");
   for (i=0; i<16; i++) {
      ir = sto.Binomial(40, 0.25);
      printf("%8i  ", ir);
   }

   // make random numbers with hypergeometric distribution
   printf("\n\nRandom numbers with hypergeometric distribution:\n");
   for (i=0; i<16; i++) {
      ir = sto.Hypergeometric(20, 20, 40);
      printf("%8i  ", ir);
   }

   // Test out normal cdf function
   for (i=0; i<8; i++) {
	   double x = -2 + (double)i/2;
	   printf("NormCdf(%f)=%f\n", x, NormCdf(x));
   }

   // Test out the inverse of the normal cdf.
   // It doesn't work well at y=0 or y=1, but works well in between.
   for (i=0; i<17; i++) {
	   double y = (double)i/16;
	   printf("InvNormCdf(%f)=%f\n", y, InvNormCdf(y));
   }

   EndOfProgram();                     // system-specific exit code
   return 0;
}
