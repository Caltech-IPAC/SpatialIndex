#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

#include <tinyhtm/htm.h>


int debug = 0;


/******************************************************************/
/*                                                                */
/* HTMCELL                                                        */
/*                                                                */
/* The purpose of this program is to compute the x,y,z HTM cell   */
/* values for a given (ra,dec).                                   */
/*                                                                */
/******************************************************************/



int main(int argc, char **argv)
{
   int     htmLevel;
   double  dtr, ra, dec, cos_dec;

   struct  htm_v3 vec;
   int64_t htmId;

   dtr = atan(1.)/45.;


   // Process command-line arguments

   if(argc < 4)
   {
      printf("[struct stat=\"ERROR\", msg=\"Usage: htmcell level ra dec\"]\n");
      fflush(stdout);
      exit(1);
   }

   htmLevel = atoi(argv[1]);

   ra  = atof(argv[2]);
   dec = atof(argv[3]);


   // Turn the coordinates into a 3-vector

   cos_dec = cos(dec*dtr);

   vec.x = cos(ra*dtr ) * cos_dec;
   vec.y = sin(ra*dtr ) * cos_dec;
   vec.z = sin(dec*dtr);

   htm_v3_normalize(&vec, &vec);

   if(debug)
   {
      printf("DEBUG> (RA,Dec) = (%.6f,%.6f) ->",      ra, dec);
      printf(" (x,y,z) = (%.8f,%.8f,%.8f) ->", vec.x, vec.y, vec.z);
      fflush(stdout);
   }


   // Find the HTM cell ID

   htmId = htm_v3_id(&vec, htmLevel);

   if(debug)
   {
      printf(" HTM cell: %" PRId64 "  (%" PRId64 ")\n", htmId, htm_idtodec(htmId));
      printf("  -----\n");
      fflush(stdout);
   }


   // Return values

   printf("[struct stat=\"OK\", ra=%.8f, dec=%.8f, x=%.17f, y=%.17f, z=%.17f, level=%d, htm=%" PRId64 "]\n", 
      ra, dec, vec.x, vec.y, vec.z, htmLevel, htmId);
   fflush(stdout);
   exit(0);
}
