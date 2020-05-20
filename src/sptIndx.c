#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

#include <tinyhtm/htm.h>


#define MAXCOL    64
#define MAXSTR 32768

static int     csvInit    (char *fname);
static int     csvReadRow ();
static char   *csvOrigStr ();
static char   *csvColVal  (int col);
static char   *csvStrip   (char *);

static int64_t sky2hpx    (int order, int64_t nside, double lon, double lat);
static int64_t xyf2nest   (int order, int64_t ix, int64_t iy, int64_t face_num);
static int64_t spread_bits(int v);

int debug = 0;


/******************************************************************/
/*                                                                */
/* SPTINDX                                                        */
/*                                                                */
/* The purpose of this program is to pre-process a table          */
/* containing astronomical coordinates (columns "ra" and "dec")   */
/* and add five new columns: x, y, z, htmindx and hpxindx.        */
/* The first three are unit sphere 3-vector coordinates           */
/* corresponding to each (RA,Dec).  The fourth is the HTM index   */
/* number for the cell covering that location and the last is the */
/* HEALPix cell index number the same location.                   */
/*                                                                */
/* The inputs are the HTM and HEALPix levels desired, the input   */
/* CSV table file, and the output CSV table file to be created.   */
/*                                                                */
/* The input is assumed to have one extra "header" line at        */
/* the top containing the column names (also in CSV format).      */
/*                                                                */
/* This program isn't foolproof.  It does not check for           */
/* pre-existing columns with the four names we create and         */
/* doesn't look for information like Equinox (J2000 is assumed).  */
/* The coordinates are assume to be decimal degrees.              */
/*                                                                */
/* While it can be used as-is, it is really meant to be more      */
/* of a template for designing code specific to local file        */
/* formats and databases.                                         */
/*                                                                */
/******************************************************************/



int main(int argc, char **argv)
{
   int     htmLevel, hpxLevel;
   int64_t nside;
   int     col, row;
   int     ncol, ira, idec;
   double  dtr, ra, dec, cos_dec;

   char    infile [MAXSTR];
   char    outfile[MAXSTR];

   FILE   *fout;

   struct  htm_v3 vec;
   int64_t htmId;

   int64_t hpxId;


   dtr = atan(1.)/45.;


   // Process command-line arguments

   if(argc < 4)
   {
      printf("Usage: sptIndx level infile outfile\n");
      fflush(stdout);
      exit(1);
   }

   htmLevel = atoi(argv[1]);

   hpxLevel = htmLevel;

   nside = 1ULL << hpxLevel;

   if(debug)
   {
      printf("hpxLevel = %d, nside = %" PRId64 "\n", hpxLevel, nside);
      fflush(stdout);
   }

   strcpy(infile,  argv[2]);
   strcpy(outfile, argv[3]);


   // Open the input CSV table

   if(csvInit(infile))
   {
      printf("ERROR opening CSV file [%s].  Exiting.\n", infile);
      fflush(stdout);
      exit(1);
   }


   // Open the output CSV file

   if((fout = fopen(outfile, "w+")) == (FILE *)NULL)
   {
      printf("ERROR opening output file [%s].  Exiting.\n", outfile);
      fflush(stdout);
      exit(1);
   }


   // Process the header

   ira  = -1;
   idec = -1;

   ncol = csvReadRow();

   if(debug)
      printf("\nLine: [%s]\n", csvOrigStr());

   for(col=0; col<ncol; ++col)
   {
      if(debug)
         printf("header col %d: [%s]\n", col, csvColVal(col));

      if(strcasecmp(csvColVal(col), "ra") == 0)
      {
         ira = col;

         if(debug)
            printf("DEBUG> RA in column %d\n\n", col);
      }

      if(strcasecmp(csvColVal(col), "dec") == 0)
      {
         idec = col;

         if(debug)
            printf("Dec in column %d\n\n", col);
      }
   }

   if(ira<0 || idec < 0)
   {
      printf("ERROR: Cannot find columns 'ra' and 'dec'.  Exiting.\n");
      fflush(stdout);
      exit(1);
   }


   // Extended header for the output file

   fprintf(fout, "%s,x,y,z,htm%d,hpx%d\n", csvOrigStr(), htmLevel, hpxLevel);
   fflush(fout);


// Loop over the input data records
   
   row = 0;

   while(1)
   {
      ncol = csvReadRow();

      if(debug)
         printf("\nDEBUG> Line: [%s]\n", csvOrigStr());

      if(ncol < 0)
         break;

      ra  = (atof(csvColVal(ira)));
      dec = (atof(csvColVal(idec)));


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


      // Find the HEALPix cell ID
   
      hpxId = sky2hpx(hpxLevel, nside, ra, dec);

      if(debug)
      {
         printf("DEBUG> HPX cell: %" PRId64 "\n", hpxId);
         printf("  -----\n");
         fflush(stdout);
      }



      // Extended output record output file

      fprintf(fout, "%s,%.17f,%.17f,%.17f,%" PRId64 ",%" PRId64 "\n", 
         csvOrigStr(), vec.x, vec.y, vec.z, htmId, hpxId);
      fflush(fout);

      ++row;
   }
   
   fflush(fout);
   fclose(fout);
    
   printf("[struct stat=\"OK\", nrow=%d]\n", row);
   fflush(stdout);
   exit(0);
}


//------------------------------------------------------------

FILE  *csvFile;

char   csvStr[MAXSTR];

int   *csvWidths;
char **csvColvals;

int    csvNcol, csvMaxcol;
int    csvNlines, csvFirst;


static int csvInit(char *fname)
{
   csvFile = fopen(fname, "r");

   if(csvFile == (FILE *)NULL)
      return 1;

   csvNcol   = 0;
   csvMaxcol = MAXCOL;

   csvWidths  = (int   *)malloc(csvMaxcol * sizeof(int));
   csvColvals = (char **)malloc(csvMaxcol * sizeof(char *));

   csvNlines  = 0;
   csvFirst   = 1;

   return 0;
}


static int csvReadRow()
{
   int    i, len, inquote;
   int    index, sublen;

   char   str     [MAXSTR];
   char   tmpstr  [MAXSTR];
   char   stripstr[MAXSTR];

   char  *end;
   char   delimChar;

   delimChar = ',';


   /* Read a line of input */

   if(fgets(str, MAXSTR, csvFile) == (char *)NULL)
   {
      fclose(csvFile);
      return -1;
   }

   while(strlen(str) > 0
         &&(   str[strlen(str) - 1] == '\n'
            || str[strlen(str) - 1] == '\r') )
      str[strlen(str) - 1]  = '\0';

   strcpy(csvStr, str);

   strcat(str, " ");

   len = strlen(str);

   if(debug)
   {
      printf("DEBUG> Input line: [%s](%d)\n", str, len);
      fflush(stdout);
   }


   /* Parse through the line, pulling out fields */

   index   = 0;
   end     = str;
   inquote = 0;

   while(1)
   {
      if(end >= str + len)
         break;

      while(end < str + len && *end == ' ')
         ++end;

      if(*end == '"')
      {
         inquote = 1;
         ++end;
      }


      // Collecting characters for the next field

      sublen = 0;

      while(1)
      {
         /* Hit the end of string */

         if(end >= str + len)
            break;


         /* Hit the next delimiter */

         if(!inquote && *end == delimChar)
         {
            ++end;
            break;
         }

         // Comma-delimited uses double double-quotes as an escaped double-quote

         if(inquote && *end == '"' && *(end+1) == '"')
         {
            tmpstr[sublen] = '"';
            ++sublen;

            end += 2;
            continue;
         }


         // End of the quoted string 

         if(inquote && *end == '"')
         {
            inquote = 0;
            ++end;

            while(end < str + len && *end != delimChar)
               ++end;

            ++end;

            break;
         }

         tmpstr[sublen] = *end;
         ++sublen;
         ++end;
      }


      /* We've collected the substring; save it */

      tmpstr[sublen] = '\0';

      strcpy(stripstr, csvStrip(tmpstr));

      if(csvFirst)
      {
         csvWidths [index] = strlen(stripstr)+32;

         csvColvals[index] = (char *)malloc(csvWidths[index] * sizeof(char));
      }
      else
      {
         if(strlen(stripstr)+32 > csvWidths[index])
         {
            csvWidths [index] = strlen(stripstr)+32;

            csvColvals[index] = (char *)realloc(csvColvals[index], csvWidths[index] * sizeof(char));
         }
      }

      strcpy(csvColvals[index], stripstr);

      if(debug)
      {
         printf("DEBUG> line %d, column %d: value = [%s], width = %d\n", 
               csvNlines, index, csvColvals[index], csvWidths[index]);
         fflush(stdout);
      }

      ++index;

      if(index > csvNcol)
         ++csvNcol;

      if(csvNcol >= csvMaxcol)
      {
         csvMaxcol += MAXCOL;

         csvWidths  = (int   *)realloc(csvWidths,  csvMaxcol * sizeof(int));
         csvColvals = (char **)realloc(csvColvals, csvMaxcol * sizeof(char *));

         for(i=csvMaxcol-MAXCOL; i<csvMaxcol; ++i)
            csvWidths[i] = 0;
      }
   }

   csvFirst = 0;

   ++csvNlines;

   return index;
}


static char *csvOrigStr()
{
   return csvStr;
}


static char *csvColVal(int col)
{
   return csvColvals[col];
}


static char *csvStrip(char *in)
{
   int   i, len;
   char *out;

   out = in;

   if(*out == '"')
      ++out;

   while(*out == ' ')
      ++out;
   
   len = strlen(out);

   for(i=len-1; i>=0; --i)
   {
      if(out[i] == ' ')
         out[i] = '\0';

      else if(out[i] == '"')
      {
         out[i] = '\0';
         break;
      }

      else
         break;
   }

      return(out);
}

//------------------------------------------------------------


// For consistence with the HEALPix library, we are using the
// same set of constants they have rather than deriving them
// through internal calculation.

static const double twothird = 2.0/3.0;

static const double pi     = 3.141592653589793238462643383279502884197;
static const double halfpi = 1.570796326794896619231321691639751442099;



// It is easiest to perform internal calculations using the 12 face numbers
// and internal "x,y" coordinates inside the face (treating each one as a
// square of nside * nside pixels).  The conversion of this to HEALPix pixel
// IDs can be done in real time but would involve a lot of unnecessary 
// repetition if done for a lot of locations so they have pre-calculated
// a bunch of constants that can just be looked up based on the 0-255 potential
// values 

/* utab[m] = (short)(
      (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
   | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7)); */

static const short utab[]={
      0,     1,     4,     5,    16,    17,    20,    21,    64,    65,    68,    69,    80,    81,    84,    85,
    256,   257,   260,   261,   272,   273,   276,   277,   320,   321,   324,   325,   336,   337,   340,   341,
   1024,  1025,  1028,  1029,  1040,  1041,  1044,  1045,  1088,  1089,  1092,  1093,  1104,  1105,  1108,  1109,
   1280,  1281,  1284,  1285,  1296,  1297,  1300,  1301,  1344,  1345,  1348,  1349,  1360,  1361,  1364,  1365,
   4096,  4097,  4100,  4101,  4112,  4113,  4116,  4117,  4160,  4161,  4164,  4165,  4176,  4177,  4180,  4181,  
   4352,  4353,  4356,  4357,  4368,  4369,  4372,  4373,  4416,  4417,  4420,  4421,  4432,  4433,  4436,  4437,
   5120,  5121,  5124,  5125,  5136,  5137,  5140,  5141,  5184,  5185,  5188,  5189,  5200,  5201,  5204,  5205,
   5376,  5377,  5380,  5381,  5392,  5393,  5396,  5397,  5440,  5441,  5444,  5445,  5456,  5457,  5460,  5461,
  16384, 16385, 16388, 16389, 16400, 16401, 16404, 16405, 16448, 16449, 16452, 16453, 16464, 16465, 16468, 16469,
  16640, 16641, 16644, 16645, 16656, 16657, 16660, 16661, 16704, 16705, 16708, 16709, 16720, 16721, 16724, 16725,
  17408, 17409, 17412, 17413, 17424, 17425, 17428, 17429, 17472, 17473, 17476, 17477, 17488, 17489, 17492, 17493,
  17664, 17665, 17668, 17669, 17680, 17681, 17684, 17685, 17728, 17729, 17732, 17733, 17744, 17745, 17748, 17749,
  20480, 20481, 20484, 20485, 20496, 20497, 20500, 20501, 20544, 20545, 20548, 20549, 20560, 20561, 20564, 20565,
  20736, 20737, 20740, 20741, 20752, 20753, 20756, 20757, 20800, 20801, 20804, 20805, 20816, 20817, 20820, 20821,
  21504, 21505, 21508, 21509, 21520, 21521, 21524, 21525, 21568, 21569, 21572, 21573, 21584, 21585, 21588, 21589, 
  21760, 21761, 21764, 21765, 21776, 21777, 21780, 21781, 21824, 21825, 21828, 21829, 21840, 21841, 21844, 21845 };


// The best way to understand this code is to analyze it in the context of
// "HEALPix: A FRAMEWORK FOR HIGH-RESOLUTION DISCRETIZATION AND FAST ANALYSIS
// OF DATA DISTRIBUTED ON THE SPHERE" (ApJ 622:759–771, 2005 April 1)
// with particular attention to Figure 4.



// Given the coordinates of a point on the sky 
// find the HEALPix "pixel" ID.

static int64_t sky2hpx(int order, int64_t nside, double lon, double lat)
{
   double dtr, z, phi;

   dtr = pi/180.;


   // z and phi are the coordinates best suited to dealing
   // with HEALPix calculations.  z (the colatitude) is by 
   // definition between -1 and +1.  phi is just longtude in
   // radians.

   z   = cos((90.-lat)*dtr);
   phi = lon * dtr;

   if(debug)
   {
      printf("DEBUG> (RA, Dec) : (%-g, %-g) -> (z, phi) (%g, %-g)\n",
         lon, lat, z, phi);
      fflush(stdout);
   }


   double zabs = fabs(z);  // Parameter used to determine if we are in
                           // equatorial or polar regime.

   // The parameter tt is phi scaled to the range 0-4

   double tt = phi / halfpi;

   if(debug)
   {
      printf("DEBUG> z = %-g, tt = %-g\n", z, tt);
      fflush(stdout);
   }


   // In the equatorial region, the faces are square boxes (diamonds)
   // in z,phi coordinates.  The "ascending" edge (on the left) therefore
   // corresponds to one cartesian axis for that face and the "descending"
   // edge (on the right) to the other axis.

   // Or sort of; if we were drawing thes "axes" they would start at 
   // the leftmost point of face zero (actually (315,0) in lon,lat) and
   // the first would run up and to the right off the diagram and the 
   // second would run down and to the right off the diagram.

   int64_t face_num, ix, iy;

   if (zabs <= twothird) /* Equatorial region */
   {
      double temp1 = nside * (0.5 + tt);  // The 0.5 here compensates for the fact that 
                                          // the first face is cut in half at longitude zero

      double temp2 = nside * (z * 0.75);  // The equatorial region runs from z=-2/3 to z=+2/3, so
                                          // the total range is 4/3 and we need to scale by 3/4


      // The above two parameters are still in z,phi (e.g. lon,lat) orientation.
      // These two are in our rotated coordinates described above running along
      // the sides of the faces.

      int64_t jp = (int64_t)(temp1-temp2);  // index of  ascending edge line
      int64_t jm = (int64_t)(temp1+temp2);  // index of descending edge line

      if(debug)
      {
         printf("DEBUG> jp = %" PRId64 ", jm = %" PRId64 "\n", jp, jm);
         fflush(stdout);
      }


      // There are 12 "faces" in HEALPix.  The equatorial region
      // includes the four along the equator and half of the four
      // around each pole.  First we determine which face we are in.

      int64_t ifp = jp/nside;  // These are the face "coordinates" using
      int64_t ifm = jm/nside;  // using the same rotate scheme.

      if(debug)
      {
         printf("DEBUG> ifp = %" PRId64 ", ifm = %" PRId64 "\n", ifp, ifm);
         fflush(stdout);
      }


      if (ifp == ifm)           // faces 4 to 7; the ones along the equator
      {                         // (the diagonal in this rotated system).
         if(ifp == 4)
            face_num = 4;       // ifp = 4
         else
            face_num = ifp + 4; // ifp = 1,2,3 -> 5,6,7
      }

      else if (ifp < ifm)       // (half-)faces 0 to 3 (the north pole)
         face_num = ifp;

      else                      // (half-)faces 8 to 11 (the south pole)
         face_num = ifm + 8;


      // The "x,y" coordinates inside a face run from
      // zero to nside-1. We can get this by masking the
      // jm, jp coordinates with nside-1 (i.e. mask off
      // the low N bits).  y has to be reversed.

      ix = jm & (nside-1);
      iy = (nside-1) - (jp & (nside-1));
   }


   // The polar regions (zabs > 2/3) use similar logic, though 
   // here the ascending side is running up to the pole and the 
   // descending side is along the next quadrant running down
   // the pole (this same edge is then the ascending side for the
   // next quadrant/face.

   else
   {
      double tp, tmp;
      int    ntt, jp, jm;

      ntt = (int)tt;

      if(ntt>=4)
         ntt=3;

      tp = tt-ntt;

      tmp = nside * sqrt(3. * (1.-zabs));  // Also, the face coordinates
                                           // are not linear in z, rather
                                           // the coordinate has to be scaled.


      jp = (int)(tp*tmp);        // increasing edge line index
      jm = (int)((1.0-tp)*tmp);  // decreasing edge line index

      if (jp>=nside) jp = nside-1;  // for points too close to the boundary
      if (jm>=nside) jm = nside-1;

      if (z >= 0)
      {
         face_num = ntt;     // North: 0-3

         ix = (nside-1) - jm;
         iy = (nside-1) - jp;
      }
      else
      {
         face_num = ntt + 8; // South: 8-11

         ix =  jp;
         iy =  jm;
      }
   }

   if(debug)
   {
      printf("DEBUG> face_num = %" PRId64 "\n",       face_num);
      printf("DEBUG> ix, iy   = %" PRId64 ", %" PRId64 "\n", ix, iy);
      fflush(stdout);
   }


   // This ix,iy "cartesian" coordinate inside a face needs to 
   // be converted to the HEALPix numbering scheme

   return xyf2nest(order, ix, iy, face_num);
}


int64_t xyf2nest (int order, int64_t ix, int64_t iy, int64_t face_num)
{
   return (face_num<<(2*order)) + spread_bits(ix) + (spread_bits(iy)<<1);
}


int64_t spread_bits (int v)
{
  return  (int64_t)(utab[ v     &0xff])      | ((int64_t)(utab[(v>> 8)&0xff])<<16)
       | ((int64_t)(utab[(v>>16)&0xff])<<32) | ((int64_t)(utab[(v>>24)&0xff])<<48);
}
