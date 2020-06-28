/*
  Copyright (c) 2020, Caltech IPAC.

  License information at
    https://github.com/Caltech-IPAC/nexsciTAP/blob/master/LICENSE
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sptQuery.h>


struct sptConstraints sptConeSearch(char *indname, int indexMode, int indexEncoding, int level, 
                                    char *xcol, char *ycol, char *zcol, double ra, double dec, double radius)
{
   int    i, ind, nranges;
   double maxpct, cos_dec, cosrad;
   double sumpct, pct, total, newtotal;

   char   tmpstr[1024];

   struct sptSkiplist *skipped;
   struct sptRanges   *list;

   struct  htm_v3 vec;
   int64_t id;

   int64_t nside, npface, skytotal;

   struct htm_ids *ids = (struct htm_ids *)NULL;

   enum   htm_errcode err;

   struct hpxRanges *hpxrng;

   struct sptConstraints constraints;


   constraints.status = 0;

   strcpy(constraints.errorMsg, "");

   strcpy(constraints.indexConstraint, "");
   strcpy(constraints.geomConstraint,  "");


   maxpct = 30.;

   if(level > 20)
   {
      sprintf(constraints.errorMsg, "Indexing level cannot be greater than 20 (which is already sub-arcsecond).");

      constraints.status = 1;

      return constraints;
   }


   // Radius

   cosrad = cos(radius * HTM_RAD_PER_DEG);


   // Get 3-vector for location

   cos_dec = cos(dec * HTM_RAD_PER_DEG);

   vec.x = cos( ra * HTM_RAD_PER_DEG) * cos_dec;
   vec.y = sin( ra * HTM_RAD_PER_DEG) * cos_dec;
   vec.z = sin(dec * HTM_RAD_PER_DEG);

   htm_v3_normalize(&vec, &vec);

   if(sptDebug)
   {
      printf("\n(RA,Dec) = (%10.6f,%10.6f)\n",      ra, dec);
      printf(  "(x,y,z)  = (%11.8f,%11.8f,%11.8f)\n", vec.x, vec.y, vec.z);
      fflush(stdout);
   }


   if(indexMode == HTM)
   {
      // Find the center cell ID

      id = htm_v3_id(&vec, level);

      if(sptDebug)
      {
         if(indexEncoding == 1)
            printf("\nCenter cell: %" PRId64 "  (%" PRId64 ")\n\n", id, htm_idtodec(id));
         else
            printf("\nCenter cell: %" PRId64 ")\n\n", id);

         printf("%llu cells, cell size ~%-g\n\n", (1ULL << (2*level)), 90./(1<<level));
         fflush(stdout);
      }


      // Get a list of the exact ranges for the region

      ids = htm_s2circle_ids(NULL, &vec, radius, level, SIZE_MAX, &err);

      // Process the HTM ranges

      nranges = ids->n;

      if(nranges <= 0)
      {
         strcpy(constraints.errorMsg, "No HTM ranges found (a physical impossibility so there is a bug in the code).");

         constraints.status = 1;

         return constraints;
      }

      if(sptDebug)
      {
         printf("\n  Radius %-g degrees  (%d Ranges)\n\n", radius, nranges);
         fflush(stdout);
      }


      // Process all this to create a linked-list of ranges
      // and info on the cell IDs that were skipped

      if(sptDebug)
      {
         printf("\n# HTM CONSTRAINT -----------------------------------------\n");
         fflush(stdout);
      }

      total = 0.;

      skipped = (struct sptSkiplist *)malloc(nranges * sizeof(struct sptSkiplist));
      list    = (struct sptRanges *)  malloc(nranges * sizeof(struct sptRanges));

      if(sptDebug)
      {
         printf("\n");
         fflush(stdout);
      }

      for(i=0; i<nranges; ++i)
      {
         list[i].removeFlag = 0;

         list[i].min = ids->range[i].min;
         list[i].max = ids->range[i].max;

         list[i].prev = -1;
         list[i].next = -1;

         if(i > 0)         list[i].prev = i-1;
         if(i < nranges-1) list[i].next = i+1;

         skipped[i].index = i;

         if(i == 0)
            skipped[i].count = 0;
         else
            skipped[i].count = ids->range[i].min - ids->range[i-1].max;

         total += ids->range[i].max - ids->range[i].min + 1.;
      }


      // Info printout

      if(sptDebug)
      {
         for(i=0; i<nranges; ++i)
         {
            if(indexEncoding == 1)
               printf("%3d: %" PRId64 " -> %" PRId64 " : %" PRId64 " (%" PRId64 " -> %" PRId64 ") skipped %" PRIu64 "\n", 
                     i, list[i].min, list[i].max, list[i].max - list[i].min + 1,
                     htm_idtodec(list[i].min), htm_idtodec(list[i].max), skipped[i].count);
            else
               printf("%3d: %" PRId64 " -> %" PRId64 " : %" PRId64 " skipped %" PRIu64 "\n", 
                     i, list[i].min, list[i].max, list[i].max - list[i].min + 1, skipped[i].count);
         }

         printf("\ntotal = %-g\nrange = %" PRId64 "\n\n", total, list[nranges-1].max - list[0].min);
      }


      // Now sort the skipped list so we can find
      // small gaps to compress out

      qsort(skipped, nranges, sizeof(struct sptSkiplist), sptSortCmpFunc);

      sumpct = 0.;

      for(i=0; i<nranges; ++i)
      {
         ind = skipped[i].index;

         pct = 100. * skipped[i].count / total;

         sumpct += pct;

         if(sumpct < maxpct && ind > 0)
         {
            if(sptDebug)
            {
               printf("Reassign index %d (%" PRIu64 " cells)\n", ind, skipped[i].count);
               fflush(stdout);
            }

            list[ind].removeFlag = 1;
         }
      }


      // Compress out until we read the max percentage addition allowed
      // We do this by going through the skipped list up to cutoff and 
      // modifying the linked list as we go (i.e. adding the range pointed
      // to to the preceding range and removing that node).

      for(i=1; i<nranges; ++i)
      {
         if(list[i].removeFlag == 1)
         {
            list[list[i].prev].max  = list[i].max;
            list[list[i].prev].next = list[i].next;

            if(list[i].next >= 0)
               list[list[i].next].prev = list[i].prev;

            list[i].min        = 0;
            list[i].max        = 0;
            list[i].removeFlag = -1;
            list[i].prev       = -1;
            list[i].next       = -1;
         }
      }


      // And print out the final ranges

      id = 0;

      i = 0;

      if(sptDebug)
      {
         printf("WHERE (\n");
         fflush(stdout);
      }

      newtotal = 0.;

      strcpy(constraints.indexConstraint, "");

      while(1)
      {
         if(list[id].min == list[id].max)
         {
            newtotal += 1.;

            if(sptDebug)
            {
               if(indexEncoding == 1)
               {
                  if(i == 0)
                     printf("      (htm%d = %" PRId64 ")                         // %" PRId64 "\n",
                           level, list[id].min, htm_idtodec(list[id].min));
                  else
                     printf("   OR (htm%d = %" PRId64 ")                         // %" PRId64 "\n",
                           level, list[id].min, htm_idtodec(list[id].min));
               }
               else
               {
                  if(i == 0)
                     printf("      (htm%d = %" PRId64 ")\n", level, list[id].min);
                  else
                     printf("   OR (htm%d = %" PRId64 ")\n", level, list[id].min);
               }
            }
            
            if(indexEncoding == 1)
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s = %" PRId64 ")", indname, htm_idtodec(list[id].min));
               else
                  sprintf(tmpstr, " OR (%s = %" PRId64 ")", indname, htm_idtodec(list[id].min));
            }
            else
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s = %" PRId64 ")", indname, list[id].min);
               else
                  sprintf(tmpstr, " OR (%s = %" PRId64 ")", indname, list[id].min);
            }
         }
         else
         {
            newtotal += list[id].max - list[id].min + 1.;

            if(sptDebug)
            {
               if(indexEncoding == 1)
               {
                  if(i == 0)
                     printf("      (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")  // %" PRId64 " to  %" PRId64 "\n",
                           level, list[id].min, list[id].max, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
                  else
                     printf("   OR (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")  // %" PRId64 " to  %" PRId64 "\n",
                           level, list[id].min, list[id].max, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
               }
               else
               {
                  if(i == 0)
                     printf("      (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")\n", level, list[id].min, list[id].max);
                  else
                     printf("   OR (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")\n", level, list[id].min, list[id].max);
               }
            }
            
            if(indexEncoding == 1)
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
               else
                  sprintf(tmpstr, " OR (%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
            }
            else
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, list[id].min, list[id].max);
               else
                  sprintf(tmpstr, " OR (%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, list[id].min, list[id].max);
            }
         }

         strcat(constraints.indexConstraint, tmpstr);

         if(list[id].next == -1)
            break;

         id = list[id].next;

         ++i;
      }
      
      if(sptDebug)
      {
         printf(")\n");
         fflush(stdout);

         printf("\nnew total = %-g\n", newtotal);
         fflush(stdout);
      }

      skytotal = 8 * (1LL << 2*level);

      if(sptDebug)
      {
         printf("\n# %-g cells in query; %" PRIu64 " cells covering whole sky\n\n",
            newtotal, skytotal);
         fflush(stdout);
      }
   }


   else  // HPX mode
   {
      hpxrng = (struct hpxRanges *)NULL;

      nranges = hpxConeSearch(level, ra, dec, radius, &hpxrng);

      if(nranges <= 0)
      {
         sprintf(constraints.errorMsg, "No HEALPix ranges found (a physical impossibility so there is a bug in the code).");

         constraints.status = 1;

         return constraints;
      }

      skipped = (struct sptSkiplist *)malloc(nranges * sizeof(struct sptSkiplist));
      list    = (struct sptRanges *)  malloc(nranges * sizeof(struct sptRanges));

      if(sptDebug)
      {
         printf("\n");
         fflush(stdout);
      }

      total = 0;

      for(i=0; i<nranges; ++i)
      {
         list[i].removeFlag = 0;

         list[i].min = hpxrng[i].min;
         list[i].max = hpxrng[i].max;

         list[i].prev = -1;
         list[i].next = -1;

         if(i > 0)         list[i].prev = i-1;
         if(i < nranges-1) list[i].next = i+1;

         skipped[i].index = i;

         if(i == 0)
            skipped[i].count = 0;
         else
            skipped[i].count = hpxrng[i].min - hpxrng[i-1].max;

         total += hpxrng[i].max - hpxrng[i].min + 1.;
      }

      free(hpxrng);


      // Info printout

      if(sptDebug)
      {
         for(i=0; i<nranges; ++i)
         {
            printf("%3d: %" PRId64 " -> %" PRId64 " : %" PRId64 " skipped %" PRIu64 "\n", 
                  i, list[i].min, list[i].max, list[i].max - list[i].min + 1, skipped[i].count);
         }

         printf("\ntotal = %-g\nrange = %" PRId64 "\n\n", total, list[nranges-1].max - list[0].min);
      }


      // Now sort the skipped list so we can find
      // small gaps to compress out

      qsort(skipped, nranges, sizeof(struct sptSkiplist), sptSortCmpFunc);

      sumpct = 0.;

      for(i=0; i<nranges; ++i)
      {
         ind = skipped[i].index;

         pct = 100. * skipped[i].count / total;

         sumpct += pct;

         if(sumpct < maxpct && ind > 0)
         {
            if(sptDebug)
            {
               printf("Reassign index %d (%" PRIu64 " cells)\n", ind, skipped[i].count);
               fflush(stdout);
            }

            list[ind].removeFlag = 1;
         }
      }


      // Compress out until we read the max percentage addition allowed
      // We do this by going through the skipped list up to cutoff and 
      // modifying the linked list as we go (i.e. adding the range pointed
      // to to the preceding range and removing that node).

      for(i=1; i<nranges; ++i)
      {
         if(list[i].removeFlag == 1)
         {
            list[list[i].prev].max  = list[i].max;
            list[list[i].prev].next = list[i].next;

            if(list[i].next >= 0)
               list[list[i].next].prev = list[i].prev;

            list[i].min        = 0;
            list[i].max        = 0;
            list[i].removeFlag = -1;
            list[i].prev       = -1;
            list[i].next       = -1;
         }
      }


      // And print out the final ranges

      id = 0;

      i = 0;

      if(sptDebug)
      {
         printf("WHERE (\n");
         fflush(stdout);
      }

      newtotal = 0.;

      strcpy(constraints.indexConstraint, "");

      while(1)
      {
         if(list[id].min == list[id].max)
         {
            newtotal += 1.;

            if(sptDebug)
            {
               if(i == 0)
                  printf("      (hpx%d = %" PRIu64 ")\n", level, list[id].min);
               else
                  printf("   OR (hpx%d = %" PRIu64 ")\n", level, list[id].min);
            }
            
            if(i == 0)
               sprintf(tmpstr, "(%s = %" PRIu64 ")", indname, list[id].min);
            else
               sprintf(tmpstr, " OR (%s = %" PRIu64 ")", indname, list[id].min);
         }
         else
         {
            newtotal += list[id].max - list[id].min + 1.;

            if(sptDebug)
            {
               if(i == 0)
                  printf("      (%s BETWEEN %" PRIu64 " AND %" PRIu64 ")\n", indname, list[id].min, list[id].max);
               else
                  printf("   OR (%s BETWEEN %" PRIu64 " AND %" PRIu64 ")\n", indname, list[id].min, list[id].max);
            }
            
            if(i == 0)
               sprintf(tmpstr, "(%s BETWEEN %" PRIu64 " AND %" PRIu64 ")", indname, list[id].min, list[id].max);
            else
               sprintf(tmpstr, " OR (%s BETWEEN %" PRIu64 " AND %" PRIu64 ")", indname, list[id].min, list[id].max);
         }

         strcat(constraints.indexConstraint, tmpstr);

         if(list[id].next == -1)
            break;

         id = list[id].next;

         ++i;
      }

      if(sptDebug)
      {
         printf(")\n");
         fflush(stdout);

         printf("\nnew total = %-g\n", newtotal);
         fflush(stdout);
      }

      nside    = 1ULL << level;
      npface   = nside << level;
      skytotal = 12 * npface;

      if(sptDebug)
      {
         printf("\n# %-g cells in query; %" PRIu64 " cells covering whole sky\n\n",
               newtotal, skytotal);
         fflush(stdout);
      }
   }

   if(sptDebug)
   {
      printf("# EXACT CONE CONSTRAINT -----------------------------------\n\n");
      fflush(stdout);

      printf("WHERE (%.12f*%s)+(%.12f*%s)+(%.12f*%s)>=%.12e\n\n", vec.x, xcol, vec.y, ycol, vec.z, zcol, cosrad);

      printf("# ---------------------------------------------------------\n\n");
      fflush(stdout);
   }

   sprintf(constraints.geomConstraint, "(%.12f*%s)+(%.12f*%s)+(%.12f*%s)>=%.12e", vec.x, xcol, vec.y, ycol, vec.z, zcol, cosrad);

   free(list);
   free(skipped);

   return constraints;
}



struct sptConstraints sptPolygonSearch(char *indname, int indexMode, int indexEncoding, int level, 
                                       char *xcol, char *ycol, char *zcol, int npoly, double *ra, double *dec)
{
   int    i, flip, ind, nranges, nreverse;

   double cos_dec, total, sumpct, pct, maxpct, newtotal, handedness;

   char   tmpstr[1024];

   struct sptSkiplist *skipped;
   struct sptRanges   *list;

   int64_t id;

   int64_t nside, npface, skytotal;

   struct htm_ids *ids = (struct htm_ids *)NULL;

   enum   htm_errcode err;

   struct htm_s2cpoly *htmPoly;

   struct htm_v3 *verts;

   Vec    *poly, *normal;

   struct sptConstraints constraints;


   constraints.status = 0;

   strcpy(constraints.errorMsg, "");

   strcpy(constraints.indexConstraint, "");
   strcpy(constraints.geomConstraint,  "");


   maxpct = 30.;

   struct hpxRanges *hpxrng;


   // Calculate the 3-vectors for the polygon vertices

   poly = (Vec *)malloc(npoly * sizeof(Vec));

   for(i=0; i<npoly; ++i)
   {
      poly[i].x = cos(ra[i]*HTM_RAD_PER_DEG) * cos(dec[i]*HTM_RAD_PER_DEG);
      poly[i].y = sin(ra[i]*HTM_RAD_PER_DEG) * cos(dec[i]*HTM_RAD_PER_DEG);
      poly[i].z = sin(dec[i]*HTM_RAD_PER_DEG);
   }


   // Find the normal vectors for each of the polygon sides
   // Correct for handedness of the polygon

   normal = (Vec *)malloc(npoly * sizeof(Vec));

   flip = 1;

   nreverse = 0;

   for(i=0; i<npoly; ++i)
   {
      vCross(&poly[i], &poly[(i+1)%npoly], &normal[i]);
      vNormalize(&normal[i]);

      handedness = vDot(&normal[i], &poly[(i+2)%npoly]);

      if(sptDebug)
      {
         printf("Normal (x,y,z)  = (%11.8f,%11.8f,%11.8f) [handedness: %11.8f]\n",
            normal[i].x, normal[i].y, normal[i].z, handedness);
      }

      if(fabs(handedness) < 1.e-10)
      {
         strcpy(constraints.errorMsg, "Degenerate polygon corner.");

         constraints.status = 1;

         return constraints;
      }

      if(i == 0 && handedness < 0.)
         flip = -1.;

      if(i > 0 && flip*handedness <= 0.)
      {
         strcpy(constraints.errorMsg, "Polygon is not convex;");

         constraints.status = 1;

         return constraints;
      }

      if(flip < 0)
      {
         if(sptDebug)
            printf("DEBUG> Reversing normal %d\n", i);

         vReverse(&normal[i]);

         ++nreverse;
      }
   }

   if(sptDebug && nreverse > 0)
      printf("\n%d reversed\n\n", nreverse);

   if(indexMode == HTM)
   {
      verts = (struct htm_v3 *)malloc(npoly * sizeof(struct htm_v3));

      for(i=0; i<npoly; ++i)
      {
         cos_dec = cos(dec[i]*HTM_RAD_PER_DEG);

         verts[i].x = cos( ra[i]*HTM_RAD_PER_DEG) * cos_dec;
         verts[i].y = sin( ra[i]*HTM_RAD_PER_DEG) * cos_dec;
         verts[i].z = sin(dec[i]*HTM_RAD_PER_DEG);

         htm_v3_normalize(&verts[i], &verts[i]);

         if(sptDebug)
         {
            printf("(RA,Dec) = (%10.6f,%10.6f), ",      ra[i], dec[i]);
            printf("(x,y,z)  = (%11.8f,%11.8f,%11.8f)\n", verts[i].x, verts[i].y, verts[i].z);
            fflush(stdout);
         }
      }

      if(sptDebug)
      {
         printf("\nFinal normals:\n");
         for(i=0; i<npoly; ++i)
         {
            printf("(x,y,z)  = (%11.8f,%11.8f,%11.8f)\n", normal[i].x, normal[i].y, normal[i].z);
            fflush(stdout);
         }
      }


      htmPoly = htm_s2cpoly_init(verts, (size_t)npoly, &err);

      ids = htm_s2cpoly_ids(ids, htmPoly, level, SIZE_MAX, &err);   // Spherical convex polygon


      // Process the HTM ranges

      nranges = ids->n;

      if(nranges <= 0)
      {
         sprintf(constraints.errorMsg, "No HTM ranges found (a physical impossibility so there is a bug in the code).");

         constraints.status = 1;

         return constraints;
      }


      // Process all this to create a linked-list of ranges
      // and info on the cell IDs that were skipped

      total = 0.;

      skipped = (struct sptSkiplist *)malloc(nranges * sizeof(struct sptSkiplist));
      list    = (struct sptRanges *)  malloc(nranges * sizeof(struct sptRanges));

      if(sptDebug)
      {
         printf("\n");
         fflush(stdout);
      }

      for(i=0; i<nranges; ++i)
      {
         list[i].removeFlag = 0;

         list[i].min = ids->range[i].min;
         list[i].max = ids->range[i].max;

         list[i].prev = -1;
         list[i].next = -1;

         if(i > 0)         list[i].prev = i-1;
         if(i < nranges-1) list[i].next = i+1;

         skipped[i].index = i;

         if(i == 0)
            skipped[i].count = 0;
         else
            skipped[i].count = ids->range[i].min - ids->range[i-1].max;

         total += ids->range[i].max - ids->range[i].min + 1.;
      }


      // Info printout

      if(sptDebug)
      {
         for(i=0; i<nranges; ++i)
         {
            if(indexEncoding == 1)
            {
               printf("%3d: %" PRId64 " -> %" PRId64 " : %" PRId64 " (%" PRId64 " -> %" PRId64 ") skipped %" PRIu64 "\n", 
                     i, list[i].min, list[i].max, list[i].max - list[i].min + 1,
                     htm_idtodec(list[i].min), htm_idtodec(list[i].max), skipped[i].count);
            }
            else
            {
               printf("%3d: %" PRId64 " -> %" PRId64 " : %" PRId64 " skipped %" PRIu64 "\n", 
                     i, list[i].min, list[i].max, list[i].max - list[i].min + 1, skipped[i].count);
            }
         }

         printf("\ntotal = %-g\nrange = %" PRId64 "\n\n", total, list[nranges-1].max - list[0].min);
      }


      // Now sort the skipped list so we can find
      // small gaps to compress out

      qsort(skipped, nranges, sizeof(struct sptSkiplist), sptSortCmpFunc);

      sumpct = 0.;

      for(i=0; i<nranges; ++i)
      {
         ind = skipped[i].index;

         pct = 100. * skipped[i].count / total;

         sumpct += pct;

         if(sumpct < maxpct && ind > 0)
         {
            if(sptDebug)
            {
               printf("Reassign index %d (%" PRIu64 " cells)\n", ind, skipped[i].count);
               fflush(stdout);
            }

            list[ind].removeFlag = 1;
         }
      }


      // Compress out until we read the max percentage addition allowed
      // We do this by going through the skipped list up to cutoff and 
      // modifying the linked list as we go (i.e. adding the range pointed
      // to to the preceding range and removing that node).

      for(i=1; i<nranges; ++i)
      {
         if(list[i].removeFlag == 1)
         {
            list[list[i].prev].max  = list[i].max;
            list[list[i].prev].next = list[i].next;

            if(list[i].next >= 0)
               list[list[i].next].prev = list[i].prev;

            list[i].min        = 0;
            list[i].max        = 0;
            list[i].removeFlag = -1;
            list[i].prev       = -1;
            list[i].next       = -1;
         }
      }


      // And print out the final ranges

      id = 0;

      i = 0;

      if(sptDebug)
      {
         printf("WHERE (\n");
         fflush(stdout);
      }

      newtotal = 0.;

      strcat(constraints.indexConstraint, "");

      while(1)
      {
         if(list[id].min == list[id].max)
         {
            newtotal += 1.;

            if(sptDebug)
            {
               if(indexEncoding == 1)
               {
                  if(i == 0)
                     printf("      (htm%d = %" PRId64 ")                         // %" PRId64 "\n",
                           level, list[id].min, htm_idtodec(list[id].min));
                  else
                     printf("   OR (htm%d = %" PRId64 ")                         // %" PRId64 "\n",
                           level, list[id].min, htm_idtodec(list[id].min));
               }
               else
               {
                  if(i == 0)
                     printf("      (htm%d = %" PRId64 ")\n", level, list[id].min);
                  else
                     printf("   OR (htm%d = %" PRId64 ")\n", level, list[id].min);
               }
            }
            
            if(indexEncoding == 1)
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s = %" PRId64 ")", indname, htm_idtodec(list[id].min));
               else
                  sprintf(tmpstr, " OR (%s = %" PRId64 ")", indname, htm_idtodec(list[id].min));
            }
            else
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s = %" PRId64 ")", indname, list[id].min);
               else
                  sprintf(tmpstr, " OR (%s = %" PRId64 ")", indname, list[id].min);
            }
         }
         else
         {
            newtotal += list[id].max - list[id].min + 1.;

            if(sptDebug)
            {
               if(indexEncoding == 1)
               {
                  if(i == 0)
                     printf("      (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")  // %" PRId64 " to  %" PRId64 "\n",
                           level, list[id].min, list[id].max, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
                  else
                     printf("   OR (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")  // %" PRId64 " to  %" PRId64 "\n",
                           level, list[id].min, list[id].max, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
               }
               else
               {
                  if(i == 0)
                     printf("      (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")\n", level, list[id].min, list[id].max);
                  else
                     printf("   OR (htm%d BETWEEN %" PRId64 " AND %" PRId64 ")\n", level, list[id].min, list[id].max);
               }
            }
            
            if(indexEncoding == 1)
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
               else
                  sprintf(tmpstr, " OR (%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, htm_idtodec(list[id].min), htm_idtodec(list[id].max));
            }
            else
            {
               if(i == 0)
                  sprintf(tmpstr, "(%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, list[id].min, list[id].max);
               else
                  sprintf(tmpstr, " OR (%s BETWEEN %" PRId64 " AND %" PRId64 ")", indname, list[id].min, list[id].max);
            }
         }

         if(list[id].next == -1)
            break;

         id = list[id].next;

         ++i;

         strcat(constraints.indexConstraint, tmpstr);
      }

      if(sptDebug)
      {
         printf(")\n");
         fflush(stdout);

         printf("\nnew total = %-g\n", newtotal);
         fflush(stdout);
      }

      skytotal = 8 * (1LL << 2*level);

      if(sptDebug)
      {
         printf("\n# %-g cells in query; %" PRIu64 " cells covering whole sky\n\n",
               newtotal, skytotal);
         fflush(stdout);
      }

      free(list);
      free(skipped);
      free(verts);
   }


   else  // HPX mode
   {
      hpxrng = (struct hpxRanges *)NULL;

      nranges = hpxPolygonSearch(level, npoly, ra, dec, &hpxrng);

      if(nranges <= 0)
      {
         sprintf(constraints.errorMsg, "No HEALPix ranges found (a physical impossibility so there is a bug in the code).");

         constraints.status = 1;

         return constraints;
      }


      skipped = (struct sptSkiplist *)malloc(nranges * sizeof(struct sptSkiplist));
      list    = (struct sptRanges *)  malloc(nranges * sizeof(struct sptRanges));

      if(sptDebug)
      {
         printf("\n");
         fflush(stdout);
      }

      total = 0;

      for(i=0; i<nranges; ++i)
      {
         list[i].removeFlag = 0;

         list[i].min = hpxrng[i].min;
         list[i].max = hpxrng[i].max;

         list[i].prev = -1;
         list[i].next = -1;

         if(i > 0)         list[i].prev = i-1;
         if(i < nranges-1) list[i].next = i+1;

         skipped[i].index = i;

         if(i == 0)
            skipped[i].count = 0;
         else
            skipped[i].count = hpxrng[i].min - hpxrng[i-1].max;

         total += hpxrng[i].max - hpxrng[i].min + 1.;
      }

      free(hpxrng);


      // Info printout

      if(sptDebug)
      {
         for(i=0; i<nranges; ++i)
         {
            printf("%3d: %" PRId64 " -> %" PRId64 " : %" PRId64 " skipped %" PRIu64 "\n", 
                  i, list[i].min, list[i].max, list[i].max - list[i].min + 1, skipped[i].count);
         }

         printf("\ntotal = %-g\nrange = %" PRId64 "\n\n", total, list[nranges-1].max - list[0].min);
      }


      // Now sort the skipped list so we can find
      // small gaps to compress out

      qsort(skipped, nranges, sizeof(struct sptSkiplist), sptSortCmpFunc);

      sumpct = 0.;

      for(i=0; i<nranges; ++i)
      {
         ind = skipped[i].index;

         pct = 100. * skipped[i].count / total;

         sumpct += pct;

         if(sumpct < maxpct && ind > 0)
         {
            if(sptDebug)
            {
               printf("Reassign index %d (%" PRIu64 " cells)\n", ind, skipped[i].count);
               fflush(stdout);
            }

            list[ind].removeFlag = 1;
         }
      }


      // Compress out until we read the max percentage addition allowed
      // We do this by going through the skipped list up to cutoff and 
      // modifying the linked list as we go (i.e. adding the range pointed
      // to to the preceding range and removing that node).

      for(i=1; i<nranges; ++i)
      {
         if(list[i].removeFlag == 1)
         {
            list[list[i].prev].max  = list[i].max;
            list[list[i].prev].next = list[i].next;

            if(list[i].next >= 0)
               list[list[i].next].prev = list[i].prev;

            list[i].min        = 0;
            list[i].max        = 0;
            list[i].removeFlag = -1;
            list[i].prev       = -1;
            list[i].next       = -1;
         }
      }


      // And print out the final ranges

      id = 0;

      i = 0;

      if(sptDebug)
      {
         printf("WHERE (\n");
         fflush(stdout);
      }

      newtotal = 0.;

      strcat(constraints.indexConstraint, "");

      while(1)
      {
         if(list[id].min == list[id].max)
         {
            newtotal += 1.;

            if(sptDebug)
            {
               if(i == 0)
                  printf("      (%s = %" PRIu64 ")\n", indname, list[id].min);
               else
                  printf("   OR (%s = %" PRIu64 ")\n", indname, list[id].min);
            }

            if(i == 0)
               sprintf(tmpstr, "      (%s = %" PRIu64 ")", indname, list[id].min);
            else
               sprintf(tmpstr, "   OR (%s = %" PRIu64 ")", indname, list[id].min);
         }
         else
         {
            newtotal += list[id].max - list[id].min + 1.;

            if(sptDebug)
            {
               if(i == 0)
                  printf("      (%s BETWEEN %" PRIu64 " AND %" PRIu64 ")\n", indname, list[id].min, list[id].max);
               else
                  printf("   OR (%s BETWEEN %" PRIu64 " AND %" PRIu64 ")\n", indname, list[id].min, list[id].max);
            }

            if(i == 0)
               sprintf(tmpstr, "(%s BETWEEN %" PRIu64 " AND %" PRIu64 ")", indname, list[id].min, list[id].max);
            else
               sprintf(tmpstr, " OR (%s BETWEEN %" PRIu64 " AND %" PRIu64 ")", indname, list[id].min, list[id].max);
         }

         if(list[id].next == -1)
            break;

         id = list[id].next;

         ++i;

         strcat(constraints.indexConstraint, tmpstr);
      }

      free(list);
      free(skipped);
   }

   if(sptDebug)
   {
      printf(")\n");
      fflush(stdout);

      printf("\nnew total = %-g\n", newtotal);
      fflush(stdout);
   }

   nside    = 1ULL << level;
   npface   = nside << level;
   skytotal = 12 * npface;

   if(sptDebug)
   {
      printf("\n# %-g cells in query; %" PRIu64 " cells covering whole sky\n\n",
         newtotal, skytotal);
      fflush(stdout);
   }


   if(sptDebug)
   {
      printf("# EXACT POLYGON CONSTRAINT ----------------------------------\n\n");
      fflush(stdout);

      printf("WHERE (%.12f*x)+(%.12f*y)+(%.12f*z)>=0.\n", normal[0].x, normal[0].y, normal[0].z);

      for(i=1; i<npoly; ++i)
         printf("AND   (%.12f*x)+(%.12f*y)+(%.12f*z)>=0.\n", normal[i].x, normal[i].y, normal[i].z);

      printf("# ---------------------------------------------------------\n\n");
      fflush(stdout);
   }


   strcpy(constraints.geomConstraint, "");

   sprintf(tmpstr, "(%.12f*%s)+(%.12f*%s)+(%.12f*%s)>=0.", normal[0].x, xcol, normal[0].y, ycol, normal[0].z, zcol);

   strcat(constraints.geomConstraint, tmpstr);

   for(i=1; i<npoly; ++i)
   {
      sprintf(tmpstr, " AND (%.12f*%s)+(%.12f*%s)+(%.12f*%s)>=0.", normal[i].x, xcol, normal[i].y, ycol, normal[i].z, zcol);

      strcat(constraints.geomConstraint, tmpstr);
   }

   free(poly);
   free(normal);

   return constraints;
}



// Structure comparison function

int sptSortCmpFunc(const void *a, const void *b)
{
   int64_t ca, cb;

   ca = ((struct sptSkiplist *)a)->count;
   cb = ((struct sptSkiplist *)b)->count;

   if(ca > cb) return  1;
   if(ca < cb) return -1;

   return 0;
}

// -----------------------------------------------------------------

#define MAXSTACK 1024
#define MAXRANGE 1024


// For consistence with the HEALPix library, we are using the
// same set of constants they have rather than deriving them
// through internal calculation.

static const double pi     = 3.141592653589793238462643383279502884197;
static const double halfpi = 1.570796326794896619231321691639751442099;


// It is easiest to perform internal calculations using the 12 face numbers
// and internal "x,y" coordinates inside the face (treating each one as a
// square of nside * nside pixels).  The conversion of this to HEALPix pixel
// IDs can be done in real time but would involve a lot of unnecessary 
// repetition if done for a lot of locations so they have pre-calculated
// a bunch of constants that can just be looked up based on the 0-255 potential
// values 

/* ctab[m] = (short)(
       (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
    | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4)); */

static const short ctab[]={
     0,    1,  256,  257,    2,    3,  258,  259,  512,  513,  768,  769,  514,  515,  770,  771,
     4,    5,  260,  261,    6,    7,  262,  263,  516,  517,  772,  773,  518,  519,  774,  775,
  1024, 1025, 1280, 1281, 1026, 1027, 1282, 1283, 1536, 1537, 1792, 1793, 1538, 1539, 1794, 1795,
  1028, 1029, 1284, 1285, 1030, 1031, 1286, 1287, 1540, 1541, 1796, 1797, 1542, 1543, 1798, 1799,
     8,    9,  264,  265,   10,   11,  266,  267,  520,  521,  776,  777,  522,  523,  778,  779,
    12,   13,  268,  269,   14,   15,  270,  271,  524,  525,  780,  781,  526,  527,  782,  783,
  1032, 1033, 1288, 1289, 1034, 1035, 1290, 1291, 1544, 1545, 1800, 1801, 1546, 1547, 1802, 1803,
  1036, 1037, 1292, 1293, 1038, 1039, 1294, 1295, 1548, 1549, 1804, 1805, 1550, 1551, 1806, 1807,
  2048, 2049, 2304, 2305, 2050, 2051, 2306, 2307, 2560, 2561, 2816, 2817, 2562, 2563, 2818, 2819,
  2052, 2053, 2308, 2309, 2054, 2055, 2310, 2311, 2564, 2565, 2820, 2821, 2566, 2567, 2822, 2823,
  3072, 3073, 3328, 3329, 3074, 3075, 3330, 3331, 3584, 3585, 3840, 3841, 3586, 3587, 3842, 3843,
  3076, 3077, 3332, 3333, 3078, 3079, 3334, 3335, 3588, 3589, 3844, 3845, 3590, 3591, 3846, 3847,
  2056, 2057, 2312, 2313, 2058, 2059, 2314, 2315, 2568, 2569, 2824, 2825, 2570, 2571, 2826, 2827,
  2060, 2061, 2316, 2317, 2062, 2063, 2318, 2319, 2572, 2573, 2828, 2829, 2574, 2575, 2830, 2831,
  3080, 3081, 3336, 3337, 3082, 3083, 3338, 3339, 3592, 3593, 3848, 3849, 3594, 3595, 3850, 3851,
  3084, 3085, 3340, 3341, 3086, 3087, 3342, 3343, 3596, 3597, 3852, 3853, 3598, 3599, 3854, 3855 };


int jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
int jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };


struct STACK 
{
   int64_t pix;
   int order;
};


int hpxConeSearch(int omax, double ra, double dec, double radius, struct hpxRanges **hpxrng)
{
   int    i, o, sdist, zone;

   int64_t pix;

   double z, phi, cangdist, cosrad, sintheta;
   double dr[256], crpdr[256], crmdr[256], dot, dist;

   struct STACK *stk;

   struct hpxRanges *rng;

   int nstack, maxstack;
   int nrange, maxrange;

   Vec pixel, user;

   user.x = cos(ra*DTR) * cos(dec*DTR);
   user.y = sin(ra*DTR) * cos(dec*DTR);
   user.z = sin(dec*DTR);

   radius = radius * DTR;

   cosrad = cos(radius);


   // Set up our stack, used to recurse down through the pixel orders
   // and pick up matching pixels or sets of pixels

   maxstack = MAXSTACK;

   nstack = 0;

   stk = (struct STACK *)malloc(maxstack * sizeof(struct STACK));


   // Set up our range set list 

   maxrange = MAXRANGE;

   nrange = 0;

   rng = (struct hpxRanges *)malloc(maxrange * sizeof(struct hpxRanges));



   // The positional comparison is an approximation.  
   // We use a bounding circle for each pixel with the 
   // radius for the circle based on the maximum pixel center
   // to corner distance for the entire set at a given order.

   for (o=0; o<=omax; ++o)
   {    
      dr[o] = hpxMaxPixRad(o); // safety distance

      if(sptDebug)
      {
         printf("DEBUG> order %2d:  nside = %llu, pixel size = %8.4f deg, R-dr = %8.4f, R+dr = %8.4f\n",
            o, (1ULL << o), dr[o]/DTR, (radius-dr[o])/DTR, (radius+dr[o])/DTR);
         fflush(stdout);
      }

      crpdr[o] = cos(radius+dr[o]);
      crmdr[o] = cos(radius-dr[o]);

      if(radius+dr[o] > pi)
         crpdr[o] = -1.;

      if(radius-dr[o] < 0.)
         crmdr[o] = 1.; 
   }    



   // Start by putting the 12 base pixels in reverse 
   // order onto the stack.  This will trigger recursive
   // processing for each one to depth.

   if(sptDebug)
      printf("\nInitialize stack\n");
   
   for (i=0; i<12; ++i) 
   {
      stk[nstack].pix   = 11-i;
      stk[nstack].order = 0;

      if(sptDebug)
      {
         printf("DEBUG> push stack[%2d]: pix = %" PRIu64 ", order = %d\n",
            nstack, stk[nstack].pix, stk[nstack].order);
         fflush(stdout);
      }

      ++nstack;


      // If necessary, increase the size of the stack

      if(nstack >= maxstack)
      {
         maxstack += MAXSTACK;

         stk = (struct STACK *)realloc(stk, maxstack * sizeof(struct STACK));
      }
   }

   while (nstack > 0) // as long as there are pixels on the stack
   {    
      // pop current pixel number and order from the stack

      pix = stk[nstack-1].pix;
      o   = stk[nstack-1].order;

      if(sptDebug)
      {
         printf("--------\nDEBUG>  POP stack[%2d]: pix = %" PRId64 ", order = %d\n",
            nstack-1, stk[nstack-1].pix, stk[nstack-1].order);
         fflush(stdout);
      }

      --nstack;

      hpxPix2Loc(o, pix, &z, &phi);

      sintheta = sqrt( (1.-z) * (1.+z) );   // Numerically the same as 1-z^2 but with less roundoff.

      if(sptDebug)
      {
         printf("DEBUG> pixel %" PRId64 " / order %d:  z,phi = (%-g, %-g), ra,dec = (%-g,%-g)\n",
            pix, o, z, phi, phi/DTR, 90.-acos(z)/DTR);
         fflush(stdout);
      }

      pixel.x = sintheta * cos(phi);
      pixel.y = sintheta * sin(phi);
      pixel.z = z;


      // cosine of angular distance between pixel center and disk center

      cangdist = user.x * pixel.x
               + user.y * pixel.y
               + user.z * pixel.z;

      if(sptDebug)
      {
         printf("DEBUG> pixel %" PRIu64 " / order %d:\n", pix, o);
         printf("DEBUG>       cangdist = %9.6f (distance = %8.4f)\n", cangdist, acos(cangdist)/DTR);
         printf("DEBUG>       cosrad   = %9.6f (R        = %8.4f)\n", cosrad,   acos(cosrad)/DTR);
         printf("DEBUG>                            (dr       = %8.4f)\n", dr[o]/DTR);
         printf("DEBUG>       crmdr    = %9.6f (R-dr     = %8.4f)\n", crmdr[o], acos(crmdr[o])/DTR);
         printf("DEBUG>       crpdr    = %9.6f (R+dr     = %8.4f)\n", crpdr[o], acos(crpdr[o])/DTR);
         fflush(stdout);
      }



      // Four cases:
      //  
      // zone = 0: pixel lies completely outside the queried shape
      //        1: pixel may overlap with the shape, pixel center is outside
      //        2: pixel center is inside the shape, but maybe not the complete pixel
      //        3: pixel lies completely inside the shape

      zone = 0;

      if(cangdist > crpdr[o])  // Otherwise Zone 0
      {    
         if(cangdist < cosrad)
            zone = 1;           

         else if(cangdist <= crmdr[o])
            zone = 2;

         else
            zone = 3;

         if(sptDebug)
         {
            if     (zone == 1) printf("DEBUG> zone 1: pixel may overlap with the shape, pixel center is outside\n");
            else if(zone == 2) printf("DEBUG> zone 2: pixel center is inside the shape, but maybe not the complete pixel\n");
            else if(zone == 3) printf("DEBUG> zone 3: pixel lies completely inside the shape\n");
            else               printf("DEBUG> zone = %d:  This can't happen.\n", zone);
            fflush(stdout);
         }

         
         // Check the pixel

         if (o < omax)  // We can go deeper (or we keep everthing)
         {    
            if (zone == 3)         // Pixel lies completely inside the shape: keep it all
            {    
               sdist = 2*(omax-o);                            // the "bit-shift distance" between map orders

               rng[nrange].min =  pix    << sdist;            // output all subpixels
               rng[nrange].max = (pix+1) << sdist;

               if(sptDebug)
               {
                  printf("DEBUG> RANGE %d (save pixels: %" PRId64 " to %" PRId64 ")\n", nrange, rng[nrange].min, rng[nrange].max);
                  fflush(stdout);
               }


               // Special check to see if we should combine this with the previous range

               if(nrange > 0 && rng[nrange].min == rng[nrange-1].max + 1)
                  rng[nrange-1].max = rng[nrange].max;
               else
                  ++nrange;


               if(nrange >= maxrange)
               {
                  maxrange += MAXRANGE;

                  rng = (struct hpxRanges *)realloc(rng, maxrange * sizeof(struct hpxRanges));
               }
            }    

            else                                                // (1<=zone<=2)
               for (i=0; i<4; ++i) 
               {
                  stk[nstack].pix   = 4*pix+3-i;                // add children to stack
                  stk[nstack].order = o+1;

                  if(sptDebug)
                  {
                     printf("DEBUG> push stack[%2d]: pix = %" PRIu64 ", order = %d\n",
                        nstack, stk[nstack].pix, stk[nstack].order);
                     fflush(stdout);
                  }

                  ++nstack;

                  if(nstack >= maxstack)
                  {
                     maxstack += MAXSTACK;

                     stk = (struct STACK *)realloc(stk, maxstack * sizeof(struct STACK));
                  }
               }
         }    

         else                     // This is the bottom level; keep the pixel
         {
            rng[nrange].min = pix;   // A "range" of one 
            rng[nrange].max = pix;

            if(sptDebug)
            {
               printf("DEBUG> RANGE %d (save pixel: %" PRId64 ")\n", nrange, pix);
               fflush(stdout);
            }


            // Special check to see if we should combine this with the previous range

            if(nrange > 0 && rng[nrange].min == rng[nrange-1].max + 1)
               rng[nrange-1].max = rng[nrange].max;
            else
               ++nrange;


            if(nrange >= maxrange)
            {
               maxrange += MAXRANGE;

               rng = (struct hpxRanges *)realloc(rng, maxrange * sizeof(struct hpxRanges));
            }
         }
      }
   }


   if(sptDebug)
   {
      printf("\n");
      printf("Center: (%8.4f,%8.4f)  Radius: %-g\n\n", ra, dec, radius/DTR);
   }

   if(sptDebug)
   {
      for(i=0; i<nrange; ++i)
      {
         printf("RANGE> %d: %" PRIu64 " - %" PRIu64 "\n", i, rng[i].min, rng[i].max);

         for(pix=rng[i].min; pix<=rng[i].max; ++pix)
         {
            hpxPix2Loc (omax, pix, &z, &phi);

            sintheta = sqrt( (1.-z) * (1.+z) );

            pixel.x = sintheta * cos(phi);
            pixel.y = sintheta * sin(phi);
            pixel.z = z;

            dot = user.x * pixel.x
                + user.y * pixel.y
                + user.z * pixel.z;

            dist = acos(dot);
            
            printf("       %" PRId64 ": z,phi (%8.5f %8.5f) -> ra,dec (%8.4f, %8.4f) [pixel size: %-g, distance: %-g]\n", 
                  pix, z, phi, phi/DTR, 90.-acos(z)/DTR, dr[omax]/DTR, dist/DTR);

            fflush(stdout);

         }
      }
   }

   *hpxrng = rng;

   free(stk);

   return nrange;
}


int hpxPolygonSearch(int omax, int npoly, double *ra, double *dec, struct hpxRanges **hpxrng)
{
   int        i, flip, nreverse, o, sdist, zone, izone;
   int64_t    pix;

   double     z, phi, cosrad, rcenter, sintheta;
   double     handedness, crad, rdist;
   double     dr[256];

   double ***crlimit;

   struct STACK *stk;

   struct hpxRanges *rng;

   int nstack, maxstack;
   int nrange, maxrange;

   Vec pixel, *poly, *normal;
   Vec center;


   // Calculate the 3-vectors for the polygon vertices

   poly = (Vec *)malloc(npoly * sizeof(Vec));

   for(i=0; i<npoly; ++i)
   {
      poly[i].x = cos(ra[i]*DTR) * cos(dec[i]*DTR);
      poly[i].y = sin(ra[i]*DTR) * cos(dec[i]*DTR);
      poly[i].z = sin(dec[i]*DTR);
   }


   // Find the normal vectors for each of the polygon sides
   // Correct for handedness of the polygon

   normal = (Vec *)malloc(npoly * sizeof(Vec));

   flip = 1;

   nreverse = 0;

   if(sptDebug)
      printf("\nIn hpxPolygonSearch()\n");
   
   for(i=0; i<npoly; ++i)
   {
      vCross(&poly[i], &poly[(i+1)%npoly], &normal[i]);
      vNormalize(&normal[i]);

      handedness = vDot(&normal[i], &poly[(i+2)%npoly]);

      if(sptDebug)
      {
         printf("Normal (x,y,z)  = (%11.8f,%11.8f,%11.8f) [handedness: %11.8f]\n",
            normal[i].x, normal[i].y, normal[i].z, handedness);
      }

      if(fabs(handedness) < 1.e-10)
         return 1;

      if(i == 0 && handedness < 0.)
         flip = -1.;
      
      if(i > 0 && flip*handedness <= 0.)
         return 1;

      if(flip < 0)
      {
         if(sptDebug)
            printf("Reversing normal %d\n", i);

         vReverse(&normal[i]);

         ++nreverse;
      }
   }

   if(sptDebug && nreverse > 0)
      printf("\n%d reversed\n\n", nreverse);


   // Find the minimum bounding circle around the polygon
   // (returns cos(radius) rather than radius because it is
   //  more useful)

   hpxBoundingCircle(npoly, poly, &center, &cosrad);

   rcenter = acos(cosrad);

   if(sptDebug)
   {
      printf("Bounding circle: (%-g, %-g) [%-g]\n", 
         atan2(center.y, center.x)/DTR, asin(center.z)/DTR, acos(cosrad)/DTR);
      fflush(stdout);
   }


   // We determine whether a HEALPix pixel is 
   //
   //    Completely outside the polygon,
   //    Completely inside the polygon, or
   //    Indeterminant (might be overlapping)
   //
   // based on where its center is with regard to the N polygon 
   // side normals.  A point on the sphere that is on one of the
   // edges will be exactly 90 degrees away from that side's 
   // normal or equivalently the dot product between the point
   // location and that normal will be exactly 0.
   //
   // If a point is dr degrees outside that edge, the dot product
   // will be cos(90-dr) (which will be greater than zero) and
   // if dr degrees inside that edge it will be cos(90+dr)
   // (which will be less than zero).
   //
   // For this algorithm, we use the maximum pixel bounding radius
   // at each HEALPix level as an upper limit for the above dr.  
   // This lets us say with certainty that when the dot product of
   // the pixel center and one of the polygon edges is less than
   // than cos(90+dr), that pixel is more than dr degrees inside
   // the polygon (as far as that edge is concerned).  When the
   // dot product is less than cos(90-dr), it is more than dr
   // degrees outside.  In between and we can't be sure and have
   // to drill down, performing the sam kind of check for the
   // four sub-pixels at the next level.
   //
   // Later we will loop over the edges to cumulative effect is,
   // but here we just calculate the cos(90-dr) and cos(90+dr
   // values for all pixel levels with appropriate failovers
   // for > 90, < 0 etc.

   crlimit = (double ***)malloc((omax+1) * sizeof(double **));

   for (o=0; o<=omax; ++o)
   {    
      crlimit[o] = (double **)malloc(npoly * sizeof(double *));

      dr[o] = hpxMaxPixRad(o); // safety distance

      for (i=0; i<npoly; ++i)
      {
         crlimit[o][i] = (double *)malloc(3 * sizeof(double));

         crlimit[o][i][0] = (halfpi+dr[o] > pi) ? -1. : cos(halfpi+dr[o]);
         crlimit[o][i][1] = (o == 0)            ?  0. : crlimit[0][i][1];
         crlimit[o][i][2] = (halfpi-dr[o] < 0.) ?  1. : cos(halfpi-dr[o]);
      }
   }


   // Set up our stack, used to recurse down through the pixel orders
   // and pick up matching pixels or sets of pixels

   maxstack = MAXSTACK;

   nstack = 0;

   stk = (struct STACK *)malloc(maxstack * sizeof(struct STACK));


   // Set up our range set list 

   maxrange = MAXRANGE;

   nrange = 0;

   rng = (struct hpxRanges *)malloc(maxrange * sizeof(struct hpxRanges));


   // Start by putting the 12 base pixels in reverse 
   // order onto the stack.  This will trigger recursive
   // processing for each one to depth.

   if(sptDebug)
      printf("\nInitialize stack\n");
   
   for (i=0; i<12; ++i) 
   {
      stk[nstack].pix   = 11-i;
      stk[nstack].order = 0;

      if(sptDebug)
      {
         printf("DEBUG> push stack[%2d]: pix = %" PRIu64 ", order = %d\n",
            nstack, stk[nstack].pix, stk[nstack].order);
         fflush(stdout);
      }

      ++nstack;


      // If necessary, increase the size of the stack

      if(nstack >= maxstack)
      {
         maxstack += MAXSTACK;

         stk = (struct STACK *)realloc(stk, maxstack * sizeof(struct STACK));
      }
   }

   while (nstack > 0) // as long as there are pixels on the stack
   {    
      // pop current pixel number and order from the stack

      pix = stk[nstack-1].pix;
      o   = stk[nstack-1].order;

      if(sptDebug)
      {
         printf("--------\nDEBUG>  POP stack[%2d]: pix = %" PRId64 ", order = %d\n",
            nstack-1, stk[nstack-1].pix, stk[nstack-1].order);

         if(nstack == 1)
            printf("DEBUG>  This is the end of the stack; we exit after this\n");

         fflush(stdout);
      }

      --nstack;

      hpxPix2Loc(o, pix, &z, &phi);

      sintheta = sqrt( (1.-z) * (1.+z) );   // Numerically the same as 1-z^2 but with less roundoff.

      if(sptDebug)
      {
         printf("DEBUG> pixel %" PRId64 " / order %d:  z,phi = (%-g, %-g), ra,dec = (%-g,%-g)\n",
            pix, o, z, phi, phi/DTR, 90.-acos(z)/DTR);
         fflush(stdout);
      }

      pixel.x = sintheta * cos(phi);
      pixel.y = sintheta * sin(phi);
      pixel.z = z;


      // As we examine each pixel from the stack, we start by assuming the pixel
      // is completely inside the polygon.  We check each polygon edge in turn
      // to see if we need to amend this assumption.  
      // 
      // If, based on our "crlimit" conditions above, we determine the pixel is
      // in the indeterminant regime for any side, that becomes the new best guess
      // state for the pixel.  If it is completely outside for any side, it is
      // outside and we can abandon it.
      // 
      // If, by the time we are done checking sides, it remains fully inside,
      // we can add all lowest-level subpixels of it to our set of pixel ranges
      // and stop checking.
      // 
      // If we reach the finest resolution with the pixel still being indeterminant,
      // we add it to the list to be safe.
      //
      // Our four cases:
      //  
      // zone = 0: pixel lies completely outside the queried shape
      //        1: pixel may overlap with the shape, pixel center is outside
      //        2: pixel center is inside the shape, but maybe not the complete pixel
      //        3: pixel lies completely inside the shape

      zone = 3;  // Assume the pixel in entirely in the polygon

      for(i=0; i<npoly; ++i)
      {
         crad = vDot(&pixel, &normal[i]);

         rdist = acos(vDot(&pixel, &center));

         if(sptDebug > 1)
         {
            printf("DEBUG> Polygon side %d of %d\n", i+1, npoly);

            printf("DEBUG> Checking side of sky: distance %-g vs. region size %-g + pixel size %-g\n",
               rdist/DTR, rcenter/DTR, dr[o]/DTR);

            fflush(stdout);
         }

         if(rdist > rcenter + dr[o])
         {
            if(sptDebug > 1)
            {
               printf("DEBUG> rdist - (rcenter + dr[o]) = %-g  =>  zone 0\n",
                  rdist/DTR - (rcenter/DTR + dr[o]/DTR));
               fflush(stdout);
            }

            zone = 0;
         }

         for(izone=0; izone<zone; ++izone)
         {
            if(sptDebug > 1)
            {
               if(izone == 0)
                  printf("DEBUG> Checking limit (outer pad):      %-g vs. %-g\n",
                     acos(crad)/DTR-90., acos(crlimit[o][i][izone])/DTR-90.);
               if(izone == 1)
                  printf("DEBUG> Checking limit (polygon exact):  %-g vs. %-g\n",
                     acos(crad)/DTR-90., acos(crlimit[o][i][izone])/DTR-90.);
               if(izone == 2)
                  printf("DEBUG> Checking limit (inner pad):      %-g vs. %-g\n",
                     acos(crad)/DTR-90., acos(crlimit[o][i][izone])/DTR-90.);
               fflush(stdout);
            }

            if(crad < crlimit[o][i][izone])  // If this side says to demote, do so
            {
               if(sptDebug > 1)
               {
                  if     (izone == 1) printf("DEBUG> zone -> 1: pixel may overlap with the shape, pixel center is outside\n");
                  else if(izone == 2) printf("DEBUG> zone -> 2: pixel center is inside the shape, but maybe not the complete pixel\n");
                  else if(izone == 3) printf("DEBUG> zone -> 3: pixel lies completely inside the shape\n");
                  fflush(stdout);
               }

               zone = izone;

               if(zone == 0)                // And if we are completely outside, get out
                  break;
            }
         }

         if(zone == 0)  // If after processing this side we are out, stop
            break;
      }

      if(sptDebug > 1 && zone == 0)
      {
         printf("DEBUG> zone 0: pixel completely outside shape\n");
         fflush(stdout);
      }

      if(zone > 0)  // We need to further process this pixel
      {    
         if(sptDebug)
         {
            if     (zone == 1) printf("DEBUG> zone 1: pixel may overlap with the shape, pixel center is outside\n");
            else if(zone == 2) printf("DEBUG> zone 2: pixel center is inside the shape, but maybe not the complete pixel\n");
            else if(zone == 3) printf("DEBUG> zone 3: pixel lies completely inside the shape\n");
            else               printf("DEBUG> zone = %d:  This can't happen.\n", zone);
            fflush(stdout);
         }
         

         // Check the pixel

         if (o < omax)  // We can go deeper (or we keep everthing)
         {    
            if (zone == 3)         // Pixel lies completely inside the shape: keep it all
            {    
               sdist = 2*(omax-o);                            // the "bit-shift distance" between map orders

               rng[nrange].min =  pix    << sdist;            // output all subpixels
               rng[nrange].max = (pix+1) << sdist;

               if(sptDebug)
               {
                  printf("DEBUG> RANGE %d (save pixels: %" PRId64 " to %" PRId64 ")\n", nrange, rng[nrange].min, rng[nrange].max);
                  fflush(stdout);
               }


               // Special check to see if we should combine this with the previous range

               if(nrange > 0 && rng[nrange].min == rng[nrange-1].max + 1)
                  rng[nrange-1].max = rng[nrange].max;
               else
                  ++nrange;


               if(nrange >= maxrange)
               {
                  maxrange += MAXRANGE;

                  rng = (struct hpxRanges *)realloc(rng, maxrange * sizeof(struct hpxRanges));
               }
            }    

            else                                                // (1<=zone<=2)
               for (i=0; i<4; ++i) 
               {
                  stk[nstack].pix   = 4*pix+3-i;                // add children to stack
                  stk[nstack].order = o+1;

                  if(sptDebug)
                  {
                     printf("DEBUG> push stack[%2d]: pix = %" PRIu64 ", order = %d\n",
                        nstack, stk[nstack].pix, stk[nstack].order);
                     fflush(stdout);
                  }

                  ++nstack;

                  if(nstack >= maxstack)
                  {
                     maxstack += MAXSTACK;

                     stk = (struct STACK *)realloc(stk, maxstack * sizeof(struct STACK));
                  }
               }
         }    

         else                     // This is the bottom level; keep the pixel
         {
            rng[nrange].min = pix;   // A "range" of one 
            rng[nrange].max = pix;

            if(sptDebug)
            {
               printf("DEBUG> RANGE %d (save pixel: %" PRId64 ")\n", nrange, pix);
               fflush(stdout);
            }


            // Special check to see if we should combine this with the previous range

            if(nrange > 0 && rng[nrange].min == rng[nrange-1].max + 1)
               rng[nrange-1].max = rng[nrange].max;
            else
               ++nrange;


            if(nrange >= maxrange)
            {
               maxrange += MAXRANGE;

               rng = (struct hpxRanges *)realloc(rng, maxrange * sizeof(struct hpxRanges));
            }
         }
      }
   }

   *hpxrng = rng;


   for (o=0; o<=omax; ++o)
   {
      for (i=0; i<npoly; ++i)
          free(crlimit[o][i]);

      free(crlimit[o]);
   }

   free(crlimit);
   free(stk);
   free(poly);
   free(normal);

   if(sptDebug)
   {
      printf("\nDEBUG>  Done with stack, returning %d ranges.\n\n", nrange);
      fflush(stdout);
   }

   return nrange;
}


// For a given order, the most "distorted" pixels are the ones
// along the transition latitude between equatorial and polar
// regions (stated without proof; see the geometry discussion
// in the HEALPix documentation).
//
// In particular, the pixel corner north of the pixel center
// (south in the southern hemisphere) gives us the max center
// to corner distance for the order.  We will use this as 
// a bounding circle for purposes of overlap calculation with
// user-supplied regions in the main code above.

double hpxMaxPixRad(int order)
{
   int64_t nside;

   double z, phi, sintheta, tmp, dot, dist;

   Vec center, corner;


   // The number of pixels along the side of a face
   // is a power of two of the order

   nside = 1ULL << order;

   if(sptDebug > 1)
   {
      printf("\nDEBUG> Pixel radius for order %d\n", order);
      printf("DEBUG> nside = %" PRId64 "\n", nside);
      fflush(stdout);
   }


   // The first vector is the center of one of the
   // transition pixels; specifically the one whose
   // edge runs along longitude zero.  
   
   // Being a transition pixel, it has a colatitude
   // value of 2/3.

   // There are 4*nside pixels total circling that
   // latitude, all the same.  So each one spans
   // a longitude angle of (2*pi)/(4*nside) and the
   // center of the first one is therefore at longitude
   // pi/(4*nside).

   z = 2./3.;

   sintheta = sqrt( (1.-z) * (1.+z) );

   phi = pi / (4.*nside);

   if(sptDebug > 1)
   {
      printf("DEBUG> center z,phi = %-g, %-g (%-g deg)\n", z, phi, phi/DTR);
      fflush(stdout);
   }

   center.x = sintheta * cos(phi);
   center.y = sintheta * sin(phi);
   center.z = z;


   // The corner of that pixel that is directly north
   // of the center is by definition at longitude zero.
   
   // For it's (co)latitude, we again refer to the 
   // HEALPix geometry discussion, where it can be
   // shown to be given a quadratic in the pixel 
   // side length

   tmp = 1. - 1./nside;

   tmp =  tmp*tmp;

   z = 1. - tmp/3.;

   phi = 0.;

   if(sptDebug > 1)
   {
      printf("DEBUG> corner z,phi = %-g, %-g (%-g deg)\n", z, phi, phi/DTR);
      fflush(stdout);
   }

   sintheta = sqrt( (1.-z) * (1.+z) );

   corner.x = sintheta * cos(phi);
   corner.y = sintheta * sin(phi);
   corner.z = z;


   // Given these two sky vectors, we can get the distance
   // between them from their dot product

   dot = center.x * corner.x
       + center.y * corner.y
       + center.z * corner.z;

   dist = acos(dot);
   
   if(sptDebug > 1)
   {
      printf("DEBUG> pixel size =  %-g (%-g deg)\n\n", dist, dist/DTR);
      fflush(stdout);
   }

   return dist;
}



// Convert a pixel number to z,phi coordinates

void hpxPix2Loc (int order, int64_t pixin, double *z, double *phi)
{
   int64_t nside, face_num, npface, npix;

   int64_t ix, iy, jr, nr, pix;
   int64_t itmp;

   double fact1, fact2;
   double tmp;

   nside  = 1ULL << order;
   npface = nside << order;
   npix   = 12 * npface;
   fact2  = 4./npix;
   fact1  = (nside << 1) * fact2;

   pix = pixin;

   face_num = pix >> (2*order);

   pix &= (npface-1);

   ix = hpxCompressBits(pix);
   iy = hpxCompressBits(pix>>1);

   jr = (jrll[face_num] << order) - ix - iy - 1; 

   if (jr < nside)
   {    
      nr = jr;

      tmp=(nr*nr)*fact2;

      *z = 1 - tmp; 
   }    
   else if (jr > 3*nside)
   {    
      nr = 4*nside - jr;

      tmp=(nr*nr) * fact2;

      *z = tmp - 1; 
   }    
   else 
   {    
      nr = nside;

      *z = (2*nside-jr) * fact1;
   }    

   itmp = jpll[face_num] * nr + ix - iy;
   
   if (itmp<0)
      itmp += 8 * nr;

   if(nr == nside)
      *phi = 0.75 * halfpi * itmp * fact1;
   else
      *phi = (0.5 * halfpi * itmp) / nr;

   if(sptDebug > 1)
   {
      printf("DEBUG> order    = %d\n", order);
      printf("DEBUG> pix      = %" PRId64 "\n", pix);
      printf("       -------------\n");
      printf("DEBUG> nside    = %" PRId64 "\n", nside);
      printf("DEBUG> npface   = %" PRId64 "\n", npface);
      printf("DEBUG> npix     = %" PRId64 "\n", npix);
      printf("DEBUG> fact2    = %-g\n", fact2);
      printf("DEBUG> fact1    = %-g\n", fact1);
      printf("DEBUG> face_num = %" PRId64 "\n", face_num);
      printf("DEBUG> pix     -> %" PRId64 "\n", pix);
      printf("DEBUG> ix       = %" PRId64 "\n", ix);
      printf("DEBUG> iy       = %" PRId64 "\n", iy);
      printf("DEBUG> jr       = %" PRId64 "\n", jr);
      printf("DEBUG> nr       = %" PRId64 "\n", nr);
      printf("DEBUG> tmp      = %0g\n", tmp);
      printf("DEBUG> itmp     = %" PRId64 "\n", itmp);
      printf("DEBUG> z       -> %-g\n", *z);
      printf("DEBUG> phi     -> %-g\n", *phi);
      fflush(stdout);
   }
}


// These two routines fine a minimum bounding circle 
// around a set of vertices.  

void hpxBoundingCircle (int np, Vec *point, Vec *center, double *cosrad)
{
   int i;

   vMidpoint(&point[0], &point[1], center);

   *cosrad = vDot(&point[0], center);

   for (i=2; i<np; ++i)
      if (vDot(&point[i], center) < *cosrad)
         hpxGetCircle(np, point, i, center, cosrad);
}


void hpxGetCircle (int np, Vec *point, int q, Vec *center, double *cosrad)
{
   int i;

   vMidpoint(&point[0], &point[1], center);

   *cosrad = vDot(&point[0], center);

   for (i=1; i<q; ++i)
      if (vDot(&point[i], center) < *cosrad)
         hpxGetCircle(np, point, i, center, cosrad);
}



int hpxCompressBits (int64_t v)
{
   int64_t raw = v&0x5555555555555555ull;

   raw|=raw>>15;

   return ctab[ raw     &0xff]      | (ctab[(raw>> 8)&0xff]<< 4)
       | (ctab[(raw>>32)&0xff]<<16) | (ctab[(raw>>40)&0xff]<<20);
}




/***************************************************/
/*                                                 */
/* vCopy()                                         */
/*                                                 */
/* Copy the contents of one vector to another      */
/*                                                 */
/***************************************************/

void vCopy(Vec *v, Vec *c)
{
   c->lon = v->lon;
   c->lat = v->lat;

   c->x = v->x;
   c->y = v->y;
   c->z = v->z;

   return;
}



/***************************************************/
/*                                                 */
/* vCalcRADec()                                    */
/*                                                 */
/* Update vector with RA and Dec based on x,y,z    */
/*                                                 */
/***************************************************/

void vCalcRADec(Vec *v)
{
   v->lon = atan2(v->y, v->x)/DTR;
   v->lat = asin(v->z)/DTR;

   while(v->lon >= 360.) v->lon -= 360.;
   while(v->lon <    0.) v->lon += 360.;

   return;
}



/***************************************************/
/*                                                 */
/* vCalcXYZ()                                      */
/*                                                 */
/* Update vector with x,y,z based on RA,Dec        */
/*                                                 */
/***************************************************/

void vCalcXYZ(Vec *v)
{
   v->x = cos(v->lat * DTR) * cos(v->lon * DTR);
   v->y = cos(v->lat * DTR) * sin(v->lon * DTR);
   v->z = sin(v->lat * DTR);

   return;
}



/***************************************************/
/*                                                 */
/* vMidpoint()                                     */
/*                                                 */
/* Finds the midpoint between two points on the    */
/* sky.                                            */
/*                                                 */
/***************************************************/

void vMidpoint(Vec *a, Vec *b, Vec *c)
{
   c->x = a->x + b->x;
   c->y = a->y + b->y;
   c->z = a->z + b->z;

   vNormalize(c);
   
   return;
}



/***************************************************/
/*                                                 */
/* vPixCenter()                                    */
/*                                                 */
/* Finds the center of a a pixel (four corners)    */
/* on the sky.                                     */
/*                                                 */
/***************************************************/

void vPixCenter(Vec *a, Vec *b, Vec *c, Vec *d, Vec *v)
{
   v->x = a->x + b->x + c->x + d->x;
   v->y = a->y + b->y + c->y + d->y;
   v->z = a->z + b->z + c->z + d->z;

   vNormalize(v);
   
   return;
}



/***************************************************/
/*                                                 */
/* vCross()                                        */
/*                                                 */
/* Vector cross product.                           */
/*                                                 */
/***************************************************/

int vCross(Vec *v1, Vec *v2, Vec *v3)
{
   v3->x =  v1->y*v2->z - v2->y*v1->z;
   v3->y = -v1->x*v2->z + v2->x*v1->z;
   v3->z =  v1->x*v2->y - v2->x*v1->y;

   if(v3->x == 0.
   && v3->y == 0.
   && v3->z == 0.)
      return 0;
   
   return 1;
}


/***************************************************/
/*                                                 */
/* vDot()                                          */
/*                                                 */
/* Vector dot product.                             */
/*                                                 */
/***************************************************/

double vDot(Vec *a, Vec *b)
{
   double sum;

   sum = a->x * b->x
       + a->y * b->y
       + a->z * b->z;

   return sum;
}


/***************************************************/
/*                                                 */
/* vNormalize()                                    */
/*                                                 */
/* Normalize the vector                            */
/*                                                 */
/***************************************************/

double vNormalize(Vec *v)
{
   double len;

   len = 0.;

   len = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);

   if(len == 0.)
      len = 1.;

   v->x = v->x / len;
   v->y = v->y / len;
   v->z = v->z / len;
   
   return len;
}


/***************************************************/
/*                                                 */
/* vReverse()                                      */
/*                                                 */
/* Reverse the vector direction                    */
/*                                                 */
/***************************************************/

void vReverse(Vec *v)
{
   v->x = -1. * v->x;
   v->y = -1. * v->y;
   v->z = -1. * v->z;
   
   return;
}


/***************************************************/
/*                                                 */
/* vPrint()                                        */
/*                                                 */
/* Print out vector (for debugging)                */
/*                                                 */
/***************************************************/

void vPrint(Vec *v, char *label)
{
   vCalcRADec(v);

   printf("VECTOR> %9.6f %9.6f %9.6f  ->  %10.6f %10.6f (%s)\n",
    v->x, v->y, v->z, v->lon, v->lat, label);
   fflush(stdout);

   return;
}
