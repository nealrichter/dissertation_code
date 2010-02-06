/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * RICE 4.0x              Copyright (C) 1993             Rene' Jager *
 *                                                                   *
 *                                                                   *
 * This toolbox is free software; you can redistribute it and/or     *
 * modify it under the terms of the GNU General Public License as    *
 * published by the Free Software Foundation; either version 2 of    *
 * the License, or (at your option) any later version.               *
 *                                                                   *
 * This toolbox is distributed in the hope that it will be useful,   *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of    *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      *
 * GNU General Public License for more details.                      *
 *                                                                   *
 * You should have received a copy of the GNU General Public License *
 * along with this toolbox; if not, write to the:                    *
 *                                                                   *
 *    Free Software Foundation, Inc.                                 *
 *    675 Mass Ave, Cambridge                                        *
 *    MA 02139, USA.                                                 *
 *                                                                   *
 * See the RICE documentation for more information on the toolbox.   *
 * The file COPYING for the complete GNU General Public License.     *
 *                                                                   *
 * You can reach me by (preferably e-mail):                          *
 *                                                                   *
 *    Rene' Jager                                                    *
 *                                                                   *
 *    Delft University of Technology                                 *
 *    Department of Electrical Engineering                           *
 *    Control Laboratory                                             *
 *                                                                   *
 *    Room ET 12.06                                                  *
 *                                                                   *
 *    Mekelweg 4                                                     *
 *    P.O.Box 5031                                                   *
 *    2600 GA  Delft                                                 *
 *    The Netherlands                                                *
 *                                                                   *
 *    e-mail: R.Jager@ET.TUDelft.NL                                  *
 *    phone:  +31-15-78 51 14                                        *
 *    fax:    +31-15-62 67 38                                        *
 *    telex:  38151 butud nl                                         *
 *                                                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
   File:   fuzzy.c
   Author: Rene' Jager
   Update: March 1, 1992
   Info:   utility source file, see file fuzzy.doc
*/


/* header file */

#include "fuzzy.h"


/* macro's */

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(a)   ((a) * (a))


/* smooth membership function(s) */

float fpi(float val, float p1, float p2, float p3, float p4)
{
   if(p2 <= val && val <= p3)
      return 1.0;
   if(val <= p1 || p4 <= val)
      return 0.0;
   if(p1 < val && val <= (p1 + p2)/2.0)
      return 2.0*SQR(val - p1)/SQR(p2 - p1);
   if((p1 + p2)/2.0 <= val && val < p2)
      return 1.0 - 2.0*SQR(val - p2)/SQR(p1 - p2);
   if(p3 < val && val <= (p3 + p4)/2.0)
      return 1.0 - 2.0*SQR(val - p3)/SQR(p3 - p4);
   if((p3 + p4)/2.0 <= val && val < p4)
      return 2.0*SQR(val - p4)/SQR(p3 - p4);

   return 0.0;
}


/* hard membership function(s) */

float ftrapezium(float val, float p1, float p2, float p3, float p4)
{
   if(p2 <= val && val <= p3)
      return 1.0;
   if(val <= p1 || p4 <= val)
      return 0.0;
   if(p1 < val && val < p2)
      return (val - p1)/(p2 - p1);
   if(p3 < val && val < p4)
      return (p4 - val)/(p4 - p3);

   return 0.0;
}


/* indexed-centre-of-gravity */

float ficog(int len, float *set, float lim)
{
   register int i;
   float up = 0.0, low = 0.0;

   for(i = 1; i <= len; set++, i++)
      if(*set >= lim)
      {
         up += *set*i;
         low += *set;
      }

   up -= 1.0;

   if(low == 0.0)
      return -1.0;

   return up/low;
}


/* centre-of-gravity */

float fcog(int len, float *set)
{
   register int i;
   float up = 0.0, low = 0.0;

   for(i = 1; i <= len; set++, i++)
   {
      up += *set*i;
      low += *set;
   }

   up -= 1.0;

   if(low == 0.0)
      return -1.0;

   return up/low;
}


/* mean-of-maximum */

float fmom(int len, float *set)
{
   int low = 0;
   float up = 0.0, top;

   top = fhgt(len, set);
   set += len;

   while(len)
   {
      if(*set == top)
      {
         up += len;
         low++;
      }
      set--;
      len--;
   }

   up -= 1.0;

   if(low == 0)
      return -1.0;

   return up/low;
}


/* general norm operators */

float *fnorm(int len, float *dest, float *src, float (*norm)(float, float))
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = (*norm)(*dest, *src);
         dest++;
      }
   else
      while(len--)
      {
         *dest = (*norm)(*dest, *src);
         dest++;
         src++;
      }

   return save;
}


float *fxnorm(int len, float *dest, float *src, 
              float (*norm)(float, float, float), float par)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = (*norm)(*dest, *src, par);
         dest++;
      }
   else
      while(len--)
      {
         *dest = (*norm)(*dest, *src, par);
         dest++;
         src++;
      }

   return save;
}


/* intersection operators */

float *fzand(int len, float *dest, float *src)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = MIN(*dest, *src);
         dest++;
      }
   else
      while(len--)
      {
         *dest = MIN(*dest, *src);
         dest++;
         src++;
      }

   return save;
}


float *fpand(int len, float *dest, float *src)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = *dest * *src;
         dest++;
      }
   else        
      while(len--)
      {
         *dest = *dest * *src;
         dest++;
         src++;
      }

   return save;
}


float *fland(int len, float *dest, float *src)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = MAX(*dest + *src - 1.0, 0.0);
         dest++;
      }
   else
      while(len--)
      {
         *dest = MAX(*dest + *src - 1.0, 0.0);
         dest++;
         src++;
      }

   return save;
}


/* union operators */

float *fzor(int len, float *dest, float *src)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = MAX(*dest, *src);
         dest++;
      }
   else
      while(len--)
      {
         *dest = MAX(*dest, *src);
         dest++;
         src++;
      }

   return save;
}


float *fpor(int len, float *dest, float *src)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = *dest + *src - *dest * *src;
         dest++;
      }
   else
      while(len--)
      {
         *dest = *dest + *src - *dest * *src;
         dest++;
         src++;
      }

   return save;
}


float *flor(int len, float *dest, float *src)
{
   float *save = dest;

   if(len < 0)
      while(len++)
      {
         *dest = MIN(*dest + *src, 1.0);
         dest++;
      }
   else
      while(len--)
      {
         *dest = MIN(*dest + *src, 1.0);
         dest++;
         src++;
      }

   return save;
}


/* negation operator */

float *fnot(int len, float *dest)
{
   float *save = dest;

   while(len--)
   {
      *dest = 1.0 - *dest;
      dest++;
   }

   return save;
}


/* height */

float fhgt(int len, float *set)
{
   float value = 0.0;

   while(len--)
   {
      value = MAX(value, *set);
      set++;
   }

   return value;
}


/* alpha-cut and strong alpha-cut */

float *fcut(int len, float *dest, float *src, float cut)
{
   float *save = dest;

   while(len--)
   {
      *dest = (*src >= cut) ? *src : 0.0;
      dest++;
      src++;
   }

   return save;
}


float *fscut(int len, float *dest, float *src, float cut)
{
   float *save = dest;

   while(len--)
   {
      *dest = (*src > cut) ? *src : 0.0;
      dest++;
      src++;
   }

   return save;
}



#ifdef TEST_FUZZY

#include <stdio.h>

int main()
{
   float val;
   val = fpi(0.1, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   val = fpi(1.1, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   val = fpi(1.6, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   val = fpi(2.1, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   val = fpi(3.1, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   val = fpi(3.6, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   val = fpi(4.1, 1.0, 2.0, 3.0, 4.0);
   printf("val = %f\n", val);
   return 0;
}

#endif
