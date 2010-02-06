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
   File:   fuzzy.h
   Author: Rene' Jager
   Update: November 16, 1992
   Info:   utility include file, see file fuzzy.doc
*/
           

/* prevent recursive inclusion */

#ifndef _FUZZY_H_
#define _FUZZY_H_

#define FALSE   0
#define TRUE    1

/* smooth membership functions */

float fpi(float val, float p1, float p2, float p3, float p4);

#define fdelta(val, p1, p2, p3)    fpi(val, p1, p2, p2, p3)
#define fgamma(val, p1, p2)        fpi(val, p1, p2, p2, p2)
#define flambda(val, p1, p2)       fpi(val, p1, p1, p1, p2)


/* straight membership functions */

float ftrapezium(float val, float p1, float p2, float p3, float p4);

#define fblock(val, p1, p2)          ftrapezium(val, p1, p1, p2, p2)
#define ftriangle(val, p1, p2, p3)   ftrapezium(val, p1, p2, p2, p3)


/* defuzzification methods */

float ficog(int len, float *set, float lim);
float fcog(int len, float *set);
float fmom(int len, float *set);


/* general norm operators */

float *fnorm(int len, float *dest, float *src, float (*norm)(float, float));
float *fxnorm(int len, float *dest, float *src,
              float (*norm)(float, float, float), float par);


/* intersection operators */

float *fzand(int len, float *dest, float *src);
float *fland(int len, float *dest, float *src);
float *fpand(int len, float *dest, float *src);


/* union operators */

float *fzor(int len, float *dest, float *src);
float *flor(int len, float *dest, float *src);
float *fpor(int len, float *dest, float *src);


/* negation operator */

float *fnot(int len, float *dest);


/* height */

float fhgt(int len, float *src);


/* alpha-cut and strong alpha-cut */

float *fcut(int len, float *dest, float *src, float cut);
float *fscut(int len, float *dest, float *src, float cut);


#endif /* _FUZZY_H_ */

                                               
