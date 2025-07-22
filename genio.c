/**************************************************************************

   This module contains FORTRAN-callable routines to read and write
   a PC-GENPLOT binary data file:

   iret = getbin( filename, xarray, yarray, npoints, comment, ncmts, maxpts )

   call putbin( filename, xarray, yarray, npoints, comment, ncmts )	


        filename  character variable containing name of datafile
	xarray	  real*8 array to receive x data (first array in file)
	y array   real*8 array to receive y data (second array in file)
	npoints   number of points found in x,y arrays (returned by getbin)
                  (note: if greater than maxpts, only maxpts points are 
		  actually returned in the x,y arrays)
	iret      returns 0 for successful read, 1 for error
	maxpts     maximum number of points to be returned by getbin

   Also included are the following utility subroutines:

   getint()        reads an integer from the input file
   getreal()       reads a floating point real number (converted from  
                   GENPLOT format to standard floating point format)
   getrecsize()    reads a Ryan-McFarlane unformatted file record marker
   getlinerec()    reads an ASCII line from the input file

***************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "fortrancall.h"

#define ERROR 1
#define OK 0
#define MAXCMTLNS 16
#define MAXLINLTH 80
#define MAXELEMENTS 1024      /* Max number of elements in an array record */

extern int errno;

FILE *inputfile, *outputfile;
double factor, offset;
char   info[81];

/*********************************************************************
  getint()        read an integer from the input file
  putint( iout )  write an integer to the output file
**********************************************************************/

unsigned int getint()
{
  int itmp;
  itmp =  fgetc(inputfile);
  itmp += (fgetc(inputfile) << 8);
  return itmp;
}

void putint( iout )
int iout;
{
  fputc( (char) (iout & 0xff), outputfile );
  fputc( (char) (iout >> 8),   outputfile );
}


/*************************************************************************
getrecsize()       Read a Ryan-McFarlane FORTRAN unformatted record 
                      marker from the input file 
putrecsize( size ) Write a record marker to the output file

Note that record size includes the two four-byte record markers
at either end of the record
***************************************************************************/

union i4_datum {
       long int number;
       unsigned int word[2];
                } reclth;

long int getrecsize()
{
  reclth.word[0] = getint();
  reclth.word[1] = getint();
  return reclth.number;
  }

void putrecsize( recsize )
int recsize;
{
   reclth.number = (long) recsize;
   putint( reclth.word[0] );
   putint( reclth.word[1] );
}


/***************************************************************
  getlinerec( line ) Read a R-M FORTRAN unformatted record into a
                     character array

  putlinerec( line ) Write a R-M FORTRAN unformatted record from 
                     a character array
 ***************************************************************/

int ill_linerec()
{
  printf("Illegal line record in getlinerec\n");
  return ERROR;
}

int getlinerec( line )
char line[];
{
int  i,  linlth;

   linlth = getrecsize()-8;
   if (linlth > MAXLINLTH) return ill_linerec();

   for ( i=0; (i<linlth) && !feof(inputfile); i++ ) 
       line[i] = (char) fgetc( inputfile );

   if ( (linlth != getrecsize()-8) || feof(inputfile)) return ill_linerec();
   else return OK;
}


void putlinerec( line )
char line[];
{
int  i;
   putrecsize( MAXLINLTH+8 );
   i = 0;
   while (line[i] && (i < MAXLINLTH) )   /* output string up to null */
     fputc( line[i++], outputfile ) ;  
   while (i++ < MAXLINLTH)               /* pad the line with blanks */
      fputc( ' ', outputfile );
   putrecsize( MAXLINLTH+8 );
}


/****************************************************************************
* getreal()       Convert 8087-format real*4 number from the disk buffer
*                 to a double precision number for an IBM 370
* putreal( rout ) Convert an IBM370 real*8 number to IEEE real*4 and
*                 write it to disk
****************************************************************************/

union fpdatum {
       float number;
       char byte[4];
              } fp;

double getreal()
{
long int i;
   for (i = 3; i >= 0; i--) 
      fp.byte[i] = (char) fgetc(inputfile);
   return (double) fp.number;
}

void putreal( double *rout )
{
int i;
   fp.number = (float) *rout;
   for (i = 3; i >= 0; i--) 
     fputc( fp.byte[i], outputfile );
}


/**********************************************************************
 getarray( array, length, maxpts )  Retrieve one or more R-M FORTRAN 
                            unformatted records containing a real*4
                            array from the disk file. Does not return
			    more than the specified maximum # of pts.
 puarray( array, length )   Write one or more R-M FORTRAN unformatted
                            records containing a real number array to
                            the disk file.
 **********************************************************************/
int ill_arrayrec()
{
   printf("Illegal array record\n");
   return ERROR;
}

int getarray( array, arrylth, maxpts )
double array[];
long int arrylth, maxpts;
{
int i, lastel, nxtarry;
   nxtarry = 0;
   do {
      lastel = (getrecsize()-8) / 4;
      for ( i=1; i<=lastel; i++ ) {
	 if (nxtarry<maxpts) 
	   array[nxtarry++] = getreal();
	 else {
          nxtarry++;
          getreal();
	 }
         if (feof(inputfile)) break;
       }
      if ( lastel != ((getrecsize()-8)/4) )  return ill_arrayrec();
    }
while ( (nxtarry < arrylth) && !feof(inputfile) );
return arrylth != nxtarry;
}

void putarray( array, arrylth )
double array[];
long int   arrylth;
{
int i, elemnts, nxtel;
   nxtel = 0;
   do {
      elemnts = arrylth - nxtel;
      if (elemnts > MAXELEMENTS) elemnts = MAXELEMENTS;
      putrecsize( (4*elemnts) + 8 );
      for ( i=0; i < elemnts; i++ ) putreal( array+(nxtel++) );
      putrecsize( (4*elemnts) + 8 );
      }
   while (nxtel < arrylth);
}



int array_err( errstring )
char *errstring;
{
 printf("%s\n",errstring);
 return ERROR;
}


long int FORTRAN(getbin)( filename, xarry, yarry, npoints, comment, ncmts, maxpts )
char filename[], comment[][MAXLINLTH];
double xarry[], yarry[];
long int *npoints, *ncmts, *maxpts;
{
int cmtlns, i, iret;
   i = 0; 
                      
   while ((filename[i++] != ' ') && (i <= MAXLINLTH)); /* trim end blanks */
   filename[i-1] = 0;
   inputfile = fopen( filename, "r" );
   if (!inputfile) {
     printf("file open error %4d for file \'%s\'\n",errno,filename);
     return ERROR;
   }

   iret = getlinerec( info );
   if ( (iret) || strncmp( info, "VERSION GENPLT 4.0", 18 ) ) 
      return array_err( "Not a GENPLOT binary file");

/* Read comment strings (including any FIR file header) */

   cmtlns = 0;
   do 
     if ( getlinerec(comment[cmtlns]) ) 
        return array_err("Bad comment line record");
   while 
     ( strncmp(comment[cmtlns++],"\001",1) && !feof(inputfile)  );
   *ncmts = (long) -- cmtlns;


/* Read number of points and number of arrays (not used) */

  i = getrecsize();
  *npoints = getint();
  (void)getint();
  if (i != getrecsize()) return array_err("Bad npoints, ncol record");

/* Read X data (Field) */

  iret = getlinerec( info );         /* comment for x array */
  if (iret) return array_err("Bad X comment");
  i = getrecsize();
  factor = getreal();           /* multiplication factor for array */
  offset = getreal();           /* offset for array                */
  if (i != getrecsize()) 
    return array_err("Bad X factor, offset record");
  if (getarray( xarry, *npoints, *maxpts ) ) 
    return array_err("Could not get X array");

/* Read Y data (Intensity) */

  iret = getlinerec( info );         /* comment for y array */
  if (iret) return array_err("Bad Y comment");
  i = getrecsize();
  factor = getreal();           /* multiplication factor for array */
  offset = getreal();           /* offset for array                */
  if (i != getrecsize()) return array_err("Bad Y factor, offset record");

  if ( getarray( yarry, *npoints, *maxpts ) ) 
    return array_err("Could not get Y array");
  fclose( inputfile );
  return OK;
}

/* ****************************************************************** */

long int putbin( filename, xarry, yarry, npoints, comment, ncmts )
char filename[], comment[][MAXLINLTH];
double xarry[], yarry[];
long int *npoints, *ncmts;
{
int i, ncol;
   i = 0; 
   while ((filename[i++] != ' ') && (i <= MAXLINLTH)); /* trim end blanks */
   filename[i-1] = 0;
                     
   outputfile = fopen( filename, "w" );
   if (!outputfile) {
     printf("file open error %4d\n",errno);
     return ERROR;
   }

   putlinerec( "VERSION GENPLT 4.0" );
   for ( i = 0; i < *ncmts; i++ ) putlinerec( comment[i] );
   putlinerec( "\001" );

   putrecsize( 12 );
   ncol = 2;
   putint( (int) *npoints );
   putint( (int) ncol );
   putrecsize( 12 );

/* Write X data (Field) */
  factor=1.0;
  offset=0.0;
  putlinerec( "Field" );
  putrecsize( 16 );
  putreal( &factor );
  putreal( &offset );
  putrecsize( 16 );
  putarray( xarry, *npoints );
  
/* Write Y data (Intensity) */

  putlinerec( "Intensity" );
  putrecsize( 16 );
  putreal( &factor );
  putreal( &offset );
  putrecsize( 16 );
  putarray( yarry, *npoints );
  fclose( outputfile );
  return OK;
}
