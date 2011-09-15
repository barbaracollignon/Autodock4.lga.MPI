/*

 $Id: readfield.cc,v 1.5 2009/05/08 23:02:17 rhuey Exp $
 $Id: readfield.cc, autodock4.lga.MPI v0.0 2010/10/10 22:55:01 collignon Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include "readfield.h"

extern FILE *logFile;
extern int num_filefld[1];
extern bool bool_filefld[1];
extern char datafld[5][1][LINE_LEN];
extern int number_of_fld;
extern bool BOOL_READ_FLD;
extern int FLD_TOT;

void readfield( GridMapSetInfo *info,
                char line[LINE_LEN],
                Clock jobStart,
                struct tms tms_jobStart,
                int num_fld_cur)  


{
    FILE *fldFile;
    char rec9[9], inputline[LINE_LEN];
    double localSpacing;
    register int i = 0;

// DEFINE HDF5 VARIABLES
    hid_t file_map;
    herr_t status;

// DEFINE OTHER VARIABLES FOR HDF5 READING
    char buf_name0[LINE_LEN];
    char buf_name1[LINE_LEN];
    char buf_name2[LINE_LEN];
    char buf_name3[LINE_LEN];
    char buf_name4[LINE_LEN];
 
    static int nb_buf0=0;
    static int nb_buf1=0;
    static int nb_buf2=0;
    static int nb_buf3=0;
    static int nb_buf4=0;
    bool BOOL_READ_BUF0=FALSE;
    bool BOOL_READ_BUF1=FALSE;
    bool BOOL_READ_BUF2=FALSE;
    bool BOOL_READ_BUF3=FALSE;
    bool BOOL_READ_BUF4=FALSE;
    bool file_map_test=FALSE;


    /**
    ** GRID_DATA_FILE
    ** Read the (AVS-format) grid data file, .fld
    ** _____________________________________________________________
    **/
    (void) sscanf( line, "%*s %s", info->FN_gdfld);
     
//    if ( openFile( info->FN_gdfld, "r", &fldFile, jobStart, tms_jobStart, TRUE )) {
//        pr( logFile, "Opening Grid Map Dimensions file:\t\t%s\n\n", info->FN_gdfld);
//    }

/* Check if the fld in *.h5 has been read */
if(num_filefld[num_fld_cur] == num_fld_cur) {
    bool_filefld[num_fld_cur]=TRUE;
} else {
    bool_filefld[num_fld_cur]=FALSE;
    num_filefld[num_fld_cur]=num_fld_cur;

    /* buf_name0,1,2...are the names of the map datasets in the *.h5 file */
    sprintf (buf_name0,"%s%s",info->FN_gdfld,"fld0");       // Set the name of the buffer containing the grid
    sprintf (buf_name1,"%s%s",info->FN_gdfld,"fld1");       // Set the name of the buffer containing the grid
    sprintf (buf_name2,"%s%s",info->FN_gdfld,"fld2");       // Set the name of the buffer containing the grid
    sprintf (buf_name3,"%s%s",info->FN_gdfld,"fld3");       // Set the name of the buffer containing the grid
    sprintf (buf_name4,"%s%s",info->FN_gdfld,"fld4");       // Set the name of the buffer containing the grid
}
    

/* Read the *.h5 file until the fld have been read */
if (BOOL_READ_FLD == FALSE){
   /* Read the map in the *.h5 file if it has not been read  */
   if (bool_filefld[num_fld_cur] == FALSE) {
       bool_filefld[num_fld_cur] = TRUE;
       /* Open the binary *.h5 file that has been previously created */
       /* The name of the *h5 file is "file3.h5". It has previously been created automatically  */
       file_map = H5Fopen("maps.h5",H5F_ACC_RDONLY,H5P_DEFAULT); // If leave the fld open can get fclose error.
       if(file_map != NULL) {
              if(BOOL_READ_BUF0 == FALSE && nb_buf0 < FLD_TOT){
                BOOL_READ_BUF0 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name0,datafld[0][num_fld_cur]);
                  number_of_fld++;
                  nb_buf0++;
              }
              if(BOOL_READ_BUF1 == FALSE && nb_buf1 < FLD_TOT){
                 BOOL_READ_BUF1 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name1,datafld[1][num_fld_cur]);
                 BOOL_READ_BUF1 = TRUE;
                  number_of_fld++;
                  nb_buf1++;
              }
              if(BOOL_READ_BUF2 == FALSE && nb_buf2 < FLD_TOT){
                 BOOL_READ_BUF2 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name2,datafld[2][num_fld_cur]);
                  number_of_fld++;
                  nb_buf2++;
              }
              if(BOOL_READ_BUF3 == FALSE && nb_buf3 < FLD_TOT){
                BOOL_READ_BUF3 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name3,datafld[3][num_fld_cur]);
                  number_of_fld++;
                  nb_buf3++;
              }
              if(BOOL_READ_BUF4 == FALSE && nb_buf4 < FLD_TOT){
                 BOOL_READ_BUF4 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name4,datafld[4][num_fld_cur]);
                 number_of_fld++;
                 nb_buf4++;
              }

        }
           status = H5Fclose(file_map); // If leave the map open can get fclose error.
   } else {
        printf ("ERROR WHEN OPENING THE *H5 FILE  \n");
   }
        if(number_of_fld == 5) {
           BOOL_READ_FLD=TRUE;
        }
}

    /** Skip over the AVS-readable .fld  header comments, until 
    **
    ** #SPACING
    */
//    while( fgets(line, LINE_LEN, fldFile) != NULL ) {
//      (void) sscanf(line,"%s", rec9);
        (void) sscanf(datafld[0][num_fld_cur],"%s", rec9);
        if (equal(rec9,"#SPACING", 8)) {
//        (void) sscanf(line,"%*s %lf", &localSpacing);
          (void) sscanf(datafld[0][num_fld_cur],"%*s %lf", &localSpacing);
            info->spacing = localSpacing;
//            break;
        }
//    } /* endwhile */
    info->inv_spacing = 1. / (info->spacing);
    pr( logFile, "Grid Point Spacing =\t\t\t\t%.3f Angstroms\n\n", info->spacing);

    /*
    ** #NELEMENTS 
    */
//    (void) fgets(inputline, LINE_LEN, fldFile);
//    (void) sscanf(inputline,"%*s %d %d %d", &info->num_points[X], &info->num_points[Y], &info->num_points[Z]);
     (void) sscanf(datafld[1][num_fld_cur],"%*s %d %d %d", &info->num_points[X], &info->num_points[Y], &info->num_points[Z]);
    /* Check that these are all even numbers... */
    for (i=0; i<SPACE; i++) {
        if ( (info->num_points[i])%2 != 0 ) {
            stop("the number of user-specified grid points must be even in the \"#NELEMENTS\" line in the \".fld\" file.");
            exit(-1);
        }
    }

  pr( logFile, "Even Number of User-specified Grid Points =\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", info->num_points[X],info->num_points[Y],info->num_points[Z]);
    for (i = 0;  i < SPACE;  i++) {
        info->num_points1[i] = info->num_points[i] + 1;
    } /* i */
// Comments removed to gain disk space and ouput lisibility
// pr( logFile, "Adding the Central Grid Point makes:\t\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", info->num_points1[X], info->num_points1[Y], info->num_points1[Z]);
//printf("Adding the Central Grid Point makes:\t\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", info->num_points1[X], info->num_points1[Y], info->num_points1[Z]);
    if ( (info->num_points[X] <= 0)||(info->num_points[Y] <= 0)||(info->num_points[Z] <= 0) ) {
 //       stop("insufficient grid points." );
 //       exit( -1 );
    } else if ((info->num_points[X] > MAX_GRID_PTS)||(info->num_points[Y] > MAX_GRID_PTS)||(info->num_points[Z] > MAX_GRID_PTS)) {
        stop("too many grid points." );
        exit( -1 );
    }

    /*
    ** #CENTER 
    */
//    (void) fgets(inputline, LINE_LEN, fldFile);
//    (void) sscanf(inputline,"%*s %lf %lf %lf", &info->center[X], &info->center[Y], &info->center[Z]);
      (void) sscanf(datafld[2][num_fld_cur],"%*s %lf %lf %lf", &info->center[X], &info->center[Y], &info->center[Z]);

  pr( logFile, "Coordinates of Central Grid Point of Maps =\t(%.3f, %.3f, %.3f)\n\n", info->center[X],  info->center[Y],  info->center[Z]);

    /*
    ** #MACROMOLECULE 
    */
//    (void) fgets(inputline, LINE_LEN, fldFile);
//    (void) sscanf(inputline,"%*s %s", info->FN_receptor);
      (void) sscanf(datafld[3][num_fld_cur],"%*s %s", info->FN_receptor); 

  pr( logFile, "Macromolecule file used to create Grid Maps =\t%s\n\n", info->FN_receptor);

    /*
    ** #GRID_PARAMETER_FILE 
    */
//    (void) fgets(inputline, LINE_LEN, fldFile);
//    (void) sscanf(inputline,"%*s %s", info->FN_gpf);
      (void) sscanf(datafld[4][num_fld_cur],"%*s %s", info->FN_gpf);
  pr( logFile, "Grid Parameter file used to create Grid Maps =\t%s\n\n", info->FN_gpf);

    /*
    ** Close Grid-dimensions data file 
    */
//    fclose(fldFile);

    /*
    ** Determine the dimensions of the grids,
    */
    for (i = 0;  i < SPACE;  i++) {
        info->lo[i] = info->center[i] - (info->num_points[i]/2) * (info->spacing);
        info->hi[i] = info->center[i] + (info->num_points[i]/2) * (info->spacing);
    }

//  Comments removed to gain disk space and ouput lisibility 
//  pr( logFile, "Minimum coordinates in grid = (%.3f, %.3f, %.3f)\n",   info->lo[X], info->lo[Y], info->lo[Z]);
//  pr( logFile, "Maximum coordinates in grid = (%.3f, %.3f, %.3f)\n", info->hi[X], info->hi[Y], info->hi[Z]);
//    fflush( logFile );
}
/*
** EOF
*/
