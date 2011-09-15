/*

 $Id: readmap.cc,v 1.11 2009/05/08 23:02:17 rhuey Exp $
 $Id: readmap.cc, autodock4.lga.MPI v0.0 2010/10/10 22:55:01 collignon Exp $

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
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <math.h>
#include "readmap.h"
//#include "hdf5.h" in readmap.h
//#include "hdf5_hl.h" in readmap.h

extern char dock_param_fn[];
extern char *programname;
extern int ignore_errors;
extern int ElecMap;
extern FILE *logFile;
extern int debug;

extern float data1[MAX_ATOM_TYPES_DB][MAX_GRID_PTS_TOT]; // the maximum number of line (- 6) in the rec.X.map files 
                                        // that is : grid points only.
                                        // has to be changed as it is working for parp only

extern char data2[6][MAX_ATOM_TYPES_DB][LINE_LEN];
extern int TAB_MAP_TOT[MAX_MAPS_DB];
extern int num_filename[MAX_ATOM_TYPES_DB];
extern bool BOOL_READ_MAP;
extern bool bool_filename[MAX_ATOM_TYPES_DB];
extern int number_of_map;
extern int MAP_TOT;

char mapf2c(Real);

Statistics readmap( char           line[LINE_LEN],
                    int            outlev,
                    Clock          jobStart,
                    struct tms     tmsJobStart,
                    Boole          B_charMap,
                    Boole          *P_B_HaveMap, 
                    int            num_maps, 
                    GridMapSetInfo *info,
                    Real           map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    // double *maps 
                    char           map_type,
                    int myrank,
                    int            num_maps_cur
                  )

{
//    FILE *map_file;

    char FileName[PATH_MAX];
    char FldFileName[PATH_MAX];
    char GpfName[PATH_MAX];
    char ExtGpfName[PATH_MAX];
    char message[LINE_LEN];
    char mmFileName[PATH_MAX];
    char xyz_str[4];
    char C_mapValue;
//    char map_line[LINE_LEN];
//    char inputline[LINE_LEN];
    char atom_type_name[MAX_CHARS];

    Real cen[SPACE];
    Real spacing = 0.;
    double map_max;
    double map_min;
    double map_total;

    int indpf = 0;
    int nel[SPACE];
    int nv=0;
    int nvExpected = 0;

    register int xyz = 0;
    register int i = 0;
    register int j = 0;
    register int k = 0;

    struct tms tms_jobEnd;
    struct tms tms_loadEnd;
    struct tms tms_loadStart;

    Clock jobEnd;
    Clock loadEnd;
    Clock loadStart;

    Statistics map_stats;

    strcpy( xyz_str, "xyz\0" );

// DEFINE HDF5 VARIABLES
    hid_t file_map;
    herr_t status;

// DEFINE OTHER VARIABLES FOR HDF5 READING
    char buf_name0[LINE_LEN];
    char buf_name1[LINE_LEN];
    char buf_name2[LINE_LEN];
    char buf_name3[LINE_LEN];
    char buf_name4[LINE_LEN];
    char buf_name5[LINE_LEN];
    char buf_name6[LINE_LEN];
    char buf_name7[LINE_LEN];
    char str[LINE_LEN];
    int i_map=0;
    bool file_map_test;

    static int nb_buf0=0;
    static int nb_buf1=0;
    static int nb_buf2=0;
    static int nb_buf3=0;
    static int nb_buf4=0;
    static int nb_buf5=0;
    static int nb_buf6=0;
    static int nb_buf7=0;
    bool BOOL_READ_BUF0=FALSE;
    bool BOOL_READ_BUF1=FALSE;
    bool BOOL_READ_BUF2=FALSE;
    bool BOOL_READ_BUF3=FALSE;
    bool BOOL_READ_BUF4=FALSE;
    bool BOOL_READ_BUF5=FALSE;
    bool BOOL_READ_BUF6=FALSE;
    bool BOOL_READ_BUF7=FALSE;

    //maps->atom_type = num_maps;

    /*
    \  ATOMIC AFFINITY or ELECTROSTATIC GRID MAP
     \  Read in active site grid map...
      \____________________________________________________________
     */

    (void) sscanf( line, "%*s %s", FileName );
//       (void) sscanf( line, "%s", FileName );
//BOOL_READ_MAP=FALSE;
//if(FileName.length() == NULL){
//  if ( openFile( FileName, "r", &map_file, jobStart,tmsJobStart,TRUE )) {
//        *P_B_HaveMap = TRUE;
//      if (debug > 0) {
//          for (i=0; i < info->num_atom_types; i++) {
//              (void) fprintf(logFile, "info->atom_type_name[%d] = \"%s\"\n", i, info->atom_type_name[i] );
//          }
//       }
//   }  
//}
    
/* Check if the map in *.h5 has been read */
if(num_filename[num_maps_cur] == num_maps_cur) {
    bool_filename[num_maps_cur]=TRUE;
} else {
    bool_filename[num_maps_cur]=FALSE;
    num_filename[num_maps_cur]=num_maps_cur;

    /* buf_name0,1,2...are the names of the map datasets in the *.h5 file */
    sprintf (buf_name0,"%s%s",FileName,"map");       // Set the name of the buffer containing the grid
    sprintf (buf_name1,"%s%s",FileName,"head1");     // Set the name of the buffer containing the first line of charcater
    sprintf (buf_name2,"%s%s",FileName,"head2");    
    sprintf (buf_name3,"%s%s",FileName,"head3");     
    sprintf (buf_name4,"%s%s",FileName,"head4");
    sprintf (buf_name5,"%s%s",FileName,"head5");
    sprintf (buf_name6,"%s%s",FileName,"head6");
    sprintf (buf_name7,"%s%s","Total","nummap");
}

/* Read the *.h5 file until all maps have been read */
if (BOOL_READ_MAP == FALSE){  
   /* Read the map in the *.h5 file if it has not been read  */
   if (bool_filename[num_maps_cur] == FALSE) {  
       bool_filename[num_maps_cur] = TRUE;       
       /* Open the binary *.h5 file that has been previously created */
       /* The name of the *h5 file is "file.h5". It has previously been created automatically  */
       file_map = H5Fopen("maps.h5",H5F_ACC_RDONLY,H5P_DEFAULT); // If leave the map open can get fclose error.
       if(file_map != NULL) {
              *P_B_HaveMap = TRUE;
              if(BOOL_READ_BUF0 == FALSE && nb_buf0 < MAP_TOT){
               BOOL_READ_BUF0 = TRUE;
               status = H5LTread_dataset_float(file_map,buf_name0,data1[num_maps_cur]);
                  number_of_map++;
                  nb_buf0++;
              }
              if(BOOL_READ_BUF1 == FALSE && nb_buf1 < MAP_TOT){
                BOOL_READ_BUF1 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name1,data2[0][num_maps_cur]);
                  number_of_map++;
                  nb_buf1++;
              }
              if(BOOL_READ_BUF2 == FALSE && nb_buf2 < MAP_TOT){
                 BOOL_READ_BUF2 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name2,data2[1][num_maps_cur]);
                 BOOL_READ_BUF2 = TRUE;
                  number_of_map++;   
                  nb_buf2++;
              }
              if(BOOL_READ_BUF3 == FALSE && nb_buf3 < MAP_TOT){
                 BOOL_READ_BUF3 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name3,data2[2][num_maps_cur]);
                  number_of_map++;
                  nb_buf3++;
              }
              if(BOOL_READ_BUF4 == FALSE && nb_buf4 < MAP_TOT){
                BOOL_READ_BUF4 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name4,data2[3][num_maps_cur]);
                  number_of_map++;
                  nb_buf4++;
              }
              if(BOOL_READ_BUF5 == FALSE && nb_buf5 < MAP_TOT){
                 BOOL_READ_BUF5 = TRUE; 
                 status = H5LTread_dataset_char(file_map,buf_name5,data2[4][num_maps_cur]);
                 number_of_map++;
                 nb_buf5++;
              }  
              if(BOOL_READ_BUF6 == FALSE && nb_buf6 < MAP_TOT){
                 BOOL_READ_BUF6 = TRUE;
                 status = H5LTread_dataset_char(file_map,buf_name6,data2[5][num_maps_cur]);
                 number_of_map++;
                 nb_buf6++;
              }
              if(BOOL_READ_BUF7 == FALSE && nb_buf7 < MAP_TOT){
                 BOOL_READ_BUF7 = TRUE;
                 status = H5LTread_dataset_int(file_map,buf_name7,&TAB_MAP_TOT[num_maps_cur]);
                 number_of_map++;
                 nb_buf7++;
              }
          
              MAP_TOT=TAB_MAP_TOT[num_maps_cur];

        status = H5Fclose(file_map); // If leave the map open can get fclose error.       
        } else {
        printf ("ERROR WHEN OPENING THE *H5 FILE  \n");
        }
   }  
        if(number_of_map == MAP_TOT*7) {
           BOOL_READ_MAP=TRUE;
        }
}
 
        file_map_test=TRUE;

        if (map_type == 'e') {
            strcpy(atom_type_name, "e\0");
        } else if (map_type == 'd') {
            strcpy(atom_type_name, "d\0");
        } else {
            strcpy(atom_type_name, info->atom_type_name[num_maps]);
//            strcpy(atom_type_name, info->atom_type_name[num_maps_cur]);
  
        }
        if (outlev >= -1) {
        pr( logFile, "Opened Grid Map %d (%s):\t\t\t\t%s\n", num_maps+1, atom_type_name, FileName );
        }
        if (!ignore_errors && outlev > -1) {
            pr( logFile, "Checking header information.\n" );
        }
         /*
         \ Check header lines of grid map... 
         /
         \ :Line 1  GRID_PARAMETER_FILE 
        */
//   } else {

//        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
          if(!file_map_test) {
            warn_bad_file( FileName,"Could not read GRID_PARAMETER_FILE line." );
        } else {
            (void) sscanf(data2[0][num_maps_cur], "%*s %s", GpfName);
            if ( strindex( dock_param_fn, ".dpf" ) == -1) {
                pr_2x( stderr, logFile,"Can't find \".dpf\" in the dock-parameter filename.\n\n" );
                pr_2x( stderr, logFile,"AutoDock needs the extension of the grid parameter file to be \".gpf\"\nand that of the docking parameter file to be \".dpf\".\n\n" );
            } else {
                /*
                \ replace ".dpf" with ".gpf".
                */
                indpf = strindex( dock_param_fn, "dpf" );
                strcpy(ExtGpfName, dock_param_fn);
                ExtGpfName[ indpf ] = 'g';
            }
        } /* endif */
         /*
         \ :Line 2  GRID_DATA_FILE 
        */
//        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
          if(!file_map_test){
            warn_bad_file( FileName,"Could not read \".fld\" GRID_DATA_FILE line." );
        } else {
            (void) sscanf(data2[1][num_maps_cur], "%*s %s", FldFileName);
            if (!ignore_errors) {
                check_header_line( FldFileName, info->FN_gdfld );
            } /* endif */
        } /* endif */
         /*
         \ :Line 3  MACROMOLECULE 
        */
        if (!file_map_test) {
            warn_bad_file( FileName,"Could not read MACROMOLECULE line." );
        } else {
            (void) sscanf(data2[2][num_maps_cur],"%*s %s", mmFileName);
            check_header_line( mmFileName, info->FN_receptor );
        } /* endif */
         /*
         \ :Line 4  SPACING 
        */
        if (!file_map_test) {
            warn_bad_file( FileName,"Could not read SPACING line." );
        } else {
            (void) sscanf(data2[3][num_maps_cur],"%*s " FDFMT, &spacing);
            check_header_float(spacing, info->spacing, "grid point spacing", FileName );
        } /* endif */
         /*
         \ :Line 5  NELEMENTS 
        */
        if (!file_map_test) {
            warn_bad_file( FileName,"Could not read NELEMENTS line." );
        } else {
            (void) sscanf(data2[4][num_maps_cur],"%*s %d %d %d", &nel[X], &nel[Y], &nel[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                //maps->num_points[xyz] = nel[xyz];
                //maps->num_points1[xyz] = nel[xyz] + 1;
                check_header_int( nel[xyz], info->num_points[xyz], xyz_str[xyz], FileName );
            } /* xyz */
        } /* endif */
         /* 
         \ :Line 6  CENTER
        */
        if (!file_map_test) {
            warn_bad_file( FileName,"Could not read CENTER line." );
        } else {
            (void) sscanf(data2[5][num_maps_cur],"%*s " FDFMT3, &cen[X], &cen[Y], &cen[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                //maps->center[xyz] = cen[xyz];
                check_header_float(cen[xyz], info->center[xyz], "grid-map center", FileName );
            } /* xyz */
        } /* endif */
//    } /* endif */
//    flushLog;

    /*
    \   Now find the extrema of the grid-map energies,
     \  While reading in the values...
      \____________________________________________________________
     */

    map_max = -BIG;
    map_min =  BIG;
    map_total = 0.;
    nvExpected = info->num_points1[X] * info->num_points1[Y] * info->num_points1[Z];
    nv = 0;
    if (outlev > 0) {
    pr( logFile, "Number of grid points expected in  x-dimension:  %d\n", info->num_points1[X] );
    pr( logFile, "Number of grid points expected in  y-dimension:  %d\n", info->num_points1[Y] );
    pr( logFile, "Number of grid points expected in  z-dimension:  %d\n", info->num_points1[Z] );
    pr( logFile, "Looking for %d energies from Grid Map %d... \n", nvExpected, num_maps+1 );
    }
//    flushLog;
    loadStart = times( &tms_loadStart );

    for ( k = 0;  k < info->num_points1[Z];  k++) {
        for ( j = 0;  j < info->num_points1[Y];  j++) {
            for ( i = 0;  i < info->num_points1[X];  i++) {
                if (B_charMap) {
//                    if (fgets(map_line, LINE_LEN, map_file) != NULL) { /*new*/
                    if (file_map_test) {
//                        if (sscanf( map_line,  "%c",  &C_mapValue ) != 1) continue;
//                        map[k][j][i][num_maps] = mapc2f(C_mapValue);

                         /* CONVERT THE GRID INTO STRING - HAS TO BE CHANGE */             
                         sprintf(str, "%f", data1[num_maps_cur][i_map]);
                         if (sscanf( str,  "%c",  &C_mapValue ) != 1) continue;
                         map[k][j][i][num_maps] = mapc2f(C_mapValue);
                         
                         ++i_map;
                         nv++;
                    }
                } else {
                  if (file_map_test) { /*new*/
//                      if (sscanf( map_line,  FDFMT,  &map[k][j][i][num_maps] ) != 1) continue;
//                      nv++;
//                  
                         sprintf(str, "%f", data1[num_maps_cur][i_map]);
                         if (sscanf( str,  FDFMT ,  &map[k][j][i][num_maps] ) != 1) continue;

                         i_map++;
                         nv++;
                    } 
                }
                map_max = max( map_max, map[k][j][i][num_maps] );
                map_min = min( map_min, map[k][j][i][num_maps] );
                map_total += map[k][j][i][num_maps];
            }

        }
    }

    i_map=0;
    map_stats.number = nv;
    map_stats.minimum = map_min;
    map_stats.maximum = map_max;
    map_stats.mean = map_total / (double) nv;

    /*
    if (map_stats.number > 1) {
        double deviation = 0.;
        double sum_squares = 0.;
        for ( k = 0;  k < info->num_points1[Z];  k++) {
            for ( j = 0;  j < info->num_points1[Y];  j++) {
                for ( i = 0;  i < info->num_points1[X];  i++) {
                    deviation = map[k][j][i][num_maps] - map_stats.mean;
                    sum_squares += deviation * deviation;
                }
            }
        }
        map_stats.standard_deviation = sqrt(sum_squares / (map_stats.number - 1));
    } else {
        map_stats.standard_deviation = 0.;
    }
    */

//    pr( logFile, "Closing file.\n" );
//    fclose( map_file );
//  status = H5Fclose(file_map);
    if (outlev >= -1) {
    pr( logFile, "%d energies found for map %d\n", nv, num_maps+1 );
    }
    if (map_type == 'e' && outlev > -1) {
          pr( logFile, "Minimum electrostatic potential = %.2f,  maximum electrostatic potential = %.2f\n\n", map_min, map_max );
    } else if (outlev > -1) {
          pr( logFile, "Minimum energy = %.2f,  maximum energy = %.2f\n\n", map_min, map_max );
    }
    if (outlev > -1) {
    pr( logFile, "Time taken (s): " );
    }

    loadEnd = times( &tms_loadEnd );
    if (outlev > -1) {
    timesys( loadEnd - loadStart, &tms_loadStart, &tms_loadEnd );
    }

    pr( logFile, "\n" );

    if (nv != nvExpected ) {
        prStr( message, "\n%s: wrong number of values read in. Check grid map!\n\n", programname  );
        pr_2x( stderr, logFile, message );

        jobEnd = times( &tms_jobEnd );
        timesys( jobEnd - jobStart, &tmsJobStart, &tms_jobEnd );
        pr_2x( logFile, stderr, UnderLine );

        /* END PROGRAM */
        exit(-1);
    } 

//    flushLog;

    return map_stats;
}

Real mapc2f(char numin)
{
    Real numout;
    if (numin == 0) {
        numout = 0.;
    } else if (numin > 0) {
        numout = numin * 10.;
    } else {
        numout = numin /10.;
    }
    return numout;
}


/*
    char mapf2c(Real numin)
    {
        char numout;
        if (numin == 0.) {
            numout = 0;
        } else if ((-12.8 < numin) && (numin < 0.)) {
            numout = numin * 10.;
        } else if ((0. < numin) && (numin < 1280.)) {
            numout = numin / 10.;
        } else if (numin >= 1280.) {
            numout = 127;
        } else {
            numout = -128;
        }
        return numout;
    }
*/

/* EOF */
