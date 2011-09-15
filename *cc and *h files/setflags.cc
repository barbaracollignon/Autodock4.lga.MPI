/*

 $Id: setflags.cc,v 1.13 2009/05/08 23:02:17 rhuey Exp $
 $Id: setflags.cc, autodock4.lga.MPI v0.0 2010/10/10 22:55:01 collignon Exp $

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "setflags.h"
#include "openfile.h"
#include "version.h"
#include "banner.h"
#define LIG_MAX 5000000
#define PAR_MAX 45

extern FILE *parFile;
extern FILE *ligFile;
extern FILE *logFile;
extern FILE *stateFile;
extern int  write_stateFile;
extern char *programname;

extern char dock_param_fn[];
extern char dock_lig_list[LIG_MAX][LINE_LEN] ;
char dock_par_list[PAR_MAX][LINE_LEN] = {0} ;
char dock_par_store[PAR_MAX][LINE_LEN] = {0} ;
extern char AutoDockHelp[];
extern int  debug;
extern int  ignore_errors;
extern int  parse_tors_mode;
extern int  keepresnum;
extern bool BOOL_LIGFILE;
extern bool BOOL_PARFILE;
extern bool BOOL_NBLIG;
extern bool BOOL_NBMAP;
//extern int int_liglist;

extern int SIZE;
extern int LIG_TOT;
extern int MAP_TOT;


int setflags( int argc, char ** argv, char * version_num, int id, int id2)

/*
** naming convention: 
** var   => var is an integer variable;
** var => var is a pointer to a pointer to a character variable
*/

/******************************************************************************/
/*      Name: setflags                                                        */
/*  Function: read flags from argv; return argindex of first non arg.   */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Matthew Morris, TSRI.                                   */
/*            (Adapted from code supplied by Bruce Duncan, TSRI.)             */
/*      Date: 02/02/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: argc,argv                                                 */
/*   Returns: argindex                                                      */
/*   Globals: *parFile;				                              */
/*            *logFile;					                      */
/*            *programname;						      */
/*            dock_param_fn[];						      */
/*            ignore_errors;						      */
/*            command_mode;						      */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/11/92 GMM     Modified for Autodock flags:                              */
/*                  -c = Command mode;                                        */
/*                  -p = Parameter filename;                                  */
/*                  -l = Log filename;                                        */
/*                  -o = Use old PDBq format (q in columns 55-61)             */
/*                  -d = Debug mode;                                          */
/* 01/22/93 GMM     -i = Ignore header-checking                               */
/* 06/14/93 GMM     -t = Parse the PDBq file to check torsions, then stop.    */
/* 02/02/94 GMM     -k = Keep original residue number in output PDB clusters. */
/******************************************************************************/

{
    float atest=0.0;
    int argindex;
    int int_liglist=0;
    int int_parlist=0; 
    int int_parlist_store=0;
    char filename[PATH_MAX];
    char lineligfile[LINE_LEN];
    char lineparfile[LINE_LEN];
    char seed_key[LINE_LEN];  
    char paralgt[LINE_LEN]="ligand_types";
    char parafld[LINE_LEN]="fld";
    char paramap[LINE_LEN]="map";
    char paraelec[LINE_LEN]="elecmap";
    char paradesolv[LINE_LEN]="desolvmap";
    char paramove[LINE_LEN]="move";
    char paraabout[LINE_LEN]="about";
    char paratorsdof[LINE_LEN]="torsdof";
/*----------------------------------------------------------------------------*/
/* Initialize                                                                 */
/*----------------------------------------------------------------------------*/
    argindex = 1;
    programname = argv[0];
//    parFile = stdin;
//    ligFile = stdin;
    logFile = stdout;
    /*
     * see autoglobal.h for initialization of debug, keepresnum and logicals...
     */
    if (argc==1) { //No arguments provided
        usage(stdout, "AutoDock");
        exit(0);
    }
/*----------------------------------------------------------------------------*/
/* Loop over arguments                                                        */
/*----------------------------------------------------------------------------*/
    while((argc > 1) && (argv[1][0] == '-')){
        if (argv[1][1] == '-') argv[1]++;

        switch(argv[1][1]){
#ifdef FOO
        case 'n':
            ncount = atoi(argv[2]);
            argv++;
            argc--;
            argindex++;
            break;
#endif
        case 'd':
            debug++;
            break;
        case 'u':
        case 'h':
            usage(stdout, "AutoDock");
            exit(0);
            break;
        case 'i':
            ignore_errors = TRUE;
            break;
        case 'k':
            keepresnum--;
            break;
        case 'C':
            //show copyright
            show_copyright(stdout);
            show_warranty(stdout);
            exit(0);
            break;
        case 'c':
            //command_mode removed with 4.1 release spring 2009, mp + rh
            fprintf(stderr, "\n%s: command mode is not supported in this version of autodock\n",programname );
            break;
        case 'N':
            if(BOOL_NBLIG == FALSE) {
            LIG_TOT = atoi(argv[2]);
            atest = (LIG_TOT / (SIZE-1)); 
                if(atest < 1.0){
                   fprintf(stderr,"\n %d LIGANDs proceed on %d CPUs\n",LIG_TOT, SIZE-1);
                   fprintf(stderr, "\n%s: The number of Ligands has to be greater that the number of CPUs\n", programname);
                   fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                   return(-1);
                 }
            BOOL_NBLIG = TRUE;
            } 
            argv++;
            argc--;
            argindex++;
            break;
//        case 'G':
//            if (BOOL_NBMAP == FALSE) {
//            MAP_TOT = atoi(argv[2]);
//            BOOL_NBMAP = TRUE;
//            }
//            argv++;
//            argc--;
//            argindex++;
//            break; 
        case 'l':
            if ( (logFile = ad_fopen(argv[2], "w")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr,"\n Log file name = %s\n",argv[2]); 
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't create log file %s\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
//            Remove setlinebuf as it could cost performance
//            setlinebuf(logFile); // to ensure output even if crash
            argv++;
            argc--;
            argindex++;
            break;
//        case 'L': /* New flag for the MPI version */
//            /* Write one output per ligand */
//            sprintf(filename,"dock%d_%d.dlg", id,id2); // id is the rank. id2 is the work id.
//            if ( (logFile = ad_fopen(filename, "w")) == NULL ){
//#ifdef DEBUG
//                fprintf(stderr,"\n Log file name = %s\n",filename);
//#endif /* DEBUG */
//                fprintf(stderr, "\n%s: can't create log file %s\n", programname,filename);
//                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
//                return(-1);
//            }
//            Remove setlinebuf as it could cost performance
//            setlinebuf(logFile); // to ensure output even if crash
//            argv++;
//            argc--;
//            argindex++;
//            break;
        case 's':
            if ( (stateFile = ad_fopen(argv[2], "w")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr,"\n State file name = %s\n",argv[2]); 
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't create state file %s\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
            else{
                fprintf(stateFile,"<?xml version=\"1.0\" ?>\n");
                fprintf(stateFile,"<autodock>\n");
                fprintf(stateFile,"\t<version>%s</version>\n", version_num);
                fprintf(stateFile,"\t<autogrid_version>%s</autogrid_version>\n", version_num);
                fprintf(stateFile,"\t<output_xml_version>%5.2f</output_xml_version>\n", OUTPUT_XML_VERSION);
                write_stateFile = TRUE;
            }
            argv++;
            argc--;
            argindex++;
            break;    
        case 'p':
            strcpy(dock_param_fn, argv[2]);
        if (BOOL_PARFILE == FALSE) {
            if ( (parFile = ad_fopen(argv[2], "r")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr,"\n Parameter file name = %s\n",argv[2]);
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't find or open parameter file %s\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            } else {
                  strcpy(dock_par_list[3],paralgt);
                  strcpy(dock_par_list[4],parafld);
                  strcpy(dock_par_list[5],paramap);
                  strcpy(dock_par_list[6],paraelec);
                  strcpy(dock_par_list[7],paradesolv);
                  strcpy(dock_par_list[8],paramove);
                  strcpy(dock_par_list[9],paraabout);  
                  strcpy(dock_par_list[17],paratorsdof);
                  while(fgets(lineparfile,LINE_LEN,parFile) !=  NULL) {
                    strcpy(dock_par_list[int_parlist],lineparfile);
                       int_parlist++;
                  }
            BOOL_PARFILE=TRUE;
            fclose(parFile);
            }
        }
            argv++;
            argc--;
            argindex++;
            break;
        case 'L':   /* New flag for the MPI version */             
        if (BOOL_LIGFILE == FALSE) { /* Do it only once */
            /* Open the list that contains the ligands and their specific parameters */ 
            if ((ligFile = ad_fopen(argv[2], "r")) == NULL ) { 
#ifdef DEBUG
                fprintf(stderr,"\n Ligand list file name = %s\n",argv[2]);
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't find or open ligand list file %s\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            } else {
                  while (fgets(lineligfile,LINE_LEN,ligFile) != NULL && int_liglist < LIG_TOT) {
                         strcpy(dock_lig_list[int_liglist],lineligfile);
                  int_liglist++;
                   }
            BOOL_LIGFILE=TRUE; 
            fclose(ligFile);
            }
        }
          /* Write one output per ligand */
            sprintf(filename,"dock%d_%d.dlg", id,id2); // id is the rank. id2 is the work id.
            if ( (logFile = ad_fopen(filename, "w")) == NULL ){
#ifdef DEBUG
                fprintf(stderr,"\n Log file name = %s\n",filename);
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't create log file %s\n", programname,filename);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
//            Remove setlinebuf as it could cost performance
//            setlinebuf(logFile); // to ensure output even if crash
            argv++;
            argc--;
            argindex++;
            break;
        case 't':
            parse_tors_mode = TRUE;
            break;
        case 'v':
            fprintf(stdout, "AutoDock %-8s\n", version_num);
            fprintf(stdout, " Copyright (C) 2009 The Scripps Research Institute.\n");
            fprintf(stdout, " License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n");
            fprintf(stdout, " This is free software: you are free to change and redistribute it.\n");
            fprintf(stdout, " There is NO WARRANTY, to the extent permitted by law.\n");
            exit(0);
            break;
        default:
            fprintf(stderr,"%s: unknown switch \"-%c\".  \n",programname,argv[1][1]);
            usage(stderr, programname);
            return(-1);
            /* break; */
        }
        argindex++;
        argc--;
        argv++;
    }
    return(argindex);
}
/* EOF */
