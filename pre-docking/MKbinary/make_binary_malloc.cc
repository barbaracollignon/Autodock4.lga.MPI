/*
$Id: make_binary_malloc.cc, last modified 2010/10/10 22:55:01 collignon Exp $
 
    Copyright (C) 2011  Author: Barbara Collignon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*

> cat map.dpf
 map myprot.A.map                       # atom-specific affinity map
 map myprot.C.map                       # atom-specific affinity map
 map myprot.Cl.map                      # atom-specific affinity map
 map myprot.F.map                       # atom-specific affinity map
 map myprot.OA.map                      # atom-specific affinity map
 map myprot.N.map                       # atom-specific affinity map
 map myprot.P.map                       # atom-specific affinity map
 map myprot.S.map                       # atom-specific affinity map
 map myprot.HD.map                      # atom-specific affinity map
 map myprot.e.map                       # electrostatic potential map
 map myprot.d.map                       # desolvation potential map
 fld myprot.maps.fld                    #field

> make_binary.exe -m map.dpf 
> ls *h5
  file.h5


*/

#include "hdf5.h"
#include "hdf5_hl.h"
#include "limits.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LINE_LEN 256

bool BOOL_H5F;
bool BOOL_H5F2;

int main(int argc, char **argv )
{
    int i_test=0;
    BOOL_H5F=false;
    BOOL_H5F2=false;

    FILE *pMAP;
    FILE *map_file_bis;

    char line[LINE_LEN];
    char FileNameMap[PATH_MAX];
    char indicMap[PATH_MAX];

//  DEFINE HDF5 VARIABLE1
    hid_t file_map,file_map2,string_dset,acc_tpl1,fid1;
    hsize_t dimsfb[1],dimsfbb[2],dimsfcc[2];
    herr_t status,ret;

//  DEFINE HDF5 VARIABLE2
    int FILE_LEN,FILE_LEN2;
    char buf_name0[LINE_LEN];
    char buf_name1[LINE_LEN];
    char buf_name2[LINE_LEN];
    char buf_name3[LINE_LEN];
    char buf_name4[LINE_LEN];
    char buf_name5[LINE_LEN];
    char buf_name6[LINE_LEN];
    char buf_name7[LINE_LEN];
    char buf_name8[LINE_LEN];
    char buf_name9[LINE_LEN];
    char buf_name10[LINE_LEN];
    char buf_name11[LINE_LEN];
    int fl_map_tmp=0;
    int fl_map_tmp2=0;
    int char_line=0;
    int count_line =0;
    float buf_tmpx=0.0;
    printf("test\n");
    char **buf_map ;
    printf("test1\n");
    char buf_head[6][LINE_LEN] = {0};
    printf("test2\n");
    char buf_fld[5][LINE_LEN] = {0};
    printf("test3\n");
    float *buf_map2 ;
    printf("test5\n");
    
//  READ MAP
    while((argc > 1) && (argv[1][0] == '-')) {
        if (argv[1][1] == '-') argv[1]++;
           switch(argv[1][1]){
               case 'm':   /* Flag to specify the file that contains the list of maps */
               /* Open the list of maps */
               if ((pMAP = fopen(argv[2], "r")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr,"\n Map list file name = %s\n",argv[2]);
#endif /* DEBUG */
                fprintf(stderr, "can't find or open %s\n", argv[2]);
                fprintf(stderr, "Unsuccessful Completion.\n\n");
                return(-1);
                }
                argv++;
                argc--;
                break;
      }
    }
    printf("Reading Map by Rank0");  
    while(fgets(line, LINE_LEN, pMAP)) {
        printf("i:%d\n",i_test);      
        (void) sscanf( line, "%*s %s", FileNameMap );
        (void) sscanf( line, "%s %*s", indicMap );
        map_file_bis=fopen( FileNameMap, "r"); 
             printf("%s\n",FileNameMap);
         
//      COUNT THE NUMBER OF LINE IN THE MAP FILE
        while ((char_line = fgetc(map_file_bis)) != EOF ) {

               if ( char_line == '\n' )
                  {
                  ++count_line; /* Bump the counter for every newline. */
                  }
        }
        rewind(map_file_bis);
        FILE_LEN=count_line;
        count_line=0;
        printf("%s %s : file length:%d\n",indicMap,FileNameMap,FILE_LEN);
        
        /* write the map in memory (buffer) */
        printf("write the map in memory (buffer)\n");

//        char buf_map[FILE_LEN][LINE_LEN] ;
        printf("free\n");
//         free(buf_map);
//         free(buf_map2);
         printf("free2\n");
         buf_map = (char**)malloc((FILE_LEN)*sizeof(char*));
         for (int i = 0; i < FILE_LEN; i++) { 
                buf_map[i] = (char*) malloc(LINE_LEN*sizeof(char));  
         }
//        char buf_head[6][LINE_LEN] = {0};
         printf("buf_map\n");
//        char buf_fld[5][LINE_LEN] = {0};
//        printf("test3\n");
//        char buf_xyz[FILE_LEN][LINE_LEN];
//        printf("test4\n");
//        float buf_map2[FILE_LEN-6] ;
         buf_map2 = (float*)malloc((FILE_LEN-6)*sizeof(float));
         printf("buf_map2\n");
//      GET EACH LINE OF THE MAP FILE , CONVERT INTO FLOAT, STORE IN buf_map 
        while ((fl_map_tmp < FILE_LEN) && (fgets(buf_map[fl_map_tmp],LINE_LEN,map_file_bis) != NULL)) {
//        printf("testxx\n");
          if(strcmp(indicMap,"map") == 0) {
              if(fl_map_tmp == 0) {
                strcpy(buf_head[0],buf_map[0]);
              }
              if(fl_map_tmp == 1){
                strcpy(buf_head[1],buf_map[1]);
              }
              if(fl_map_tmp == 2){
                strcpy(buf_head[2],buf_map[2]);
              }
              if(fl_map_tmp == 3){
                strcpy(buf_head[3],buf_map[3]);
              }
              if(fl_map_tmp == 4){
                strcpy(buf_head[4],buf_map[4]);
              }
              if(fl_map_tmp == 5){
                strcpy(buf_head[5],buf_map[5]);
              }
              if(fl_map_tmp >= 6) {
                  buf_tmpx = atof(buf_map[fl_map_tmp]);
                  buf_map2[fl_map_tmp2]= buf_tmpx ;
                  ++fl_map_tmp2;
              }
//          printf("write fld in memory buffer\n");
          } else if(strcmp(indicMap,"fld") == 0) {
              if(fl_map_tmp == 6) {
                strcpy(buf_fld[0],buf_map[6]);
              }
              if(fl_map_tmp == 7){
                strcpy(buf_fld[1],buf_map[7]);
              }
              if(fl_map_tmp == 8){
                strcpy(buf_fld[2],buf_map[8]);
              }
              if(fl_map_tmp == 9){
                strcpy(buf_fld[3],buf_map[9]);
              }
              if(fl_map_tmp == 10){
                strcpy(buf_fld[4],buf_map[10]);
              }
          } else if(strcmp(indicMap,"xyz") == 0) {
//                  strcpy(buf_xyz[fl_map_tmp],buf_map[fl_map_tmp]);
          }

        ++fl_map_tmp;
        }
        FILE_LEN2=fl_map_tmp2;
        fl_map_tmp2=0;
        fl_map_tmp=0;
        fclose(map_file_bis);

        printf("end writing map in memory (buffer)\n");

//  OPEN FILE - INDEPENDENT CALL
    if(BOOL_H5F == false){
       file_map = H5Fcreate("maps.h5",H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);
       printf("create file map1:%p",&file_map);
       BOOL_H5F = true;
    } else {
//    file_map = H5Fopen("maps.h5",H5F_ACC_RDWR,H5P_DEFAULT);
      file_map2=H5Freopen(file_map);
      printf("open file map1:%p",&file_map);
    }
    printf("TEST1");
//OPEN FILE - COLLECTIVE CALL
//      if(BOOL_H5F == false){
//      /* setup file access template with parallel IO access. */
//      acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
//        /* set Parallel access with communicator */
//        ret = H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);
//      /* open the file collectively */
//      fid1=H5Fopen("file2.h5",H5F_ACC_RDWR,acc_tpl1);
//      /* Release file-access template */
//      ret=H5Pclose(acc_tpl1);
//      printf("create file map1:%p",&fid1);
//      BOOL_H5F = true;
//      } else {

        /* To create the dataset we: */
        dimsfb[0] = FILE_LEN2 ;
        dimsfbb[0] = 1 ;
        dimsfbb[1] = LINE_LEN;
        dimsfcc[0] = 1;
        dimsfcc[1] = LINE_LEN;
        printf("TEST2");
        if(strcmp(indicMap,"map") == 0) {
            sprintf (buf_name0,"%s%s",FileNameMap,"map");
            printf("before create dset0:%s\n",buf_name0);
            sprintf (buf_name1,"%s%s",FileNameMap,"head1");
            printf("before create dset1:%s\n",buf_name1);
            sprintf (buf_name2,"%s%s",FileNameMap,"head2");
            printf("before create dset2:%s\n",buf_name2);
            sprintf (buf_name3,"%s%s",FileNameMap,"head3");
            printf("before create dset3:%s\n",buf_name3);
            sprintf (buf_name4,"%s%s",FileNameMap,"head4");
            printf("before create dset4:%s\n",buf_name4);
            sprintf (buf_name5,"%s%s",FileNameMap,"head5");
            printf("before create dset5:%s\n",buf_name5);
            sprintf (buf_name6,"%s%s",FileNameMap,"head6");
            printf("before create dset6:%s\n",buf_name6);
        } else if (strcmp(indicMap,"fld") == 0) {
          sprintf (buf_name7,"%s%s",FileNameMap,"fld0");
          sprintf (buf_name8,"%s%s",FileNameMap,"fld1");
          sprintf (buf_name9,"%s%s",FileNameMap,"fld2");
          sprintf (buf_name10,"%s%s",FileNameMap,"fld3");
          sprintf (buf_name11,"%s%s",FileNameMap,"fld4");
        } else if (strcmp(indicMap,"xyz") == 0) {
//          sprintf (buf_name8,"%s%s",FileNameMap,"xyz");  
        }

        if(strcmp(indicMap,"map") == 0) {
        string_dset = H5LTmake_dataset(file_map,buf_name0,1,dimsfb,H5T_NATIVE_FLOAT,buf_map2);
        printf("after create datas0 %p\n",&string_dset);
        string_dset = H5LTmake_dataset(file_map,buf_name1,2,dimsfbb,H5T_NATIVE_CHAR,buf_head[0]);
        printf("after create datas1 %p\n",&string_dset);
        string_dset = H5LTmake_dataset(file_map,buf_name2,2,dimsfbb,H5T_NATIVE_CHAR,buf_head[1]);
        printf("after create datas1 %p\n",&string_dset);
        string_dset = H5LTmake_dataset(file_map,buf_name3,2,dimsfbb,H5T_NATIVE_CHAR,buf_head[2]);
        printf("after create datas1 %p\n",&string_dset);
        string_dset = H5LTmake_dataset(file_map,buf_name4,2,dimsfbb,H5T_NATIVE_CHAR,buf_head[3]);
        printf("after create datas1 %p\n",&string_dset);
        string_dset = H5LTmake_dataset(file_map,buf_name5,2,dimsfbb,H5T_NATIVE_CHAR,buf_head[4]);
        printf("after create datas1 %p\n",&string_dset);
        string_dset = H5LTmake_dataset(file_map,buf_name6,2,dimsfbb,H5T_NATIVE_CHAR,buf_head[5]);
        printf("after create datas1 %p\n",&string_dset);
        } else if (strcmp(indicMap,"fld") == 0) {
        string_dset = H5LTmake_dataset(file_map,buf_name7,2,dimsfcc,H5T_NATIVE_CHAR,buf_fld[0]);
        string_dset = H5LTmake_dataset(file_map,buf_name8,2,dimsfcc,H5T_NATIVE_CHAR,buf_fld[1]);
        string_dset = H5LTmake_dataset(file_map,buf_name9,2,dimsfcc,H5T_NATIVE_CHAR,buf_fld[2]);
        string_dset = H5LTmake_dataset(file_map,buf_name10,2,dimsfcc,H5T_NATIVE_CHAR,buf_fld[3]);
        string_dset = H5LTmake_dataset(file_map,buf_name11,2,dimsfcc,H5T_NATIVE_CHAR,buf_fld[4]);
        
        } else if (strcmp(indicMap,"xyz") == 0) {
//        string_dset = H5LTmake_dataset(file_map,buf_name8,2,dimsfcc,H5T_NATIVE_CHAR,buf_xyz);
        }

        printf("H5FCLOSE1\n");
//      status = H5Fclose(file_map);
    ++i_test;
    }
    status = H5Fclose(file_map2);
    status = H5Fclose(file_map); 
    fclose(pMAP);
    
return 0;
}


