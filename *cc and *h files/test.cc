#include "hdf5.h"
#include "hdf5_hl.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
// DEFINE HDF5 VARIABLE1
    hid_t file_map,string_dset0;
    hsize_t dimsfb[1],dimsfb2[1];
    herr_t status;

   int MAX=256;
   int LINE_LEN=256;
   FILE *map_file;

// DEFINE HDF5 VARIABLE2
    int FILE_LEN;
    char buf_name0[MAX];
    char buf_name1[MAX];
    int fl_map_tmp=0;
    int fl_map_tmp2=0;
    int fl_map=0;
    int char_line=0; 
    int count_line=0;
    float fl_buf_map=0.0;
    float buf_tmpx=0.0;
    int n_values,nrow,i;

    map_file = fopen("rec.A.map", "r"); 
         
//      COUNT THE NUMBER OF LINE IN THE MAP FILE
        while ((char_line = fgetc(map_file)) != EOF ) {
//               printf("character:%c\n",char_line);               
               if ( char_line == '\n' )
                  {
                  ++count_line; /* Bump the counter for every newline. */
                  }
        }
        rewind(map_file); 
        FILE_LEN=count_line;
        count_line=0;
        printf("map %s : file length:%d\n","FileName",FILE_LEN);
        /* write the map in memory (buffer) */
        printf("write the map in memory (buffer)\n");

        char buf_map[FILE_LEN][LINE_LEN];
        float buf_map2[FILE_LEN];

//      GET EACH LINE OF THE MAP FILE , CONVERT INTO FLOAT, STORE IN buf_map2
        while ((fl_map_tmp < FILE_LEN) && (fgets(buf_map[fl_map_tmp],LINE_LEN,map_file) != NULL)) {
               if(fl_map_tmp >= 6){
               buf_tmpx = atof(buf_map[fl_map_tmp]);
               buf_map2[fl_map_tmp2]=buf_tmpx;
//               printf("line ori:%s\n",buf_map[fl_map_tmp]);
//               printf("line float:%f\n",buf_map2[fl_map_tmp]*2);
               ++fl_map_tmp2;
               }
               ++fl_map_tmp;
               
        }
        fl_map=fl_map_tmp;
        fl_map_tmp=0;
        fclose(map_file);
        
        printf("end writing map in memory (buffer)\n");

//      OPEN THE .h5 FILE
        file_map = H5Fcreate("file.h5", H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);
        printf("after create mapfile %d:\n",file_map);
        /* To create the dataset we: */

        dimsfb[0] = FILE_LEN ;

        sprintf (buf_name0,"%s%s","FileName","ha");
        printf("before create dset0:%s\n",buf_name0);
       
           string_dset0 = H5LTmake_dataset(file_map,buf_name0,1,dimsfb,H5T_NATIVE_FLOAT,buf_map2);
           printf("after create datas0 %p\n",&string_dset0);

// call close functions on all references 
         status = H5Fclose(file_map);
         printf("after close 1\n");
     
          float  data[FILE_LEN];
          char str[LINE_LEN];
          printf("before re-open \n");
            
         /* open file from ex_lite1.c */
         file_map = H5Fopen("file.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

         printf("before read/get DB\n");
         status = H5LTread_dataset_float(file_map,buf_name0,data);
         status = H5LTget_dataset_info(file_map,buf_name0,dimsfb2,NULL,NULL);

         /* print it by rows */
         nrow = (size_t)dimsfb2[0];
         for (i=0; i< nrow; i++ )
         {
          printf ("%f\n", data[i]);
          sprintf(str, "%f", data[i]);
          printf ("%s\n", str);
         }

         /* close file */
         status = H5Fclose (file_map);

         return 0;
}
