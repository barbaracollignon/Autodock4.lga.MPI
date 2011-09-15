#include <stdio.h>
       
int main( void )
{
    char *src = "test 2 2.9999 ligand_types X C D E T";
    char word[12];
    int n;
    int b=0;
    while ((sscanf ( src, "%s%n", word, &n ) == 1)) {
       if (b > 2) {
         puts ( word );
       }
      src += n;
      b+=1;
    }
    return 0;
}
