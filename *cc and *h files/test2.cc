#include <stdio.h>
       
int main( void )
{
    const char *src = "this is a test";
    char word[5];
    int n;
    while ( sscanf ( src, "%6s%n", word, &n ) == 1 ) {
       puts ( word );
       src += n;
    }
    return 0;
}
