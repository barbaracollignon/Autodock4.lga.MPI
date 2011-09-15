#include <vector>
#include <string>
#include <sstream>
#include <iostream>


//using namespace std;

int main()
{
    string str("Split me by whitespaces");
    string buf; // Have a buffer string
    char buf1[56]="Split me by whitespaces";
    char buf2[56];
    stringstream ss(str); // Insert the string into a stream

    vector<string> tokens; // Create vector to hold our words

    while (ss >> buf)
        tokens.push_back(buf);
 
    int vsize = tokens.size();
    int n;
    
    for(n=0;n<vsize; n++)
        cout << tokens[n] << "\n" ;

    sprintf(buf2,"%%s\n",buf1);
    printf("buf2:%s\n",buf2);

}
