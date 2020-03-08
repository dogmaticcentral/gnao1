

// compile with  g++ -o  tinyx tinyxml2.cpp test2.cpp  -O3

#include <stdio.h>
#include "tinyxml2.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
using namespace tinyxml2;
using namespace std;



int main( int argc, const char ** argv ) {

    if (argc<3) {
        fprintf(stderr, "Usage: %s <xml file>  <max_depth>\n", argv[0]);
        exit(1);
    }

 	XMLDocument* doc = new XMLDocument();
	clock_t startTime = clock();
	doc->LoadFile( argv[1] );
 	clock_t loadTime = clock();

	int errorID = doc->ErrorID();
	delete doc; doc = 0;
 	clock_t deleteTime = clock();

	printf( "Test file '%s' loaded. ErrorID=%d\n", argv[1], errorID );
	if ( !errorID ) {
		printf( "Load time   = %u s\n", (unsigned)((loadTime - startTime)/CLOCKS_PER_SEC));
		printf( "Delete time = %u s\n", (unsigned)((deleteTime - loadTime)/CLOCKS_PER_SEC));
		printf( "Total time  = %u s\n", (unsigned)((deleteTime - startTime)/CLOCKS_PER_SEC));
	}
	exit(0);


}
