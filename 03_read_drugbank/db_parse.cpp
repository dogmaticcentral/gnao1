
// compile with  g++ -o  tinyx tinyxml2.cpp db_parse.cpp  -O3

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

// not we have elements and nodes here
// (nodes have parents and siblings and such
// and elements have attributes, text itd)
void parseDrug(XMLNode *node) {
    XMLElement *e = node->ToElement();
    if(strcmp(e->Attribute("type"), "small molecule")) return; // otherwise we have all sorts of shit here, including fish
    printf("%s  %s  %s\n", e->Value(), e->Attribute("type"), node->FirstChildElement("name")->GetText());
    for(XMLNode *sib = node->FirstChildElement("synonyms")->FirstChildElement("synonym"); sib; sib = sib->NextSibling()){
          printf("     synonym %s\n", sib->ToElement()->GetText());
    }
    for(XMLNode *sib = node->FirstChildElement("products")->FirstChildElement("product"); sib; sib = sib->NextSibling()){
          printf("     product %s\n", sib->FirstChildElement("name")->ToElement()->GetText());
    }
    for(XMLNode *sib = node->FirstChildElement("targets")->FirstChildElement("target"); sib; sib = sib->NextSibling()){
          printf("     target %s\n",  sib->FirstChildElement("name")->ToElement()->GetText());
          printf("            id %s\n",  sib->FirstChildElement("id")->ToElement()->GetText());
          for(XMLNode *nephew = sib->FirstChildElement("actions")->FirstChildElement("action");
                        nephew; nephew = nephew->NextSibling()){
                printf("            action %s\n", nephew->ToElement()->GetText());
          }
    }

    //exit(1);
}


int main( int argc, const char ** argv ) {


    const char * filename = "/storage/databases/drugbank/drugbank.xml";

 	XMLDocument* doc = new XMLDocument();
	clock_t startTime = clock();
	doc->LoadFile( filename );
 	clock_t loadTime = clock();
	int errorID = doc->ErrorID();


	printf( "Test file '%s' loaded. ErrorID=%d.\n", filename, errorID );
	if ( !errorID ) {
		printf( "Load time   = %u s\n", (unsigned)((loadTime - startTime)/CLOCKS_PER_SEC));
	} else {
	    exit(1);
    }

    XMLHandle docHandle(doc);
    XMLElement *entry = docHandle.FirstChildElement().ToElement();
    if(entry){
        for(XMLNode *node = entry->FirstChildElement(); node; node = node->NextSibling()){
            parseDrug(node);
        }
    } else {
        printf("niente palente\n");
    }

    exit(0);
}
