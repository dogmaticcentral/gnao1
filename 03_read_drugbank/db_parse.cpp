
// compile with  g++ -o  readb tinyxml2.cpp db_parse.cpp  -O3

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
void checkProdrug(XMLNode *node) {
   // is this a prodrug? this way of looking for prodrugs does not guarantee anything,
    // but drugbank itslef has only 14 names listed as prodrugs - this was for example
    // you cannot find that levodopa is prodrug of dopamine
    int prodrug = 0;
    for(XMLNode *sib = node->FirstChildElement("categories")->FirstChildElement("category"); sib; sib = sib->NextSibling()){
          for(XMLNode *nephew = sib->FirstChildElement(); nephew; nephew = nephew->NextSibling()){
                if ( nephew->ToElement()->GetText() && !strcmp(nephew->ToElement()->GetText(), "Prodrugs")) {
                    prodrug  = 1;
                }
          }
    }
    if (! prodrug) {
        const char * text = node->FirstChildElement("description")->GetText();
        if (text) {
            if (strstr(text, "prodrug") || strstr(text, "Prodrug") || strstr(text, "pro-drug"))  prodrug  = 1;
        }
    }
    if (! prodrug) {
        const char * text = node->FirstChildElement("pharmacodynamics")->GetText();
        if (text) {
            if (strstr(text, "prodrug") || strstr(text, "Prodrug") || strstr(text, "pro-drug"))  prodrug  = 1;
        }
    }
     if (! prodrug && node->FirstChildElement("mechanism-of-action") && node->FirstChildElement("mechanism-of-action")->GetText()) {
        const char * text = node->FirstChildElement("mechanism-of-action")->GetText();
        if (text) {
            if (strstr(text, "prodrug") || strstr(text, "Prodrug") || strstr(text, "pro-drug"))  prodrug  = 1;
        }
    }

    if (prodrug)  printf("     prodrug\tmaybe\n");

}

void parseDrug(XMLNode *node) {
    XMLElement *e = node->ToElement();
    // if(strcmp(e->Attribute("type"), "small molecule")) return; // otherwise we have all sorts of shit here, including fish
    //if ( strcmp(node->FirstChildElement("name")->GetText(), "Levodopa")) return;
    // printf("%s\t%s  %d\n", e->Value(),  node->FirstChildElement("name")->GetText(), node->GetLineNum());

    printf("%s\t%s\n", e->Value(),  node->FirstChildElement("name")->GetText());

    checkProdrug(node);

    for(XMLNode *sib = node->FirstChildElement("external-identifiers")->FirstChildElement("external-identifier"); sib; sib = sib->NextSibling()){
          XMLNode *nephew = sib->FirstChildElement("resource");
          if (!nephew) continue;
          if ( strcmp(nephew->ToElement()->GetText(), "PubChem Compound")) continue;
          printf("     pubchem\t%s\n", nephew->NextSibling()->ToElement()->GetText());
          break;
    }


    for(XMLNode *sib = node->FirstChildElement("synonyms")->FirstChildElement("synonym"); sib; sib = sib->NextSibling()){
          printf("     synonym\t%s\n", sib->ToElement()->GetText());
    }
    for(XMLNode *sib = node->FirstChildElement("products")->FirstChildElement("product"); sib; sib = sib->NextSibling()){
          printf("     product\t%s\n", sib->FirstChildElement("name")->ToElement()->GetText());
    }
    for(XMLNode *sib = node->FirstChildElement("international-brands")->FirstChildElement("international-brand"); sib; sib = sib->NextSibling()){
          printf("     brand\t%s\n", sib->FirstChildElement("name")->ToElement()->GetText());
    }
    for(XMLNode *sib = node->FirstChildElement("targets")->FirstChildElement("target"); sib; sib = sib->NextSibling()){
          printf("     target\t%s\n",  sib->FirstChildElement("name")->ToElement()->GetText());
          printf("            identifier\t%s\n",  sib->FirstChildElement("id")->ToElement()->GetText());
          if (  sib->FirstChildElement("polypeptide")
                && sib->FirstChildElement("polypeptide")->FirstChildElement("gene-name")
                && sib->FirstChildElement("polypeptide")->FirstChildElement("gene-name")->GetText() ) {
                printf("            gene-name\t%s\n",   sib->FirstChildElement("polypeptide")->FirstChildElement("gene-name")->GetText() );
          }
          for(XMLNode *nephew = sib->FirstChildElement("actions")->FirstChildElement("action");
                        nephew; nephew = nephew->NextSibling()){
                if ( nephew->ToElement()->GetText()) {
                    printf("            action\t%s\n", nephew->ToElement()->GetText());
                }
          }
    }


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
