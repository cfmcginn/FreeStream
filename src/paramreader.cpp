#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>


// external vars defined in UVH2+1.cpp which are loaded here
extern int    NUMT,UPDATE,SNAPUPDATE;
extern long int STEPS;
extern double AT,EPS,TINIT,ETAOS,TF,TSTART,COEFF;
extern double B,L1COEF,L2COEF;
extern double IC;
extern int PTASIZE,PHIPASIZE;
extern int FREEZE;
extern char EOSNAME[255];
extern char ETANAME[255];
extern char ZETANAME[255];
extern char BETANAME[255];
extern char LAMBDANAME[255];
extern double PTMAX, TRIEPS, TRIANGLE, QUADEPS, QUADANGLE;
extern double QUINTEPS, QUINTANGLE, SEXEPS, SEXANGLE, SEPTEPS, SEPTANGLE;
extern double BIEPS, BIANGLE, MONOEPS, MONOANGLE;
extern int FULL;
extern int PCE, NS, IFLOW;
extern double RNUC, ANUC, SIGMANN, TANORM;
extern double SMOOTH;
extern int SMOOTHING;
extern long int B3DEVENTS;
extern int preeqflow;
extern double SCAL;


using namespace std;

// this workhorse examines a key to see if corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(char *key, char *value) {
	// integer params
	if (strcmp(key,"NUMT")==0) NUMT=atoi(value);
	if (strcmp(key,"B")==0) B=atof(value);
	if (strcmp(key,"L1COEF")==0) L1COEF=atof(value);
	if (strcmp(key,"L2COEF")==0) L2COEF=atof(value);
	if (strcmp(key,"TINIT")==0) TINIT=atof(value);
	if (strcmp(key,"STEPS")==0) STEPS=atoi(value);
	if (strcmp(key,"B3DEVENTS")==0) B3DEVENTS=atoi(value);
	if (strcmp(key,"UPDATE")==0) UPDATE=atoi(value);
	if (strcmp(key,"SNAPUPDATE")==0) SNAPUPDATE=atoi(value);
	if (strcmp(key,"AT")==0) AT=atof(value);
	if (strcmp(key,"EPS")==0) EPS=atof(value);
	if (strcmp(key,"ETAOS")==0) ETAOS=atof(value);
	if (strcmp(key,"COEFF")==0) COEFF=atof(value);
	if (strcmp(key,"TF")==0) TF=atof(value);
	if (strcmp(key,"TSTART")==0) TSTART=atof(value);
	if (strcmp(key,"IC")==0) IC=atof(value);
	if (strcmp(key,"PTASIZE")==0) PTASIZE=atoi(value);
	if (strcmp(key,"PHIPASIZE")==0) PHIPASIZE=atoi(value);
	if (strcmp(key,"FREEZE")==0) FREEZE=atoi(value);
	if (strcmp(key,"EOSNAME")==0) strcpy(EOSNAME,value);
	if (strcmp(key,"ETANAME")==0) strcpy(ETANAME,value);
	if (strcmp(key,"ZETANAME")==0) strcpy(ZETANAME,value);
	if (strcmp(key,"BETANAME")==0) strcpy(BETANAME,value);
	if (strcmp(key,"LAMBDANAME")==0) strcpy(LAMBDANAME,value);
	if (strcmp(key,"PTMAX")==0) PTMAX=atof(value);
	if (strcmp(key,"MONOEPS")==0) MONOEPS=atof(value);
	if (strcmp(key,"MONOANGLE")==0) MONOANGLE=atof(value);
	if (strcmp(key,"BIEPS")==0) BIEPS=atof(value);
	if (strcmp(key,"BIANGLE")==0) BIANGLE=atof(value);
	if (strcmp(key,"TRIEPS")==0) TRIEPS=atof(value);
	if (strcmp(key,"TRIANGLE")==0) TRIANGLE=atof(value);
	if (strcmp(key,"QUADEPS")==0) QUADEPS=atof(value);
	if (strcmp(key,"QUADANGLE")==0) QUADANGLE=atof(value);
	if (strcmp(key,"QUINTEPS")==0) QUINTEPS=atof(value);
	if (strcmp(key,"QUINTANGLE")==0) QUINTANGLE=atof(value);
	if (strcmp(key,"SEXEPS")==0) SEXEPS=atof(value);
	if (strcmp(key,"SEXANGLE")==0) SEXANGLE=atof(value);
	if (strcmp(key,"SEPTEPS")==0) SEPTEPS=atof(value);
	if (strcmp(key,"SEPTANGLE")==0) SEPTANGLE=atof(value);
	if (strcmp(key,"FULL")==0) FULL=atoi(value);
	if (strcmp(key,"PCE")==0) PCE=atoi(value);
	if (strcmp(key,"NS")==0) NS=atoi(value);
	if (strcmp(key,"IFLOW")==0) IFLOW=atoi(value);
	if (strcmp(key,"RNUC")==0) RNUC=atof(value);
	if (strcmp(key,"ANUC")==0) ANUC=atof(value);
	if (strcmp(key,"SIGMANN")==0) SIGMANN=atof(value);
	if (strcmp(key,"TANORM")==0) TANORM=atof(value);
	if (strcmp(key,"SMOOTHING")==0) SMOOTHING=atoi(value);
	if (strcmp(key,"SMOOTH")==0) SMOOTH=atof(value);
	if (strcmp(key,"SCAL")==0) SCAL=atof(value);
	if (strcmp(key,"preeqflow")==0) preeqflow=atoi(value);
	return;
}

//
// This routine assumes that paramters are in text file with
// each parameter on a new line in the format 
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
// 
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParameters(const char *filename) {
		
	string commentmarker = "//"; 
	char space = ' '; 
	char tab = '\t';

	int maxline = 128; // maximum line length used in the buffer for reading
	char buffer[maxline];
	ifstream paramFile(filename);
	
	while(!paramFile.eof()) {
		paramFile.getline(buffer,maxline,'\n');
		string line = buffer; int length = strlen(buffer);
		if (line.substr(0,commentmarker.length())!=commentmarker && line.length()>0) {
			char key[32]="",value[32]=""; int founddelim=0;
			for (int i=0;i<length;i++) {
				if (buffer[i]==space || buffer[i]==tab) founddelim=1;
				else {
					if (founddelim==0) key[strlen(key)] = buffer[i];
					else value[strlen(value)] = buffer[i];
				}
			}
			if (strlen(key)>0 && strlen(value)>0) {
				setParameter(key,value);
				cout << key << " = " << value << endl;
			}
		}
	}
	
	return;	
}

