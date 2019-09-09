//Original Author: Paul Romatschke
//Updates (mostly personal readibility and command line configurations): Chris McGinn

//c and cpp dependencies
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <omp.h>

//Local dependencies
#include <include/paramreader.h>
#include <include/linear_int.h>
#include <include/FS.h>

//GSL dependencies
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//#include <unistd.h>

using namespace std;

int reachedTf=0;
//effective degrees of freedom for Nc=Nf=3 ideal gas
double dof=15.6269;

double **Tref;

//MORE Local dependencies
#include <include/fromvh2.h>
#include <include/fromrelid.h>

//return temperature*lattice spacing
double T(int sx,int sy)
{
  double temp=0;
  temp=e[sx][sy]/dof;
  return pow(temp,0.25);
}

double Tlast(int sx,int sy)
{
  double temp=0;
  temp=elast[sx][sy]/dof;
  return pow(temp,0.25);
}

double eos(double mye)
{
  double temp=0;
  temp=mye/3.;
  return temp;
}

int PBC(int j,int DIM)
{
  int i=j;
  while(i<1){i+=DIM;}
  while(i>DIM){i-=DIM;}
  return i;
}

double f(double px, double py, double pz, double tn, double t0, int sx, int sy, int mu, int nu)
{  
  if(sx == 4 && sy == 5 && mu == 0 && nu == 0 && false){ 
    std::cout << "ParamSet: ";
    std::cout << px << ", " << py << ", " << pz << ", " << tn << ", " << t0 << ", " << sx << ", " << sy << ", " << mu << ", " << nu << std::endl;
  }

  double pperp2 = px*px+py*py;
  double pt = sqrt(pperp2+pz*pz);
  double pt0 = sqrt(pperp2+pz*pz*tn*tn/(t0*t0));
  double pvec[4] = {pt, px, py, pz};  
  
  double xx = sx - px*(tn*pt-t0*pt0)/(pperp2);
  double yy = sy - py*(tn*pt-t0*pt0)/(pperp2);  

  int x0 = PBC(int(floor(xx)), NUMT);
  int x1 = PBC(int(ceil(xx)), NUMT);
  int y0 = PBC(int(floor(yy)), NUMT);
  int y1 = PBC(int(ceil(yy)), NUMT);
  
  double Lambda = bilinear_int(x0,x1, y0, y1, Tref[x0][y0], Tref[x1][y0], Tref[x1][y1], Tref[x0][y1], xx, yy);
  double ff = exp(-pt0/Lambda);
  if(fabs(Lambda)<1.e-16) ff=0;

  if(!isfinite(ff)) printf("problem here %i %i %f\n",sx,sy,ff);

  double res= ff*pvec[mu]*pvec[nu]/pvec[0];
  if(!isfinite(res)) printf("problem here %i %i %f %g\n", sx, sy, res, Lambda);
  if(sx == 4 && sy == 5 && mu == 0 && nu == 0 && false){std::cout << " ANSWER: " << res << std::endl;}
  
  return res;
}

double fshell(double *k, size_t dim, void *params)
{
  struct my_f_params {int sx; int sy; int mu; int nu; double tt; double t0;};
  struct my_f_params *fp = (struct my_f_params *)params;
  
  double A = f(k[0], k[1], k[2], fp->tt, fp->t0, fp->sx, fp->sy, fp->mu, fp->nu);
  return A;
}

double intf(int sx, int sy, double *tab, double tt, double t0, int mcalls, bool doPrint = false)
{
  doPrint=false;
  if(doPrint){
    std::cout << "PRINTING INTF" << std::endl;
    std::cout << "sx, sy, tt, t0, mcalls: " << sx << ", " << sy << ", " << tt << ", " << t0 << ", " << mcalls << std::endl;
    std::cout << "TABULAR: " << std::endl;
    int maxVal = 0;
    for(int bIX = 0; bIX < 4; ++bIX){
      for(int bIY = 0; bIY < 4; ++bIY){
	std::string tempStr = std::to_string(tab[bIX + bIY*4]);
	if(maxVal > tempStr.size()) maxVal = tempStr.size();
      }
    }

    for(int bIX = 0; bIX < 4; ++bIX){
      for(int bIY = 0; bIY < 4; ++bIY){
	std::string tempStr = std::to_string(tab[bIX + bIY*4]);

	std::cout << " " << tempStr;
      }
      std::cout << std::endl;
    }
  }
  
  double res, err;
  double relerr=0.05;

  //ROMATSCHKE ORIG
  double xl[3] = {-5, -5, -5};
  double xu[3] = {5, 5, 5};
  //CFM EDIT 2019.09.04
  //  double xl[3] = { -15, -15, -15 };
  //  double xu[3] = { 15, 15, 15};
  
  const gsl_rng_type *T;
  gsl_rng *r;

  struct my_f_params{ int sx; int sy; int mu; int nu; double tt; double t0;};
  struct my_f_params params;
  params.sx=sx;
  params.sy=sy;
  params.tt=tt;
  params.t0=t0;
  gsl_monte_function G = { &fshell, 3, &params };
  size_t calls = mcalls;

  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_monte_vegas_state *s3 = gsl_monte_vegas_alloc (3);

  //CFM EDIT t03,t13,t23 unused
  //  double t00,t01,t02,t03;
  double t00,t01,t02;
  //  double t11,t12,t13;
  double t11,t12;
  //  double t22,t23;
  double t22;
  double t33;
  
  //size_t totcalls=calls;
  params.mu=0;
  params.nu=0;
  //fairly hard cut-off at around 100 MeV energy scale
  //the regulator can be brought down when simultaneously
  //increasing the number of function evaluations
  double regulator=pow(AT,4)*2.e-2;

  //while (!isfinite(res))

  //  std::cout << "START INTEGRATE (Calls=" << calls << "): " << std::endl;
  gsl_monte_vegas_integrate (&G, xl, xu, 3, calls, r, s3, &res, &err);  
  //  std::cout << "STOP INTEGRATE: " << std::endl;
  double tres=res;
  int mbreak=0;
  
  do{
    gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
    if(std::fabs((res-tres)/res) < relerr) mbreak=1;
    if(res < regulator) mbreak=1;
    tres=res;
  } 
  while (!mbreak);
  
  t00=res;
  if((isfinite(res)) && (res>regulator)){    
    params.mu=0;
    params.nu=1;
    gsl_monte_vegas_init(s3);      
    gsl_monte_vegas_integrate(&G, xl, xu, 3,calls, r, s3, &res, &err);
    mbreak=0;
    tres=res;

    do{
      gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
      if(fabs((res-tres)/res)<relerr) mbreak=1;	 
      tres=res;
    } 
    while (!mbreak);
      
    t01 = res;
    if(!isfinite(res)) t01=0;
      
    params.mu=0;
    params.nu=2;
    gsl_monte_vegas_init(s3);
    gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err);
    mbreak=0;
    tres=res;

    do{
      gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
      if(std::fabs((res-tres)/res) < relerr) mbreak=1;	 
      tres=res;
    } 
    while (!mbreak);

    t02=res;
    if(!isfinite(res)) t02=0;
      
    params.mu=1;
    params.nu=1;
    gsl_monte_vegas_init(s3);
    gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
    mbreak=0;
    tres=res;

    do{
      gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
      if(std::fabs((res-tres)/res) < relerr) mbreak=1;	 
      tres=res;
    } 
    while (!mbreak);

    t11=res;
    if (!isfinite(res))
      t11=0;
    
    params.mu=1;
    params.nu=2;
    gsl_monte_vegas_init(s3);
    gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err);
    mbreak=0;
    tres=res;

    do{
      gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
      if(fabs((res-tres)/res)<relerr) mbreak=1;	 
      tres=res;
    } 
    while (!mbreak);

    t12=res;
    if(!isfinite(res)) t12=0;

    
    params.mu=2;
    params.nu=2;
    gsl_monte_vegas_init(s3);
    gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err);
    mbreak=0;
    tres=res;

    do{
      gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
      if(fabs((res-tres)/res)<relerr) mbreak=1;	 
      tres=res;
    } 
    while(!mbreak);

    t22=res;
    if(!isfinite(res)) t22=0;
      
    params.mu=3;
    params.nu=3;
    gsl_monte_vegas_init(s3);
    gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err);
    mbreak=0;
    tres=res;
    do{
      gsl_monte_vegas_integrate (&G, xl, xu, 3,calls, r, s3, &res, &err); 
      if (fabs((res-tres)/res)<relerr) mbreak=1;	 
      tres=res;
    } 
    while (!mbreak);


    t33=res;
    if (!isfinite(res)) t33=0;
  }
  else{
    //assume this is a low energy state;
      
    t00=0;
    t11=0;
    t22=0;
    t33=0;
    t01=0;
    t02=0;
    t12=0;
  }
  gsl_monte_vegas_free(s3); 
  gsl_rng_free(r);

  //  std::cout << "TOO: " << t00 << std::endl;
  if(t00 < regulator){
    //    std::cout << "SETTING TAB FROM BELOW REGULATOR" << std::endl;
    
    tab[0]=regulator;
    tab[1]=0;
    tab[2]=0;
    tab[3]=0.;
    tab[4]=0;
    tab[5]=-regulator/3.;
    tab[6]=0;
    tab[7]=0.;
    tab[8]=0;
    tab[9]=0;
    tab[10]=-regulator/3.;
    tab[11]=0;
    tab[12]=0;
    tab[13]=0;
    tab[14]=0;
    tab[15]=-regulator/3.;
  }
  else{
    tab[0]=t00;
    tab[1]=-t01;
    tab[2]=-t02;
    tab[3]=0.;
    tab[4]=t01;
    tab[5]=-t11;
    tab[6]=-t12;
    tab[7]=0.;
    tab[8]=t02;
    tab[9]=-t12;
    tab[10]=-t22;
    tab[11]=0;
    tab[12]=0;
    tab[13]=0;
    tab[14]=0;
    tab[15]=-t33;

    if(doPrint){
      std::cout << "ABOVE REGULATOR" << std::endl;
      std::cout << " NEW MATRIX: " << std::endl;
      
      int maxVal = 0;
      for(int bIX = 0; bIX < 4; ++bIX){
	for(int bIY = 0; bIY < 4; ++bIY){
	  std::string tempStr = std::to_string(tab[bIX + bIY*4]);
	  if(maxVal > tempStr.size()) maxVal = tempStr.size();
	}
      }
      
      for(int bIX = 0; bIX < 4; ++bIX){
	for(int bIY = 0; bIY < 4; ++bIY){
	  std::string tempStr = std::to_string(tab[bIX + bIY*4]);	  
	  std::cout << " " << tempStr;
	}
	std::cout << std::endl;
      }
    }
    
  }

  //  double tab[16];
  //*tab={t00,t01,t02,0.,t01,t11,t12,0,t02,t12,t22,0,0,0,0,t33};
  //printNxNmatrix(tab,4,1);

  //CFM EDIT ADDED return type
  return -1;
}

void comppi(double *tab,int sx, int sy)
{
  pib[sx][sy]=tab[0]+tab[5]+tab[10]+tab[15]; //trace of pi minus trace anomaly (which is 0 here)

  pixx[sx][sy]=-tab[5]-(e[sx][sy]*4./3.-pib[sx][sy])*u[0][sx][sy]*u[0][sx][sy]-(e[sx][sy]/3.-pib[sx][sy]);
  piyy[sx][sy]=-tab[10]-(e[sx][sy]*4./3.-pib[sx][sy])*u[1][sx][sy]*u[1][sx][sy]-(e[sx][sy]/3.-pib[sx][sy]);
  pixy[sx][sy]=-tab[6]-(e[sx][sy]*4./3.-pib[sx][sy])*u[0][sx][sy]*u[1][sx][sy];
}

void updatevars(double tt)
{
  long int position;
  int sx,sy;
  double tab[16];
  double stat[4];
  int nthreads,tid;
  double gamma;

  /*
  sx=Middle;
  sy=Middle;
  f(0.5,0.5,0.5,0.25,0.25,sx,sy,0,0);

  //usleep(10000000);

  intf(sx,sy,tab,tt,tt);
  printNxNmatrix(tab,4,1);
  calcedvquick(tab,stat,0);
  e[sx][sy]=fabs(stat[0]/(24*M_PI)*dof);
  u[0][sx][sy]=stat[1];
  u[1][sx][sy]=stat[2];
  
  printf("eps=%f T=%f ux=%f uy=%f\n",e[sx][sy],T(sx,sy),u[0][sx][sy],u[1][sx][sy]);*/
  /*
  sx=90;
  sy=100;
  intf(sx,sy,tab,7.495912,2.495912,2000);
  printNxNmatrix(tab,4,1);
  intf(sx,sy,tab,7.495912,2.495912,2000);
  printNxNmatrix(tab,4,1);
  intf(sx,sy,tab,7.495912,2.495912,8000);
  printNxNmatrix(tab,4,1);
  */


#pragma omp parallel shared(nthreads,tid,u,e,position) private(sx,sy,gamma)
  {
    tid = omp_get_thread_num();
    if (tid == 0)
      nthreads = omp_get_num_threads();

    //    std::cout << "NUMBER OF THREADS: " << tid << ", " << nthreads << std::endl;

    bool notPrinted = true;
    //    std::cout << "PRE UPDATE: " << std::endl;
#pragma omp for schedule(dynamic,1)
    for(position=0;position<NUMT*NUMT;position++)
    {
      sx=position%NUMT+1;
      sy=position/NUMT+1;
      //      std::cout << "INTF Call " << sx << ", " << sy << " (tt, TINIT, fmtoGeV, AT, TINIT*fmtoGeV/AT): " << tt << ", " << TINIT << ", " << fmtoGeV << ", " << AT << ", " << TINIT*fmtoGeV/AT << std::endl;      
      intf(sx, sy, tab, tt, TINIT*fmtoGeV/AT, 40000, sx == 4 && sy == 5);
      calcedvquick(tab,stat,0, sx == 4 && sy == 5);

      if(e[sx][sy] > 0 && notPrinted){
	//	std::cout << "FIRST NON-ZERO: " << e[sx][sy] << " AT SX,SY: " << sx << ", " << sy << std::endl;
	notPrinted = false;
      }
      
      //      std::cout << e[sx][sy] << " ";
      //      std::cout << "REPLACING " << e[sx][sy] << " WITH " << fabs(stat[0]) << std::endl;
      e[sx][sy]=fabs(stat[0]);
      gamma=1/sqrt(1-(stat[1]*stat[1]+stat[2]*stat[2]));
      u[0][sx][sy]=gamma*stat[1];
      u[1][sx][sy]=gamma*stat[2];
      comppi(tab,sx,sy);

      e[sx][sy]*=dof*pow(AT,4)/(24*M_PI);
      pib[sx][sy]*=dof*pow(AT,4)/(24*M_PI);
      pixx[sx][sy]*=dof*pow(AT,4)/(24*M_PI);
      pixy[sx][sy]*=dof*pow(AT,4)/(24*M_PI);
      piyy[sx][sy]*=dof*pow(AT,4)/(24*M_PI);
    }
  }//end of parallel session


  double tmax=0;
  long int mpos=0;
  for(position=0;position<NUMT*NUMT;position++)
    {
      sx=position%NUMT+1;
      sy=position/NUMT+1;
      if (T(sx,sy)>tmax)
	{
	  tmax=T(sx,sy);
	  mpos=position;
	}
    }
  /*
  sx=mpos%NUMT+1;
  sy=mpos/NUMT+1;
  printf("max T=%f at sx=%i sy=%i tt=%f t0=%f\n",tmax/AT,sx,sy,tt,TINIT*fmtoGeV/AT);
  intf(sx,sy,tab,tt,TINIT*fmtoGeV/AT,2000);
  printNxNmatrix(tab,4,1);
  intf(sx,sy,tab,tt,TINIT*fmtoGeV/AT,4000);
  printNxNmatrix(tab,4,1);
   intf(sx,sy,tab,tt,TINIT*fmtoGeV/AT,8000);
   printNxNmatrix(tab,4,1);*/
}

//put free-streamed energy density etc into files
void outputnewinitialconditions(double tt, std::string appendToName)
{
  while(appendToName.find(".") != std::string::npos){appendToName.replace(appendToName.find("."), 1, "p");}

  char cbuffer[255];

  fstream inited,initux,inituy,initpixx,initpixy,initpiyy,itime,initpi,initxyed,initeduxuy;
  sprintf(cbuffer,"output/time-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  itime.open(cbuffer,ios::out);
  sprintf(cbuffer,"output/inited-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  inited.open(cbuffer, ios::out);
  sprintf(cbuffer,"output/initux-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initux.open(cbuffer, ios::out);
  sprintf(cbuffer,"output/inituy-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  inituy.open(cbuffer, ios::out);
  sprintf(cbuffer,"output/initpixx-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initpixx.open(cbuffer, ios::out);
  sprintf(cbuffer,"output/initpixy-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initpixy.open(cbuffer, ios::out);
  sprintf(cbuffer,"output/initpiyy-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initpiyy.open(cbuffer, ios::out);
  sprintf(cbuffer,"output/initpi-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initpi.open(cbuffer, ios::out);

  //CFM MOD OUTPUT
  sprintf(cbuffer,"output/initxyed-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initxyed.open(cbuffer, ios::out);

  //REMOVE LATTICE N FOR READABLE OUTPUT
  //  initxyed << "LatticeN\tDummy\tDummy\ted\n";
  initxyed << "Dummy\tDummy\ted\n";
  
  sprintf(cbuffer,"output/initeduxuy-%.5f_%s.dat",t/fmtoGeV*AT, appendToName.c_str());
  initeduxuy.open(cbuffer, ios::out);
  //  initeduxuy << "LatticeN\tDummy\tux\tuy\tDummy\tDummy\tDummy\n";

  
  double spacingFM = AT/fmtoGeV;
  double start = spacingFM/2. - 1*(spacingFM*NUMT)/2.;
  
  //  std::cout << "ENERGY CHECK (SCAL " << SCAL << "): " << std::endl;
  unsigned int counter=0;
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	//	std::cout << e[sx][sy]/SCAL <<"\t";

	inited << e[sx][sy]/SCAL <<"\t";
	initux << u[0][sx][sy] <<"\t";
	inituy << u[1][sx][sy] <<"\t";
	initpixx << pixx[sx][sy] <<"\t";
	initpixy << pixy[sx][sy] <<"\t";
	initpiyy << piyy[sx][sy] <<"\t";
	initpi << pib[sx][sy] <<"\t";

	double xfm = start + spacingFM*(sx - 1);
	double yfm = start + spacingFM*(sy - 1);
	//Lattice sites
	initxyed << xfm << "\t" << yfm << "\t" << e[sx][sy]/SCAL << "\n";
	initeduxuy << e[sx][sy]/SCAL << "\t" << u[0][sx][sy] << "\t" << u[1][sx][sy] << "\t0\t0\t0\n";
	++counter;
      }

  //  std::cout << std::endl;
  
  inited<<endl;
  initux<<endl;
  inituy<<endl;
  initpixx<<endl;
  initpixy<<endl;
  initpi<<endl;
  initpiyy<<endl;
  
  itime << t*AT/fmtoGeV;

  inited.close();
  initux.close();
  inituy.close();
  initpixx.close();
  initpixy.close();
  initpiyy.close();
  initpi.close();

  initxyed.close();
  initeduxuy.close();
}


int main(int argc, char* argv[])
{
  char *paramsfile;
  char cbuffer[1255];
  char outdir[255];
  int ABORT = 0;
  
  if(argc>1) paramsfile=argv[1];
  else{
    printf("");
    printf("ERROR please specify a params file name");
    return 1;
  }

  fstream parameters;
  parameters.open(paramsfile, ios::in);
  if(parameters.is_open()){      
    extern void readParameters(const char*); 
    printf("This is FS 1.0 \n");
    readParameters(paramsfile);
  }
  else{
    printf("\nERROR: params file %s could not be opened\n",paramsfile);
    ABORT=1;
  }
  parameters.close();

  int outputIter = 0;
  
  if(!ABORT){
    // sprintf(cbuffer,"%i-%i-%i-A%.5f-tini%.2f-dof%.1f-dt%.4f-ETA%.4f-%s-O%i",NX,NY,NZ,AT,TINIT,dof,dt,ETAOS,desc,ORDER);
    //outdir=cbuffer;

    printf("===> Info: Setting data directory to data ...");
    sprintf(outdir,"data");
    printf("done.\n");
    
    //outdir="data_0000";
    printf("===> Info: allocating memory...");
    allocateMemory();
    Tref = new double*[NUMT+2];
    for(int i = 0; i < NUMT+2; ++i){ 
      Tref[i] = new double[NUMT+2];
    }
    printf("done.\n");
    
    printf("===> Info: setting initial conditions...");
    setInitialConditions();
    //make ed physical
    //      std::cout << "POST INIT MOD, divide all E by pow(AT,4), AT=" << AT << ": " << std::endl;
    for(int sx = 0; sx <= NUMT+1; ++sx){
      for(int sy = 0; sy <= NUMT+1; ++sy){
	e[sx][sy]/=pow(AT,4);
      }
    }
      
      
    for(int sx = 0;sx <= NUMT+1; ++sx){
      for(int sy = 0;sy <= NUMT+1; ++sy){
	Tref[sx][sy]=T(sx,sy);
      }
    }

    printf("done.\n");
      
    //printf("test: %g\n",e[Middle][Middle]/pow(AT,4));
    //printf("Test: central temperature %g\n",T(Middle,Middle)/AT);
    
    //generate paramter file for hadronic afterburner
    //printf("===> Info: generatin hadronic parameter files...");
    //int suc=generatehadronparameters();
    //if (suc)
    //printf("done.\n");
    //else
    //	printf("FAILED!!!\n");
    
    sprintf(cbuffer,"%s/freezeout_bulk.dat",outdir);
    freeze_out.open(cbuffer, ios::out);
    sprintf(cbuffer,"%s/meta.dat",outdir);
    meta.open(cbuffer, ios::out);
    meta << "# tau [fm]\t" << "T [GeV]\t" << "epsilon [GeV4]\t" << "Phi[GeV4]\t" << "\n";
    sprintf(cbuffer,"%s/ecc.dat",outdir);
    ecces.open(cbuffer, ios::out);
    ecces << "#tau\t e_x\t e_p\n";
    
    printf("     Time[fm/c]   T_cent [GeV] \t Ecc(spatial)\t Ecc(momentum)\n");
    
    double eps = EPS;
    long int i=0;
      
    //main loop here
    bool noPrint = true;
    while((reachedTf==0)&&(!ABORT)) {
      ++i;
      //      printf("i=%i done\n",i);
      if ((i-1)%UPDATE==0) 
	{

	  {
	    //Copy fields to memory to compare with next UPDATE time step
	    copyUPDATE();
	  }

	  updatevars(t);
	      //printf("checking ed=%f\n",e[Middle][Middle]);
	  outputMeasurements(t, outputIter);
	  ++outputIter;
	  outputnewinitialconditions(t, argv[2]);
	  //	  noPrint=false;
	  if(!noPrint){
	    std::cout << "TERMINATING ON SECOND STEP" << std::endl;
	    return 1;
	  }
	  
	  //CFM EDIT
	  double time = t/fmtoGeV*AT;
	  const double cutoff = 100;
	  if(time >= cutoff){
	    std::cout << "We are terminating at " << cutoff << " fm/c. return 1" << std::endl;
	    std::cout << "Done." << std::endl;
	    return 1;
	  }	      
	}
      
      if ((i-1)%SNAPUPDATE==0) 
	{      
	  snapshot(t); 
	}
      
      t += eps;
      // if (i==UPDATE+1)
      //  ABORT=1;
    }
    
    freeze_out.close();
    meta.close();
    ecces.close();
    
    if (!ABORT) cout << "Done.\n";
    else cout << "Aborted because encountered nan.\n";
    
    //if necessary, show system vh2 is done
    int stat=system("cp data/params.txt logdir/vh2-is-done.log");
  }
}
