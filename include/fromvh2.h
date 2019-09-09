const double fmtoGeV=5.0677;



// these hold the current values
double ***u,**e,**pixy,**pixx,**piyy,**pib;

// these hold the updated values
double ***U,**E,**Pixy,**Pixx,**Piyy,**Pib;

// these hold the values from the last UPDATE time step
double ***ulast,**elast,**pixylast,**pixxlast,**piyylast,**pilast;

//overall time
double t = 0;

//these are global for convenience; used in doInc
double ****dtmat;

double ****vec;/*bei Pauli rhs*/

//center of lattice 
int Middle;

//for convenience 
//global definition of ut(sx,sy) -- in order not to have it calculated
//a gazillion times
double **globut;

double ***thf;
double ***dtpixx,***dtpixy,***dtpiyy,***dtpi;
double **mypixt,**mypiyt,**mypitt,**mypiee;

// output files
fstream freeze_out;
fstream meta;
fstream ecces;

// initialize global arrays
void allocateMemory() {

  //cout << "==> Allocating memory\n";


 u = new double**[2];

 for (int i=0;i<2;i++) 
   u[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     u[i][j] = new double[NUMT+2];

	 
 e = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
    e[i] = new double[NUMT+2];
 
 pixy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixy[i] = new double[NUMT+2];


 pixx = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixx[i] = new double[NUMT+2];

 piyy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   piyy[i] = new double[NUMT+2];

 pib = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pib[i] = new double[NUMT+2];

//////////////////////
 ulast = new double**[2];

 for (int i=0;i<2;i++) 
   ulast[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     ulast[i][j] = new double[NUMT+2];

	 
 elast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
    elast[i] = new double[NUMT+2];
 
 pixylast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixylast[i] = new double[NUMT+2];


 pixxlast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixxlast[i] = new double[NUMT+2];

 piyylast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   piyylast[i] = new double[NUMT+2];

 pilast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pilast[i] = new double[NUMT+2];
//////////////////////

 U = new double**[2];

 for (int i=0;i<2;i++) 
   U[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     U[i][j] = new double[NUMT+2];

	 
 E = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   E[i] = new double[NUMT+2];


 Pixy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pixy[i] = new double[NUMT+2];


 Pixx = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pixx[i] = new double[NUMT+2];


 Piyy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Piyy[i] = new double[NUMT+2];

 Pib = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pib[i] = new double[NUMT+2];

 dtmat= new double***[NUMT+2];
 vec=new double***[NUMT+2];

 for (int i=0;i<NUMT+2;i++)
   { 
     dtmat[i]=new double**[NUMT+2];
     vec[i]=new double**[NUMT+2];
     for (int j=0;j<NUMT+2;j++) 
       {
	 dtmat[i][j]=new double*[3];
	 vec[i][j]=new double*[3];
	 
	 for (int k=0;k<3;k++) 
	   {
	     dtmat[i][j][k] = new double[3];
	     vec[i][j][k] = new double[3];
	   }
	 //dtmat = new double*[3];
	 //for (int i=0;i<3;i++) dtmat[i] = new double[3];
	 // vec = new double*[3];
	 //for (int i=0;i<3;i++) vec[i] = new double[1];
       }}
  
 globut = new double*[NUMT+2];
 thf = new double**[NUMT+2];
 dtpixx = new double**[NUMT+2];
 dtpixy = new double**[NUMT+2];
 dtpiyy = new double**[NUMT+2];
 dtpi   = new double**[NUMT+2];
 mypixt = new double*[NUMT+2];
 mypiyt = new double*[NUMT+2];
 mypitt = new double*[NUMT+2];
 mypiee = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++)
   { 
     globut[i] = new double[NUMT+2];
     mypixt[i] = new double[NUMT+2];
     mypiyt[i] = new double[NUMT+2];
     mypitt[i] = new double[NUMT+2];
     mypiee[i] = new double[NUMT+2];
     thf[i] = new double*[NUMT+2];
     dtpixx[i] = new double*[NUMT+2];
     dtpixy[i] = new double*[NUMT+2];
     dtpi[i]   = new double*[NUMT+2];
     dtpiyy[i] = new double*[NUMT+2];
     for (int j=0;j<NUMT+2;j++) 
       {
	 thf[i][j] = new double[4];
	 dtpixx[i][j] = new double[4];
	 dtpixy[i][j] = new double[4];
	 dtpi[i][j] = new double[4];
	 dtpiyy[i][j] = new double[4];
       }
   }

}


void smeare(double **pe)
{
  // smooth small densities -MPM
   //======================================================
      // j.nagle - 11/01/2013 - copied in McCumber's smoothing
      // for lumpy condition hydro to run 
      // Set to -1.0 to turn off ....
      // p.romatschke - 10.9.2014 modified to run larger systems
      //double SMOOTH = 0.001;  
      // turn smooth on and off via the params file. -MPM

  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	if(pe[sx][sy]/pow(AT,4)<SMOOTH)
	  {
	    int sx_up = sx+1;
	    int sx_dn = sx-1;
	    int sy_up = sy+1;
	    int sy_dn = sy-1;
	    
	    if (sx_up > NUMT) sx_up = NUMT;
	    if (sx_dn < 1)   sx_dn = 1;
	    if (sy_up > NUMT) sy_up = NUMT;
	    if (sy_dn < 1)   sy_dn = 1;
	    
	    pe[sx][sy] = (pe[sx_up][sy_up]+pe[sx_up][sy]+pe[sx_up][sy_dn]+pe[sx][sy_up]+pe[sx][sy]+pe[sx][sy_dn]+pe[sx_dn][sy_up]+pe[sx_dn][sy]+pe[sx_dn][sy_dn])/9.0;
	  }
      }
}

void smearu(double ***pu, double **pe)
{
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	if(pe[sx][sy]/pow(AT,4)<SMOOTH)
	  {
	    int sx_up = sx+1;
	    int sx_dn = sx-1;
	    int sy_up = sy+1;
	    int sy_dn = sy-1;
	    if (sx_up > NUMT) sx_up = NUMT;
	    if (sx_dn < 1)   sx_dn = 1;
	    if (sy_up > NUMT) sy_up = NUMT;
	    if (sy_dn < 1)   sy_dn = 1;
	    
	    pu[0][sx][sy] = (pu[0][sx_up][sy_up]+pu[0][sx_up][sy]+pu[0][sx_up][sy_dn]+pu[0][sx][sy_up]+pu[0][sx][sy]+pu[0][sx][sy_dn]+pu[0][sx_dn][sy_up]+pu[0][sx_dn][sy]+pu[0][sx_dn][sy_dn])/9.0;
	    pu[1][sx][sy] = (pu[1][sx_up][sy_up]+pu[0][sx_up][sy]+pu[1][sx_up][sy_dn]+pu[1][sx][sy_up]+pu[1][sx][sy]+pu[1][sx][sy_dn]+pu[1][sx_dn][sy_up]+pu[1][sx_dn][sy]+pu[1][sx_dn][sy_dn])/9.0;
	  }
      }
}

//enforce periodic BoundaryConditions

void enforcePBCs()
{
  
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[0][0][sy]=u[0][NUMT][sy];
      u[0][NUMT+1][sy]=u[0][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[0][sx][0]=u[0][sx][NUMT];
      u[0][sx][NUMT+1]=u[0][sx][1];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[1][0][sy]=u[1][NUMT][sy];
      u[1][NUMT+1][sy]=u[1][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[1][sx][0]=u[1][sx][NUMT];
      u[1][sx][NUMT+1]=u[1][sx][1];
   }
  for(int sy=1;sy<=NUMT;sy++)
    {
      e[0][sy]=e[NUMT][sy];
      e[NUMT+1][sy]=e[1][sy]; 
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      e[sx][0]=e[sx][NUMT];
      e[sx][NUMT+1]=e[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      pixx[0][sy]=pixx[NUMT][sy];
      pixx[NUMT+1][sy]=pixx[1][sy];     
    }

  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixx[sx][0]=pixx[sx][NUMT];
      pixx[sx][NUMT+1]=pixx[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      pixy[0][sy]=pixy[NUMT][sy];
      pixy[NUMT+1][sy]=pixy[1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixy[sx][0]=pixy[sx][NUMT];
      pixy[sx][NUMT+1]=pixy[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      piyy[0][sy]=piyy[NUMT][sy];
      piyy[NUMT+1][sy]=piyy[1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      piyy[sx][0]=piyy[sx][NUMT];
      piyy[sx][NUMT+1]=piyy[sx][1];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      pib[0][sy]=pib[NUMT][sy];
      pib[NUMT+1][sy]=pib[1][sy];     
    }
}


//reads in initial energy density profile
//from inited.dat (generate by either initE.cpp or your
//favorite routine) as well as flow velocities and shear tensor
void setInitialConditions()
{

  //CFM EDIT sig unused
  //  double sig;
  //extern double randGauss(double); // generates a gaussian random number with mean 0

  //Middle
  Middle=(NUMT-1)/2+1;
  //Middle=1;

  //printf("mi %i\n",Middle);

  //freeze-out temperature
  TF=TF*AT;

  //cout << "TF=" << TF/AT << endl;

  //cout << "TSTART=" << TSTART << endl;

  //load equation of state

  //loadeos();
  //loadeta();
  //loadzeta();
  //loadbeta2();
  //loadlambda1();
  

  fstream inited,initux,inituy,initpixx,initpixy,initpiyy,itime,initpi;
  itime.open("input/time.dat",ios::in);
  inited.open("input/inited.dat", ios::in);
  initux.open("input/initux.dat", ios::in);
  inituy.open("input/inituy.dat", ios::in);  
  initpixx.open("input/initpixx.dat", ios::in);
  initpixy.open("input/initpixy.dat", ios::in);
  initpiyy.open("input/initpiyy.dat", ios::in);
  initpi.open("input/initpi.dat", ios::in);
  

  std::cout << "Initialize CHECK: " << std::endl;
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {

	
	inited >> e[sx][sy];
	e[sx][sy]*=SCAL;
	std::cout << e[sx][sy] << " ";
	initux >> u[0][sx][sy];
	inituy >> u[1][sx][sy];
	initpixx >> pixx[sx][sy];
	initpixy >> pixy[sx][sy];
	initpiyy >> piyy[sx][sy];
	initpi >> pib[sx][sy];
      }

  std::cout << std::endl;
 
  itime >> TINIT;

  //convert fm/c to lattice units
  t=TINIT*fmtoGeV/AT;

  


  inited.close();
  initux.close();
  inituy.close();
  initpixx.close();
  initpixy.close();
  initpiyy.close();
  initpi.close();

  enforcePBCs();

  //if this option is chosen, pre-flow is being generated
  //from ed distribution
  if (preeqflow==2)
    {
      printf("===> Info: preeqflow=2 option chosen, ignoring velocity input files\n");
      


      /*   for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    u[0][sx][sy]=0.;
	    u[1][sx][sy]=0.;
	    }*/

      
      
      if (SMOOTHING)
	for (int a=0;a<40;a++)
	smeare(e);
      
      for (int sx=2;sx<=NUMT-1;sx++)
	for (int sy=2;sy<=NUMT-1;sy++)
	  {
	    u[0][sx][sy]=-0.5*(e[sx+1][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	    u[1][sx][sy]=-0.5*(e[sx][sy+1]-e[sx][sy-1])/e[sx][sy]/3.*t;
	  }

      /*
      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=1;
	  u[0][sx][sy]=-0.5*(e[sx+1][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-(e[sx][2]-e[sx][1])/e[sx][sy]/3.*t;
	}
      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=NUMT;
	  u[0][sx][sy]=-0.5*(e[sx+1][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-(e[sx][NUMT]-e[sx][NUMT-1])/e[sx][sy]/3.*t;
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=1;
	  u[0][sx][sy]=-(e[sx+1][sy]-e[sx][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-0.5*(e[sx][sy+1]-e[sx][sy-1])/e[sx][sy]/3.*t;
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=NUMT;
	  u[0][sx][sy]=-(e[sx][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-0.5*(e[sx][sy+1]-e[sx][sy-1])/e[sx][sy]/3.*t;
	  }*/

      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=1;
	  u[0][sx][sy]=u[0][sx][sy+1];
	  u[1][sx][sy]=u[1][sx][sy+1];
	}
      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=NUMT;
	  u[0][sx][sy]=u[0][sx][sy-1];
	  u[1][sx][sy]=u[1][sx][sy-1];
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=1;
	  u[0][sx][sy]=u[0][sx+1][sy];
	  u[1][sx][sy]=u[1][sx+1][sy];
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=NUMT;
	  u[0][sx][sy]=u[0][sx-1][sy];
	  u[1][sx][sy]=u[1][sx-1][sy];
	}

      
      double stopval=1.;
      for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    if (fabs(u[0][sx][sy])>stopval)
	      u[0][sx][sy]*=stopval/fabs(u[0][sx][sy]);
	    if (fabs(u[1][sx][sy])>stopval)
	      u[1][sx][sy]*=stopval/fabs(u[1][sx][sy]);	   
	      } 
      
      if (SMOOTHING)
	for (int a=0;a<10;a++)
	smearu(u,e);
    }

  //printf("ED at center %f\n",e[Middle][Middle]);
}

int generatehadronparameters()
{
  FILE * pFile;
  
  int worked=1;
  pFile = fopen ("parameters/default/fixed.param","w");
  
  if (pFile!=NULL)
    {
      fprintf(pFile,"###################\n");
      fprintf(pFile,"bool B3D_PRHYDRO true\n");
      fprintf(pFile,"double B3D_PR_OMEGA_XYSPACING %f\n",AT/fmtoGeV);
      fprintf(pFile,"double B3D_PR_OMEGA_TAUSPACING %f\n",AT*EPS*UPDATE/fmtoGeV);
      fprintf(pFile,"int B3D_PR_NPRCELLSMAX 10000000\n");
      fprintf(pFile,"int B3D_NEVENTSMAX %li\n",B3DEVENTS);
      fprintf(pFile,"int B3D_NPARTSMAX 20000\n");
      fprintf(pFile,"double HYDRO_FOTEMP %f\n",TF/AT*1000);
      fprintf(pFile,"int B3D_NACTIONSMAX 100000\n");
      fprintf(pFile,"###################\n");
      
      fprintf(pFile,"bool B3D_OSUHYDRO false\n");
      fprintf(pFile,"double B3D_TAUCOLLMAX 50.0\n");
      fprintf(pFile,"double B3D_ETAMAX 1.0\n");
      fprintf(pFile,"double B3D_XYMAX 20.0\n");
      //fprintf(pFile,"double B3D_XYMAX 40.0\n");
      fprintf(pFile,"int B3D_NETA 4\n");
      //fprintf(pFile,"int B3D_NXY 20\n");
      fprintf(pFile,"int B3D_NXY 10\n");
      fprintf(pFile,"int B3D_NSAMPLE 1\n");
      fprintf(pFile,"string B3D_RESONANCES_DECAYS_FILE progdata/resinfo/decays_pdg.dat\n");
      fprintf(pFile,"string B3D_RESONANCES_INFO_FILE progdata/resinfo/resonances_pdg.dat\n");
      fprintf(pFile,"string B3D_INPUT_DATAROOT output\n");
      fprintf(pFile,"string B3D_OUTPUT_DATAROOT output\n");
      fprintf(pFile,"bool B3D_BJORKEN true\n");
      fprintf(pFile,"bool B3D_ERROR_PRINT false\n");
      
      fprintf(pFile,"double B3D_SIGMAMAX 30.0\n");
      fprintf(pFile,"double B3D_SIGMADEFAULT 1.0\n");
      fprintf(pFile,"bool B3D_BJMAKER_GAUSSIAN false\n");
      fprintf(pFile,"bool B3D_BJMAKER_BALANCE false\n");
      fprintf(pFile,"bool B3D_STAR_ACCEPTANCE false\n");
      fprintf(pFile,"bool B3D_VIZWRITE false\n");
      fprintf(pFile,"double B3D_SIGMAINELASTIC 0.0\n");
      fprintf(pFile,"double B3D_INELASTIC_Q0 200\n");
      fprintf(pFile,"bool B3D_INELASTIC false\n");
      fprintf(pFile,"string B3D_INELASTIC_INFO_FILE inelastic.tmpdouble GLAUBER_PP_ENERGY_FRAC 1.0\n");
      fprintf(pFile,"#######################\n");
      
      
      fclose (pFile);
    }
  else
    {
      //printf("==> Info: opening b3d params file failed!\n");
      worked=0;
    }
  return worked;

  //dparams.close();
}

//freeze-out that doesn't assume a monotonically 
//decreasing temperature from the center out
void blockfreeze()
{


  //All grid points that are nearest neighbors in x, y, or tau are checked pairwise.
  //A rectangular piece of the freezeout surface is defined to be half-way between
  //any pair whose temperatures straddle TF.
  
  reachedTf = 1; //hydro evolution is finished if everywhere in the time slice T < TF
  for (int sx=1;sx<=NUMT;sx++) {
    for (int sy=1;sy<=NUMT;sy++)
    {
      
      if (T(sx,sy) > TF) reachedTf = 0; //if one cell has T>TF, freezeout not reached
      
      if (sy != NUMT)
      {
	if (((T(sx,sy) > TF) && (T(sx,sy+1) <= TF)) || ((T(sx,sy) <= TF) && (T(sx,sy+1) > TF)))
	{ 
	  reachedTf = 0;
	  int direction = 2;
	  if ((T(sx,sy) <= TF) && (T(sx,sy+1) > TF)) direction = -2;
	  freeze_out << (sx-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << (sy-Middle+0.5)/fmtoGeV*AT << "\t";
	  freeze_out << t/fmtoGeV*AT << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + u[0][sx][sy+1]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + u[1][sx][sy+1]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixx[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1]))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixy[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1])) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + piyy[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1])))<< "\t";
	  freeze_out << 0.5 * (pib[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pib[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1])) ) << "\t";
	  freeze_out << 0.5 * (T(sx,sy) + T(sx,sy+1))/AT << "\n";
	}
      }
      if (sx != NUMT)
      {
      	if (((T(sx,sy) > TF) && (T(sx+1,sy) <= TF)) || ((T(sx,sy) <= TF) && (T(sx+1,sy) > TF)))
	{ 
	  reachedTf = 0;
	  int direction = 1;
	  if ((T(sx,sy) <= TF) && (T(sx+1,sy) > TF)) direction = -1;
	  freeze_out << (sx-Middle+0.5)/fmtoGeV*AT << "\t";
	  freeze_out << (sy-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << t/fmtoGeV*AT << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + u[0][sx+1][sy]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + u[1][sx+1][sy]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixx[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy]))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixy[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy])) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + piyy[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy])))<< "\t";
	  freeze_out << 0.5 * (pib[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + pib[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy])))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy) + T(sx+1,sy))/AT << "\n";
	}
      }

      if (t>TINIT*fmtoGeV/AT)
      {
	if ((T(sx,sy) <= TF && (Tlast(sx,sy) > TF)) || ((T(sx,sy) > TF) && (Tlast(sx,sy) <= TF)))
	{
	  int direction = 3;
	  if ((T(sx,sy) > TF) && (Tlast(sx,sy) <= TF)) direction = -3;
	  freeze_out << (sx-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << (sy-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << t/fmtoGeV*AT - 0.5 * UPDATE*EPS*AT/fmtoGeV << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + ulast[0][sx][sy]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + ulast[1][sx][sy]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 	
			      + pixxlast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy]))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixylast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy])) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + piyylast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy])))<< "\t";
	  freeze_out << 0.5 * (pib[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + pilast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy])))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy) + Tlast(sx,sy))/AT << "\n";
	}
      }
    }
  }
}


double ut(int sx,int sy)
{
  double temp=1.;
  temp+=u[0][sx][sy]*u[0][sx][sy];
  temp+=u[1][sx][sy]*u[1][sx][sy];
  return sqrt(temp);
}

//returns correct pi
double pi (int delta,int beta,int sx,int sy)
{

  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta

  int phi=0;
  int deltap,betap;

  if (beta==2)
    betap=0;
  if (beta==3)
    betap=3;
  if (beta==1)
    betap=2;
  if (beta==0)
    betap=1;
  if (delta==2)
    deltap=0;
  if (delta==3)
    deltap=3;
  if (delta==1)
    deltap=2;
  if (delta==0)
    deltap=1;
  

  if (deltap>betap)
    {
      phi=betap;
      betap=deltap;
      deltap=phi;
    }
  if (deltap<betap)
    {
      if (betap==1 && deltap==0)
	{
	  return u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]+u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy];
	}
      if (betap==2)
	{
	  if (deltap==0)
	    return u[0][sx][sy]/ut(sx,sy)*pixy[sx][sy]+u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy];
	  if (deltap==1)
	    return pixy[sx][sy];
	}
      else return 0;
    }

  if (deltap==betap)
    {
      if (deltap==0)
	return u[0][sx][sy]/ut(sx,sy)*u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]+2*u[0][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy]+u[1][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy];
      if (deltap==1)
	return pixx[sx][sy];
      if (deltap==2)
	return piyy[sx][sy];
      if (deltap==3)
	return (-1)/(t*t)*((-1)*u[0][sx][sy]/ut(sx,sy)*u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]-2*u[0][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy]-u[1][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy]+pixx[sx][sy]+piyy[sx][sy]);
    }

  //CFM EDIT return -1 for double function at end
  return -1;
}


void outputMeasurements(double t, int outStep = 0)
{
  cout.precision(5);
  int dwidth = 13;

  std::string outStr = "Step" + std::to_string(outStep);
  std::string outStr2 = std::to_string(t/fmtoGeV*AT);
  if(outStr2.find(".") != std::string::npos){
    while(outStr2.substr(outStr2.size()-1, 1).find("0") != std::string::npos){
      outStr2 = outStr2.substr(0, outStr2.size()-1);
    }
  }
  
  std::string outStr3 = "";
  while(outStr.size() + outStr2.size() + outStr3.size() < dwidth){outStr3 = outStr3 + " ";}

  std::cout << outStr << outStr3 << outStr2;
  //  cout.width(dwidth); cout << outStr << t/fmtoGeV*AT;
  //globali=geti(e[Middle][Middle]);
  //globalx=getx(globali,e[Middle][Middle]);
//   double ex=anisospace();
//   double ep=anisomomentum();
//   double px = totalmomentumx();
//   double py = totalmomentumy();

  //CFM EDIT totpx, totpy, e1c, e1s, e2c, e2s, e3c, e3s unused
  //  double ex, ep, totpx, totpy, e1c, e1s, e2c, e2s, e3c, e3s;
  double ex, ep;
  //CFM EDIT eps1c, eps1s, eps2c, eps2s, eps3c, eps3s, eps4c, eps4s unused;
  //  double eps1c, eps1s, eps2c, eps2s, eps3c, eps3s, eps4c, eps4s;

  //CFM EDIT eps5c, eps5s, eps6c, eps6s, eps7c, eps7s unused;
  //  double eps5c, eps5s, eps6c, eps6s, eps7c, eps7s;
  /*allanisomomentum(ex,ep, totpx, totpy, e1c, e1s, e2c, e2s, e3c, e3s,
		  eps1c, eps1s, eps2c, eps2s, eps3c, eps3s, eps4c, 
		  eps4s, eps5c, eps5s, eps6c, eps6s, eps7c, eps7s,pacc,cs2acc);*/

  cout.width(dwidth); cout << T(Middle,Middle)/AT;
  cout << "\t" << ex << "\t" << ep;
  //cout << "\t" << overlapS()/AT/AT;
  //cout << "\t" << l1coeff(T(Middle,Middle));  
  cout << endl;

  meta << t/fmtoGeV*AT <<"\t";
  meta << T(Middle,Middle)/AT << "\t";
  meta << e[Middle][Middle]/AT/AT/AT/AT << "\t";
  meta << -pi(3,3,Middle,Middle)*t*t/AT/AT/AT/AT << "\t";
  meta << pib[Middle][Middle]/AT/AT/AT/AT << "\t";
  meta << "\n";
  
  ecces << t/fmtoGeV*AT <<"\t";
  ecces << ex << "\t";
  ecces << ep << "\t";
  

  switch (FREEZE){
    case 0:
      //stupidfreeze(pacc,cs2acc);
      break;
      /*    case 1:
      fancyfreeze(pacc,cs2acc);
      break;*/
    case 2:
      blockfreeze();
      break;
    default:
      blockfreeze();
  }
}

void snapTcontour(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tcontour_%.3f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);

   out << "#x [fm] \t y [fm] \t T [GeV] \n";
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  //globali=geti(e[sx][sy]);
	  //globalx=getx(globali,e[sx][sy]);
	  out << (sx-Middle)/fmtoGeV*AT;
	  out << "\t";
	  out << (sy-Middle)/fmtoGeV*AT;
	  out << "\t";
	  out << T(sx,sy)/AT << "\n";
	}
      out << endl;
    }
  out.close();
}

void snapVcontour(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vcontour_%.3f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  out << "#x [fm] \t y [fm] \t ux \t uy \t gamma \n";
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  //globali=geti(e[sx][sy]);
	  //globalx=getx(globali,e[sx][sy]);
	  out << (sx-Middle)/fmtoGeV*AT << "\t";
	  out << (sy-Middle)/fmtoGeV*AT << "\t";
	  out << u[0][sx][sy] << "\t";
	  out << u[1][sx][sy] << "\t";
	  out << sqrt(u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]) <<"\n";
	}
    }
  out.close();
}

void snapTprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    //globali=geti(e[Middle][s]);
    //globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << T(Middle,s)/AT << endl;
  }
  out.close();
}

void snapEDprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/EDprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << e[Middle][s]/AT/AT/AT/AT << endl;
  }
  out.close();
}

void snapVprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    //globali=geti(e[Middle][s]);
    //globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << u[1][Middle][s]/ut(Middle,s) << endl;
  }
  out.close();
}

void snappieeprofile(double time)
{
  fstream out;
  char fname[255];

  //CFM EDIT vr unused
  //double vr=0;
  sprintf(fname,"data/snapshot/Piprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    //globali=geti(e[Middle][s]);
    //globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << -pi(3,3,Middle,s)*t*t/(e[Middle][s]*4./3.)<< endl;
  }
  out.close();
}

void snappirrprofile(double time)
{
  fstream out;
  char fname[255];
  //CFM EDIT vr unused
  //double vr=0;
  sprintf(fname,"data/snapshot/PiRprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {   
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << -pi(1,1,Middle,s)/(4./3.*e[Middle][s])<< endl;
  }
  out.close();
}

void snappibulkprofile(double time)
{
  fstream out;
  char fname[255];
  //CFM EDIT vr unused
  //double vr=0;
  sprintf(fname,"data/snapshot/PiBulkprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {    
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << pib[Middle][s]/pow(AT,4)<< endl;
  }
  out.close();
}



void snapshot(double tt)
{
  //CFM EDIT hel2, dt, tb, td unused
  //  double hel2[4];
  //  double dt[4];
  //  double tb[4],td[4];
 
  //fordebug(10,10);
  
  snapTcontour(tt);
  snapVcontour(tt);
  //snapFOdata(tt,pacc,cs2acc);

  snapTprofile(tt);   
  snapEDprofile(tt);
  snapVprofile(tt);
  //snapVxprofile(tt);
  //snapV2profile(tt);
  //snapV2profile(tt);
  //snapuphiprofile(tt);
  snappieeprofile(tt);
  snappirrprofile(tt);
  snappibulkprofile(tt);
}

//Copy fields at each UPDATE time step
void copyUPDATE()
{
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	for (int i=0;i<2;i++) 
	  ulast[i][sx][sy]=u[i][sx][sy];
	elast[sx][sy]=e[sx][sy];
	pixylast[sx][sy]=pixy[sx][sy];
	pixxlast[sx][sy]=pixx[sx][sy];
	piyylast[sx][sy]=piyy[sx][sy];
	pilast[sx][sy]=pib[sx][sy];
      }
}
