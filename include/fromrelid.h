void printNxNmatrix(double *mat,int DIM,int flag)
{
  //standard
  if (flag==0)
    {
      for (int a=0;a<DIM;a++)
	{
	  printf("(");
	  for (int b=0;b<DIM;b++)
	    printf(" %.16g ",mat[a*DIM+b]);
	  printf(")\n");
	}
    }

  //Mathematica
  if (flag==1)
    {
      printf("{");
      for (int a=0;a<DIM;a++)
	{
	  printf("{");
	  for (int b=0;b<DIM;b++)
	    {
	      printf("%f",mat[a*DIM+b]);
	      if (b<DIM-1)
		printf(",");
	    }
	  printf("}");
	  if(a<DIM-1)
	    printf(",");
	  else
	    printf("}");
	  printf("\n");
	}
    }
  

}

int calcedvquick(double *tab, double *statevec,int mdebug)
{
  int WARN=0;
  double uguess[]={1,0,0,0};
  double uguessnew[]={0,0,0,0};
  double eu0=-1.0;
  int stepper=0;

  if (mdebug)
    printNxNmatrix(tab,4,1);

  while((fabs((uguessnew[0]-eu0)/eu0)>1.e-5)&&(stepper<100))
    {
      eu0=uguessnew[0];      
      for (int i=0;i<4;i++)
	uguessnew[i]=0.0;
      
      for (int i=0;i<4;i++)
	for (int j=0;j<4;j++)
	  uguessnew[i]+=tab[4*i+j]*uguess[j];

      for (int i=0;i<4;i++)
	{ 
	  if (mdebug)
	    printf("i=%i u=%g\n",i,uguessnew[i]);
	  uguess[i]=uguessnew[i]/uguessnew[0];
	}      
      stepper++;
      if (mdebug)
	printf("guess %f\n",uguessnew[0]);
    }
  
  if (isnan(uguessnew[0]))
    WARN=1;


  for (int a=1;a<4;a++)
    {
      if (isnan(uguess[a]))
	WARN=1;
      else
	statevec[a]=uguess[a];
    }

  double u2=statevec[1]*statevec[1]+statevec[2]*statevec[2]+statevec[3]*statevec[3];

  statevec[0]=uguessnew[0];

  
  

  if (mdebug)
    printf("u2=%f 1=%f 2=%f 3=%f\n",u2,statevec[1],statevec[2],statevec[3]);
  return WARN;
}
