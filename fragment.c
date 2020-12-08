 //FILE *cubic;
  //cubic = fopen( "cubic.out", "w" ); 
  //fprintf(cubic,"# 1: a \n");
  //fprintf(cubic,"# 2: cs2 \n");
  //int i;
  //for (i=0; i<NBINS; i++)
    //{

  // fprintf(cubic,"%12lg %12lg\n",x,cs2);
   
     //}


  //fclose(cubic);













#ifdef Cubic_Galileon

  double M,c2,c3,ksi;
  double H,dHdt,dphidt,ddphiddt;                     

  double a2,a3,x1,x3,epsilon,h1,Qs,cs2,lambda2,zeta;
  double e1,e2,e3,b2,b3;
  
  //Omega0=0.25
  //OmegaLambda=0.75
  //H_over_c=100./299792.458                 //c is not dimensionless,and H0 is 100,so the dimension of 100 can get from H0=100=100h km/s/Mpc=v/d, 
  // Mpl=1.                                  //so h/Mpc=1/d=k
  
  //H_over_c = 100. / SPEEDOFLIGHT;        
  M=pow(100.*params.Hubble100,2./3.);         
  c2=-1.;
  c3=1./6./sqrt(6.*params.OmegaLambda);
  ksi=sqrt(6.*params.OmegaLambda);


  H=100.*params.Hubble100*sqrt((params.Omega0/pow(x,3.)+OmegaRad/pow(x,4.)+sqrt(pow(params.Omega0/pow(x,3.)+OmegaRad/pow(x,4.),2.)+4.*params.OmegaLambda))/2.);
  dHdt=(pow(100.*params.Hubble100,4.)*params.OmegaLambda/H/H-H*H-100.*params.Hubble100*100.*params.Hubble100*(params.Omega0/pow(x,3.)+2.*OmegaRad/pow(x,4.))/2.)/(1.+pow(100.*params.Hubble100/H,4.)*params.OmegaLambda);
 

  dphidt=ksi*100.*params.Hubble100*100.*params.Hubble100/H;                                             
  ddphiddt=-ksi*100.*params.Hubble100*100.*params.Hubble100*dHdt/H/H;

  
  //parameters for Cubic Galileon
  a2=-c2/2.;
  a3=c3/M/M/M/3.;
  x1=-a2*dphidt*dphidt/H/H/3.;
  x3=6.*a3*dphidt*dphidt*dphidt/H;
  epsilon=ddphiddt/H/dphidt;
  h1=dHdt/H/H;
 
  Qs=3.*(4.*x1+4.*x3+x3*x3)/(2.-x3)/(2.-x3);
  cs2=2.*(1.+3.*epsilon)*x3-x3*x3-4.*h1-6.*params.Omega0*100.*params.Hubble100*100.*params.Hubble100/H/H/x/x/x;   
  lambda2=12.*a3*dphidt*dphidt/(H*H*cs2*Qs)/(2.-x3)/(2.-x3);
  zeta=lambda2*dphidt*dphidt;


  //parameters of growth functions for Cubic Galileon
  e1=-(3./x + 0.25*100.*params.Hubble100*100.*params.Hubble100*(-3.*params.Omega0/pow(x,4.0)-
3.*params.Omega0*params.Omega0/pow(x,7.0)/sqrt(params.Omega0*params.Omega0/pow(x,6.)+4.*params.OmegaLambda))/H/H);
  
  b2=1.5 *params.Omega0*100.*params.Hubble100*100.*params.Hubble100/pow(x,5.0)/H/H;
  e2=b2*(1.-zeta*3.*a3*dphidt*dphidt/x/x);
  
  b3=3.*a3*dphidt*dphidt*b2*b2*H*H*zeta*zeta*lambda2;            
  e3=0.5*e2-b3;


/* cosmic time */
  dydx[0] = 1.0*100.*params.Hubble100/x/H;


  for (int j=0; j<NkBINS; j++)
    {      

      dydx[1+j*8] = e1*y[1+j*8] + e2*y[2+j*8];                                           /* d^2 D1/ da^2 */         
      dydx[2+j*8] = y[1+j*8];                                                        /* d D1/ da */
      dydx[3+j*8] = e1*y[3+j*8] + e2*y[4+j*8] + e3*y[2+j*8]*y[2+j*8];                            /* d^2 D2/ da^2 */
      dydx[4+j*8] = y[3+j*8];                                                        /* d D2/ dt */

    /* third-order growth is not used, it is set to the LCDM one */
      dydx[5+j*8] = e1*y[5+j*8] + b2*y[6+j*8] - 2.*b2*y[2+j*8]*y[2+j*8]*y[2+j*8]; /* d^2 D31/ da^2 */
      dydx[6+j*8] = y[5+j*8];                                                      /* d D31/ da */
      dydx[7+j*8] = e1*y[7+j*8] + b2*y[8+j*8] - 2.*b2*y[2+j*8]*y[4+j*8] 
	+ 2.*b2*y[2+j*8]*y[2+j*8]*y[2+j*8];                                     /* d^2 D32/ da^2 */                                      
      dydx[8+j*8] = y[7+j*8];                                                /* d D32/ da */


    }

  return GSL_SUCCESS;
#endif








#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".scale_factor.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# 1: scale factor\n");
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  fprintf(fd," %12lg",pow(10.,SPLINE[SP_TIME]->x[i]));
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif


#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".D1.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: linear growth rate\n",                     2,  NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW1+ik]->y[i]));
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".D2.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: 2nd-order growth rate\n",           NkBINS+2,2*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW2+ik]->y[i]));
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".D31.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: first 3rd-order growth rate\n",   2*NkBINS+2,3*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW31+ik]->y[i]));
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".D32.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: second 3rd-order growth rate\n",  3*NkBINS+2,4*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW32+ik]->y[i]));
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".f1.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: linear d ln D/d ln a\n",          4*NkBINS+2,5*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",SPLINE[SP_FOMEGA1+ik]->y[i]);
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".f2.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: 2nd-order d ln D/d ln a\n",       5*NkBINS+2,6*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",SPLINE[SP_FOMEGA2+ik]->y[i]);
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".f31.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: first 3rd-order d ln D/d ln a\n", 6*NkBINS+2,7*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",SPLINE[SP_FOMEGA31+ik]->y[i]);
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif

#ifdef SCALE_DEPENDENT
      /* writes scale-dependent growth rates on a file */
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".f32.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
#if defined(MOD_GRAV_FR) || defined(Cubic_Galileon)
	  /* with modified gravity the first wavenumber is set to zero */
	  if (!ik)
	    k=0.0;
	  else
#endif
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# %d-%d: second 3rd-order d ln D/d ln a\n",7*NkBINS+2,8*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",SPLINE[SP_FOMEGA32+ik]->y[i]);
	  
	  fprintf(fd,"\n");
	}
      fclose(fd);
