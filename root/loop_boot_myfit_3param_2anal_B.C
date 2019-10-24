

//
//  Information on fitting
//    http://root.cern.ch/root/html/tutorials/fit/Ifit.C.html
//
//

enum fitmodel {  A_SQ , A_ALPHA_S , A_ONLY  } ;


int ndata , ndata_total   ;
int ncopy ;
double *a_fm, *hyper,*hyper_err, *meta_s, *m_pi , *alpha_s ;
double *m_pi_sq, *meta_s_sq , *hyper_copy , *Q_prob  ;
double *m_pi_err, *meta_s_err, *m_pi_sq_err, *meta_s_sq_err ;

double *hyper_copy , *Q_prob_copy ;

FILE *fp ;

const double m_etas_phys = 0.6858 ;
const double m_pi_phys = 0.135 ;

#include "TMatrixT.h"
TMatrixD *all_data ;
TMatrixD *all_prob ;

const int noboot = 1000 ;

char input_tag[100] ;


void load_boot(char fff[],  int noboot, 
	       double & quantity, double & quantity_err, 
	       double *boot_quantity)
{
  string fff_name(fff) ;

  read_quantity_from_disk(fff_name,noboot,quantity, boot_quantity) ;
  quantity_err = boot_variance(boot_quantity, noboot) ;

}


double calc_alpha(double a_fm)
{
  double alpha_s_init = 0.286474  ;
  double scale_init = 2.131 ; 
  double scale_want   ;
  double alpha_want ;

  if( a_fm == 0.0 )
    {
      cout << "ERROR: zero lattice spacing\n" ;
      alpha_want = 0.0 ;
      scale_want = 0.0 ;
    }
  else
    {
      scale_want  = 0.1973 /a_fm  ;
      int nf=3;
      double beta[4] ;
      set_beta(beta, nf) ;
      alpha_want = evolve_alphas(scale_want,scale_init, 
				 alpha_s_init,
				 beta,nf, 2000) ;
    }

  //      cout << "DEBUG " << a_fm << " " <<  scale_want << " " << alpha_want << "\n" ;
      return alpha_want ;
}


void load_data_init()
{
  cout << "Reading header information from \n"  ;
  //  char filename[] = "./loop_data.txt"  ;
  cout << "Reading data from  " << fileN << "\n" ;
  if( (fp = fopen(fileN,"r") ) == NULL )
    {
      cerr << "Error reading " << fileN << endl ;
      exit(0) ; 
    }

  int ncols ;

  ncols = fscanf(fp,"%s",input_tag);
  if( ncols == 1 )
    {
      cout << "Input tag = " << input_tag << "\n" ;
    }
  else
    {
      cout << "ERROR in reading the tag\n" ;
      exit(1) ;
    }



  ncols = fscanf(fp,"%d",&ndata_total);
  if( ncols == 1 )
    {
      cout << "Number of data items found = " << ndata_total << "\n" ;
    }
  else
    {
      cout << "ERROR in reading number of values\n" ;
      exit(1) ;
    }
  ndata = ndata_total ;

  ncols = fscanf(fp,"%d",&ncopy);
  if( ncols == 1 )
    {
      cout << "Number of copies = " << ncopy << "\n" ;
    }
  else
    {
      cout << "ERROR in reading number of copies\n" ;
      exit(1) ;
    }

  cout << "The IO has been initialized and header read\n" ;

}

void reserve_memory()
{
  a_fm       = new double[ndata] ;
  alpha_s       = new double[ndata] ;

  hyper     = new double[ndata] ;
  hyper_err = new double[ndata] ;

  hyper_copy     = new double[ncopy] ;
  Q_prob_copy = new double[ncopy] ;
  
  meta_s     = new double[ndata] ;
  meta_s_err = new double[ndata] ;

  m_pi       = new double[ndata] ;
  m_pi_err   = new double[ndata] ;
  
  meta_s_sq     = new double[ndata] ;
  meta_s_sq_err = new double[ndata] ;

  m_pi_sq       = new double[ndata] ;
  m_pi_sq_err   = new double[ndata] ;

  Q_prob        = new double[ndata] ;
  
  all_data = new TMatrixD(ndata, noboot) ;
  all_prob = new TMatrixD(ndata, noboot) ;

  cout << "The memory has been reserved\n" ;
}

void write_residual(Double_t *par)
{
   Int_t i;

   cout << "List of residuals \n" ;
//calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   Double_t er ;
   Double_t tmp ;


   for (i=0;i<ndata; i++) {
     tmp = func_err(meta_s_sq[i],meta_s_sq_err[i],
		    m_pi_sq[i],m_pi_sq_err[i],
     		    a_fm[i],alpha_s[i],par) ;
       
     er = sqrt(tmp*tmp + hyper_err[i]*hyper_err[i] ) ;

     delta  = (hyper[i] - func(a_fm[i], m_pi_sq[i],
			       meta_s_sq[i],par))/er;

     cout << "Residual[ " << i << " ]= " << delta << " pion mass " <<  m_pi[i]  << " GeV\n" ;
     chisq += delta*delta;
   }

   cout << "chi**2 from residuals = " << chisq << "\n" ;

}



void plot_func(double par[], int w_fit)
{
  const int n_plot_pts = 200 ;
  const double a_start = 0.0 ;
  const double a_end = 0.1 ;

  double a_delta = (a_end - a_start) / n_plot_pts;
  double a  =  a_start ; 
  double mpi_sq_phys = m_pi_phys * m_pi_phys ;
  double mpi_etaS_phys = m_etas_phys * m_etas_phys ;

   char filename[200] ;

  if( w_fit == A_SQ )
    sprintf(filename,"plots/%s_plotfit_asq.dat",input_tag) ;
  else if( w_fit == A_ALPHA_S )
    {
      sprintf(filename,"plots/%s_plotfit_a_alpha.dat",input_tag) ;
    }
  else if( w_fit == A_ONLY )
    sprintf(filename,"plots/%s_plotfit_a_alpha.dat",input_tag) ;
  else
    {
      cout << "plot_func::ERROR w_fit= " << w_fit << " out of range\n" ;
      exit(0);
    }


  cout << "SUM Writing plotted fit to the file  " << filename << "\n" ;

  FILE *f_plot ;
  f_plot = fopen(filename,"w");

  for(int i = 0 ; i < n_plot_pts ; ++i)
    {
      double ff ;

      if( w_fit == A_SQ )
	ff = func(a, mpi_sq_phys,mpi_etaS_phys, par);
      else if( w_fit == A_ALPHA_S )
	{
	  // org	ff = func_alpha(a, mpi_sq_phys,mpi_etaS_phys, par);
	  double alpha_want = calc_alpha(a) ;
	  ff = func_alpha(a, alpha_want,mpi_sq_phys,mpi_etaS_phys, par);
	}
      else if( w_fit == A_ONLY )
	ff = func_A(a, mpi_sq_phys,mpi_etaS_phys, par);

      double a_sq = a * a ;
      fprintf(f_plot,"%g %g\n",a_sq,ff); 

      a += a_delta ;
    }

  fclose(f_plot) ; 
}




void load_data()
{
  cout << "Reading one chunk from input file\n"  ;
  int ncols ;
  double  quantity,  quantity_err ;
  double *boot_quantity ;

  cout << "Number of bootstrap samples " << noboot  << "\n" ;
  boot_quantity = new double[noboot] ;

  char mass[300] ;
  char q_value[300] ;

  int ii = 0 ;

  for(int i = 0 ; i < ndata_total ; ++i)
    {

      ncols = fscanf(fp,"%lf %lf %s %s %lf %lf",
		     &a_fm[ii],&alpha_s[ii],
		     mass,q_value, 
		     &m_pi[ii],&m_pi_err[ii]);
      if( ncols != 6 )
	{
	  cout << "Error reading data \n" ;
	  exit(1) ;
	}


      ncols = fscanf(fp,"%lf %lf",&meta_s[ii], &meta_s_err[ii]) ;
      if( ncols != 2 )
	{
	  cout << "Error reading data \n" ;
	  exit(1) ;
	}

      double a_inv_gev = 0.1973 / a_fm[ii] ;
      double m_pi_gev = m_pi[ii] * a_inv_gev ;


      //      cout << "DEBUG a_fm_cut = " << a_fm_cut << " a_fm[ii]= " << a_fm[ii] << endl ;
      if( a_fm[ii] < a_fm_cut && m_pi_gev  < m_pi_gev_cut )
	{


	  load_boot(mass,noboot,quantity,quantity_err,boot_quantity) ; 
	       
	  hyper[ii]     = quantity ;
	  hyper_err[ii] = quantity_err ;

	  for(int ib=0 ; ib < noboot ; ++ ib)
	    all_data(ii,ib) = boot_quantity[ib] ; 

	  load_boot(q_value,noboot,quantity,quantity_err,boot_quantity) ; 
	  Q_prob[ii] =  quantity ;

	  for(int ib=0 ; ib < noboot ; ++ib)
	{
	  all_prob(ii,ib) = boot_quantity[ib] ; 
	  //	cout << "DEBUG all_prob = " << all_prob(i,ib) << "\n" ;
	}
	  // exit(0) ;

	  cout << "q_value = " << q_value << " Q_prob= " << Q_prob[ii] << "\n" ;


	  meta_s[ii] *= a_inv_gev ;
	  m_pi[ii]   *= a_inv_gev;

	  meta_s_err[ii] *= a_inv_gev ;
	  m_pi_err[ii]   *= a_inv_gev;

	  meta_s_sq[ii] = meta_s[ii] * meta_s[ii] ;
	  meta_s_sq_err[ii] = 2.0 * meta_s[ii] * meta_s_err[ii] ;

	  m_pi_sq[ii] = m_pi[ii] *  m_pi[ii] ;
	  m_pi_sq_err[ii] = 2.0 * m_pi_err[ii] * m_pi[ii] ;

	  ++ii; 
	}


    }

  cout << "Number of ensembles included after cut = " << ii << "\n" ; 
  ndata = ii ;

  for(int i = 0 ; i < ndata ; ++i)
    {
      cout << "PION_MASS[" ;
      cout << i << "] = " ;
      cout <<  m_pi[i] << " +/- " << m_pi_err[i] << " " ;
      cout << " ETA_S= " << meta_s[i] << " " ;
      cout << meta_s_err[i] << "\n" ;

    }




  for(int i = 0 ; i < ndata ; ++i)
    {
      cout << "mc/ms[" ;
      cout << i << "] = " ;
      cout <<  hyper[i] << " +/- " << hyper_err[i] << " " ;
      cout << "\n" ;

    }


  cout << "Square of masses \n"  ;
  for(int i = 0 ; i < ndata ; ++i)
    {
      cout << "SQUARE " ;
      cout << i << " " ;
      cout << " " << meta_s_sq[i] << " " ;
      cout << " " << meta_s_sq_err[i] << " " ;

      cout << " " << m_pi_sq[i]  ;
      cout << " " << m_pi_sq_err[i]  ;

      cout << "\n" ;
    }

      free [] boot_quantity ;

      cout << "The data has been loaded\n" ;
}


#include<math.h>
#include "TMinuit.h"

//______________________________________________________________________________


double func(double afm, double mpi_sq,double metas_sq, double *par)
{
  double mpi_sq_phys = m_pi_phys * m_pi_phys ;
  double mpi_etaS_phys = m_etas_phys * m_etas_phys ;

  double xx_ss = (metas_sq - mpi_etaS_phys) / mpi_etaS_phys ;

  double value= par[0] * ( 1.0 + par[1]*afm*afm)/( 1.0 + par[2] * xx_ss) ;

 return value;
}


double func_err(double meta_s_sq, double meta_s_sq_err, 
		    double mpi_sq, double mpi_sq_err,
		    double a, double alpha, double *par)
{
  double tmp1 = par[0] * (1.0 + par[1] *a*a  ) ;

  double tmp2 = (m_etas_phys*m_etas_phys)/ meta_s_sq ;
  double tmp2_err  = meta_s_sq_err * (m_etas_phys*m_etas_phys)/ (meta_s_sq*meta_s_sq) ;

  double tmp3 = (1.0 + par[2] * (mpi_sq - m_pi_phys*m_pi_phys ) / m_pi_phys*m_pi_phys )  ; 
  double tmp3_err = mpi_sq_err * par[2]  / (m_pi_phys*m_pi_phys)  ;

  double a_err ; 
  double b_err ; 

  //  a_err = fabs(tmp1 * tmp2_err * tmp3) ; 
  a_err = (tmp1 * tmp2_err * tmp3) ; 

  b_err = tmp1 * tmp2 * tmp3_err ;

  double value=  sqrt(a_err*a_err + b_err*b_err ) ;
  //  cout << "DEBUG value = " << a_err << " " << b_err  << " " << value << "\n" ;

 return value;
}



double func_A(double afm, double mpi_sq,double metas_sq, double *par)
{
  double mpi_sq_phys = m_pi_phys * m_pi_phys ;
  double mpi_etaS_phys = m_etas_phys * m_etas_phys ;

  double value= par[0] * ( 1.0 + par[1]*afm)*( 1.0 + par[2]*(mpi_sq - mpi_sq_phys)/mpi_sq_phys ) *  (mpi_etaS_phys / metas_sq)    ;

 return value;
}


double func_A_err(double meta_s_sq, double meta_s_sq_err, 
		    double mpi_sq, double mpi_sq_err,
		    double a, double alpha, double *par)
{
  double tmp1 = par[0] * (1.0 + par[1] *a ) ;

  double tmp2 = (m_etas_phys*m_etas_phys)/ meta_s_sq ;
  double tmp2_err  = meta_s_sq_err * (m_etas_phys*m_etas_phys)/ (meta_s_sq*meta_s_sq) ;

  double tmp3 = (1.0 + par[2] * (mpi_sq - m_pi_phys*m_pi_phys ) / m_pi_phys*m_pi_phys )  ; 
  double tmp3_err = mpi_sq_err * par[2]  / (m_pi_phys*m_pi_phys)  ;

  double a_err ; 
  double b_err ; 

  //  a_err = fabs(tmp1 * tmp2_err * tmp3) ; 
  a_err = (tmp1 * tmp2_err * tmp3) ; 

  b_err = tmp1 * tmp2 * tmp3_err ;

  double value=  sqrt(a_err*a_err + b_err*b_err ) ;
  //  cout << "DEBUG value = " << a_err << " " << b_err  << " " << value << "\n" ;

 return value;
}




double func_alpha(double afm, double alpha, double mpi_sq,
double metas_sq, double *par)
{
  double mpi_sq_phys = m_pi_phys * m_pi_phys ;
  double mpi_etaS_phys = m_etas_phys * m_etas_phys ;

  double xx_ss = (metas_sq - mpi_etaS_phys) / mpi_etaS_phys ;

  double value= par[0] * ( 1.0 + par[1]*afm*alpha) /( 1.0 + par[2] * xx_ss) ;

 return value;
}


double func_alpha_err(double meta_s_sq, double meta_s_sq_err, 
		      double mpi_sq, double mpi_sq_err,
		      double a, double alpha, double *par)
{
  double tmp1 = par[0] * (1.0 + par[1] *a ) ;

  double tmp2 = (m_etas_phys*m_etas_phys)/ meta_s_sq ;
  double tmp2_err  = meta_s_sq_err * (m_etas_phys*m_etas_phys)/ (meta_s_sq*meta_s_sq) ;

  double tmp3 = (1.0 + par[2] * (mpi_sq - m_pi_phys*m_pi_phys ) / m_pi_phys*m_pi_phys )  ; 
  double tmp3_err = mpi_sq_err * par[2]  / (m_pi_phys*m_pi_phys)  ;

  double a_err ; 
  double b_err ; 

  //  a_err = fabs(tmp1 * tmp2_err * tmp3) ; 
  a_err = (tmp1 * tmp2_err * tmp3) ; 

  b_err = tmp1 * tmp2 * tmp3_err ;

  double value=  sqrt(a_err*a_err + b_err*b_err ) ;
  //  cout << "DEBUG value = " << a_err << " " << b_err  << " " << value << "\n" ;

 return value;
}




//____________________________________________________________________
//
//  Chi**2 definition
//
//____________________________________________________________________


void fcn(Int_t &npar, double *gin, double &f, double *par, Int_t iflag)
{
   Int_t i;

//calculate chisquare
   double chisq = 0;
   double delta;
   double er ;
   double tmp ;

   for (i=0;i<ndata; i++) {
       
     er =hyper_err[i] ;

     delta  = (hyper[i] - func(a_fm[i],m_pi_sq[i],meta_s_sq[i], par))/ er ;
     chisq += delta*delta;
   }
   f = chisq;
}



void fcn_A(Int_t &npar, double *gin, double &f, double *par, Int_t iflag)
{
   Int_t i;

//calculate chisquare
   double chisq = 0;
   double delta;
   double er ;
   double tmp ;

   for (i=0;i<ndata; i++) {
     tmp = func_A_err(meta_s_sq[i],meta_s_sq_err[i],
		    m_pi_sq[i],m_pi_sq_err[i],
     		    a_fm[i],alpha_s[i],par) ;
       
     er = sqrt(tmp*tmp + hyper_err[i]*hyper_err[i] ) ;
     delta  = (hyper[i] - func_A(a_fm[i],m_pi_sq[i],meta_s_sq[i], par))/er ;
     chisq += delta*delta;
   }
   f = chisq;
}



void fcn_alpha(Int_t &npar, double *gin, double &f, double *par, Int_t iflag)
{
   Int_t i;

//calculate chisquare
   double chisq = 0;
   double delta;
   double er ;
   double tmp ;

   for (i=0;i<ndata; i++) {
       
     er = hyper_err[i]  ;

     delta  = (hyper[i] - func_alpha(a_fm[i],alpha_s[i],
				     m_pi_sq[i],
				     meta_s_sq[i], par))/er;


     chisq += delta*delta;
   }
   f = chisq;
}


//_____________________________________________________________________
//

void single_anal(int w_fit, double result[], double prob[], 
		 TMatrixD *l_boot_result, TMatrixD *l_boot_weight,
		 double & chisq_, int & dof_)
{
  double * latt_coeff      = new double[ncopy] ;
  double * mpi_sq_coeff    = new double[ncopy] ;

  //initialize TMinuit with a maximum of 5 params
  TMinuit *gMinuit = new TMinuit(5);  

  if( w_fit == A_SQ )
    gMinuit->SetFCN(fcn);
  else if( w_fit == A_ALPHA_S )
    gMinuit->SetFCN(fcn_alpha);
  else if( w_fit == A_ONLY )
    gMinuit->SetFCN(fcn_A);
  else
    {
      cout << "single_anal::ERROR w_fit= " << w_fit << " out of range\n" ;
      exit(0);
    }


  double arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  // Set starting values and step sizes for parameters
  static double vstart[4] = {10, 10 , 0.1 , 1};
  static double step[4] = {0.1 , 0.1 , 0.01 , 0.1};
  gMinuit->mnparm(0, "con", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "slope", vstart[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "slope_mpisq", vstart[2], step[2], 0,0,ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;

  // Print results
  double amin,edm,errdef;
  Int_t nvpar = 3 , nparx,icstat;

  double  currentValue,  currentError ; 
  double  slopeValue,  slopeError ; 

  
  Int_t dof = ndata - nvpar ;
  double ans = amin / dof ; 
  double amin_tot = 0.0 ;

  for(int ic=0 ; ic < ncopy ; ++ic)
    {
       cout << "Working on copy number " << ic << "\n" ;
       load_data() ;

       gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
       gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

       double  currentValue_b,  currentError_b ;
       gMinuit->GetParameter(0, currentValue_b, currentError_b) ;
       hyper_copy[ic] = currentValue_b ;
       double tnp_err ; 
       gMinuit->GetParameter(1, latt_coeff[ic], tnp_err) ;
       gMinuit->GetParameter(2, mpi_sq_coeff[ic], tnp_err) ;


       write_residual(&gMinuit->fX[0]) ;

       amin_tot += amin ;
       ans = amin / dof ; 
       cout << "Chi**2 / dof = " << amin << " / " << dof  << " = " << ans <<  "\n" ;
       cout << "mc/ms " << currentValue_b << " +/-" << currentError_b << "\n" ;

       Q_prob_copy[ic] = 0 ; 
       for(int ii = 0 ; ii < ndata ; ++ii)
	 {
	 Q_prob_copy[ic] += Q_prob[ii] ;
	 }
       cout << "DEBUG ndata = " << ndata << endl ;
       Q_prob_copy[ic] /= ndata ;

       cout << "Starting bootstrap analysis\n" ;
       for(int ib=0 ; ib < noboot; ++ib)
	 {
	   for(int i = 0 ; i < ndata ; ++i)
	     hyper[i] = all_data(i,ib)  ;

	   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	   double  b_currentValue_b,  b_currentError_b ;
	   gMinuit->GetParameter(0, b_currentValue_b, b_currentError_b) ;
	   l_boot_result(ic,ib) =  b_currentValue_b ;

	   l_boot_weight(ic,ib) = 0.0 ;
	   for(int i = 0 ; i < ndata ; ++i)
	     l_boot_weight(ic,ib) += all_prob(i,ib) ;
	   l_boot_weight(ic,ib) /= ndata ;


	   //	   cout << "DEBUG l_boot_weight(ic,ib) = " << l_boot_weight(ic,ib) << "\n" ;

	 } // bootstrap loop


    } // loop over copies

   cout << "Results ------- \n" ;
   cout << "i    prob  hyperfine \n"; 
   double q_mean = 0.0 ;
   double q_total = 0.0 ;
   double latt_mean = 0.0 ;
   double mpi_sq_coeff_mean = 0.0 ;


   for(int ic=0 ; ic < ncopy ; ++ic)
     {
       result[ic] = hyper_copy[ic] ;
       prob[ic]   = Q_prob_copy[ic] ;

       q_mean  += Q_prob_copy[ic] * hyper_copy[ic] ;
       q_total += Q_prob_copy[ic] ;
       latt_mean           += Q_prob_copy[ic] * latt_coeff[ic] ;
       mpi_sq_coeff_mean   += Q_prob_copy[ic] * mpi_sq_coeff[ic] ;


       cout << ic << " " << Q_prob_copy[ic]  << " " << hyper_copy[ic] << " GeV\n" ; 
     }

   q_mean /= q_total ;
   latt_mean  /= q_total ;
   mpi_sq_coeff_mean  /= q_total ;


   double q_err = weight_variance(hyper_copy, Q_prob_copy,
				  q_mean, ncopy) ;
   double latt_err = weight_variance(latt_coeff, Q_prob_copy,
				     latt_mean, ncopy) ;

   double mpi_sq_err = weight_variance(mpi_sq_coeff, Q_prob_copy,
				       mpi_sq_coeff_mean, ncopy) ;



  if( w_fit == A_SQ )
    {
      cout << "SUM Lattice spacing fit model = 1 + A * a**2\n" ;
    }
  else if( w_fit == A_ALPHA_S ) 
    {
      cout << "SUM Lattice spacing fit model = 1 + A * a * alpha\n" ;
    }
  else if( w_fit == A_ONLY ) 
    {
      cout << "SUM Lattice spacing fit model = 1 + A * a \n" ;
    }

  amin_tot /= ncopy ;

  dof_ = dof ;
  chisq_ = amin_tot ;

  ans = amin_tot / dof ; 
  cout << "SUM Chi**2 / dof = " << amin_tot << " / " << dof  << " = " << ans <<  "\n" ;

   cout << "SUM mc/ms = " << q_mean << " +/- " <<  q_err << " GeV\n" ;
   cout << "SUM latt_coeff = " << latt_mean << " +/- " << latt_err  << "\n" ;
   cout << "SUM mpi**2     = " << mpi_sq_coeff_mean << " +/- " << mpi_sq_err  << "\n" ;

   cout << "SUM double par[4] = { " << q_mean << " , " ;
   cout << latt_mean << " , " << mpi_sq_coeff_mean << " } ; \n" ;

   double ppp[3] = { q_mean , latt_mean , mpi_sq_coeff_mean }; 
  plot_func(ppp, w_fit) ; 

   cout <<  "SUM \n" ;
}


//
//  Main driver routine
//


void loop_boot_myfit_3param_2anal_B()
{
  int offset ; 

  double chisq ; 
  int dof ;
  double prob_cont ;


  load_data_init() ;
  const int n_anal = 2 ;
  const int dim = n_anal * ncopy ;

  //
  //  reserve memory
  //

  double *l_result = new double[ncopy]  ;
  double *l_weight   = new double[ncopy] ;
  TMatrixD *l_boot_result = new TMatrixD(ncopy, noboot) ;
  TMatrixD *l_boot_weight = new TMatrixD(ncopy, noboot) ;

  double *all_result = new double[dim]  ;
  double *all_weight   = new double[dim] ;
  TMatrixD *all_boot_result = new TMatrixD(dim, noboot) ;
  TMatrixD *all_boot_weight = new TMatrixD(dim, noboot) ;

  reserve_memory() ;

  //
  //  loop over analysis
  //


  single_anal(A_SQ,l_result, l_weight,
	      l_boot_result,l_boot_weight, chisq,dof) ;
  fclose(fp); 
  prob_cont = compute_Q(dof,chisq) ;

  for(int ii=0 ; ii < ncopy ; ++ii)
    {
      all_result[ii] = l_result[ii] ; 
      all_weight[ii]   = l_weight[ii] ; 
      for(int ib=0 ; ib < noboot ; ++ib)
	{
	  all_boot_result(ii,ib) = l_boot_result(ii,ib) ;
	  all_boot_weight(ii,ib) = prob_cont*l_boot_weight(ii,ib) ;

#if 0
	  cout << "DEBUG-all_boot_weight =   "  << all_boot_weight(ii,ib) << "\n" ;
	  cout << "DEBUG-l_boot_weight = " <<  l_boot_weight(ii,ib) << "\n" ;
#endif

	}
    }


  load_data_init() ;
  single_anal(A_ALPHA_S,l_result, l_weight,
	      l_boot_result,l_boot_weight,chisq,dof) ;
  fclose(fp); 
  prob_cont = compute_Q(dof,chisq) ;

  offset = ncopy ;
  for(int ii=0 ; ii < ncopy ; ++ii)
    {
      all_result[ii + offset] = l_result[ii] ; 
      all_weight[ii + offset] = l_weight[ii] ; 
      for(int ib=0 ; ib < noboot ; ++ib)
	{
	  all_boot_result(ii + offset,ib) = l_boot_result(ii,ib) ;
	  all_boot_weight(ii + offset,ib) = prob_cont*l_boot_weight(ii,ib) ;
	}
    }




  //
  //  save data for final combined analysis
  //
  //  TFile *MyFile = new TFile("boot.root","RECREATE") ;

  TFile *MyFile = new TFile("boot.root","UPDATE") ;
  string m_tag(input_tag) ;

  string all_boot_result_tag = "all_boot_result_" + m_tag ;
  MyFile->WriteTObject(all_boot_result,all_boot_result_tag.c_str()) ;

  string all_boot_weight_tag = "all_boot_weight_" + m_tag ;
  MyFile->WriteTObject(all_boot_weight,all_boot_weight_tag.c_str()) ;

  TVectorD  all_result_t(dim) ;
  TVectorD  all_weight_t(dim) ;
  for(int ii=0 ; ii < dim ; ++ii)
    {
      all_result_t(ii) = all_result[ii] ; 
      all_weight_t(ii) = all_weight[ii] ;
    }

  string all_result_tag = "all_result_" + m_tag ;
  string all_weight_tag = "all_weight_" + m_tag ;
  
  MyFile->WriteTObject(&all_result_t,all_result_tag.c_str()) ;
  MyFile->WriteTObject(&all_weight_t,all_weight_tag.c_str()) ;
  
  delete MyFile;

  //
  //  final summary
  //

  double ans     = weight_mean(all_result,all_weight,dim) ;
  double ans_err = weight_variance(all_result,all_weight,
				   ans, dim) ;

  cout << "\n\nSUM Final over all analysis\n" ;
  cout << "SUM Total number of ensembles " << ndata_total << "\n" ;
  cout << "SUM Cut on pion mass = " << m_pi_gev_cut << " GeV\n" ;
  cout << "SUM Cut on lattice spacing = " << a_fm_cut << " fm\n" ;
  cout << "SUM Number of ensembles in analysis " << ndata << "\n" ;
  cout << "SUM Number of copies = " << ncopy << "\n" ;
  cout << "SUM Number of analysis = " << n_anal << "\n" ;
  cout << "SUM Number of bootstrap samples = " << noboot << "\n" ;

  cout << "SUM mc/ms = " << ans << " +/- " << ans_err << "\n" ;



  //
  //   bootstrap analysis
  //
  double *all_result_tmp = new double[dim]  ;
  double *all_weight_tmp   = new double[dim] ;

  double * boot_tmp  = new double[noboot] ;

  for(int ib=0 ; ib < noboot ; ++ ib)
    {
      for(int ii = 0 ; ii < dim ; ++ii)
	{
	  all_result_tmp[ii] = all_boot_result(ii,ib) ;
	  all_weight_tmp[ii] = all_boot_weight(ii,ib) ;
	  //	  cout << "all_weight_tmp[" << ii << "] = " << all_boot_weight(ii,ib) << "\n" ;
	}
      boot_tmp[ib] = weight_mean(all_result_tmp,all_weight_tmp,dim) ;
    }

  double boot_err = boot_variance(boot_tmp, noboot) ;

  cout << "SUM Statistical errors = " << boot_err << "\n" ;

  //
  //  free up memory
  //

  delete [] all_result_tmp ;
  delete [] all_weight_tmp ;

  delete [] boot_tmp ;


  delete [] l_result ;
  delete [] l_weight ;

  delete [] all_result ;
  delete [] all_weight ;
    
  delete all_boot_result ;
  delete l_boot_result ;

  delete l_boot_weight ;
  delete all_boot_weight ;
    

  exit(0) ;
}



