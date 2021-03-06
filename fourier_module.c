#include "allvars.h"

/*************************************************************************
 NAME:       Dec_Alia_appr
 FUNCTION:   Deconvolves and eliminates Aliasing using Jing(2005) approx. 
 INPUTS:     magnitude k, parameters prm   
 RETURN:     SN
**************************************************************************/ 

double Dec_Alia_appr( double k ){

    double kn, SN;
    //Nyquist Frequency 
    kn = M_PI*PRM.Nc/(1.0*PRM.Lbox);

    #ifdef CIC
      SN = 1-2*sin(M_PI*k/(2.*kn))*sin(M_PI*k/(2.*kn))/3.;                      
    #endif  
    #ifdef TSC
	    SN = 1- sin(M_PI*k/(2.*kn))*sin(M_PI*k/(2.*kn))+ 2*pow(sin(M_PI*k/(.2*kn)),4)/15.;        
    #endif 

    return SN;
}

/*************************************************************************
 NAME:       DEC_WINDOW
 FUNCTION:   Deconvolution function, p=2 for CIC and p=3 for TSC
 INPUTS:     p, vector k, value obtained to deconvolve deconv, parameter prm
 RETURN:     0
**************************************************************************/ 

int Dec_window(  double kx, double ky, double kz, double *deconv ){

    int p;
    double kn, wx,wy,wz;
    //Nyquist Frequency 
    kn = M_PI*PRM.Nc/PRM.Lbox;
    
    #ifdef CIC
       p=2;
    #endif     
    #ifdef TSC
       p=3;
    #endif 

    //Finding the window function in FS 
    if (kx!=0) wx = sin(0.5*M_PI*kx/kn)/(0.5*M_PI*kx/kn);
    else wx = 1.;
    if (ky!=0) wy = sin(0.5*M_PI*ky/kn)/(0.5*M_PI*ky/kn);
    else wy = 1.;
    if (kz!=0) wz = sin(0.5*M_PI*kz/kn)/(0.5*M_PI*kz/kn);
    else wz = 1.;
    
    *deconv = pow( wx*wy*wz,p );

    return 0;
}

/*************************************************************************                                                                  
 NAME:      Deconv_Alias                                                                                                                    
 FUNCTION:   Deconvolution and aliasing                                                                                                     
 INPUTS:     magnitude k, parameters prm                                                                                                    
 RETURN:     ERR                                                                                                                            
**************************************************************************/

double Deconv_Alias( int N_cell, double kx, double ky, double kz ){

  double kn, ERR, fsin;
  //Nyquist Frequency                                                                                                                        
  kn = M_PI*PRM.Nc/PRM.Lbox;
  
  #ifdef NGP
    ERR = 1.0;
  #endif
  #ifdef CIC
    fsin = sin(0.5*M_PI*kx/kn);
    ERR = 1.0 - (2.0/3.0)*fsin*fsin;

    fsin = sin(0.5*M_PI*ky/kn);
    ERR = ERR*(1.0 - (2.0/3.0)*fsin*fsin);
  
    fsin = sin(0.5*M_PI*kz/kn);
    ERR = ERR*(1.0 - (2.0/3.0)*fsin*fsin);
  #endif 

  #ifdef TSC  
    fsin = sin(0.5*M_PI*kx/kn);                                                  
    ERR = 1- fsin*fsin+ 2.*pow(fsin,4)/15.;

    fsin = sin(0.5*M_PI*ky/kn);
    ERR = ERR*( 1- fsin*fsin+ 2.*pow(fsin,4)/15.); 

    fsin = sin(0.5*M_PI*kz/kn);       
    ERR = ERR*( 1- fsin*fsin+ 2.*pow(fsin,4)/15.);        
  #endif                                                                                                                                   

  return ERR;
}
/*************************************************************************
 NAME:       RTC
 FUNCTION:   A fourier transformed of a real entry array is done using FFTW
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/ 
		
int RTC( fftw_complex *FT_cd ){

    //Dimensions                                                                                                                              
    int rank = 3;
    //Size fftw array                                                                                                                         
    int dumb = floor(PRM.Nc/2)+1;
    //Dimensions Array                                                                                                                        
    int const n[3]= {PRM.Nc,PRM.Nc,PRM.Nc};
    //Declaring the plan                                                                                                                      
    fftw_plan plan;
    //Density contrast array                                                                                                                  
    double *in;
    in = (double *)calloc((size_t)PRM.NcTot,sizeof(double));
    int pos, i, j, k;

    //Creating plan                                                                                                                           
    plan = fftw_plan_dft_r2c( rank,n,in,FT_cd,FFTW_ESTIMATE );

    //Initializing In array                                                                                                                     
    for( i=0; i<PRM.Nc; i++){
	   for( j=0; j<PRM.Nc; j++){
	     for( k=0; k<PRM.Nc; k++){

    		pos = k + PRM.Nc*(j + i*PRM.Nc );
    		in[pos] = cells[pos].den_con;                                                                                                        

    		}
	    }
	}

    //Executing plan                                                                                                                          
    fftw_execute(plan);

    //Destroying plan
    fftw_destroy_plan(plan);
    //Releasing memory
    fftw_free(in);

    return 0;
}

/*************************************************************************
 NAME:       CTC
 FUNCTION:   A fourier transformed of a complex entry array is done using FFTW
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/ 
        
int CTC( fftw_complex *FT_cd ){

    //Dimensions                                                                                                                              
    int rank = 3;
    //Dimensions Array                                                                                                                        
    int const n[3]= {PRM.Nc,PRM.Nc,PRM.Nc};
    //Declaring the plan                                                                                                                      
    fftw_plan plan;
    //Density contrast array                                                                                                                  
    fftw_complex *in;
    in = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*PRM.NcTot );
    int pos,i,j,k;

    //Creating plan                                                                                                                         
    plan = fftw_plan_dft( rank,n,in,FT_cd,FFTW_FORWARD,FFTW_ESTIMATE );

    //Initializing In array                                                                                                                 
    for(i=0;i<PRM.Nc;i++){
     for(j=0;j<PRM.Nc;j++){
       for(k=0;k<PRM.Nc;k++){

            pos = k + PRM.Nc*(j + i*PRM.Nc);
            in[pos][0] = cells[pos].den_con; 
            in[pos][1] = 0;   
                                                                                                                                     }
        }
    }

    //Executing plan                                                                                                                       
    fftw_execute(plan);

    //Destroying plan
    fftw_destroy_plan(plan);
    //Releasing memory
    fftw_free(in);
                                                                  
    return 0;
}
/*************************************************************************
 NAME:       CTR
 FUNCTION:   Inverse Fourier transform of a complex array is done using FFTW.
             Test for RTC 
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/ 
    
int CTR( fftw_complex *FT_cd ){ 

  //Dimensions                                                                                                                              
  int rank = 3;
  //Dimensions Array                                                                                                                        
  int const n[3]= {PRM.Nc,PRM.Nc,PRM.Nc};
  //Declaring the plan                                                                                                                      
  fftw_plan pinv;
  //Density contrast array                                                                                                                  
  fftw_complex *out, *in;
  out = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*PRM.NcTot );
  in = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*PRM.NcTot );
  int pos,i,j,k;
  double sum,temp1;


  //Creating plan 
  pinv = fftw_plan_dft( rank,n,in,out,FFTW_BACKWARD,FFTW_ESTIMATE );
    
  //Initializing In array                                                                                                                
  for(i=0;i<PRM.Nc;i++){
    for(j=0;j<PRM.Nc;j++){
      for(k=0;k<PRM.Nc;k++){

         pos = k + PRM.Nc*(j + i*PRM.Nc);
         in[pos][0] = FT_cd[pos][0];
         in[pos][1] = FT_cd[pos][1];   

      }
    }
  }
  //Executing plan 
  fftw_execute(pinv);
    
  sum = 0;
  //Test info 
  for( k=0; k<PRM.Nc*PRM.Nc*PRM.Nc; k++){
    temp1 = fabs(out[k][0]/(PRM.NcTot)- cells[k].den_con);
    sum = sum + temp1; }            
  printf("\nTest FFTW: Squared difference %e\n\n",sum);  
  

  //Destroying plan
  fftw_destroy_plan(pinv);
  //Releasing memory 
 fftw_free(in);
 fftw_free(out);

  return 0;
}

/*************************************************************************
 NAME:       FS_Grid
 FUNCTION:   Filling cells with k modes and its contrast density
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/ 

int FS_Grid( fftw_complex *FT_cd ){

  int i,j,k,N_cell;
  double kx, ky, kz;
  
  //int dumb = (floor(PRM.Nc/2)+1)*PRM.Nc*PRM.Nc;

  g_k_max = 0;
  //Filling Structure with kx, ky, kz and density contrast k 
  for (i = 0; i < PRM.Nc; i++){
        for ( j = 0; j < PRM.Nc; j++){
             for( k= 0; k< PRM.Nc ;k++){    

                //Position of cell 
		N_cell = k+PRM.Nc*(j+PRM.Nc*i);
	        //Calculating kx
                if (i<PRM.Nc/2)
                    kx = 2*i*M_PI/(1.0*PRM.Lbox);  
                else
                    kx = 2*(-PRM.Nc+i)*M_PI/(1.0*PRM.Lbox);
                //Calculating ky
                if (j<PRM.Nc/2) 
                    ky = 2*j*M_PI/(1.0*PRM.Lbox);
                else
                    ky = 2*(-PRM.Nc+j)*M_PI/(1.0*PRM.Lbox);
                //Calculating kz
                if (k<PRM.Nc/2) 
                    kz = 2*k*M_PI/(1.0*PRM.Lbox); 
                else
                    kz = 2*(-PRM.Nc+k)*M_PI/(1.0*PRM.Lbox); 
                
		//Magnitud of k vector
                fcells[N_cell].k = sqrt( kx*kx + ky*ky + kz*kz );
 
		if ( fcells[N_cell].k > g_k_max ){ g_k_max =  fcells[N_cell].k;}

                //fcells[N_cell].r_dconk = FT_cd[N_cell][0];  
                //fcells[N_cell].i_dconk = FT_cd[N_cell][1]; 
                //fcells[N_cell].pk = (fcells[N_cell].r_dconk*fcells[N_cell].r_dconk+ fcells[N_cell].i_dconk*fcells[N_cell].i_dconk);
            	fcells[N_cell].pk = ( FT_cd[N_cell][0]*FT_cd[N_cell][0]+ FT_cd[N_cell][1]*FT_cd[N_cell][1] );
		fcells[N_cell].pk = fcells[N_cell].pk/(1.0*Deconv_Alias(N_cell,kx,ky,kz)); 
	   
	       }  
	    } 
      }
 
return 0;  
}

/*************************************************************************
 NAME:          Modes_FS
 FUNCTION:   Calculating the Power Spectrum 
 INPUTS:        array to save mean
 RETURN:       0
**************************************************************************/ 

int Modes_FS( double *mean ){

  int i, k, counter;
  double ps_norm = (1.0*PRM.Lbox*PRM.Lbox*PRM.Lbox)/pow(1.0*PRM.Nc,6);
  double shotnoise = PRM.vcell/(1.0*PRM.Npart);
      
  double kmax, kmin;
  kmin = PRM.kmin;
 
  for( k=0; k<PRM.Nbins; k++  ){
     
      kmax = kmin+PRM.deltak;
     
      counter=0;
      mean[k] = 0;
    
      for(i=0; i<PRM.NcTot; i++){
          
      
      if( (fcells[i].k >= kmin) && (fcells[i].k < kmax) ) {
          
          mean[k]  += fcells[i].pk;
          counter++;
        }
      }

      if (counter == 0) mean[k] = 0;
      else mean[k] = ps_norm*mean[k]/(1.0*counter) - shotnoise;
      
      kmin = kmax;
       
    }
 
  return 0; 
}

