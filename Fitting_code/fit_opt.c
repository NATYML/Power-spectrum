#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

/* number of fit coefficients */
//#define NCOEFFS  25
#define NCOEFFS  14

/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK   (NCOEFFS - 2)

/*  Bigger k  */
#define KBREAK 0.1

int main (  int argc, char *argv[]  ){
    
    int N, value, cont, cont2;
    const size_t ncoeffs = NCOEFFS;
    const size_t nbreak = NBREAK;
    size_t i, j;
    gsl_bspline_workspace *bw;
    gsl_vector *B;
    double dy, Bj;
    gsl_rng *r;
    gsl_vector *c, *w;
    gsl_vector *x, *y;
    gsl_matrix *X, *cov;
    gsl_multifit_linear_workspace *mw;
    double chisq, Rsq, dof, tss;
    char *filename, cmd[100], path[100];
    
    filename = argv[1];
    
    // Obtaining number of lines or file
    sprintf( cmd," wc -l < ./%s", filename );
    FILE *f = popen( cmd,"r" );
    
    while (  fgets( path, sizeof(path)-1, f ) != NULL  ) {
            //printf( "\n\nLines in ascci input file: %s\n", path );
             N = atoi( path );
    }
    
    pclose(f);
    
    FILE *fp;
    double k[N], Pk[N];
    
    //Open the data file  	
    if( !( fp = fopen( filename, "r" )) ){
        printf( "can't open file `%s`\n", filename );
        exit(0);
    }

    //Reading data from file
    for (   cont = 0; cont < N; cont++  ){
            fscanf(  fp, "%lf %lf", &k[cont] ,&Pk[cont] );
            //printf(  "%lf %lf\n", k[cont] ,Pk[cont] );
    } 

    fclose(fp);
    
    #ifdef SUBSAMP
    
             int step, n, pos;
             double L, kN;
             
             L = atof( argv[2] );
             step = atoi( argv[3] );
                
            kN = 2*M_PI/L ; 
            // n = (k0- kmin)/deltak
            n = ( k[0] - 0.5*kN )/( step*kN ) ;
             
            if(n%2==0) pos = 2;
            else pos=1;
                 
    #endif
   
                            // Performing fitting 
    
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* allocate a cubic bspline workspace (k = 4) */
    bw = gsl_bspline_alloc( 4, nbreak );
    B = gsl_vector_alloc( ncoeffs );
    /* number of data points to fit */
    const size_t n = N;

    x = gsl_vector_alloc( n );
    y = gsl_vector_alloc( n );
    X = gsl_matrix_alloc( n, ncoeffs );
    c = gsl_vector_alloc( ncoeffs );
    w = gsl_vector_alloc( n );
    cov = gsl_matrix_alloc( ncoeffs, ncoeffs );
    mw = gsl_multifit_linear_alloc( n, ncoeffs );
            
                
    int LEN = 50;
    double xmin, xmax;
    xmin = 1.0, xmax = 100;
    double deltax = (xmax-xmin)/(1.0*LEN);
    double porc[LEN], chi2[LEN];
    
    for ( cont = 0; cont< LEN; cont++){
        
                //printf( "ITER %d\n\n",cont);   
                porc[cont] = deltax*( cont+1.0 );

                /* this is the data to be fitted */
                
                for ( i = 0; i < n; i++ ){
                    
                    gsl_vector_set( x, i, k[i] );
                    gsl_vector_set( y, i, Pk[i] );
                    
                    if ( k[i] > KBREAK )  value = xmin + porc[cont];
                    else {
                        
                       // #ifdef DEF SUBSAMP
                                
                        
                        //#else  
                                value = 1; 
                        //#endif
                    }    
                    gsl_vector_set(  w, i, value  );

                    //printf( "%f %f %d\n", k[i], Pk[i] ,value);
                    
                }
             
                /* use uniform breakpoints on [kmin, kmax] */
                gsl_bspline_knots_uniform( k[0], k[N-1], bw );

                /* construct the fit matrix X */
                for ( i = 0; i < n; i++ ){
                    
                    double xi = gsl_vector_get( x, i );

                    /* compute B_j( xi) for all j */
                    gsl_bspline_eval( k[i], B, bw );

                    /* fill in row i of X */
                    for ( j = 0; j < ncoeffs; j++ ){
                        
                        Bj = gsl_vector_get( B, j );
                        //printf("%lf\n",Bj);
                        gsl_matrix_set( X, i, j, Bj );
                    }
                }

                /* do the fit */
                gsl_multifit_wlinear( X, w, y, c, cov, &chisq, mw );

                dof = n - ncoeffs;
                tss = gsl_stats_wtss( w->data, 1, y->data, 1, y->size );
                Rsq = 1.0 - chisq / tss;
                
                //fprintf( stderr, "chisq/dof = %e, Rsq = %f\n", chisq / dof, Rsq );
                
                FILE *pt;
                char filename[80];

                //Create document to write out
                sprintf( filename, "./pk_fit/pk_w%1.2f.txt", (porc[cont]+1) );
                pt = fopen(filename,"w"); 

                /* output the smoothed curve */
                double xi, yi, yerr;

                //printf("#m=1,S=0\n");
                /*for (xi = k[0]; xi < k[N-1]; xi += 0.001){
                    gsl_bspline_eval(  xi, B, bw  );
                    gsl_multifit_linear_est(  B, c, cov, &yi, &yerr  );
                    fprintf( pt , "%f %f\n", xi, yi  );
                }*/
                
                chi2[cont] = 0;
                for ( cont2=0; cont2 < N; cont2++ ){
                    
                    gsl_bspline_eval(  k[cont2], B, bw  );
                    gsl_multifit_linear_est(  B, c, cov, &yi, &yerr  );
                    fprintf( pt , "%f %f\n", k[cont2], yi  );
                    if ( k[cont2]>KBREAK ){
                          chi2[cont] = chi2[cont] + ( yi - Pk[cont2] )*( yi - Pk[cont2] ); 
                    }
                }
                printf("\n w %lf xi2 %lf  \n",  (porc[cont]+1), chi2[cont]);
                
                fclose(pt);
    }          
    
    double minvalue = chi2[0];
    
    for( cont=0; cont< LEN; cont++ ){
            if ( chi2[cont] < minvalue  )
                minvalue = chi2[cont];
    }
    printf("\nMin value %lf\n",minvalue);

    
    gsl_rng_free(r);  
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(w);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);
    
    return 0;
    
} 

