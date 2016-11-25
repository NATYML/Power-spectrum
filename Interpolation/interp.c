#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


int main ( int argc,  char **argv ) {
    
    int i,size;
    double ki, pki, kmin, kmax, deltak, Niter,Lbox;
    char *file, cmd[100];
    FILE *pt, *fl;
    
    file = argv[1];
    kmin = atof( argv[2] );
    Niter = atoi( argv[3] );
    Lbox = atof( argv[4] );
    deltak = M_PI/Lbox;
    kmax = kmin+deltak*Niter;
    
    pt = fopen( file,"r" );
    
    sprintf( cmd," wc -l < ./%s", file );
    FILE *f = popen( cmd,"r" );
    char path[100];

    while ( fgets(path, sizeof(path)-1, f ) != NULL) {
            printf("\nLines in ascci input file: %s\n", path);
            size = atoi(path)-1;
    }
    
    pclose(f);
    
    double k[size], pk[size];
    
    for (i=0; i<size;i++){
      
        fscanf(  pt, "%lf %lf", &k[i],&pk[i]  );
        
    }
    
    if ( kmax> k[size-1] ) {
        
        printf("%lf %lf %lf\n",deltak,kmax, k[size-1]);
        printf("Value of k out of range of interpolation\n");
        exit(0);
        
    }   
        
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc ( gsl_interp_cspline, size );

    gsl_spline_init ( spline, k, pk, size );

    fl = fopen( "pk.dat","w" );
    
    for ( i = 0; i < Niter; i++ ){
    
        ki = kmin+deltak*i;
        pki = gsl_spline_eval ( spline, ki, acc );
        fprintf ( fl,"%lf %lf\n", ki, pki );
        
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    fclose( fl );
    
    return 0;
    }