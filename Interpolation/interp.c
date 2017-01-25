#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>


int main ( int argc,  char **argv ) {
    
    int i, size, nk;
    double pki, dumb;
    char *file, *file2, cmd[100];
    FILE *pt, *fl;
    
    file = argv[1];
    file2 = argv[2];
 
    /* Obtaining numer of lines of file 1,
     * file to interpolate*/
    
    sprintf( cmd," wc -l < ./%s", file );
    FILE *f = popen( cmd,"r" );
    char path[100];

    while ( fgets(path, sizeof(path)-1, f ) != NULL) {
            printf("\nLines in ascci input file: %s\n", path);
            size = atoi(path);
    }
    
    pclose(f);
    
    double k[size], pk[size];
    
    /*Reading file to perform the interpolation*/
    pt = fopen( file,"r" );
    
    for (i=0; i<size;i++){
      
        fscanf(  pt, "%lf %lf", &k[i],&pk[i]  );
        
    }
    pclose(pt);
    
    /* Obtaining numer of lines of 
     * file containing k values*/
  
    sprintf( cmd," wc -l < ./%s", file2 );
    f = popen( cmd,"r" );
    char path2[100];
    
    while ( fgets(path2, sizeof(path2)-1, f) != NULL) {
        
            printf("\nLines in k input file: %s\n", path2);
            nk = atoi(path2);
    }
    pclose(f);
    
    /*Reading file to obtain k values to perform interpolation*/
    double ki[nk];
    FILE *pt2 = fopen( file2,"r" );
    
    for (i=0; i<nk;i++){
      
        fscanf(  pt2, "%lf %lf", &ki[i],&dumb  );
        
    }
    pclose(pt2);
    
    /* Performing interpolation*/
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc ( gsl_interp_cspline, size );

    gsl_spline_init ( spline, k, pk, size );

    /* Writing interpolated values*/
    fl = fopen( "pk_interp.dat","w" );
    
    for ( i = 0; i < nk; i++ ){
    
        pki = gsl_spline_eval ( spline, ki[i], acc );
        fprintf ( fl,"%lf %lf\n", ki[i], pki );
        
    }
    printf("Interpolation performed\n\n");
    
    // Realising memory 
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return 0;
    }
