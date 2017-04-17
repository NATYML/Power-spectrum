#ifndef ALLVARS_H
#include "allvars.h"
#endif

// IO Module
int conf2dump( char *);
int read_parameters( float *, char *);
int read_ascci( char * );
int read_binary( char * );
int read_FOF_PART( char * );
int read_HDF5( char * );
int write_DField( double );
int Write_PS( double * );

// Density map Module
int Grid( );
int locate_cell( double, double, double, int * );
int CIC_Wfunction( double , double * );
int TSC_Wfunction( double , double * );
double Mass_assignment( );
int density_contrast(double );

// Fourier Module 
double Dec_Alia_appr( double );
int Dec_window(  double, double, double, double * );
double Deconv_Alias( int , double, double ,double );
int RTC( fftw_complex * );
int CTC( fftw_complex * );
int CTR( fftw_complex * );
int FS_Grid( fftw_complex * );
int Modes_FS( double * );
