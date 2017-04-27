#ifndef ALLVARS_H
#define ALLVARS_H

/**********************Macros**************************/

#define NMAX 40000000 
#define N_PARAMS 100
/*-------------------------*/ 
#define N 0
#define L 1
#define NP 2
#define NF 3
#define MASS 4

#define CENTER 1 
#define I 1
#define J 2
#define K 3
/*-------------------------*/
extern double g_k_max;

/********************Structures************************/

extern struct Particle{
	
	double xp;
	double yp;
	double zp;
	double mp;

}*parts;    
  
extern struct Cell{

	double xc;
	double yc;
	double zc;
	//double mc;
        //double mngp;
	double den_con;

}*cells;

extern struct FCell{
	
	//double kx;
	//double ky;
	//double kz;
	double k;
	//double r_dconk;
	//double i_dconk;
        double pk;

}*fcells;

extern struct param{

	double Lbox;
        int Nfiles;
	long int Npart;
	int Nc;
	int Nbins;
	int NcTot;
	double deltax;
	double vcell;
	double deltak;
	double kmin;
        double global_mass;

}PRM;

extern struct gadget_header{

	  int npart[6];
	  double mass[6];
	  double time;
	  double redshift;
	  int flag_sfr;
	  int flag_feedback;
	  int npartTotal[6];
	  int flag_cooling;
	  int num_files;
	  double BoxSize;
	  double Omega0;
	  double OmegaLambda;
	  double HubbleParam;
	  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */

}header;


/********************Libraries Used********************/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include <string.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <time.h>
#include <algorithm>
//#include "hdf5.h"
// Asus
///#include "/home/nataly/Programs/FOFReaderLib-master/FOFReaderLib/FOFReaderLib.h"
// Store    
#include "/scratch/nmateusl/local/src/FOFReaderLib-master/FOFReaderLib/FOFReaderLib.h"
#include "proto.h"

#endif 
