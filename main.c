#include "allvars.h"

/*Usage:	
 	   ./main.out  
	   *	<Gadget snapshot file or ascci file> 
	   *	<parameter_file> 
./main.out /datastore/simulations/illutrisz0/snap_135 ./params.dat 
*/
   
int main( int argc, char *argv[] ){
  	 
    	//Defining variables
  	float prm[N_PARAMS];
  	char *filename, *param_filename; 
	char cmd[100];
	clock_t t;
        double time_taken; 
	t = clock();
    
  	//Reading snapshot and param filename 
	filename = argv[1];
	param_filename = argv[2]; 
	 	
	read_parameters( prm, param_filename );		
	
	//Reading particle positions  
	#ifdef HDF5
		read_HDF5( filename );
	#endif	
	#ifdef GADGET_BINARY 
		read_binary( filename ); 	 
	#endif	
	#ifdef ASCCI
		PRM.Lbox = prm[L];
		sprintf(cmd," wc -l < ./%s", filename);
		FILE *f = popen(cmd,"r");
		char path[100];

		while (fgets(path, sizeof(path)-1, f) != NULL) {
		     printf("\n\nLines in ascci input file: %s", path);
		     PRM.Npart = atoi(path);
		}
		
		pclose(f);
		read_ascci( filename );
	#endif  
        #ifdef FOF_PART
             PRM.Lbox = prm[L];
             PRM.Npart = prm[NP]*prm[NP]*prm[NP];
	     PRM.Nfiles = prm[NF];
	     PRM.global_mass = prm[MASS];
             read_FOF_PART( filename );   
        #endif
        
		
	//Parameters 		
	PRM.Nc = prm[N]  ;
	PRM.deltax = PRM.Lbox/(1.0*PRM.Nc);
	PRM.vcell = PRM.deltax*PRM.deltax*PRM.deltax;
	PRM.NcTot = PRM.Nc*PRM.Nc*PRM.Nc;
               
	//Array of Structures Cells, size N^3
	cells = (struct Cell *)calloc((size_t)PRM.NcTot, sizeof( struct Cell) );	
	if(cells==NULL){ 
            printf("Cells structure not allocated\n");
            exit(0);            
        }
		
	printf("\n%s read\n",filename);
	printf("\n Lbox %0.1f Npart %ld Ncells %d \n",PRM.Lbox, PRM.Npart, PRM.Nc);
	
	printf( "\n\n******************************** Mass Assignment  *********************************\n\n" ); 

	//The cells are created 		
	Grid( );    
	
	//Assignment of mass 
 	double mt_part = Mass_assignment(  );
        //Free memory
        free( parts );
 	
	//Density Contrast 
	density_contrast(mt_part);
	
	//Associated field stored per cell
	//write_DField( mt_part );
       
        //FT contrast density
	printf( "\n\n************************* Fourier transform ***************************\n\n" ); 

        //Size Out array FFTW  
	//int dumb = floor(PRM.Nc/2)+1;
	// FFTW array   
	fftw_complex *FT_cd;
	//(N/2 + 1)  positions of real FFT out array
	//FT_cd = fftw_malloc( sizeof(fftw_complex)*PRM.Nc*PRM.Nc*dumb );
	FT_cd = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*PRM.NcTot );

	//RTC( FT_cd );	
	CTC( FT_cd );
	//Free memory
        free( cells );

        // Test RFFT: 
	/* Comparing contrast density with the 
	   one obtained after FFT and inverse FFT */
  	#ifdef TEST_FT
	CTR( FT_cd );  
	#endif
       	    
                    //Array of Structures FCells, size N^3 
	fcells = ( struct FCell * )calloc( (size_t) PRM.NcTot, sizeof( struct FCell) );
	if(fcells==NULL){ 
            printf("Fcells structure not allocated\n");
            exit(0);
        }
	//Power Spectrum
	//Filling the k-space cells
	FS_Grid( FT_cd );  
		   
        //Free memory
        fftw_free( FT_cd );
        
        //Mean of contrast density module  
	double *mean;
        double kf = 2.0*M_PI/PRM.Lbox;
        PRM.deltak = kf;
        PRM.kmin = kf;//PRM.deltak;
	PRM.Nbins = (int) ceil( g_k_max/PRM.deltak );
	mean = (double *)calloc((size_t) PRM.Nbins,sizeof(double));

	printf("\n kmin %lf deltak %lf Nbins %d\n\n",PRM.kmin,PRM.deltak,PRM.Nbins);

	//Calculating the power spectrum
	Modes_FS( mean );
        //Free memory
        free( fcells );
	
	//Out file
	Write_PS( mean ); 
	
	//Free memory
	free( mean );

	//Time taken 
	t = clock() - t;
        time_taken = ((double)t)/CLOCKS_PER_SEC;
        printf("Time %f\n",time_taken/60);
        
  
return 0;
}       

 
