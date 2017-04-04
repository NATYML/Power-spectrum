#include "allvars.h"

/*************************************************************************
 NAME:       conf2dump
 FUNCTION:   Convert data file text in plain text 
 INPUTS:     Filename of parameter file 
 RETURN:     0
**************************************************************************/ 
  
int conf2dump( char filename[] ){
    char cmd[100];
    sprintf( cmd, 
    "grep -v \"#\" %s | grep -v \"^$\" | awk -F\"=\" '{print $2}' > %s.dump", 
	filename, filename );
    system( cmd );

    return 0;
}

/*************************************************************************
 NAME:       read_parameters
 FUNCTION:   Read parameter file and saves the information in param array
 INPUTS:     Parameter array 
 			 Filename of parameter file 
 RETURN:     0
**************************************************************************/ 

int read_parameters( float parameters[], char filename[] ) {

    char cmd[100], filenamedump[100];
    int i=0;
    FILE *file;

    //Load of File
    file = fopen( filename, "r" );
    if( file==NULL ){
		printf( "  * The file '%s' don't exist!\n", filename );
		exit(0);
		}
    fclose(file);
    
    //Converting to plain text
    conf2dump( filename );
    sprintf( filenamedump, "%s.dump", filename );
    file = fopen( filenamedump, "r" );
    
    //Reading
    while( getc( file ) != EOF ){
	fscanf( file, "%f", &parameters[i] );
	i++;}
    
    fclose( file );
    
    //printf( "  * The file '%s' has been loaded!\n", filename );

    sprintf( cmd, "rm -rf %s.dump", filename );
    system( cmd );
    
    return 0;
}

/*************************************************************************
 NAME:       read_binary
 FUNCTION:   Read a Gadget binary snapshot 
 INPUTS:     Filename of parameter file 
 			 Parameter array 
 RETURN:     0
**************************************************************************/ 

int read_binary( char *filename ){

	FILE *file = NULL;
	int i, j;
	int N_tot, N_min, N_max, dummy, nread=0;
	float Maux, faux[3];
	double dumb;
	unsigned int uintaux;
	
	//Load of File
	if(!(file = fopen( filename, "r" ))){
	  printf("can't open file `%s`\n", filename);
	  exit(0);
	}

	//printf( "\n\n***********************************************************************\n" ); 
	//printf("Reading snapshot %s\n", filename );	

	nread = fread( &dummy, sizeof(dummy), 1, file );
	nread = fread( &header, sizeof( struct gadget_header ), 1, file );
	nread = fread( &dummy, sizeof(dummy), 1, file );

	N_tot = 0;
	
	N_tot = header.npartTotal[1];
	printf("Type %d particles have Npart=%12d NpartTotal=%12d with mass %16.8lf\n", 1,
	    	    header.npart[1], header.npartTotal[1], header.mass[1]);

	PRM.Npart = N_tot;

	printf(" * Redshift... %16.8f\n",header.redshift);
	printf(" * Boxsize... %16.8f\n",header.BoxSize);
	printf(" * Omega0... %16.8f\n",header.Omega0);
	printf(" * OmageLa... %16.8f\n",header.OmegaLambda);
	printf(" * Hubbleparam... %16.8f\n",header.HubbleParam);
	
	PRM.Lbox = header.BoxSize;
      
		// Memory allocation for part variable
	parts = (struct Particle *) calloc( (size_t) N_tot,sizeof(struct Particle ) );

	if(parts == NULL){
	   printf("Structure particles could not be allocated\n");
	   exit(0);
	}
			// Positions 
	nread = fread( &dummy, sizeof(dummy), 1, file );

	for( i=0; i<N_tot; i++ ){

		nread = fread( &faux[0], sizeof(float), 3, file );
		parts[i].xp = faux[0];
		parts[i].yp = faux[1];
		parts[i].zp = faux[2];
	}

	nread = fread( &dummy, sizeof(dummy), 1, file );

	if( dummy != (3*N_tot*sizeof(float)) ){
	printf(" Can not read properly ids %d %lu\n",dummy,3*N_tot*sizeof(float));
	exit(0);
	}
				// Velocities 

	nread = fread( &dummy,sizeof(dummy),1,file );

	for( i=0; i<N_tot; i++ ){

	nread = fread( &faux[0],sizeof(float),3,file );
	dumb = faux[0];
	dumb = faux[1];
	dumb = faux[2];
	}

	nread = fread( &dummy, sizeof(dummy), 1, file );

	if(dummy != (3*N_tot*sizeof(float))){
		printf(" Can not read properly ids %d %lu\n", dummy, 3*N_tot*sizeof(float));
		exit(0);
	}

			// ID's 

	nread = fread(&dummy, sizeof(dummy), 1, file);

	for(i=0; i<N_tot; i++){
		nread=fread(&uintaux, sizeof(unsigned int), 1, file);
		dumb = uintaux;
	}

	nread = fread( &dummy, sizeof(dummy), 1, file );

	if(dummy != (N_tot*sizeof(unsigned int))){
		printf(" Can not read properly ids %d %lu\n", dummy, N_tot*sizeof(unsigned int));
		exit(0);
	}

			// Masses 

	nread = fread( &dummy, sizeof(dummy),1,file );

	N_min = N_max=0;

	for(j=0; j<=5; j++){

		N_max = N_max+header.npartTotal[j];

		if( (header.mass[j]==0) && (header.npartTotal[j]!=0)){
	  		for(i=N_min;i<N_max;i++){
				nread=fread(&Maux,sizeof(float),1,file);
			parts[i].mp = Maux;
	  		}
		}
		if((header.mass[j]!=0) && (header.npartTotal[j]!=0)){
		  for(i=N_min;i<N_max;i++){
			parts[i].mp = header.mass[j];
		  }
		}

	N_min=N_max;
	}

	nread = fread(&dummy,sizeof(dummy),1,file);

	fclose(file);
  
  return 0;
}

/*************************************************************************
 NAME:       read_ascci
 FUNCTION:   Read ascci file
 INPUTS:     Filename of parameter file 
 			 Number of particles variable 
 RETURN:     0
**************************************************************************/ 

int read_ascci( char *filename ){
  
  int idumb, ihalo;
  long int i;
  double dumb;
  FILE *pt=NULL;
  char string[100];
  
  //Open the data file  	
  if(!(pt = fopen( filename, "r" ))){
    printf("can't open file `%s`\n", filename);
    exit(0);
  }
  
  printf( "\n\n************************** Reading Ascci  ****************************\n" );
   
  parts = (struct Particle *)calloc( (size_t) PRM.Npart, sizeof( struct Particle) ); 
  if( parts==NULL ){printf("Particles structure could not be allocated \n");exit(0);}
  /*
  for(i=0; i<PRM.Npart; i++ ){
    fscanf( pt,"%lf %lf %lf %lf",&parts[i].xp, &parts[i].yp, &parts[i].zp, &parts[i].mp );  
  }*/
  
  //fscanf(pt,"%s",string);  
  for( i=0; i<PRM.Npart; i++ ){   
    fscanf( pt,"%d%*[,] %lf%*[,] %lf%*[,] %lf%*[,] %lf%*[,] %lf",&ihalo, &parts[i].xp, &parts[i].yp, &parts[i].zp, &parts[i].mp, &dumb);  
    parts[i].mp = parts[i].mp/(1.e10);
  
  if (parts[i].xp<0 ||parts[i].yp<0||parts[i].zp<0) printf("NEGATIVE VALUE\n");    
  }
  
  /*  for(i=PRM.Npart-10; i<PRM.Npart; i++ ){
    printf("%d %lf %lf %lf %lf %lf\n",ihalo, parts[i].xp, parts[i].yp, parts[i].zp, parts[i].mp, dumb);
  }
  getchar();
  */
  //printf("%d %lf %lf %lf %lf %lf\n",ihalo, parts[i-1].xp, parts[i-1].yp, parts[i-1].zp, parts[i-1].mp, dumb);
  
  /*
  parts = (struct Particle *)calloc( (size_t) PRM.Npart, sizeof( struct Particle) ); 
  if(parts==NULL){printf("Particles structure could not be allocated \n");exit(0);}
    
  for(i=0; i<PRM.Npart; i++ ){		
    
    fscanf(pt,"%d %lf %lf %lf %lf %lf %lf %d %lf %lf\n",
	      &ihalo, &parts[i].xp, &parts[i].yp, &parts[i].zp, &dumb,  &dumb,  &dumb, &idumb,&parts[i].mp, &dumb);
    
    //if(i==42611) {printf("%lf %lf %lf\n",parts[i].xp, parts[i].yp, parts[i].zp);
      
  }}   */ 
//printf("%e %e %e %e\n", parts[i].xp, parts[i].yp, parts[i].zp, parts[i].mp);
  
  
  fclose(pt);
  
  return 0;
}

/*************************************************************************
 NAME:       read_HDF5
 FUNCTION:   Read a HDF5 snapshot 
 INPUTS:     Filename of parameter file 
 			 Parameter array 
 RETURN:     0
**************************************************************************/ 

int read_HDF5( char *filename ){

	FILE *file = NULL;
	int i,j,counter,rank=2;
	double sdata[1][3],MT[6];
	char box_filename[80];
	int dumb,dumb1[6],dumb2[6],dumb0;

	hsize_t dim[2],count[2],offset[2],stride[2],block[2];
	hsize_t memspace_id, dataspace_id, dset_id;
	hid_t att, root, file_id, dataset_id;
	herr_t status;

	dim[0]=1,dim[1]=3;
	//Determines how many blocks to select from the 
	//dataspace, in each dimension
	count[0]=1,count[1]=3;
	//Steps beteween two elements
	stride[0]=1,stride[1]=1;
	block[0]=1,block[1]=1;

	printf( "\n\n***********************************************************************\n" ); 
	printf("Reading HDF5 snapshot %s\n", filename );	

	//Open an existing file
	sprintf(box_filename, "%s.%d.hdf5", filename,0);
	//Load of File
	if(!(file = fopen( box_filename, "r" ))){
	  printf("can't open file `%s`\n\n", box_filename);
	  exit(0);
	}
	file_id = H5Fopen(box_filename, H5F_ACC_RDONLY, H5P_DEFAULT);	

	root = H5Gopen(file_id, "Header",H5P_DEFAULT);

	//Reading total particle number
	att = H5Aopen(root,"NumPart_Total", H5P_DEFAULT); 
	status = H5Aread( att,H5T_NATIVE_INT,dumb1 );
	PRM.Npart = dumb1[1];
	//printf("%d\n",PRM.Npart);
	
	//Memory allocation for part variable
	parts = (struct Particle *) calloc( (size_t) PRM.Npart,sizeof(struct Particle ) );

	//Reading Lbox 
	att = H5Aopen_name(root, "BoxSize"); 
	status = H5Aread(att,H5T_NATIVE_INT, &dumb);
	PRM.Lbox = dumb;

	//Reading DM mass 
	att = H5Aopen_name(root, "MassTable"); 
	status = H5Aread(att,H5T_NATIVE_DOUBLE, MT);

	//Reading Num Files Per Snapshot
	att = H5Aopen_name(root, "NumFilesPerSnapshot");
	status = H5Aread(att,H5T_NATIVE_INT,&dumb0 );

	PRM.Nfiles = dumb0;
	
	j = 0;

	for ( counter = 0; counter < PRM.Nfiles; counter++ ){	

		if (counter != 0){
		   sprintf(box_filename, "%s.%d.hdf5", filename,counter);
		   file_id = H5Fopen(box_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		   root = H5Gopen(file_id, "Header",H5P_DEFAULT);
		   }
		
		//Reading Particles in this file
		att = H5Aopen(root,"NumPart_ThisFile", H5P_DEFAULT); 
		status = H5Aread( att,H5T_NATIVE_INT,dumb2 );
		//printf("%d\n",dumb2[1]);
		//parts = (struct Particle *) calloc( (size_t) PRM.Npart,sizeof(struct Particle ) );

		// Close the attribute
		status = H5Aclose( att );
		// Close the Group
		status = H5Gclose( root );	 

		//Open an existing dataset	
		dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
					   
		memspace_id = H5Screate_simple(rank, dim, NULL); 
		dataspace_id = H5Dget_space (dataset_id);
		
		
		for ( i = 0; i < dumb2[1]; i++ ){	
			
		  //Specifies the offset of the starting element 
		  //of the specified hyperslab.
		  offset[0] = i;
		  offset[1] = 0;
		
		  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, 
	    								 stride, dim, block);

		  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, 
	    								 H5P_DEFAULT, sdata);
		  parts[j].xp = sdata[0][0];
		  parts[j].yp = sdata[0][1];
		  parts[j].zp = sdata[0][2];
		  parts[j].mp = MT[1];
	   	  j++;
			
		}
		
		status = H5Sclose (memspace_id);
		status = H5Sclose (dataspace_id);
		//Close dataset
		status = H5Dclose( dataset_id );
		// Close the file
		status = H5Fclose(file_id);	
		
   } 

return 0;	
}


/*************************************************************************
 NAME:       write_DField
 FUNCTION:   Writes density field obtained from CIC algorithm 
 INPUTS:     Parameter array prm
 RETURN:     0
**************************************************************************/ 

int write_DField( ){


	int n,m,l,N_cell; 
	FILE *fp;
	char filename[80];

	//Create document to write out
	sprintf(filename, "../output_files/DF_%d.txt", PRM.Nc);
	fp = fopen(filename,"w"); 

	fprintf( fp, "#5 first raws: Divisions Grid N, Box length L(kpc), Size cells,\
			     #Number Particles, density_contrast\n");
	fprintf( fp, "#mass($10^{12}M_{sun}$) \n ");

	//Some parameters are stored
	fprintf( fp, "%d %lf\n%lf %d\n", 
		     PRM.Nc, PRM.Lbox, PRM.deltax, PRM.Npart );

	//Position and mass assigned is stored per cell
	for ( n = 0; n < PRM.Nc; n++){
		for ( m = 0; m < PRM.Nc; m++){
			for ( l= 0; l< PRM.Nc; l++){
				 N_cell = l+PRM.Nc*(m+PRM.Nc*n);
				 //fprintf( fp, "%lf \n",cells[N_cell].den_con);
				 fprintf( fp, "%lf %lf\n",cells[N_cell].mc,cells[N_cell].mngp);
			}
		}
	} 


	fclose(fp);

	return 0;
}
/*************************************************************************/
int Write_FT(  ){
  
	FILE *fp;
	int n, m, l, N_celda;

	//Create document to write out
	fp = fopen("../output_files/Test_FT/Fourier_transf.dat","w"); 

	//Some parameters are stored
	fprintf( fp, "%d %d %lf %lf %f\n", 
		     PRM.Nc, PRM.Nbins, PRM.deltak, PRM.kmin, 1.0);

	//FT is stored per shell
	for ( n = 0; n < PRM.Nc; n++){
		for ( m = 0; m < PRM.Nc; m++){
			for ( l= 0; l< PRM.Nc; l++){
				 N_celda = l+PRM.Nc*(m+PRM.Nc*n);
				 fprintf( fp, "%lf %lf %lf %lf %lf\n",
				 	fcells[N_celda].kx, fcells[N_celda].ky, fcells[N_celda].kz,
					fcells[N_celda].r_dconk, fcells[N_celda].i_dconk );	 		

			}
		}
	} 
	
	fclose(fp);

	return 0;
}

/*************************************************************************
 NAME:       Write_PS
 FUNCTION:   Write output file with the Power Spectrum calculated
 INPUTS:     Mean array
 			 Variance array 
 			 Parameter array prm
 RETURN:     0
**************************************************************************/ 
				
int Write_PS( double *mean ){
  
	int N_celda; 
	FILE *fp;
	int n;
	char file[80];
        double k;

	//Create document to write out
	sprintf(file, "../output_files/PS_%.2lf_kf.dat",PRM.deltak*PRM.Lbox/(2.*M_PI));
	//fp = fopen("../output_files/PS.dat","w");
	fp = fopen( file,"w" );

	//Media is stored per shell
	for ( n = 0; n < PRM.Nbins; n++){

	        //k = PRM.kmin + (n+0.5)*PRM.deltak;
	        k = PRM.kmin + n*PRM.deltak;
                fprintf( fp, "%lf %lf\n", k ,mean[n] );
		 
	} 

	fclose(fp);

	return 0;
}
