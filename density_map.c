#include "allvars.h"
 
/*************************************************************************
 NAME:       Grid
 FUNCTION:   Construct the grid to get cells coordinates
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/	

int Grid( ){

int i, j, k, N_cell;
int Niter = PRM.Nc;

for ( i = 0; i < Niter; i++){

	for ( j = 0; j < Niter; j++){
		for ( k = 0; k< Niter;k++){

			//index
			N_cell = k+Niter*(j+Niter*i);

			cells[N_cell].xc = i*PRM.deltax;
			cells[N_cell].yc = j*PRM.deltax;
			cells[N_cell].zc = k*PRM.deltax;
			}
		}
	} 

return 0;
}

/*************************************************************************
 NAME:       locate_cell
 FUNCTION:   Find Cell where particle is located
 INPUTS:     position of the particle to be located on the cell
 			 id of particle
 			 variable to save the id of cell 
 			 parameters prm 
 RETURN:     0
**************************************************************************/	

int locate_cell(double xp, double yp, double zp, int *index ){

       //Numbers in order to find the cell
       //Counters in x,y y z
       double tol = 1e-4;

	// i, j, k
       if ( fabs(xp-PRM.Lbox) <=tol )
           index[I] = PRM.Nc-1;
       else 
	  index[I] = floor( (xp / PRM.Lbox) * PRM.Nc);
       
       if (fabs(yp-PRM.Lbox)<=tol)
	  index[J] = PRM.Nc -1;
       else
          index[J] = floor( (yp / PRM.Lbox) * PRM.Nc);

       if (fabs(zp-PRM.Lbox)<=tol)
	  index[K] = PRM.Nc -1;
       else
	 index[K] = floor( (zp / PRM.Lbox) * PRM.Nc );

	//Index of cell where particle falls
	//k + N( j + N i )
	index[0] = index[K]+PRM.Nc*( index[J]+PRM.Nc*index[I] );

	return 0;
	} 

/*************************************************************************
 NAME:       CIC_Wfunction 
 FUNCTION:   Cloud in cell window
 INPUTS:     
 RETURN:     0
**************************************************************************/	

int CIC_Wfunction( double x, double *W ){

	//Cell where particle falls
	W[CENTER] = 1 - fabs(x)/PRM.deltax;

	//Left Cell
	if ( fabs(x + PRM.deltax) < PRM.deltax )
	   W[0] = 1 - fabs(x+ PRM.deltax)/PRM.deltax;	

	//Right Cell
	if ( fabs( - x + PRM.deltax)< PRM.deltax )
	   W[2] = 1 - fabs(-x + PRM.deltax)/PRM.deltax;		

return 0;
}

/*************************************************************************
 NAME:       TSC_Wfunction 
 FUNCTION:   Triangular shape cloud
 INPUTS:     
 RETURN:     0
**************************************************************************/	

int TSC_Wfunction( double x, double *W ){

	//Cell where particle falls
	W[CENTER] = 0.75 - (x*x)/(PRM.deltax*PRM.deltax);
	
	//Left cell
	W[0] = 0.5*(1.5 - fabs(x + PRM.deltax)/PRM.deltax)*\
		   (1.5 - fabs(x + PRM.deltax)/PRM.deltax);		
	//Right Cell
	W[2] = 0.5*(1.5 - fabs(x - PRM.deltax)/PRM.deltax)*\
		   (1.5 - fabs(x - PRM.deltax)/PRM.deltax);		

return 0;
}
/*************************************************************************
 NAME:       Mass_assignment 
 FUNCTION:   
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/	

double Mass_assignment(  ){
	
	int i, l, k, p, q, r, pos, index[4];
	int i0, j0, k0, counter[3] = {-1,0,1};
	double Xc, Yc, Zc;
	double Wx[3],Wy[3],Wz[3],mt_cells,mt_part= 0;
        
        #ifdef CIC 
		 printf("CIC scheme chost \n\n");   
        #endif 
        #ifdef TSC
                 printf("TSC scheme chost \n\n");   
        #endif           
	
	double temp;
 	for ( i = 0; i < PRM.Npart; i++){					

	        mt_part = parts[i].mp + mt_part ;	

	        //locates cell where particle i falls in 
	        locate_cell( parts[i].xp, parts[i].yp, parts[i].zp, index );
	
	        //Nearest grid point scheme
	        //cells[index[0]].mngp = cells[index[0]].mngp + parts[i].mp;	
		  
	       	
	        //Initializing window functions	
	        for ( l = 0; l < 3; l++){ Wx[l]=0; Wy[l]=0; Wz[l]=0; }	
		
		/* Cell center */
		Xc = cells[index[0]].xc + PRM.deltax*0.5;
		Yc = cells[index[0]].yc + PRM.deltax*0.5;
		Zc = cells[index[0]].zc + PRM.deltax*0.5;
      	
		//Finding contribution of particle i to the cell's mass
		//index[0]: Index of cell where particle falls in 

			//Reading particle positions  
	
		#ifdef CIC 
		{
		 
		CIC_Wfunction( parts[i].xp - Xc, Wx );
		CIC_Wfunction( parts[i].yp - Yc, Wy ); 
		CIC_Wfunction( parts[i].zp - Zc, Wz );

		}
		#endif
		#ifdef TSC 
		{
                    
		TSC_Wfunction( parts[i].xp - Xc, Wx );
		TSC_Wfunction( parts[i].yp - Yc, Wy ); 
		TSC_Wfunction( parts[i].zp - Zc, Wz );
                
		}
		#endif 
	
		//Contribution to the cells mass to due to particle that falls in 
		cells[index[0]].mc = cells[index[0]].mc + Wx[CENTER]*Wy[CENTER]*Wz[CENTER]*parts[i].mp;		
                //printf("%lf\n",cells[index[0]].mc);
		
		temp = 0;
		//Neighbour cells	
		for ( p = 0; p <3; p++ ){
				for ( q = 0; q <3; q++ ){
						for ( r = 0; r <3; r++ ){					
							if ( (p==1 && q==1 && r==1) ) continue;

							//Index of cell 
							i0 = index[I]+ counter[p];
							j0 = index[J]+ counter[q];
							k0 = index[K]+ counter[r];

							//Periodical boundary conditions 
							//Neighbor in end of box
							if ( i0 >= PRM.Nc ) i0 = 0;
							if ( j0 >= PRM.Nc ) j0 = 0;
							if ( k0 >= PRM.Nc ) k0 = 0;

							//Neighbor in the beginning of box
							if ( i0 < 0 ) i0 = PRM.Nc-1;
							if ( j0 < 0 ) j0 = PRM.Nc-1;							
							if ( k0 < 0 ) k0 = PRM.Nc-1;
							
							//Calculating position
							pos = ( i0*PRM.Nc + j0 )*PRM.Nc + k0;
							cells[pos].mc = Wx[p]*Wy[q]*Wz[r]*parts[i].mp + cells[pos].mc;   					                                    
							temp = temp +  Wx[p]*Wy[q]*Wz[r]*parts[i].mp;
							//if (Wx[p]*Wy[q]*Wz[r]>0 ) printf(" %lf  \n",Wx[p]*Wy[q]*Wz[r] );   
				   		}
				}
					
		}
			
                //if( (parts[i].mp - temp) < 0 ){ printf("Problem: part id %d Mp %lf Mcells %lf Mass cc%lf\n",i,parts[i].mp,temp,cells[index[0]].mc);getchar(); }
                double borr;
                borr = Wx[CENTER]*Wy[CENTER]*Wz[CENTER]*parts[i].mp + temp;
                if( fabs(parts[i].mp - borr) > 1e-4 ){ printf("More mass! \n"); getchar(); }
		
	}
	
	//Test for cloud in cell
	#ifdef TEST_MASS 
	mt_cells=0;
        for( k=0; k<PRM.Nc*PRM.Nc*PRM.Nc;k++ ) mt_cells = mt_cells + cells[k].mc; 
 		printf("Cells total mass %lf \t Part total mass %lf \n\n", mt_cells, mt_part );
		printf("Relative Error Mass cells %e\n",fabs(mt_cells-mt_part)/mt_part);	
	#endif
		
return mt_part;
}

/*************************************************************************
 NAME:       density_contrast
 FUNCTION:   Using the mass assignated to every cell to calculate 
 			 the contrast density
 INPUTS:     parameters prm 
 RETURN:     0
**************************************************************************/	 

int density_contrast( double mt_part ){
  
  int i; 
  double bck_den;	
  
  bck_den = mt_part/(PRM.Lbox*PRM.Lbox*PRM.Lbox);
    
  //Density contrast 
  for(i= 0; i<PRM.NcTot; i++ ){
   
    //#ifdef NGP
    //cells[i].den_con = (cells[i].mngp/PRM.vcell)/bck_den - 1.0;	
    // #endif 
    #ifdef CIC
      cells[i].den_con = (cells[i].mc/PRM.vcell)/bck_den - 1.0;	 
    #endif  
    #ifdef TSC
      cells[i].den_con = (cells[i].mc/PRM.vcell)/bck_den - 1.0;	 
    #endif  

  } 
  
  return 0;	
}
