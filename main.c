/***********************************************************************
 *
 *
 *
 *	cGRILL. A simple affinity map generator written in C
 *
 *
 *	
 *	(c) 2011-2024. Alvaro Cortes Cabrera. <alvarocortesc@gmail.com>
 *
 *	 This code is licensed under GPL v3 terms
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

/* This is for reading molecules and other stuff */
#include <reader.c>
#include <waste.c>
#include <metals.c>

void print_header();
void print_usage();

int main(argc, argv)
int argc;
char *argv[];
{

	/* Parameters block */
	static int verbose_flag = 0;
	float center_x=0.0f, center_y = 0.0f, center_z = 0.0f;    /* Coordinates of the center of the grid */
	int size_x = 0, size_y = 0, size_z = 0;                   /* _Total_ size of the grid */
	float spacing = 0.5; 					                  /* Grid spacing in A */
	char *protein_name = NULL;
	char c;
	int pqr_type = 0;
	/* End parameters block */


	/* Molecule processing block */
        double mytot = 0;
        int ringers = 0;
        MOL2 *mymol = NULL;
	float vertex[8][3];
        /* End molecule processing block */


	/* Grid calcs and processing block     */

	float max_grid[3];  /* Maximum coordinates of the grid */
	float min_grid[3];  /* Minimum coordinates of the grid */
	int i = 0, j = 0, k = 0, l = 0;
	int x_points = 0, y_points = 0, z_points = 0;
	float point[3];
	float ****grids; /* 4D dynamic array for n grids */
	float dx = 0.0f, dy = 0.0f, dz = 0.0f;
	float dist = 0.0f;
	float dist2 = 0.0f, distr = 1.0, tmpdist = 0.0f;
	float tmpE = 0.0f, sigma = 0.0f, epsilon = 0.0f, r6 = 0.0f, r12 = 0.0f, sigma6 = 0.0f;
	float *all_points = NULL;
	int counter_points = 0;
	FILE *output[7];
	static char *output_names[7] = {"NH4.dx","O.dx","OH.dx","CH3.dx","H2O.dx","Hydrophobic.dx","S.dx"};
	int neig_index = 0, neig_j = 0;
	float heavy_nb[3];
	float alfa = 0.0f;
	float heavy_hb2[3];
	float vector1[3], vector2[3];
	float modulus1 = 0.0f, modulus2 = 0.0f;

        /* End grid calcs and processing block */

        /* "New" hydrogen bond term */

        float h_dist[3] = {3.2f,3.0f,2.8f};
        float h_ener[3] = {-2.0f,-2.8f,-4.0f};
	float C_hbond[3];
	float D_hbond[3];

	float univec[3]; /* Normal vector of the plane */
        float projection[3]; /* Probe projection onto the lone pairs plane */
	float tmpvec1[3]; /* Temporal vector 1 */
        float tmpvec2[3]; /* Temporal vector 2 */
        float tmpvec3[3]; /* Temporal vector 3 */
	float dtoplane = 0.0f; /* Distance to the plane */
	float norm = 0.0f;
	float tmpd1 = 0.0f;
	float tmpd2 = 0.0f;
        float tmpd3 = 0.0f;
        float tmpd4 = 0.0f;
	float tmpd5 = 0.0f;
	float k1 = 0.9f/ (powf( cos(1.9198),6.0f));
	float k2 = powf( cos(1.9198),2.0f);

	/* End of the hydrogen bond block */

	/* New filtering for hotspots */

	float **sel_points = NULL;
/*        float thresholds[7] = { -12.00, -4.50, -7.00, -999.9, -999.9, -1.7f, -999.9 };*/
        float thresholds[7] = { -12.00, -6.00, -7.00, -999.9, -999.9, -1.7f, -999.9 };
        int current_point = 0, total_points = 0;
        int best_grid = 0;
	float best_energy = 0.0f;
	FILE *hotspots = NULL;
	char *hots[7] = { "N", "O", "F", "X", "X", "S", "X" };
	char *hots_res[7] = { "POS", "HBA", "HBD", "X", "X", "HPH", "X" };
	int ch = 0;
	

	/* End of hotspots filtering */

        /************************************ //
 	*                                     //
	*                                     //
	*   PART A. PARAM PROCESSING A CHECK  //
	*                                     //
	*                                     //
	*                                     //
	**************************************/

	for( i = 0; i < 3; i++)
	{
		C_hbond[i] = -3.0f*h_ener[i]*powf(h_dist[i],8.0f);
                D_hbond[i] = -4.0f*h_ener[i]*powf(h_dist[i],6.0f);
	}
	


	print_header();

         while (1)
         {
           static struct option long_options[] =
             {
               {"verbose", no_argument,       &verbose_flag, 1},
               {"brief",   no_argument,       &verbose_flag, 0},
               {"spacing",   required_argument, 0, 's'},
               {"protein",   required_argument, 0, 'p'},
               {"center_x",  required_argument, 0, 'x'},
               {"center_y",  required_argument, 0, 'y'},
               {"center_z",  required_argument, 0, 'z'},
               {"size_x",    required_argument, 0, 'h'},
               {"size_y",    required_argument, 0, 'w'},
               {"size_z",    required_argument, 0, 't'},
               {"pqr",       no_argument,       0, 'q'},
               {0, 0, 0, 0}
             };
           int option_index = 0;
     
           c = getopt_long (argc, argv, "p:x:y:z:s:h:w:t:q",
                            long_options, &option_index);
           if (c == -1)
             break;
     
           switch (c)
             {
             case 0:
               /* If this option set a flag, do nothing else by now. */
               if (long_options[option_index].flag != 0)
                 break;
               printf ("option %s", long_options[option_index].name);
               if (optarg)
                 printf (" with arg %s", optarg);
               printf ("\n");
               break;
     
             case 'q':
                pqr_type = 1;
                break;
             case 'p':
		protein_name = optarg;
                break;
     
             case 'x':
		center_x = atof(optarg);
               break;
     
             case 'y':
                center_y = atof(optarg);
               break;
     
             case 'z':
                center_z = atof(optarg);
               break;
     
             case 's':
                spacing = atof(optarg);
               break;

             case 't':
                size_z = atoi(optarg);
               break;

             case 'w':
                size_y = atoi(optarg);
               break;

             case 'h':
                size_x = atoi(optarg);
               break;

             case '?':
               /* getopt_long already printed an error message. */
               break;
     
             default:
               abort ();
             }
         }


	if( size_x <= 0 || size_y <= 0 | size_z <= 0 || protein_name == NULL || spacing <= 0.0)
	{
		print_usage();
		fprintf(stderr,"ERROR-MAIN%i: Problem with parameters.\n\n",__LINE__);
		fprintf(stderr,"\tParameters: Box size: %ix%ix%i:%f for protein %s.\n\n",size_x,size_y,size_z,spacing,protein_name);
		fprintf(stderr," *** END ERROR *** \n");
		fflush(stderr);
		return -1;
	}

                fprintf(stderr,"\tParameters: Box size: %ix%ix%i:%f for protein %s.\n\n",size_x,size_y,size_z,spacing,protein_name);

        /************************************ //
 	*                                     //
	*                                     //
	*   PART B. PROTEIN PROCESSING        //
	*                                     //
	*                                     //
	*                                     //
	**************************************/

	
       /* Reading protein and typing */
       fprintf(stderr,"Reading protein ...");
       fflush(stderr);
       PDB_reader(&mymol, protein_name,0,pqr_type);
       fprintf(stderr,"done\n");


        /************************************ //
 	*                                     //
	*                                     //
	*   PART C. CALCULATE GRID SPACE      //
	*                                     //
	*                                     //
	*                                     //
	**************************************/


	fprintf(stderr,"Calculating grid space...");
	fflush(stderr);


        min_grid[0] = center_x - (size_x / 2);
        min_grid[1] = center_y - (size_y / 2);
        min_grid[2] = center_z - (size_z / 2);

	
        max_grid[0] = center_x + (size_x / 2);
        max_grid[1] = center_y + (size_y / 2);
        max_grid[2] = center_z + (size_z / 2);


	x_points = size_x / spacing;
        y_points = size_y / spacing;
        z_points = size_z / spacing;

	fprintf(stderr,"done\n");
	fprintf(stderr,"Box center: C(%8.3f,%8.3f,%8.3f)\n",center_x,center_y,center_z);
        fprintf(stderr,"From: v(%8.3f,%8.3f,%8.3f) to V(%8.3f,%8.3f,%8.3f)\n",min_grid[0],min_grid[1],min_grid[2],max_grid[0],max_grid[1],max_grid[2]);
	fprintf(stderr,"Grid points: %i\n", x_points * y_points * z_points);
	

	/* Allocating grids stuff */


        grids = (float ****)malloc(7*sizeof(float ***));  /* Allocate 7 probes (only 5 will be calculated) */

	for( i = 0 ; i < 7; ++i)
	{
		grids[i] = (float ***) malloc(x_points*sizeof(float **));
		for( j = 0; j < x_points; ++j)
		{
			grids[i][j] = (float **) malloc(y_points*sizeof(float *));
			for( k = 0; k < y_points; ++k)
			{
				grids[i][j][k] = (float *) malloc(z_points*sizeof(float));
				for( l = 0; l < z_points; ++l)
				{
					grids[i][j][k][l] = 0.0f;
				}
			}

		}
	}



        /************************************ //
 	*                                     //
	*                                     //
	*   PART D. 3D GRID CALCULATION       //
	*                                     //
	*                                     //
	*                                     //
	**************************************/

	fprintf(stderr,"Calculating energy at the grid ...");
	fflush(stderr);

        for( i = 0 ; i < x_points; ++i)
	{
		fprintf(stderr,".");
		fflush(stderr);
                for( j = 0; j < y_points; ++j)
		{
                        for( k = 0; k < z_points; ++k)
			{
				point[0] = min_grid[0] + (i * spacing);
                                point[1] = min_grid[1] + (j * spacing);
                                point[2] = min_grid[2] + (k * spacing);
//                                printf("Point %f,%f,%f\n",point[0],point[1],point[2]);

		                for( l=0; l < mymol->n_atoms; ++l)
				{
		                        dx = mymol->x[l] - point[0];
		                        dy = mymol->y[l] - point[1];
		                        dz = mymol->z[l] - point[2];
		                        dist = dx*dx +dy*dy+dz*dz;
		                        if( dist <= 100 && dist > 0)  /* 10 A cut-off */
					{

			                      dist2 = dist;
			                      dist = sqrt(dist);
			                      distr = 1 / dist;


						/* 
						   0 - NH4+
						   1 - =O
                  				   2 - OH 
                                                   3 - CH3
                                                */
					      if( mymol->gaff_types[l] > 0 )
					      {
					      grids[0][i][j][k] = grids[0][i][j][k] + (332.0 * mymol->pcharges[l]  *  0.33 * 0.5f * distr * distr );
                                              grids[1][i][j][k] = grids[1][i][j][k] + (332.0 * mymol->pcharges[l]  * -0.37 * 0.5f * distr * distr);
                                              grids[2][i][j][k] = grids[2][i][j][k] + (332.0 * mymol->pcharges[l]  * -0.12 * 0.5f * distr * distr);

					      tmpE = 0.0f;
			                      sigma = prevdw[ mymol->gaff_types[l]-1][ 39 ][0]; /* N4 */
			                      epsilon = prevdw[ mymol->gaff_types[l]-1][ 39 ][1];
			                      sigma6 = prevdw[ mymol->gaff_types[l]-1][ 39 ][2];
			                      r6 = sigma6 / (dist2*dist2*dist2);
			                      r12 = r6 * r6;
			                      tmpE = epsilon * (r12 - 2.0*r6);
					      grids[0][i][j][k] = grids[0][i][j][k] + tmpE;


                                              sigma = prevdw[ mymol->gaff_types[l]-1][ 13 ][0]; /* O */
                                              epsilon = prevdw[ mymol->gaff_types[l]-1][ 13 ][1];
                                              sigma6 = prevdw[ mymol->gaff_types[l]-1][ 13 ][2];
                                              r6 = sigma6 / (dist2*dist2*dist2);
                                              r12 = r6 * r6;
                                              tmpE = epsilon * (r12 - 2.0*r6);
                                              grids[1][i][j][k] = grids[1][i][j][k] + tmpE;


                                              sigma = prevdw[ mymol->gaff_types[l]-1][ 14 ][0]; /* OH */
                                              epsilon = prevdw[ mymol->gaff_types[l]-1][ 14 ][1];
                                              sigma6 = prevdw[ mymol->gaff_types[l]-1][ 14 ][2];
                                              r6 = sigma6 / (dist2*dist2*dist2);
                                              r12 = r6 * r6;
                                              tmpE = epsilon * (r12 - 2.0*r6);
                                              grids[2][i][j][k] = grids[2][i][j][k] + tmpE;

                                              sigma = prevdw[ mymol->gaff_types[l]-1][ 20 ][0]; /* C3 */
                                              epsilon = prevdw[ mymol->gaff_types[l]-1][ 20 ][1];
                                              sigma6 = prevdw[ mymol->gaff_types[l]-1][ 20 ][2];
                                              r6 = sigma6 / (dist2*dist2*dist2);
                                              r12 = r6 * r6;
                                              tmpE = epsilon * (r12 - 2.0*r6);
                                              grids[3][i][j][k] = grids[3][i][j][k] + tmpE;
					      grids[5][i][j][k] = grids[5][i][j][k] + tmpE;


					      /* If acceptor, assume 2pi radians */

                                              if( ( dist <= 5.0 && dist >= 1.0 ) && ((mymol->atoms[l] == 3 && has_hydrogens(mymol,l) == 0 &&  mymol->gaff_types[l] != N4) || mymol->gaff_types[l] == OH || mymol->gaff_types[l] == O ))
					      {

                                                 for(neig_j=0; neig_j < mymol->n_bonds; ++neig_j)
                                                 {
                                                   if( mymol->bond_a1[neig_j] == (l+1))
                                                   {
                                                        neig_index = mymol->bond_a2[neig_j]-1;
                                              if ( mymol->atoms[neig_index] != 4 )
                                              {
                                                       heavy_nb[0] = mymol->x[neig_index];
                                                       heavy_nb[1] = mymol->y[neig_index];
                                                       heavy_nb[2] = mymol->z[neig_index];

                                              }
                                                   }else if( mymol->bond_a2[neig_j] == (l+1)){
                                                       neig_index = mymol->bond_a1[neig_j]-1;
                                              if ( mymol->atoms[neig_index] != 4 )
                                              {
                                                       heavy_nb[0] = mymol->x[neig_index];
                                                       heavy_nb[1] = mymol->y[neig_index];
                                                       heavy_nb[2] = mymol->z[neig_index];


                                              }

                                                   }
                                                 }

                                                 for(neig_j=0; neig_j < mymol->n_bonds; ++neig_j)
                                                 {
                                                   if( mymol->bond_a2[neig_j] == (neig_index+1) && mymol->bond_a1[neig_j] != (l+1) )
						   {
							heavy_hb2[0] = mymol->x[mymol->bond_a1[neig_j]-1];
                                                        heavy_hb2[1] = mymol->y[mymol->bond_a1[neig_j]-1];
                                                        heavy_hb2[2] = mymol->z[mymol->bond_a1[neig_j]-1];
						   }else if( mymol->bond_a1[neig_j] == (neig_index+1) &&  mymol->bond_a2[neig_j] != (l+1)){
                                                        heavy_hb2[0] = mymol->x[mymol->bond_a2[neig_j]-1];
                                                        heavy_hb2[1] = mymol->y[mymol->bond_a2[neig_j]-1];
                                                        heavy_hb2[2] = mymol->z[mymol->bond_a2[neig_j]-1];
                                                   }
						 }

							/* Plane vector where the lone pairs are*/
							tmpvec2[0] = mymol->x[l] - heavy_nb[0];
                                                        tmpvec2[1] = mymol->y[l] - heavy_nb[1];
                                                        tmpvec2[2] = mymol->z[l] - heavy_nb[2];

							tmpvec3[0] = heavy_hb2[0] - heavy_nb[0];
                                                        tmpvec3[1] = heavy_hb2[1] - heavy_nb[1];
                                                        tmpvec3[2] = heavy_hb2[2] - heavy_nb[2];
							
                                                        tmpvec1[0] = (tmpvec2[1]*tmpvec3[2]) - (tmpvec3[1]*tmpvec2[2]);
                                                        tmpvec1[1] = (tmpvec2[2]*tmpvec3[0]) - (tmpvec3[2]*tmpvec2[0]);
							tmpvec1[2] = (tmpvec2[0]*tmpvec3[1]) - (tmpvec3[0]*tmpvec2[1]);

							tmpd1 = sqrt( (tmpvec1[0]*tmpvec1[0]) + (tmpvec1[1]*tmpvec1[1]) + (tmpvec1[2]*tmpvec1[2]));
							univec[0] = tmpvec1[0]/tmpd1;
                                                        univec[1] = tmpvec1[1]/tmpd1;
                                                        univec[2] = tmpvec1[2]/tmpd1;
							
                                                        /* H...Acceptor vector */
                                                        tmpvec1[0] = point[0] - mymol->x[l];
                                                        tmpvec1[1] = point[1] - mymol->y[l];
                                                        tmpvec1[2] = point[2] - mymol->z[l];

							dtoplane = (tmpvec1[0]*univec[0])+(tmpvec1[1]*univec[1])+(tmpvec1[2]*univec[2]);
							projection[0] = point[0] - (dtoplane*univec[0]);
                                                        projection[1] = point[1] - (dtoplane*univec[1]);
                                                        projection[2] = point[2] - (dtoplane*univec[2]);

							/* projection now is the projection onto the plane of lone pairs */

							/* Angle t0 is acos( projection-l dot product point-l  */
							tmpvec2[0] = projection[0] - mymol->x[l];
                                                        tmpvec2[1] = projection[1] - mymol->y[l];
                                                        tmpvec2[2] = projection[2] - mymol->z[l];
                                                        norm = sqrt( (tmpvec2[0]*tmpvec2[0])+(tmpvec2[1]*tmpvec2[1])+(tmpvec2[2]*tmpvec2[2]));
                                                        tmpvec2[0] /= norm;
                                                        tmpvec2[1] /= norm;
                                                        tmpvec2[2] /= norm;

							/* Normalize *after* projection*/
                                                        norm = sqrt( (tmpvec1[0]*tmpvec1[0])+(tmpvec1[1]*tmpvec1[1])+(tmpvec1[2]*tmpvec1[2]));
                                                        tmpvec1[0] /= norm;
                                                        tmpvec1[1] /= norm;
                                                        tmpvec1[2] /= norm;



							tmpd2 = (tmpvec2[0]*tmpvec1[0])+(tmpvec2[1]*tmpvec1[1])+(tmpvec2[2]*tmpvec1[2]); /* t0 */
							
							/* Angle ti is acos (l-heavy_nb dot product projection-l) */
							/* C=O axis */
							tmpvec1[0] = mymol->x[l] - heavy_nb[0];
                                                        tmpvec1[1] = mymol->y[l] - heavy_nb[1];
                                                        tmpvec1[2] = mymol->z[l] - heavy_nb[2];
							norm = sqrt( (tmpvec1[0]*tmpvec1[0])+(tmpvec1[1]*tmpvec1[1])+(tmpvec1[2]*tmpvec1[2]));
							tmpvec1[0] /= norm;
                                                        tmpvec1[1] /= norm;
                                                        tmpvec1[2] /= norm;


							tmpd1 = (tmpvec1[0]*tmpvec2[0])+(tmpvec1[1]*tmpvec2[1])+(tmpvec1[2]*tmpvec2[2]); /* t1. No need to divide since norm is 1 in both cases */
							tmpd5 = tmpd1;
                                                        tmpd1 = acos(tmpd1);

                                                if(  mymol->atoms[l] == 2)
                                                {
                                                        tmpE = (C_hbond[2]/powf(dist,8.0f))- (D_hbond[2]/powf(dist,6.0f));
                                                        if( tmpE < -0.25 || 1 == 1 )
							{
								if( mymol->gaff_types[l] == O)
								{
									if((fabs(tmpd1) < 1.9198 && fabs(tmpd1) > 1.5707))
									{
									  tmpd4 = tmpd5;
									  tmpd4 *= tmpd4;
									  tmpd4 = k2 - tmpd4;
									  tmpd4 = powf(tmpd4,3.0f);
									  tmpd4 *= k1;
									  tmpd4 *= tmpd2;
									  
		                                                          tmpE *= (tmpd4);
		                                                          grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                          grids[5][i][j][k] = grids[5][i][j][k] - tmpE;
                                                                          grids[0][i][j][k] = grids[0][i][j][k] + tmpE;

	
									}else if((fabs(tmpd1) < 1.5707 && fabs(tmpd1) > 0)){
										tmpE *= ((0.9f+(0.1f*sin(2.0f*tmpd1)))*tmpd2);
		                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;
                                                                          	grids[0][i][j][k] = grids[0][i][j][k] + tmpE;
									}
								}else if( mymol->gaff_types[l] == OH){

			                                                for(neig_j=0; neig_j < mymol->n_bonds; ++neig_j)
			                                                {
                       			                                 if( mymol->bond_a1[neig_j] == (l+1))
                                       			                 {
			                                                        neig_index = mymol->bond_a2[neig_j]-1;
			                                              		if ( mymol->atoms[neig_index] == 4 )
			                                              		{
			                                                       		heavy_nb[0] = mymol->x[neig_index];
			                                                       		heavy_nb[1] = mymol->y[neig_index];
			                                                       		heavy_nb[2] = mymol->z[neig_index];
			
			                                              		}
			                                                 }else if( mymol->bond_a2[neig_j] == (l+1)){
			                                                       neig_index = mymol->bond_a1[neig_j]-1;
			                                              	       if ( mymol->atoms[neig_index] == 4 )
			                                                       {
			                                                       		heavy_nb[0] = mymol->x[neig_index];
			                                                       		heavy_nb[1] = mymol->y[neig_index];
			                                                       		heavy_nb[2] = mymol->z[neig_index];
			                                              	       }
			
			                                                 }
                                                 			}

			                                                vector1[0] = heavy_nb[0] - mymol->x[l];
			                                                vector1[1] = heavy_nb[1] - mymol->y[l];
			                                                vector1[2] = heavy_nb[2] - mymol->z[l];
			                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);
			                                                vector2[0] = heavy_nb[0] - point[0];
			                                                vector2[1] = heavy_nb[1] - point[1];
			                                                vector2[2] = heavy_nb[2] - point[2];
			                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
			                                                alfa = ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );


									if((fabs(tmpd1) < 0.9599))
									{
										tmpE *= tmpd2 * tmpd2;
                                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;
	                                                                        grids[0][i][j][k] = grids[0][i][j][k] + tmpE;
									}else if( (fabs(tmpd1) > 0.9599)){
                                                                                tmpE *= alfa * alfa;
                                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;
                                                                          	grids[0][i][j][k] = grids[0][i][j][k] + tmpE;
									}
								}
							}
                                                }else if( mymol->atoms[l] == 3){
                                                        tmpE = (C_hbond[1]/powf(dist,8.0f))- (D_hbond[1]/powf(dist,6.0f));
                                                        if( tmpE < -0.25)
							{	
								if((fabs(tmpd1) < 1.5707 && fabs(tmpd1) > 0)){
                                                                                tmpE *= ((0.9f+(0.1f*sin(2.0f*tmpd1)))*tmpd2);
                                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;
										grids[0][i][j][k] = grids[0][i][j][k] + tmpE;
								}
							}

                                                       if( tmpE < 0)
                                                        {
                                                                if((fabs(tmpd1) < 1.5707 && fabs(tmpd1) > 0)){
                                                                tmpE *= ((0.9f+(0.1f*sin(2.0f*tmpd1)))*tmpd2);
                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;
								}
                                                        }


                                                }



					      }
					      /* If donor, where is the heavy atom to calc angles */
					      if( (dist <= 5.0 && dist >= 1.0 ) && ((mymol->atoms[l] == 3 && has_hydrogens(mymol,l) > 0) || mymol->gaff_types[l] == OH ))
					      {

						 heavy_nb[0] = mymol->x[l];
                                                 heavy_nb[1] = mymol->y[l];
                                                 heavy_nb[2] = mymol->z[l];

				                 for(neig_j=0; neig_j < mymol->n_bonds; ++neig_j)
				                 {
					   	   if( mymol->bond_a1[neig_j] == (l+1))
				                   {
							neig_index = mymol->bond_a2[neig_j]-1;
			                      if ( mymol->atoms[neig_index] == 4 )
					      {
						       heavy_nb[0] = mymol->x[neig_index];
                                                       heavy_nb[1] = mymol->y[neig_index];
                                                       heavy_nb[2] = mymol->z[neig_index];


                                                /* [N,O]-H */
                                                vector1[0] = (heavy_nb[0] - mymol->x[l]);
                                                vector1[1] = (heavy_nb[1] - mymol->y[l]);
                                                vector1[2] = (heavy_nb[2] - mymol->z[l]);

                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                /* H...[N,O] */
                                                vector2[0] = (point[0] - heavy_nb[0]);
                                                vector2[1] = (point[1] - heavy_nb[1]);
                                                vector2[2] = (point[2] - heavy_nb[2]);

                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);


                                                alfa = ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                /* OH acceptor */
                                                if( mymol->atoms[l] == 3){
                                                        tmpE = (C_hbond[1]/powf(dist,8.0f))- (D_hbond[1]/powf(dist,6.0f));
                                                        tmpE *=  alfa*alfa;
                                                        if( tmpE < -0.25 && alfa >= 0 ){
                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                grids[1][i][j][k] = grids[1][i][j][k] + tmpE;
                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;

                                                        }
                                                }else if ( mymol->gaff_types[l] == OH){
                                                        tmpE = (C_hbond[2]/powf(dist,8.0f))- (D_hbond[2]/powf(dist,6.0f));
                                                        tmpE *=  (alfa*alfa*alfa*alfa);
                                                        if( tmpE < -0.25 && alfa >= 0){

                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                grids[1][i][j][k] = grids[1][i][j][k] + tmpE;
                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;

                                                        }
                                                }


					      }
				                   }else if( mymol->bond_a2[neig_j] == (l+1)){
                                                       neig_index = mymol->bond_a1[neig_j]-1;
                                              if ( mymol->atoms[neig_index] == 4 )
                                              {
                                                       heavy_nb[0] = mymol->x[neig_index];
                                                       heavy_nb[1] = mymol->y[neig_index];
                                                       heavy_nb[2] = mymol->z[neig_index];


                                                /* [N,O]-H */
                                                vector1[0] = (heavy_nb[0] - mymol->x[l]);
                                                vector1[1] = (heavy_nb[1] - mymol->y[l]);
                                                vector1[2] = (heavy_nb[2] - mymol->z[l]);

                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                /* H...[N,O] */
                                                vector2[0] = (point[0] - heavy_nb[0]);
                                                vector2[1] = (point[1] - heavy_nb[1]);
                                                vector2[2] = (point[2] - heavy_nb[2]);

                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);


                                                alfa = ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                /* OH acceptor */
                                                if( mymol->atoms[l] == 3){
                                                        tmpE = (C_hbond[1]/powf(dist,8.0f))- (D_hbond[1]/powf(dist,6.0f));
                                                        tmpE *=  alfa*alfa;
                                                        if( tmpE < -0.25 && alfa >= 0 ){
                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                grids[1][i][j][k] = grids[1][i][j][k] + tmpE;
                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;

                                                        }
                                                }else if ( mymol->gaff_types[l] == OH){
                                                        tmpE = (C_hbond[2]/powf(dist,8.0f))- (D_hbond[2]/powf(dist,6.0f));
                                                        tmpE *=  (alfa*alfa*alfa*alfa);
                                                        if( tmpE < -0.25 && alfa >= 0){

                                                                grids[2][i][j][k] = grids[2][i][j][k] + tmpE;
                                                                grids[1][i][j][k] = grids[1][i][j][k] + tmpE;
                                                                grids[5][i][j][k] = grids[5][i][j][k] - tmpE;

                                                        }
                                                }

                                              }

				                   }
				                 }


/***/

					      }


					 } /* GAFF type */

                                             if( mymol->atoms[l] >= 100)
                                              {
/*                                                  if( get_metal_vdw_energy(mymol->atoms[l],40,dist) > 0)
                                                  fprintf(stderr,"%f\n",get_metal_vdw_energy(mymol->atoms[l],40,dist));*/
                                                  grids[0][i][j][k] = grids[0][i][j][k] + get_metal_vdw_energy(mymol->atoms[l],40,dist);
                                                  grids[1][i][j][k] = grids[1][i][j][k] + get_metal_vdw_energy(mymol->atoms[l],14,dist);
                                                  grids[2][i][j][k] = grids[2][i][j][k] + get_metal_vdw_energy(mymol->atoms[l],15,dist);
                                                  grids[3][i][j][k] = grids[3][i][j][k] + get_metal_vdw_energy(mymol->atoms[l],21,dist);
                                                  grids[5][i][j][k] = grids[5][i][j][k] + get_metal_vdw_energy(mymol->atoms[l],21,dist);
                                              }


					} /* Distance check */

				} /* Atoms */

			} /* Z */
	
		} /* Y */
	} /* X */


	fprintf(stderr,"done\n");
	fflush(stderr);
       cleanup(&mymol);


        /************************************ //
 	*                                     //
	*                                     //
	*   PART E. DUMP THE GRID TO THE WORLD//
	*                                     //
	*                                     //
	*                                     //
	**************************************/


	fprintf(stderr,"Dumping grid to files ...");
	fflush(stderr);


for(l = 0; l < 7; ++l)
{

	if( (output[l] = fopen(output_names[l],"w")) == NULL)
	{
                fprintf(stderr,"ERROR-MAIN%i: Cannot open grid file to write.\n\n",__LINE__);
                fprintf(stderr," *** END ERROR *** \n");
                fflush(stderr);
                return -3;
	}

	if( (all_points = (float *) calloc(x_points*y_points*z_points,sizeof(float))) == NULL)
	{
                fprintf(stderr,"ERROR-MAIN%i: Cannot allocate memory.\n\n",__LINE__);
                fprintf(stderr," *** END ERROR *** \n");
                fflush(stderr);
                return -2;
	}


	counter_points = 0;
        for( i = 0 ; i < x_points; ++i)
        {
                for( j = 0; j < y_points; ++j)
                {
                        for( k = 0; k < z_points; ++k)
                        {
				if( grids[l][i][j][k] > 200)
				 grids[l][i][j][k] = 200;

				all_points[counter_points] = grids[l][i][j][k];
				++counter_points;
			}

		}
	}

	fprintf(stderr,"Grid points: %i.\n",counter_points);

	fprintf(output[l],"object 1 class gridpositions counts %i %i %i\n", x_points, y_points, z_points);

	fprintf(output[l],"origin %8.3f %8.3f %8.3f\n",min_grid[0],min_grid[1], min_grid[2]);
	fprintf(output[l],"delta %8.3f 0.0 0.0\n",spacing);
	fprintf(output[l],"delta 0.0 %8.3f 0.0\n", spacing);
        fprintf(output[l],"delta 0.0 0.0 %8.3f\n", spacing);
	fprintf(output[l],"object 2 class gridconnections count %i %i %i\n",x_points, y_points, z_points);
	fprintf(output[l],"object 3 class array type double rank 0 items %i data follows\n",x_points*y_points*z_points);
	for( i = 0; i < counter_points; ++i)
	{
		 fprintf(output[l],"%18.3f",all_points[i]);

		 if( (i+1)%3  == 0)
		 {
			fprintf(output[l],"\n");
		 }
	}

	if( counter_points%3 != 0)
              fprintf(output[l],"\n");


	fprintf(stderr,"Dumped %i points.\n",i);

	fprintf(output[l],"attribute \"dep\" string \"positions\"\n");
        fprintf(output[l],"object \"regular positions regular connections\" class field\n");
	fprintf(output[l],"component \"positions\" value 1\n");
	fprintf(output[l],"component \"connections\" value 2\n");
	fprintf(output[l],"component \"data\" value 3\n");

	fclose(output[l]);
	free(all_points);
}
	fprintf(stderr,"done\n");

	/************************************ //
	*
	*
	*   PART F. Hotspots identification
	*
	*
	*
	*
	**************************************/
	
	sel_points = (float**)calloc(sizeof(float **), x_points * y_points * z_points);
	for ( i = 0; i < x_points * y_points * z_points; ++i)
	    sel_points[i] = (float*)calloc(sizeof(float*), 8);


        /* Thresholding and filtering */

        current_point = 0;
        best_grid = 0;
        for ( i = 0; i < x_points; ++i) {
                for ( j = 0; j < y_points; ++j) {
                        for ( k = 0; k < z_points; ++k) {

                                for( l = 0 ; l < 6; l++)
                                {
                                	if ( grids[l][i][j][k] > thresholds[l] ) {
					 grids[l][i][j][k] = 9999.9f;
	                                }
				}

				best_energy = grids[0][i][j][k];
				best_grid = 0;
				for( l = 0 ; l < 6; l++)
				{
					if( grids[l][i][j][k] < best_energy)
					{
					   best_grid = l;
					   best_energy = grids[l][i][j][k];
					}			
				}

                                if ( best_energy < 0.0f ) {
                                        sel_points[current_point][0] = min_grid[0] + (i * spacing);
                                        sel_points[current_point][1] = min_grid[1] + (j * spacing);
                                        sel_points[current_point][2] = min_grid[2] + (k * spacing);
                                        sel_points[current_point][4] = grids[best_grid][i][j][k];
                                        sel_points[current_point][5] = 1.0f;    /* Active */
                                        sel_points[current_point][6] = best_grid;
                                        sel_points[current_point][7] = 0.0f;    /* For sticky spots */
                                        ++current_point;
                                }

                        }

                }
        }

        total_points = current_point;
        current_point = 0;


        current_point = 0;
        /* Remove redundant points */
        for ( i = 0; i < total_points; i++) {
                if ( sel_points[i][5] > 0) {
                                dist2 = 4.0f;

                        for (j = 0; j < total_points; j++) {
                                if ( i != j && sel_points[j][5] > 0) {
                                        dx = sel_points[i][0] - sel_points[j][0];
                                        dy = sel_points[i][1] - sel_points[j][1];
                                        dz = sel_points[i][2] - sel_points[j][2];
                                        dist = dx * dx + dy * dy + dz * dz;
                                        if ( sel_points[j][4] <= sel_points[i][4] && dist <= dist2  ){ /*&& sel_points[j][6] == sel_points[i][6] ) { && sel_points[j][3] == sel_points[i][3])  && sel_points[j][6] == sel_points[i][6])*/
                                                sel_points[i][5] = 0.0f;
                                                break;
                                        }

                                }
                        }
                }
        }





	hotspots = fopen("hotspots.pdb","wb");
	ch = 1;
	for( current_point = 0; current_point < total_points; current_point++)
	{
		if( sel_points[current_point][5] == 1.0f)
		{
		fprintf(hotspots,"ATOM %6i  %2s  %3s  %4i    %8.3f%8.3f%8.3f  1.2 %7.4f\n",current_point+1,hots[  (int) sel_points[current_point][6] ], hots_res[  (int) sel_points[current_point][6] ],ch,sel_points[current_point][0],sel_points[current_point][1],sel_points[current_point][2],sel_points[current_point][4]);
			++ch;
		}
	}

	fclose(hotspots);



        /************************************ //
 	*                                     //	
	*                                     //
	*   PART G. CLEAN UP                  //
	*                                     //
	*                                     //
	*                                     //
	**************************************/


        for( i = 0 ; i < 7; ++i)
        {
                for( j = 0; j < x_points; ++j)
                {
                        for( k = 0; k < y_points; ++k)
                        {
                                free(grids[i][j][k]);
                        }
			free(grids[i][j]);
                }
		free(grids[i]);
        }



	fprintf(stderr," *** END NORMAL *** \n");
	fflush(stderr);


	return 0;


}

void print_header()
{

	fprintf(stderr,"cGRILL - A simple affinity map generator written in C.\n");
	fprintf(stderr,"Alvaro Cortes Cabrera <alvarocortesc@gmail.com>\n");
	fflush(stderr);

}

void print_usage()
{
	fprintf(stderr,"Usage: cGRILL.exe <params>\n\n");
	fprintf(stderr,"Valid params:\n");
	fprintf(stderr,"\t-p or --protein <file>. Set protein file.\n");
	fprintf(stderr,"\t-x or --center_x <value>. Set center of the box.\n");
        fprintf(stderr,"\t-y or --center_y <value>. Set center of the box.\n");
        fprintf(stderr,"\t-z or --center_z <value>. Set center of the box.\n");
        fprintf(stderr,"\t-h or --size_x <value>. Set size of the box in Angstrom.\n");
        fprintf(stderr,"\t-w or --size_y <value>. Set size of the box in Angstrom.\n");
        fprintf(stderr,"\t-t or --size_z <value>. Set size of the box in Angstrom.\n");
        fprintf(stderr,"\t-s or --spacing <value>. Set spacing of the box in Angstrom.\n");
        fprintf(stderr,"\t-q or --pqr. Input file is PQR (e.g. from AMBER) instead of ATM-like.\n\n");
	fflush(stderr);
}
