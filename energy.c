/***********************************************************************
 *
 *
 *	(c) 2011-2024. Alvaro Cortes Cabrera. <alvarocortesc@gmail.com>
 *
 *       This code is licensed under GPL v3 terms
 *
 *  Internal energy module
 *  
 *
 *
 */

/*#include <stdio.h>
#include <math.h>
#include <libmol2.h>
#include <gaff.h>

#define G_PI 3.14159265358979323846*/


int get_energy(MOL2 **mymol, double *tenergy, int **tortype1)
{
        MOL2 *mols;
        float dx,dy,dz;
        float df;
        float ddf1;
        float df1;
        float e1;
        float dgx,dgy,dgz;
        float dtxi,dtxj,dtyi,dtyj,dtzi,dtzj;
        float dfx,dfy,dfz;
        float gaa, ra2, rb2, gbb;
        float fg,hg,fga,hgb;
        float r2,r;
        float RJR,RIR,cst;
        int *tor_type;

        int i = 0;
        int j = 0;
        int k = 0;
        int l = 0;

        float idi,pn,pk,phi0;
        int cur_angle = 0;

        int t1,t2,t3,t4;

        float A,B;
        float sigma;
        float sigma6;
        float epsilon;
        float tmpE;
        float vdwE;
        float vdwE14;
        float  r6,r12;

        float bondE = 0;
        float angleE = 0;
        float torE = 0;
        double totalE = 0;
        double vec1[3],vec2[3];
        double uvec1[3],uvec2[3];
        double *lvec1,*lvec2;
        double tor1;
        double tor2;
        double tor3;
        float thetha;
        float dist = 0.0f;
        float dist2 = 0.0f;
        float distr = 0.0f;
        int ttor = 0;


        /* FLAGS */
        int bondflag;
        int angleflag;
        int torflag;
        int gradflag;
        int vdwflag;
	int verboseflag;

        mols = *mymol;
        tor_type = *tortype1;

        bondflag = 1;
        torflag = 1;
        vdwflag = 1;
        gradflag = 1;
        angleflag = 1;
	verboseflag = 0;


/*        printf("ENERGY ROUTINE\n");
        printf("--------------\n");
        printf("Bonds: %i.\n",mols->n_bonds);
        printf("Angles: %i.\n",mols->n_angles);
        printf("Dihedrals: %i.\n",mols->n_torsionals);*/



/******************************************************************************
 *
 *                                 BONDS STUFF
 *
 *****************************************************************************/

       if( bondflag == 1)
       {
        for(j=0; j < mols->n_bonds; ++j)
        {
             
            if( mols->gaff_types[ mols->bond_a1[j] - 1 ] < mols->gaff_types[ mols->bond_a2[j] - 1 ])
            {
                 t1 = mols->gaff_types[ mols->bond_a1[j] - 1 ]-1;
                 t2 = mols->gaff_types[ mols->bond_a2[j] - 1 ]-1;
            }else{
                 t1 = mols->gaff_types[ mols->bond_a2[j] - 1 ]-1;
                 t2 = mols->gaff_types[ mols->bond_a1[j] - 1 ]-1;
            }

#ifdef DEBUG
            printf("%i-%i | %i-%i %f %f.\n",mols->bond_a1[j] - 1,mols->bond_a2[j] - 1, t1,t2,mols->bond_dist[j], bonds[t1][t2][1]);
            printf("%f * %f = %f.\n",bonds[t1][t2][0], pow( mols->bond_dist[j] - bonds[t1][t2][1],2), bonds[t1][t2][0] * pow( mols->bond_dist[j] - bonds[t1][t2][1],2));
#endif

                 dx = mols->x[ mols->bond_a1[j] - 1 ] - mols->x[ mols->bond_a2[j] - 1 ];
                 dy = mols->y[ mols->bond_a1[j] - 1 ] - mols->y[ mols->bond_a2[j] - 1 ];
                 dz = mols->z[ mols->bond_a1[j] - 1 ] - mols->z[ mols->bond_a2[j] - 1 ];
                 mols->bond_dist[j] = sqrt(dx*dx+dy*dy+dz*dz);
                 df = bonds[t1][t2][0] * ( mols->bond_dist[j] - bonds[t1][t2][1]);

            bondE = bondE +  ( df * ( mols->bond_dist[j] - bonds[t1][t2][1]) );


            if( gradflag == 1)
            {
                 df *= 2 / mols->bond_dist[j];
                 dx *= df;
                 dy *= df;
                 dz *= df;
                 mols->grads_X[ mols->bond_a1[j] - 1] += dx; 
                 mols->grads_X[ mols->bond_a2[j] - 1] -= dx;

                 mols->grads_Y[ mols->bond_a1[j] - 1] += dy;
                 mols->grads_Y[ mols->bond_a2[j] - 1] -= dy;

                 mols->grads_Z[ mols->bond_a1[j] - 1] += dz;
                 mols->grads_Z[ mols->bond_a2[j] - 1] -= dz;

            }

        }

/*        printf("BOND STRETCHING ENERGY: %f kcal/mol.\n",bondE);*/
       }

/******************************************************************************
 *
 *                               ANGLES STUFF
 *
 ******************************************************************************/

       if( angleflag == 1)
       {
        for(j=0; j < mols->n_angles; ++j)
        {


                 t1 = mols->gaff_types[ mols->ia[j] ]-1;
                 t2 = mols->gaff_types[ mols->ja[j] ]-1;
                 t3 = mols->gaff_types[ mols->ka[j] ]-1;

               vec1[0] = mols->x[ mols->ia[j] ] - mols->x[ mols->ja[j] ];
               vec1[1] = mols->y[ mols->ia[j] ] - mols->y[ mols->ja[j] ];
               vec1[2] = mols->z[ mols->ia[j] ] - mols->z[ mols->ja[j] ];

               vec2[0] = mols->x[ mols->ka[j] ] - mols->x[ mols->ja[j] ];
               vec2[1] = mols->y[ mols->ka[j] ] - mols->y[ mols->ja[j] ];
               vec2[2] = mols->z[ mols->ka[j] ] - mols->z[ mols->ja[j] ];

               RIR = sqrt((vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]));
               RJR = sqrt((vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]));
               cst = (vec1[0]*vec2[0] + vec1[1] * vec2[1] + vec1[2]*vec2[2])  / ( RJR*RIR) ;
               thetha = acos( cst ) ;
/*                 printf("%s-%s-%s. %f.\n",let[mols->gaff_types[mols->ia[j]]-1],let[mols->gaff_types[mols->ja[j]]-1],let[mols->gaff_types[mols->ka[j]]-1],thetha);*/
/*               printf("%i-%i-%i. %f.\n",t1,t2,t3,thetha);*/
/*               printf("Diff theta: %f - %f = %f. K=%f.\n",thetha, (angl[t1][t2][t3][1]*G_PI/180), thetha - (angl[t1][t2][t3][1]*G_PI/180), angl[t1][t2][t3][0]);*/
               dx = ( thetha - (angl[t1][t2][t3][1]*G_PI/180));
               df =  angl[t1][t2][t3][0] * dx;
               angleE = angleE +  ( df * dx);
            if( gradflag == 1)
            {
                df *= 2;

/*                printf("DF: %f.\n",df);*/
                dtxi = (cst * RJR * vec1[0] + RIR* (0-vec2[0])) / ( sin(thetha) * RIR * RIR * RJR);
                dtyi = (cst * RJR * vec1[1] + RIR* (0-vec2[1])) / ( sin(thetha) * RIR * RIR * RJR);
                dtzi = (cst * RJR * vec1[2] + RIR* (0-vec2[2])) / ( sin(thetha) * RIR * RIR * RJR);
                mols->grads_X[ mols->ia[j] ] += ( df  *dtxi);
                mols->grads_Y[ mols->ia[j] ] += ( df  *dtyi);
                mols->grads_Z[ mols->ia[j] ] += ( df  *dtzi);

                dtxi = ((cst * RIR * vec2[0]) + RJR*(0-vec1[0])) / (sin(thetha) * RIR * RJR * RJR);
                dtyi = ((cst * RIR * vec2[1]) + RJR*(0-vec1[1]))  / (sin(thetha) * RIR * RJR * RJR);
                dtzi = ((cst * RIR * vec2[2]) + RJR*(0-vec1[2]))  / (sin(thetha) * RIR * RJR * RJR);

                mols->grads_X[ mols->ka[j] ] += ( df  *dtxi);
                mols->grads_Y[ mols->ka[j] ] += ( df  *dtyi);
                mols->grads_Z[ mols->ka[j] ] += ( df  *dtzi);

                dtxi = ((cst * ( RIR*RIR*(0-vec2[0]) + RJR*RJR*(0-vec1[0]))) / (sin(thetha) * RIR * RIR * RJR * RJR)) + (( vec1[0]+vec2[0] ) / (sin(thetha)*RJR*RIR));
                dtyi = ((cst * ( RIR*RIR*(0-vec2[1]) + RJR*RJR*(0-vec1[1]))) / (sin(thetha) * RIR * RIR * RJR * RJR)) + (( vec1[1]+vec2[1] ) / (sin(thetha)*RJR*RIR));
                dtzi = ((cst * ( RIR*RIR*(0-vec2[2]) + RJR*RJR*(0-vec1[2]))) / (sin(thetha) * RIR * RIR * RJR * RJR)) + (( vec1[2]+vec2[2] ) / (sin(thetha)*RJR*RIR));

                mols->grads_X[ mols->ja[j] ] += ( df  *dtxi);
                mols->grads_Y[ mols->ja[j] ] += ( df  *dtyi);
                mols->grads_Z[ mols->ja[j] ] += ( df  *dtzi);


                
/*                dtxi = (1/RIR)*( vec2[0] - cst*vec1[0]);
                dtxj = (1/RJR)*( vec1[0] - cst*vec2[0]);

                dtyi = (1/RIR)*( vec2[1] - cst*vec1[1]);
                dtyj = (1/RJR)*( vec1[1] - cst*vec2[1]);

                dtzi = (1/RIR)*( vec2[2] - cst*vec1[2]);
                dtzj = (1/RJR)*( vec1[2] - cst*vec2[2]);
                dx = df * dtxi;
                dgx = df * dtxj;
                dy = df * dtyi;
                dgy = df * dtyj;
                dz = df * dtzi;
                dgy = df * dtzj;

                 mols->grads_X[ ia[j] ] += dx;
                 mols->grads_Y[ ia[j] ] += dy;
                 mols->grads_Z[ ia[j] ] += dz;

                 mols->grads_X[ ja[j] ] -= (dgx + dx);
                 mols->grads_Y[ ja[j] ] -= (dgy + dy);
                 mols->grads_Z[ ja[j] ] -= (dgz + dz);

                 mols->grads_X[ ka[j] ] += dgx;
                 mols->grads_Y[ ka[j] ] += dgy;
                 mols->grads_Z[ ka[j] ] += dgz;*/



            }


        }
/*        free(ia);
        free(ja);
        free(ka);*/

/*        printf("ANGLE BENDING ENERGY: %f kcal/mol.\n",angleE);*/
       }



/******************************************************************************
 * 
 *                            TORSIONAL STUFF
 *
 ******************************************************************************/

        torE = 0.0;

       if( torflag == 1)
       {
/*        ik = calloc(sizeof(int),1000);
        jk = calloc(sizeof(int),1000);
        kk = calloc(sizeof(int),1000);
        lk = calloc(sizeof(int),1000);
        ttor = get_number_of_rings(mols->&ik, &jk, &kk, &lk);
        printf("Dihedrals %i.\n",ttor);*/

        for( i = 0; i < mols->n_torsionals; ++i)
            tor_type[i] = 0;  /* Not assigned */


        for( i = 0; i < mols->n_torsionals; ++i)
        {

        vec1[0] = (mols->x[mols->ik[i]] - mols->x[mols->jk[i]]);
        vec1[1] = (mols->y[mols->ik[i]] - mols->y[mols->jk[i]]);
        vec1[2] = (mols->z[mols->ik[i]] - mols->z[mols->jk[i]]);

        vec2[0] = (mols->x[mols->jk[i]] - mols->x[mols->kk[i]]);
        vec2[1] = (mols->y[mols->jk[i]] - mols->y[mols->kk[i]]);
        vec2[2] = (mols->z[mols->jk[i]] - mols->z[mols->kk[i]]);

        uvec1[0] = vec1[1]*vec2[2] - vec2[1]*vec1[2];
        uvec1[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
        uvec1[2] = vec1[0]*vec2[1] - vec2[0]*vec1[1];

        if( gradflag == 1)
        {
        ra2 = sqrt(pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2));
        gaa = 0 - (1 / (pow(uvec1[0],2) + pow(uvec1[1],2) + pow(uvec1[2],2))) * ra2;
        ra2 = 1 / ra2;
        fg = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2] * vec2[2];
        }

        vec1[0] = (mols->x[mols->jk[i]] - mols->x[mols->kk[i]]);
        vec1[1] = (mols->y[mols->jk[i]] - mols->y[mols->kk[i]]);
        vec1[2] = (mols->z[mols->jk[i]] - mols->z[mols->kk[i]]);

        vec2[0] = (mols->x[mols->kk[i]] - mols->x[mols->lk[i]]);
        vec2[1] = (mols->y[mols->kk[i]] - mols->y[mols->lk[i]]);
        vec2[2] = (mols->z[mols->kk[i]] - mols->z[mols->lk[i]]);

        uvec2[0] = vec1[1]*vec2[2] - vec2[1]*vec1[2];
        uvec2[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
        uvec2[2] = vec1[0]*vec2[1] - vec2[0]*vec1[1];

        if( gradflag == 1)
        {
        hg = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2] * vec2[2];

        rb2 = sqrt(pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2));
        gbb = (1 / (pow(uvec2[0],2) + pow(uvec2[1],2) + pow(uvec2[2],2))) * rb2;
        rb2 = 1 / rb2 ;
        }
        tor1 = uvec1[0]*uvec2[0] + uvec1[1]*uvec2[1] + uvec1[2]*uvec2[2];
        tor2 = sqrt(pow(uvec1[0],2) + pow(uvec1[1],2) + pow(uvec1[2],2)) * sqrt( pow(uvec2[0],2) + pow(uvec2[1],2)+pow(uvec2[2],2));
        tor3 = acos(tor1/tor2);
        if (tor1 > 0.0)
          tor3 = 0-tor3;

        /* Check improper */

        /* Check specific quartet */

        for( j=0; j < 90; ++j)
        {
            if( ((torsionals2[j][1] == mols->gaff_types[mols->jk[i]] && torsionals2[j][2] == mols->gaff_types[mols->kk[i]]) || (torsionals2[j][2] == mols->gaff_types[mols->jk[i]] && torsionals2[j][1] == mols->gaff_types[mols->kk[i]])) && ( (torsionals2[j][0] == mols->gaff_types[mols->ik[i]] && torsionals2[j][3] == mols->gaff_types[mols->lk[i]]) || (torsionals2[j][3] == mols->gaff_types[mols->ik[i]] && torsionals2[j][0] == mols->gaff_types[mols->lk[i]])) )
            {
              tor_type[i] = 1;
              idi = torsionals2[j][7];
              pn = fabs(torsionals2[j][4]);
              pk = torsionals2[j][5];
              phi0 = (torsionals2[j][6]*G_PI/180);
         
/*printf("Special Dihedral %i.  %s - %s - %s - %s. %f.\n",i,let[mols->gaff_types[mols->ik[i]]-1],let[mols->gaff_types[mols->jk[i]]-1],let[mols->gaff_types[mols->kk[i]]-1],let[mols->gaff_types[mols->lk[i]]-1],tor3*180/G_PI);*/

torE = torE + ( (torsionals2[j][5] / torsionals2[j][7]) * (1.0 + cos( (fabs(torsionals2[j][4])*tor3) - ( (torsionals2[j][6]*G_PI/180)))));
/*printf("%f.\n", (torsionals2[j][5] / torsionals2[j][7])* (1.0 + cos( (fabs(torsionals2[j][4])*tor3) - ( (torsionals2[j][6]*G_PI/180)))));*/


            }
        }
      if (tor_type[i] == 0)  /* Normal torsional */
      {
/*        printf("Dihedral %i.  %s - %s - %s - %s. %f.\n",i,let[mols->gaff_types[mols->ik[i]]-1],let[mols->gaff_types[mols->jk[i]]-1],let[mols->gaff_types[mols->kk[i]]-1],let[mols->gaff_types[mols->lk[i]]-1],tor3*180/G_PI);*/

        if( mols->gaff_types[mols->jk[i]]  > mols->gaff_types[mols->kk[i]])
        {
         t1 = mols->gaff_types[mols->kk[i]]-1;
         t2 = mols->gaff_types[mols->jk[i]]-1;
        }else{
         t2 = mols->gaff_types[mols->kk[i]]-1;
         t1 = mols->gaff_types[mols->jk[i]]-1;
        }

/*         printf("%i-%i. IDivf: %f. Pn: %f. K: %f. phi0: %f.\n",mols->gaff_types[mols->jk[i]],mols->gaff_types[mols->kk[i]],torsionals[t1][t2][0],torsionals[t1][t2][2],fabs(torsionals[t1][t2][1]),torsionals[t1][t2][3]);*/
          
         if(  torsionals[t1][t2][0] != 0.0)
         {
         torE = torE + ( (torsionals[t1][t2][2]/ torsionals[t1][t2][0]) * (1.0 + cos( ( fabs(torsionals[t1][t2][1])*tor3) - ( (torsionals[t1][t2][3]*G_PI/180)))));
/*         printf("%f.\n", ( torsionals[t1][t2][2] / torsionals[t1][t2][0])* (1.0 + cos( (fabs(torsionals[t1][t2][1])*tor3) - ( (torsionals[t1][t2][3]*G_PI/180)))));*/
         }
         if( torsionals[t1][t2][4] != 0.0)
         {
         torE = torE + ( (torsionals[t1][t2][6]/ torsionals[t1][t2][4]) * (1.0 + cos( ( fabs(torsionals[t1][t2][5])*tor3) - ( (torsionals[t1][t2][7]*G_PI/180)))));
/*         printf("%f.\n", ( torsionals[t1][t2][6] / torsionals[t1][t2][4])* (1.0 + cos( (fabs(torsionals[t1][t2][5])*tor3) - ( (torsionals[t1][t2][7]*G_PI/180)))));*/
         }
         if( torsionals[t1][t2][8] != 0.0)
         {
         torE = torE + ( (torsionals[t1][t2][10]/ torsionals[t1][t2][8]) * (1.0 + cos( ( fabs(torsionals[t1][t2][9])*tor3) - ( (torsionals[t1][t2][11]*G_PI/180)))));
/*         printf("%f.\n", ( torsionals[t1][t2][10] / torsionals[t1][t2][8])* (1.0 + cos( (fabs(torsionals[t1][t2][9])*tor3) - ( (torsionals[t1][t2][11]*G_PI/180)))));*/
         }
         if( torsionals[t1][t2][12] != 0.0)
         {
         torE = torE + ( (torsionals[t1][t2][14]/ torsionals[t1][t2][12]) * (1.0 + cos( ( fabs(torsionals[t1][t2][13])*tor3) - ( (torsionals[t1][t2][15]*G_PI/180)))));
/*         printf("%f.\n", ( torsionals[t1][t2][14] / torsionals[t1][t2][12])* (1.0 + cos( (fabs(torsionals[t1][t2][13])*tor3) - ( (torsionals[t1][t2][15]*G_PI/180)))));*/

         }

      }

         if( gradflag == 1)
         {
              e1 = 1.0;
              df1 = 0.0;
              ddf1 = 0.0;
           if( tor_type[i] == 0)
           {
              for(j = 0; j < torsionals[t1][t2][1]; ++j)
              {
                    df1 = sin(tor3) + df1 * cos(tor3);
                    ddf1 = e1 * cos(tor3) - df1 *  sin(tor3);
                    e1 = ddf1;
              }

              e1 = e1 * cos( (torsionals[t1][t2][3]*G_PI/180)) + df1 * sin( (torsionals[t1][t2][3]*G_PI/180));
              df1 = df1 * cos( (torsionals[t1][t2][3]*G_PI/180)) - ddf1 * sin( (torsionals[t1][t2][3]*G_PI/180));
              df1 = (0 - torsionals[t1][t2][1]) * df1;
              ddf1 = ( 0 - (torsionals[t1][t2][1] * torsionals[t1][t2][1])) * e1;
              e1 = e1 + 1.0;
/*              printf("Energy: %f.\n", e1 * torsionals[t1][t2][2]);
              printf("DF: %f.\n", df1 * torsionals[t1][t2][2]);*/
              dfx = gaa * uvec1[0] * df1 * torsionals[t1][t2][2];
              dfy = gaa * uvec1[1] * df1 * torsionals[t1][t2][2];
              dfz = gaa * uvec1[2] * df1 * torsionals[t1][t2][2];

              mols->grads_X[mols->ik[i]] += dfx;
              mols->grads_Y[mols->ik[i]] += dfy;
              mols->grads_Z[mols->ik[i]] += dfz;

              mols->grads_X[mols->jk[i]] -= dfx;
              mols->grads_Y[mols->jk[i]] -= dfy;
              mols->grads_Z[mols->jk[i]] -= dfz;


              dfx = gbb * uvec2[0] * df1 * torsionals[t1][t2][2];
              dfy = gbb * uvec2[0] * df1 * torsionals[t1][t2][2];
              dfz = gbb * uvec2[0] * df1 * torsionals[t1][t2][2];

              mols->grads_X[mols->lk[i]] += dfx;
              mols->grads_Y[mols->lk[i]] += dfy;
              mols->grads_Z[mols->lk[i]] += dfz;

              mols->grads_X[mols->kk[i]] -= dfx;
              mols->grads_Y[mols->kk[i]] -= dfy;
              mols->grads_Z[mols->kk[i]] -= dfz;


              fga = fg * (1 / (pow(uvec1[0],2) + pow(uvec1[1],2) + pow(uvec1[2],2))) * sqrt(pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2));
              hgb = hg * (1 / (pow(uvec2[0],2) + pow(uvec2[1],2) + pow(uvec2[2],2))) * sqrt(pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2));
              
              dfx = (fga*uvec1[0] - hgb*uvec2[0]) * df1 * torsionals[t1][t2][2];
              dfy = (fga*uvec1[1] - hgb*uvec2[1]) * df1 * torsionals[t1][t2][2];
              dfz = (fga*uvec1[2] - hgb*uvec2[2]) * df1 * torsionals[t1][t2][2];

              mols->grads_X[mols->jk[i]] += dfx;
              mols->grads_Y[mols->jk[i]] += dfy;
              mols->grads_Z[mols->jk[i]] += dfz;
              mols->grads_X[mols->kk[i]] -= dfx;
              mols->grads_Y[mols->kk[i]] -= dfy;
              mols->grads_Z[mols->kk[i]] -= dfz;
              }



        }

/*        printf("TOTAL TORSIONAL ENERGY = %f kcal/mol.\n",torE);
        fflush(stdout);*/

        }
        }



/*****************************************************************************
 *
 *          
 *                             VAN DER WAALS STUFF
 *
 *
 ****************************************************************************/

        vdwE = 0.0;
        vdwE14 = 0.0;
       if( vdwflag == 1)
       { 
        for( j=0; j < mols->n_pairs; ++j)
        {

              dx = (mols->x[ mols->vdwpairs_a[j]] - mols->x[ mols->vdwpairs_b[j]]);
              dy = (mols->y[ mols->vdwpairs_a[j]] - mols->y[ mols->vdwpairs_b[j]]);
              dz = (mols->z[ mols->vdwpairs_a[j]] - mols->z[ mols->vdwpairs_b[j]]);

                      dist = dx*dx + dy*dy + dz*dz;
                      dist2 = dist;
                      dist = sqrt(dist);
                      distr = 1 / dist;

                      tmpE = (332.0 * mols->pcharges[ mols->vdwpairs_a[j] ] * mols->pcharges[ mols->vdwpairs_b[j] ]) * distr;
                      vdwE += tmpE;

/*                      sigma = vdw_r[ mols->gaff_types[ mols->vdwpairs_a[j]] -1] + vdw_r[mols->gaff_types[ mols->vdwpairs_b[j]]-1 ];
                      epsilon = sqrt (  vdw_epsilon[ mols->gaff_types[mols->vdwpairs_a[j]] -1] * vdw_epsilon[ mols->gaff_types[mols->vdwpairs_b[j]] -1]);*/
                      sigma = prevdw[ mols->gaff_types[ mols->vdwpairs_a[j]] -1][mols->gaff_types[ mols->vdwpairs_b[j]]-1 ][0];
                      epsilon = prevdw[ mols->gaff_types[ mols->vdwpairs_a[j]] -1][mols->gaff_types[ mols->vdwpairs_b[j]]-1 ][1];
                      sigma6 = prevdw[ mols->gaff_types[ mols->vdwpairs_a[j]] -1][mols->gaff_types[ mols->vdwpairs_b[j]]-1 ][2];

                      r6 = sigma6 / (dist2*dist2*dist2);
/*                      r6 = r6 * r6 * r6;
                      r6 = r6 * r6;*/
                      r12 = r6 * r6;

                      epsilon /= mols->vdw_type[j];

                      tmpE = epsilon * (r12 - 2.0*r6);

                      vdwE = vdwE + tmpE;

                      if( gradflag == 1 )
                      {
                          r6 = r6 * (distr); /* 7*/
                          r12 = r12 * (distr); /* 13*/
                          df = (12.0*epsilon) * ( r6 - r12);
/*                          printf("DF: %f.\n",df);*/
                          df += 0-((332.0 * mols->pcharges[ mols->vdwpairs_a[j] ] * mols->pcharges[ mols->vdwpairs_b[j] ]) / (dist2));

	                  dx *= df;
	                  dy *= df;
	       	          dz *= df;
	       	          mols->grads_X[ mols->vdwpairs_a[j] ] += dx;
	                  mols->grads_X[ mols->vdwpairs_b[j] ] -= dx;

	                  mols->grads_Y[ mols->vdwpairs_a[j] ] += dy;
	                  mols->grads_Y[ mols->vdwpairs_b[j] ] -= dy;
	
	                  mols->grads_Z[ mols->vdwpairs_a[j] ] += dz;
	                  mols->grads_Z[ mols->vdwpairs_b[j] ] -= dz;
                     

                      }

        }
/*        printf("TOTAL VAN DER WAALS ENERGY = %f kcal/mol.\n",vdwE);*/

        }

        totalE = vdwE+torE+bondE+angleE;

                   if( verboseflag == 1)
        printf("ENERGY ROUTINE> TERMS. VDW: %f. TOR: %f. BOND: %f. ANGLE: %f.\n",vdwE,torE,bondE,angleE);

/*        printf("ENERGY ROUTINE> ENERGY: %f.\n",totalE);*/
        *tenergy = totalE;
        *mymol = mols;
        return 0;
}
