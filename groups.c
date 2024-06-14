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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_RING_MEMBERS 300  /* Puff! but simple linked lists = heck! */
#define MAX_RINGS 30         /* Little hack to prevent infinite recursive problem */

char *elems[10] = {"C ","O ","N ","H ","P ","S ","I ","Br","Cl","F "};


int has_hydrogens( MOL2 *mol, int atom )
{
        int res = 0, j = 0;
        for (j = 0; j < mol->n_bonds; ++j)
        {
                if ( mol->bond_a1[j] == (atom + 1)){
                        if( mol->atoms[mol->bond_a2[j]-1] == 4)
                                res = 1;
                }else if( mol->bond_a2[j] == (atom + 1)){
                        if( mol->atoms[mol->bond_a1[j]-1] == 4)
                                res = 1;
                }
        }
        return res;
}


int get_romaters(int type, int max,ROTAMER **tmp_rota)
{
          int i = 0;;
          ROTAMER *rot = NULL;
          FILE *lib = NULL;
          char *lib_name;
          char *line = NULL;
          char myx[12],myy[12],myz[12];
          rot = *tmp_rota;
          rot->n_atoms = aasize[ type - 1];
          rot->n_rotamers = max;

          lib_name = (char *) calloc(sizeof(char),100);

          switch(type)
          {
             case 1: strcpy(lib_name,  "lib/thr.pdb"); break;
             case 2: strcpy(lib_name,  "lib/val.pdb"); break;
             case 3: strcpy(lib_name,  "lib/ser.pdb"); break;
             case 4: strcpy(lib_name,  "lib/his.pdb"); break;
             case 5: strcpy(lib_name,  "lib/lys.pdb"); break;
             case 6: strcpy(lib_name,  "lib/arg.pdb"); break;
             case 7: strcpy(lib_name,  "lib/met.pdb"); break;
             case 8: strcpy(lib_name,  "lib/cys.pdb"); break;
             case 9: strcpy(lib_name,  "lib/glu.pdb"); break;
             case 10: strcpy(lib_name, "lib/gln.pdb"); break;
             case 11: strcpy(lib_name,  "lib/asp.pdb"); break;
             case 12: strcpy(lib_name,  "lib/asn.pdb"); break;
             case 13: strcpy(lib_name,  "lib/phe.pdb"); break;
             case 14: strcpy(lib_name,  "lib/tyr.pdb"); break;
             case 15: strcpy(lib_name,  "lib/trp.pdb"); break;
             case 16: strcpy(lib_name,  "lib/pro.pdb"); break;
             case 17: strcpy(lib_name,  "lib/gly.pdb"); break;
             case 18: strcpy(lib_name,  "lib/ala.pdb"); break;
             case 19: strcpy(lib_name,  "lib/leu.pdb"); break;
             case 20: strcpy(lib_name,  "lib/ile.pdb"); break;
          }

           printf("Library is: %s.\n",lib_name);

          lib = fopen(lib_name,"rb");
          if( lib == NULL)
          {
               printf("Error. Cannot open rotamer library.\n");
               return -1;
          }

          rot->x = (float *) calloc(sizeof(float),(max*rot->n_atoms)+10);
          rot->y = (float *) calloc(sizeof(float),(max*rot->n_atoms)+10);
          rot->z = (float *) calloc(sizeof(float),(max*rot->n_atoms)+10);
          line = (char *) calloc(sizeof(char),1024);

           printf("Max is: %i\n",max*rot->n_atoms);
          for( i = 0; i < max*rot->n_atoms; ++i)
          {
                 fgets(line,1024,lib);
                 strncpy(myx,&line[29],10);
                 strncpy(myy,&line[38],10);
                 strncpy(myz,&line[46],10);

                 rot->x[i] = atof(myx);
                 rot->y[i] = atof(myy);
                 rot->z[i] = atof(myz);

/*                 printf("%i - %f,%f,%f.\n",i,rot->x[i],rot->y[i],rot->z[i]);*/
          }
          fclose(lib);
          free(line);
          free(lib_name);


}



int super_transform(MOL2 *mol,RESIDUE *tmp_r, float target1[3], float target2[3],int ca2, int cb)
{
        RESIDUE *res;
        float ca;
        float ref[2][3];
        float from[3],to[3];
        float v[3],vs[3],vt[3];
        float norm = 0.0;
        float dx,dy,dz;
        res = tmp_r;

        /* CA ? Maybe. Depends on nitrogen protonation state. Do not mind at all.*/
        ref[0][0] = mol->x[ca2];
        ref[0][1] = mol->y[ca2];
        ref[0][2] = mol->z[ca2];
        /* CB  Same shity sutff */
        ref[1][0] = mol->x[cb];
        ref[1][1] = mol->y[cb];
        ref[1][2] = mol->z[cb];

        dx= ref[0][0];
        dy= ref[0][1];
        dz= ref[0][2];

        to[0] = ref[0][0] - ref[1][0];
        to[1] = ref[0][1] - ref[1][1];
        to[2] = ref[0][2] - ref[1][2];

        norm = sqrt(to[0]*to[0]+ to[1]*to[1] + to[2]*to[2]);
        to[0] /= norm;
        to[1] /= norm;
        to[2] /= norm;


        from[0] = target1[0] - target2[0];
        from[1] = target1[1] - target2[1];
        from[2] = target1[2] - target2[2];

        norm = sqrt(from[0]*from[0]+ from[1]*from[1] + from[2]*from[2]);

/*        printf("Norma: %f.\n",norm);*/
        from[0] /= norm;
        from[1] /= norm;
        from[2] /= norm;


        vs[0] = from[1]*to[2] - from[2]*to[1];
        vs[1] = (-from[0]*to[2] + from[2]*to[0]);
        vs[2] = from[0]*to[1] - from[1]*to[0];

        v[0] = vs[0];
        v[1] = vs[1];
        v[2] = vs[2];

        norm = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;

        ca = from[0]*to[0] + from[1]*to[1] + from[2] * to[2];

/*        printf("Angle %f.\n",acos(ca)*180/3.14159);*/


        vt[0] = v[0] * (1.0f - ca);
        vt[1] = v[1] * (1.0f - ca);
        vt[2] = v[2] * (1.0f - ca);

        res->rotM[0][0] = vt[0] * v[0] + ca;
        res->rotM[1][1] = vt[1] * v[1] + ca;
        res->rotM[2][2] = vt[2] * v[2] + ca;

        vt[0] *= v[1];
        vt[2] *= v[0];
        vt[1] *= v[2];

        res->rotM[0][1] = vt[0] - vs[2];
        res->rotM[0][2] = vt[2] + vs[1];
        res->rotM[1][0] = vt[0] + vs[2];
        res->rotM[1][2] = vt[1] - vs[0];
        res->rotM[2][0] = vt[2] - vs[1];
        res->rotM[2][1] = vt[1] + vs[0];

        res->transvector[0] = dx;
        res->transvector[1] = dy;
        res->transvector[2] = dz;


        return 0;
}


int name_to_number_residue(char tmp_resn[10])
{

          int i;

          for(i = 0; i < 20; ++i)
          {
                if( residues[i][0]  == tmp_resn[0] && residues[i][1]  == tmp_resn[1] && residues[i][2]  == tmp_resn[2] )
                   return i+1;
          }

        return -1;
}

int check_residue(int *resn,int total_res,int tmp_resn_i)
{
         int i;

         for( i = 0; i < total_res; ++i)
         {
               if( resn[i] == tmp_resn_i)
               return i;
         }

         return -1;
}


int check_angle(MOL2 **mmol, int tor, int i, int j, int k)
{
          MOL2 *mol;
          mol = *mmol;
          int step = 0;

          for( step=0; step < tor; ++step)
          {
                if( (mol->ia[step] == i && mol->ja[step] == j && mol->ka[step] == k) || (mol->ia[step] == k && mol->ja[step] == j && mol->ka[step] == i) )
                  return 1;

          }
        
return 0;
}

int check_dihedral(MOL2 **mmol, int tor, int i, int j, int k, int l)
{
          MOL2 *mol;
          mol = *mmol;
          int step = 0;

          for( step=0; step < tor; ++step)
          {
                if( (mol->ik[step] == i && mol->jk[step] == j && mol->kk[step] == k && mol->lk[step] == l) || (mol->ik[step] == l && mol->jk[step] == k && mol->kk[step] == j && mol->lk[step] == i) )
                  return 1;

          }
          return 0;
}

int import_topology(MOL2 **mymol)
{
       MOL2 *mols;
       FILE *in1 = NULL;
       char *line = NULL;
       int n_atoms = 0;
       int n_angles = 0;
       int n_dihedrals = 0;
       int n_bonds = 0;
       int n_pairs = 0; 
       int i = 0;
       int j = 0;
       int k = 0;
       int l = 0;

       mols = *mymol;
       line = (char *) calloc(sizeof(char),1024);

       in1 = fopen("types.top","rb");
       if( in1 == NULL)
       {
         free(line);
         return -1;
       }
       while( fgets(line,1024,in1))
       {
            ++n_atoms;
       }
       rewind(in1);

       printf("Old atoms: %i. New atoms: %i.\n",mols->n_atoms,n_atoms);
/*       if( (mols->n_atoms) != n_atoms);
       {
            printf("Number of atoms mismatch.\n");
            return -1;
       }*/

       n_atoms = 0;
       while( fgets(line,1024,in1))
       {
            sscanf(line,"%*i %i",&i);
            mols->gaff_types[n_atoms] = i;
/*            printf("Tipos: %i %i.\n",n_atoms,mols->gaff_types[n_atoms]);*/

            n_atoms++;
       }
       fclose(in1);

       n_atoms = 0;
       in1 = fopen("bonds.top","rb");
       if( in1 == NULL)
       {
         free(line);
         return -1;
       }
       while( fgets(line,1024,in1))
       {
            ++n_atoms;
       }
       rewind(in1);

       mols->n_bonds = n_atoms;
       mols->bonds = (int *) calloc(sizeof(int),mols->n_bonds);
       mols->bond_a1 = (int *) calloc(sizeof(int),mols->n_bonds);
       mols->bond_a2 = (int *) calloc(sizeof(int),mols->n_bonds);

       n_atoms = 0;
       while( fgets(line,1024,in1))
       {
            sscanf(line,"%i %i %i",&i,&j,&k);
            mols->bond_a1[n_atoms] = i;
            mols->bond_a2[n_atoms] = j;
            mols->bonds[n_atoms] = k;

/*            sscanf(line,"%i %i %i",mols->bond_a1[n_atoms],mols->bond_a2[n_atoms],mols->bonds[n_atoms]);*/
            n_atoms++;
       }
       fclose(in1);
       printf("Bonds in.\n");
       n_atoms = 0;
       in1 = fopen("angles.top","rb");
       if( in1 == NULL)
       {
         free(line);
         return -1;
       }
       while( fgets(line,1024,in1))
       {
            ++n_atoms;
       }
       rewind(in1);

       mols->n_angles = n_atoms;
       mols->ia = (int *) calloc( sizeof(int), mols->n_angles);
       mols->ja = (int *) calloc( sizeof(int), mols->n_angles);
       mols->ka = (int *) calloc( sizeof(int), mols->n_angles);

       n_atoms = 0;
       while( fgets(line,1024,in1))
       {
            sscanf(line,"%i %i %i",&i,&j,&k);
            mols->ia[n_atoms] = i;
            mols->ja[n_atoms] = j;
	    mols->ka[n_atoms] = k;

/*            sscanf(line,"%i %i %i",mols->ia[n_atoms],mols->ja[n_atoms],mols->ka[n_atoms]);*/
            n_atoms++;
       }
       fclose(in1);
       printf("Angles in.\n");
       n_atoms = 0;
       in1 = fopen("dihedrals.top","rb");
       if( in1 == NULL)
       {
         free(line);
         return -1;
       }
       while( fgets(line,1024,in1))
       {
            ++n_atoms;
       }
       rewind(in1);

       mols->n_torsionals = n_atoms;
       mols->ik = (int *) calloc( sizeof(int), mols->n_torsionals);
       mols->jk = (int *) calloc( sizeof(int), mols->n_torsionals);
       mols->kk = (int *) calloc( sizeof(int), mols->n_torsionals);
       mols->lk = (int *) calloc( sizeof(int), mols->n_torsionals);

       n_atoms = 0;
       while( fgets(line,1024,in1))
       {
            sscanf(line,"%i %i %i %i",&i,&j,&k,&l);
            mols->ik[n_atoms] = i;
            mols->jk[n_atoms] = j;
            mols->kk[n_atoms] = k;
            mols->lk[n_atoms] = l;


/*            sscanf(line,"%i %i %i %i",mols->ik[n_atoms],mols->jk[n_atoms],mols->kk[n_atoms],mols->lk[n_atoms]);*/
            n_atoms++;
       }
       fclose(in1);
       printf("Torsionals in.\n");

       n_atoms = 0;


       n_atoms = 0;
       in1 = fopen("vdwpairs.top","rb");
       if( in1 == NULL)
       {
         free(line);
         return -1;
       }
       while( fgets(line,1024,in1))
       {
            ++n_atoms;
       }
       rewind(in1);

       mols->n_pairs = n_atoms;
       mols->vdwpairs_a = (int *) calloc( sizeof(int), mols->n_pairs);
       mols->vdwpairs_b = (int *) calloc( sizeof(int), mols->n_pairs);
       mols->vdw_type = (int *) calloc( sizeof(int), mols->n_pairs);

       n_atoms = 0;
       while( fgets(line,1024,in1))
       {
            sscanf(line,"%i %i %i",&i,&j,&k);
            mols->vdwpairs_a[n_atoms] = i;
            mols->vdwpairs_b[n_atoms] = j;
            mols->vdw_type[n_atoms] = k;
/*            sscanf(line,"%i %i %i",mols->vdwpairs_a[n_atoms],mols->vdwpairs_b[n_atoms],mols->vdw_type[n_atoms]);*/
            n_atoms++;
       }
       fclose(in1);
       printf("Pairs in.\n");

       printf("Import: \n");
       printf("Atoms: %i.\n",mols->n_atoms);
       printf("Bonds: %i.\n",mols->n_bonds);
       printf("Angles: %i.\n",mols->n_angles);
       printf("Tortsionals: %i.\n",mols->n_torsionals);
       printf("VDW Pairs: %i.\n",mols->n_pairs);

       free(line);
       return 0;
}

int dump_topology(MOL2 mol)
{
       FILE *types = NULL;
       FILE *bonds = NULL;
       FILE *angles = NULL;
       FILE *dihedrals = NULL;
       FILE *pairs = NULL;
       int i = 0;

       if( (bonds = fopen("bonds.top","w")) == NULL)
       {
           printf("Error writing bonds.\n");
       }else{
           for(i = 0; i < mol.n_bonds; ++i)
           {
                 fprintf(bonds,"%i %i %i\n",mol.bond_a1[i],mol.bond_a2[i],mol.bonds[i]);
           }
           fclose(bonds);
       }

       if( (types = fopen("types.top","w")) == NULL)
       {
           printf("Error writing types.\n");
       }else{
           for(i = 0; i < mol.n_atoms; ++i)
           {
                 fprintf(types,"%i %i\n",i,mol.gaff_types[i]);
           }
           fclose(types);
       }

       if( (angles = fopen("angles.top","w")) == NULL)
       {
           printf("Error writing angles.\n");
       }else{
           for(i = 0; i < mol.n_angles; ++i)
           {
                 fprintf(angles,"%i %i %i\n",mol.ia[i],mol.ja[i],mol.ka[i]);
           }
           fclose(angles);
       }

       if( (dihedrals = fopen("dihedrals.top","w")) == NULL)
       {
           printf("Error writing dihedrals.\n");
       }else{
           for(i = 0; i < mol.n_torsionals; ++i)
           {
                 fprintf(dihedrals,"%i %i %i %i\n",mol.ik[i],mol.jk[i],mol.kk[i],mol.lk[i]);
           }
           fclose(dihedrals);
       }

       if( (pairs = fopen("vdwpairs.top","w")) == NULL)
       {
           printf("Error writing VDW pairs.\n");
       }else{
           for(i = 0; i < mol.n_pairs; ++i)
           {
                 fprintf(pairs,"%i %i %i\n",mol.vdwpairs_a[i],mol.vdwpairs_b[i],mol.vdw_type[i]);
           }
           fclose(pairs);
       }

       return 0;
}



int dump_pdb(MOL2 mol)
{
	int i;
        
        for(i = 0; i < mol.n_atoms; ++i)
        {
          printf("ATOM %6i  %2s  UNK     0    %8.3f%8.3f%8.3f\n",i+1,elems[ mol.atoms[i]-1],mol.x[i],mol.y[i],mol.z[i]);

        }

}

int q14_bond( MOL2 *mol, int atom, int second_atom)
{
          int j,k;
          int flag = 0;
          int vecinos[4];

           for(j=0; j < 4; ++j)
           {
             vecinos[j] = 0;
           }
           k=0;
           for(j=0; j < mol->n_bonds; ++j)
           {
                   if( mol->bond_a1[j] == (atom+1))
                   {
                       vecinos[k] = mol->bond_a2[j];
                       k++;
                   }else if( mol->bond_a2[j] == (atom+1)){
                       vecinos[k] = mol->bond_a1[j];
                       k++;
                   }
           }

           for(j=0; j < k; ++j)
           {
              if( q13_bond(mol, vecinos[j] -1 , second_atom) > 0)
              {
                  flag = 1;
              }
           }

           return flag;


}

int bonded( MOL2 mol, int atom, int second_atom)
{
           int j;

           for(j=0; j < mol.n_bonds; ++j)
           {
                  if( (mol.bond_a1[j] == (atom+1) && mol.bond_a2[j] == (second_atom+1)) || (mol.bond_a2[j] == (atom+1) && mol.bond_a1[j] == (second_atom+1)))
                       return 1;
           }
           return 0;

}

int q13_bond( MOL2 *mol, int atom, int second_atom)
{
           int j,k;
           int vecinos[4];
           int flag = 0;

           for(j=0; j < 4; ++j)
           {
             vecinos[j] = 0;
           }
           k=0;
           for(j=0; j < mol->n_bonds; ++j)
           {
                   if( mol->bond_a1[j] == (atom+1))
                   {
                       vecinos[k] = mol->bond_a2[j];
                       k++;
                   }else if( mol->bond_a2[j] == (atom+1)){
                       vecinos[k] = mol->bond_a1[j];
                       k++;
                   }
           }

           for(j=0; j < k; ++j)
           {
              if( get_number_any_bond_p(mol, vecinos[j] , second_atom+1) > 0)
              {
/*                  printf("  Excluidos: 1,3 %i - %i - %i.\n",atom,vecinos[j],second_atom);*/
                  flag = 1;
              }
           }

           return flag;

}

int get_bonds( MOL2 mol, int atom, int bond_type)
{
        int res = 0;
        int i = 0;
        for( i = 0; i < mol.n_bonds; ++i)
        {
                if( (mol.bond_a1[i] == atom || mol.bond_a2[i] == atom ) && (bond_type == 0 || mol.bonds[i] == bond_type))
                    ++res;
        }
        return res;

}


int get_bonds_p( MOL2 *mol, int atom, int bond_type)
{
        int res = 0;
        int i = 0;
        for( i = 0; i < mol->n_bonds; ++i)
        {
                if( (mol->bond_a1[i] == atom || mol->bond_a2[i] == atom ) && (bond_type == 0 || mol->bonds[i] == bond_type))
                    ++res;
        }
        return res;

}

int get_number_any_bond_p( MOL2 *mol, int atom, int second_atom) /* Second atom = 0 -> Any */
{
        int res = 0;
        int i = 0;
        for( i = 0; i < mol->n_bonds; ++i)
        {
                if( (mol->bond_a1[i] == atom &&  (mol->bond_a2[i] == second_atom || second_atom == 0) ) || (mol->bond_a2[i] == atom && (mol->bond_a1[i] == second_atom || second_atom == 0)) )
                {
                        ++res;
                }
        }
        return res;
}


int get_number_any_bond( MOL2 mol, int atom, int second_atom) /* Second atom = 0 -> Any */
{
        int res = 0;
        int i = 0;
        for( i = 0; i < mol.n_bonds; ++i)
        {
                if( (mol.bond_a1[i] == atom &&  (mol.bond_a2[i] == second_atom || second_atom == 0) ) || (mol.bond_a2[i] == atom && (mol.bond_a1[i] == second_atom || second_atom == 0)) )
                {
                        ++res;
                }
        }
        return res;
}

int get_number_bond_by_atom(MOL2 mol, int atom, int bond_type, int atom_type)
{
        int res = 0;
        int i = 0;
        for( i = 0; i < mol.n_bonds; ++i)
        {
            if( mol.bonds[i] == bond_type || bond_type == 0)
            {
              if( mol.bond_a1[i] == atom){
                  if( mol.atoms[(mol.bond_a2[i]-1)] == atom_type)
                    ++res;
              }else if( mol.bond_a2[i] == atom){
                  if( mol.atoms[(mol.bond_a1[i]-1)] == atom_type)
                    ++res;
              }
            }
        }
        return res;
}


int get_number_type_bonds( MOL2 mol, int atom, int type, int second_atom) /* Second atom = 0 -> Any */
{
        int res = 0;
        int i = 0;
        for( i = 0; i < mol.n_bonds; ++i)
        {
                if( (mol.bond_a1[i] == atom && ( mol.atoms[(mol.bond_a2[i]-1)] == second_atom || second_atom == 0) && ( mol.bonds[i] == type || type == 0 )) || (mol.bond_a2[i] == atom && ( mol.atoms[(mol.bond_a1[i]-1)] == second_atom || second_atom == 0) && ( mol.bonds[i] == type || type == 0) ))
                {
                        ++res;
                }
        }
        return res;
}


int last_link(int *path)
{
        int i = 0;
        for( i = 0; i < MAX_RING_MEMBERS; ++i)
        {
                if( path[i] < 0)
                {
                        return i;
                }
        }

        return 0;
}


/* Only for debugging */
int show_path(int *path)
{
        int i = 0;
        for(i=0 ; i < MAX_RING_MEMBERS; ++i)
        {
                printf("%i,",path[i]);
        }
        printf("\n");
}



int find_rings(MOL2 mol, int *path, int n_atom, int last_n_atom, int num_rings,int **myik,int **myjk,int **mykk,int **mylk,
               int *elem)
{
        int *new_path = NULL;
        int children = 0;
        int atom_number = 0;
        int i = 0;
        int j = 0;
        int k = 0;
        int i2 = 0;
        int j2 = 0;
	int k2 = 0;
	int tmp_num_rings = 0;
	float com_x = 0.0f;	
	float com_y = 0.0f;
	float com_z = 0.0f;
        int *ik = NULL;
        int *jk = NULL;
        int *kk = NULL;
        int *lk = NULL;
        int myelem;

	int flag = 0;
        int flag2 = 0;
        int dc = 0;

	tmp_num_rings = num_rings;
        ik = *myik;
        jk = *myjk;
        kk = *mykk;
        lk = *mylk;

	/* Append new member */
        k = last_link(path);
        path[k] = n_atom;
        myelem = *elem;
        myelem++;

	if( tmp_num_rings > MAX_RINGS)
	    return 0;
	/* Check ring closure */
/*        show_path(path);*/
        for( i2 = 0; i2 < (MAX_RING_MEMBERS ); ++i2)
        {
           for( j2 = 0; j2 < MAX_RING_MEMBERS; ++j2)
           {
                        if( path[i2] == path[j2] && path[i2] != -1 && i2!=j2)  /* Ring closure */
                        {

				tmp_num_rings = tmp_num_rings + 1;
				return 1;
                        }
            }

        }

        if( (children = get_number_type_bonds(mol,n_atom,0,0)) >= 0 ) /* Prune */
        {
                for(i = 0; i < mol.n_bonds; ++i)
                {
                        if( mol.bond_a1[i] == n_atom || mol.bond_a2[i] == n_atom )
                        {
                            if( mol.bond_a1[i] != last_n_atom && mol.bond_a2[i] != last_n_atom)
                            {
                                if( (new_path = calloc(MAX_RING_MEMBERS,sizeof(int))) == NULL)
                                {
                                         fprintf(stderr,"Error. Cant allocate memory.\n");
                                         fflush(stderr);
                                         return -2;
                                }
                                for(k = 0; k < MAX_RING_MEMBERS; ++k)
                                {
                                        new_path[k] = path[k];
                                }

                                if( mol.bond_a1[i] == n_atom)
                                                atom_number = mol.bond_a2[i];
                                else
                                                atom_number = mol.bond_a1[i];

				tmp_num_rings++;
                                if( *elem >= 3)
                                {
                                    flag = flag2 = -1;
                                    for( dc = 0; dc < 1000; ++dc)
                                    {
                                      if(ik[dc] == -1)
                                      {
                                        flag2 = dc;
                                        break;
                                      }
                                      if( jk[dc] == (new_path[myelem-2]-1) && kk[dc] == (new_path[myelem-1]-1))
                                       flag = 1;
                                    }
                                    if( flag != 1)
                                    {
                                          ik[flag2] = new_path[myelem-3]-1;
                                          jk[flag2] = new_path[myelem-2]-1;
                                          kk[flag2] = new_path[myelem-1]-1;
                                          lk[flag2] = new_path[myelem]-1;
/*                                    printf("Dihedral: %i - %i - %i - %i \n",new_path[myelem-3],new_path[myelem-2],new_path[myelem-1],new_path[myelem]);*/

                                    }

                                }
                                if( find_rings(mol,new_path, atom_number ,n_atom, tmp_num_rings,&ik,&jk,&kk,&lk,&myelem) < 0)
                                        {
                                                fprintf(stderr,"Error. Rescursive fail.\n");
                                                fflush(stderr);
                                                return -5;
                                        }
                                        atom_number = 0;
                                free(new_path);

                            }else{}

                        }
                }
        }else{
                printf("Atomo %i tiene %i enlaces.\n",n_atom,children);
                return 0;
        }

        return 0;
}


int get_number_of_rings(MOL2 mol,int **myik, int **myjk, int **mykk, int **mylk)
{
        int *new_path = NULL;
        int i = 0;
        int *ik, *jk, *kk, *lk;
        ik = *myik;
        jk = *myjk;
        kk = *mykk;
        lk = *mylk;
        int elem;
        elem = 0;

/*        printf("Iniciando...\n");*/

        for( i=0; i< 1000; i++)
        {
              ik[i] = -1;
              jk[i] = -1;
              kk[i] = -1;
              lk[i] = -1;
        }
        if( (new_path = calloc(MAX_RING_MEMBERS,sizeof(int))) == NULL)
        {
                fprintf(stderr,"Error. Cant allocate memory.\n");
                fflush(stderr);
                return -2;
        }
	
        for(i = 0 ; i < MAX_RING_MEMBERS; ++i)
                new_path[i] = -1;
        new_path[0] = 0;
        if( find_rings(mol,new_path,1,0,0,&ik,&jk,&kk,&lk,&elem) < 0)
        {
                fprintf(stderr,"Error. Recursive fail.\n");
                fflush(stderr);
                return -3;
        }
	free(new_path);

        for(i = 0; i< 1000; ++i)
        {
              if( ik[i] == -1)
                break;
        }

        return i;
}

