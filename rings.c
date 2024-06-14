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

#define MAX_RING_MEMBERS 300  /* Puff! but simple linked lists = hell! */
#define MAX_RINGS 30         /* Little hack? to prevent infinite recursive problem */



int find_rings2(MOL2 mol, int *path, int n_atom, int last_n_atom, int num_rings, int **ringer, int **aromatic)
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
        int *myringer = NULL;
        int *myaromatic = NULL;
        int elec = 0;
        int aroflag = 0;

        PPP *my_res = NULL;
        PPP *tmp_ppp = NULL;
        PPP *current_ppp = NULL;
	int flag = 0;
	PPP *sort_ppp = NULL;

	tmp_num_rings = num_rings;

        myringer = *ringer;
        myaromatic = *aromatic;
	/* Append new member */
        k = last_link(path);
        path[k] = n_atom;

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
				/* Aromatic ?*/
                                elec = 0;
                                aroflag = 0;
/*                                printf("Ring of %i elements.\n",j2-i2);*/
				for( k2 = i2 ; k2 < j2; ++k2)
				{
/*					if( mol.atoms[ path[k2]-1] != 4) 
					return 1;*/
                                        if(  j2 -i2 > myringer[path[k2]-1])
                                        myringer[path[k2]-1] = j2 -i2;
                                        if( get_bonds(mol,path[k2],0) > 3)
                                        aroflag = 1;

                                        if( mol.atoms[ path[k2]-1 ] == 1 )
                                          elec++;
                                        else if( mol.atoms[ path[k2]-1 ] == 6  )
                                          elec = elec +4;
                                        else if( mol.atoms[ path[k2]-1 ] == 3 && get_bonds(mol,path[k2],0) == 3 )
                                          elec=elec+2;
                                        else if( mol.atoms[ path[k2]-1 ] == 3 && get_bonds(mol,path[k2],0) < 3 )
                                          elec++;
                                        else if( mol.atoms[ path[k2]-1 ] == 6  )
                                          elec = elec +4;


				}

                                if( aroflag == 0 &&  (elec -2)%4 == 0 )
                                {
/*                                printf("AROMATIC RING.\n");*/
	                                for( k2 = i2 ; k2 < j2; ++k2)
	                                {
                                             myaromatic[path[k2]-1] = 1;
                                        }
                                }

/*                                        show_path(path);*/
/*                                printf("Ring found.\n");*/
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
                                if( find_rings2(mol,new_path, atom_number ,n_atom, tmp_num_rings,&myringer,&myaromatic) < 0)
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


int get_number_of_rings2(MOL2 mol,int **myringer, int **myaromatic)
{
        int *new_path = NULL;
        int i = 0;
	int *ringer = NULL;
        int *aromatic = NULL;
        ringer = *myringer;
        aromatic = *myaromatic;

/*        printf("Iniciando...\n");*/
        if( (new_path = calloc(MAX_RING_MEMBERS,sizeof(int))) == NULL)
        {
                fprintf(stderr,"Error. Cant allocate memory.\n");
                fflush(stderr);
                return -2;
        }
	
        for(i = 0 ; i < MAX_RING_MEMBERS; ++i)
                new_path[i] = -1;
        new_path[0] = 0;
        if( find_rings2(mol,new_path,1,0,0,&ringer,&aromatic) < 0)
        {
                fprintf(stderr,"Error. Recursive fail.\n");
                fflush(stderr);
                return -3;
        }
	free(new_path);

        return 0;
}

