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

int cleanup(MOL2 **mymols)
{
                int i = 0;
                MOL2 *mols;
                mols = *mymols;
                free(mols[i].x);
                free(mols[i].y);
                free(mols[i].z);
                free(mols[i].pcharges);
                free(mols[i].atoms);
                free(mols[i].ringer);
                free(mols[i].aromatic);
                free(mols[i].bonds);
                free(mols[i].bond_a1);
                free(mols[i].bond_a2);
                free(mols[i].bond_dist);
                free(mols[i].grads_X);
                free(mols[i].grads_Y);
                free(mols[i].grads_Z);
                free(mols[i].ik);
                free(mols[i].jk);
                free(mols[i].kk);
                free(mols[i].lk);
                free(mols[i].ia);
                free(mols[i].ja);
                free(mols[i].ka);
                free(mols[i].gaff_types);
                free(mols[i].backbone);
                free(mols[i].selection);
                free(mols[i].vdw_type);
                free(mols[i].vdwpairs_a);
                free(mols[i].vdwpairs_b);
                free(mols);
}
