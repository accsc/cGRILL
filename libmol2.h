/*****************************************************************
 *
 /***********************************************************************
 *
 *
 *	cGRILL. A simple affinity map generator written in C
 *
 *
 *	
 *	(C) 2011-2024. Alvaro Cortes Cabrera. <alvarocortesc@gmail.com>
 *
 *	 This code is licensed under GPL v3 terms
 *
 *
 */

#define SINGLE_BOND 1
#define DOUBLE_BOND 2
#define TRIPLE_BOND 3
#define AROMATIC 5


typedef struct{
int n_atoms;
int *atoms;
int type;
int ca;
int cb;
float transvector[3];
float rotM[3][3];
} RESIDUE;


typedef struct{
int n_atoms;
int n_rotamers;
float *x;
float *y;
float *z;
} ROTAMER;


typedef struct{

int n_atoms;
int n_bonds;
char *comment;  /* Ni p*** caso */
int charges;
int small_flag; /* Para que?*/

int *bonds;  
int *bond_a1;     
int *bond_a2;
float *bond_dist;
int *atoms;   
int *ringer;
int *aromatic;

int *ia;
int *ja;
int *ka;
int n_angles;

int *ik;
int *jk;
int *kk;
int *lk;
int n_torsionals;

int *vdwpairs_a;
int *vdwpairs_b;
int *vdw_type;
int n_pairs;

int *gaff_types;

float *pcharges;

float *x;
float *y;
float *z;

float *grads_X;
float *grads_Y;
float *grads_Z;

int *backbone;
int *selection;


} MOL2;

typedef struct{
int index;
float x;
float y;
float z;
struct PPP *next;
} PPP;

typedef struct _mol_link {
        int n_atom;
        struct _mol_link *next;
} MOL_LINK;


typedef struct{
        MOL_LINK *first_link;
        int members;
        int saturated;
        struct RING *next;
} RING;


char *residues[20] = {"THR","VAL","SER","HIS","LYS","ARG","MET","CYS","GLU","GLN","ASP","ASN","PHE","TYR","TRP","PRO","GLY","ALA","LEU","ILE"};
int Nrotamers[20] =  {  415,   38,  228,  545,  878,  812,  350,  222,  292,  496,  292,  543,  410, 1230,    0,   47,    0,    0,  256,  293};
int aasize[ 20 ] =   {   14,   16,   11,   18,   22,   24,   17,   11,   15,   17,   12,   14,   20,   21,   24,   14,    0,    0,   19,   19};

#define THR 1
#define VAL 2
#define SER 3
#define HIS 4
#define LYS 5
#define ARG 6
#define MET 7
#define CYS 8
#define GLU 9
#define GLN 10
#define ASP 11
#define ASN 12
#define PHE 13
#define TYR 14
#define TRP 15
#define PRO 16
#define GLY 17
#define ALA 18
#define LEU 19
#define ILE 20
