/****************************************************************
 *
 *
/***********************************************************************
 *
 *	cGRILL. A simple affinity map generator written in C
 *
 *
 *	
 *	(c) 2011-2024. Alvaro Cortes Cabrera. <alvarocortesc@gmail.com>
 *
 *	 This code is licensed under GPL v3 terms
 *
 *
 ***************************************************************/


#include <metals.h>

float get_metal_vdw_energy(int type_metal, int gaff_type, float dist)
{

	float sigma = 0.0f;
	float epsilon = 0.0f;
	float r2 = 0.0f;
	float r6 = 0.0f;
	float r12 = 0.0f;
	float s6 = 0.0f, s12 = 0.0f;
        
	sigma = met_vdw_r[type_metal-100] + vdw_r[ gaff_type -1 ];
	epsilon = sqrt(met_vdw_epsilon[type_metal-100])*sqrt(vdw_epsilon[gaff_type-1]);

	sigma = sigma * sigma; /* **2 */
	sigma = sigma * sigma * sigma; /* **6 */
	
	dist = dist * dist;
	dist = dist * dist * dist; /* **6 */

	s6 = sigma/dist;
	s12 = s6 * s6;

	return (epsilon * ( s12 - (2.0f* s6)));
}

