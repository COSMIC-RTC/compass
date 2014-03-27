#include "A1_2vec.h"
void A1_2vec(kp_vector& A1_00, kp_vector& A1_01, const kp_smatrix& A1)
{
	int col, row;
	real val;
	int nb_az = A1_00.size();

	for (int i=0;i<A1.nnz;i++)
	{
		col = A1.colind[i];
		row = A1.rowind[i];
		val = A1.values[i];
		
		if((col <= nb_az) && (row <= nb_az))
		{
			if (col==row) A1_00.d[row-1]=val;
			else A1_00.d[row-1]=0;
		}
		else if (row <= nb_az)
		{
			if((col - nb_az) == row) A1_01.d[row-1]=val;
			else A1_01.d[row-1]=0;
		}
	}

}



