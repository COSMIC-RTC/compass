//kp_smatrix.cpp


#ifndef __SEGER__KP_SMATRIX_HPP__
#define __SEGER__KP_SMATRIX_HPP__

template <typename real>
kp_smatrix<real>::kp_smatrix()
{
   _create(0,0,0);
}
//                                                                                            
template <typename real>
template <typename T>
kp_smatrix<real>::kp_smatrix(const kp_smatrix<T>& M)
{
   _create(0,0,0);
   this->operator=(M);
}
template <typename real>
kp_smatrix<real>::kp_smatrix(const kp_smatrix<real>& M)
{
   _create(0,0,0);
   this->operator=(M);
}


//                                                                                            
template <typename real>
kp_smatrix<real>::~kp_smatrix()
{
   _clear();
}
//                                                                                           
template <typename real>
void kp_smatrix<real>::resize(int nnz_, int dim1_, int dim2_)
{
   if (nnz != nnz_)
     {
	_clear();
	_create(nnz_, dim1_, dim2_);
     }
   else
     {
	dim1 = dim1_;
	dim2 = dim2_;
	majorDim = 'U';
     }
}
//                                                                                      
template <typename real>
template <typename T>
void kp_smatrix<real>::operator=(const kp_smatrix<T>& M)
{
   resize(M.nnz, M.dim1, M.dim2);
   memcpy(rowind, M.rowind, sizeof(int)  * nnz);   
   memcpy(colind, M.colind, sizeof(int)  * nnz);
   for (int i=0 ; i<nnz ; i++)
      values[i] = (real) M.values[i];
        
   majorDim = M.get_majorDim();
   
}
template <typename real>
void kp_smatrix<real>::operator=(const kp_smatrix<real>& M)
{
   resize(M.nnz, M.dim1, M.dim2);
   memcpy(values, M.values, sizeof(real) * nnz);
   memcpy(rowind, M.rowind, sizeof(int)  * nnz);   
   memcpy(colind, M.colind, sizeof(int)  * nnz);
   
   majorDim = M.majorDim;
}

//                                                                                    
template <typename real>
void kp_smatrix<real>::init_from_rowidx(const kp_smatrix<real>& M, const vector<int>& idx)
{
   if (this == &M)
     {
	cerr<<"Error | kp_smatrix::init_from_rowidx | the same matrix"<<endl;
	exit(EXIT_FAILURE);
     }
   for (size_t i = 0 ; i < idx.size() ; i++)
     if (idx[i] < 0 || idx[i] >= M.dim1)
       {
	  cerr<<"Error | kp_smatrix::init_from_rowidx | index error"<<endl;
	  cout<<idx[i]<<" "<<M.dim1<<endl;
	  exit(EXIT_FAILURE);
       }
   
   for (size_t i = 1 ; i < idx.size() ; i++)
     if (idx[i] <= idx[i - 1])
       {
	  cerr<<"Error | kp_smatrix::init_from_rowidx | idx haven't been sorted"<<endl;
	  exit(EXIT_FAILURE);
       }
   
   //sort rows
   vector< pair<int, pair<int, double> > > sM(M.nnz);
   //sM[i].first         -> rowind
   //sM[i].second.first  -> colind
   //sM[i].second.second -> value
   for (int i = 0 ; i < M.nnz ; i++)
     {
	sM[i].first         = M.rowind[i];
	sM[i].second.first  = M.colind[i];
	sM[i].second.second = M.values[i];
     }
   sort(sM.begin(), sM.end());
   
   
   //make first loop to calculate number of non zero elements
   size_t i_idx = 0;
   int i_this = 0;
   //idx in zero-based indexing but matrix in one-based indexing
   for (int i = 0 ; i < M.nnz ; i++)
     {
	while (i_idx < idx.size() && (idx[i_idx] + 1) < sM[i].first)
	  i_idx++;
	if (i_idx < idx.size() && sM[i].first == idx[i_idx] + 1)	
	  {
	     i_this++;
	  }
     }   
   resize(i_this, idx.size(), M.dim2);
      
   //make second loop to initialize non zero elements
   
   i_idx  = 0;
   i_this = 0;
   for (int i = 0 ; i < M.nnz ; i++)
     {
	while (i_idx < idx.size() && idx[i_idx] + 1 < sM[i].first)
	  i_idx++;
	if (i_idx < idx.size() && sM[i].first == idx[i_idx] + 1)	
	  {
	     rowind[i_this] = i_idx + 1;
	     colind[i_this] = sM[i].second.first;
	     values[i_this] = sM[i].second.second;
	     i_this++;
	  }
     }   
   
   majorDim = 'R'; // *this est une matrice row-major a la fin de cette methode
   
}
//                                                                                            
template <typename real>
void kp_smatrix<real>::init_from_transpose(const kp_smatrix<real>&M)
{
   resize(M.nnz, M.dim2, M.dim1);
   memcpy(values, M.values, sizeof(real) * nnz);
   memcpy(colind, M.rowind, sizeof(int) * nnz);
   memcpy(rowind, M.colind, sizeof(int) * nnz);
   
   if (M.majorDim == 'C') 
     majorDim = 'R';
   else if (M.majorDim == 'R') 
     majorDim = 'C';
   else 
     majorDim = 'U';
   
   
   
   
}
//                                                                                           

template <typename real>
void kp_smatrix<real>::init_from_matrix(const kp_matrix<real>& B, double epsilon /* = 0.0 */)
{
	int nb_el = B.getDim1()*B.getDim2();
	int*  row = new int[nb_el]; 
	int*  col = new int[nb_el]; 
	real* val = new real[nb_el];
	int cmpt = 0;
	for(int j=0 ; j <B.getDim2() ; j++)
	{
		for (int i=0 ; i<B.getDim1() ; i++)
		{
			if (fabs(B(i,j)) > epsilon)
			{
				row[cmpt] = i+1;
				col[cmpt] = j+1;
				val[cmpt] = B(i,j);
				cmpt++;
			}
		}
	}
	resize(cmpt,B.getDim1(),B.getDim2());
	for (int i=0 ; i<nnz ; i++)
	{
		rowind[i] = row[i];
		colind[i] = col[i];
		values[i] = val[i];
	}

	delete[] row;
	delete[] col;
	delete[] val;
}

 
template <typename real>
void kp_smatrix<real>::check()
{
   cerr<<"CHECK smatrix (DEBUG)"<<endl;
   for (int i = 0 ; i < nnz ; i++)
     {
	if (rowind[i] < 1 || rowind[i] > dim1 ||
	    colind[i] < 1 || colind[i] > dim2)
	  {	    
	     cerr<<"Error | kp_smatrix::check | bad index"<<endl;
	     cerr<<dim1<<"x"<<dim2<<" "<<rowind[i]<<"x"<<colind[i]<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
}
//                                                                                           
template <typename real>
void kp_smatrix<real>::_create(int nnz_, int dim1_, int dim2_)
{
   nnz  = nnz_;
   dim1 = dim1_;
   dim2 = dim2_;
   values = new real[nnz];
   rowind = new int[nnz];
   colind = new int[nnz];
   majorDim = 'U' ;
}
//                                                                                             
template <typename real>
void kp_smatrix<real>::_clear()
{
   if (values == NULL || rowind == NULL || colind == NULL)
     {
	cerr<<"Error | kp_smatrix::_clear | double clear"<<endl;
	exit(EXIT_FAILURE);
     }
   delete [] values;
   delete [] rowind;
   delete [] colind;      
   values   = NULL;
   rowind   = NULL;
   colind   = NULL;
   nnz      = 0;
   dim1     = 0;
   dim2     = 0;
   majorDim = 'U';
}


template <typename real>
void kp_gemv(real alpha, kp_smatrix<real>& A, kp_vector<real>& x, real beta, kp_vector<real>& y)
{

char* NotTrans = new char[2];
NotTrans[0] = 'N';
NotTrans[1] = '\0';

char* matDescra = new char[8];
matDescra[0] = 'G';
for (int i=1 ; i<=6 ;i++)
	matDescra[i] = 'F';
matDescra[7] = '\0';


   if (A.dim2 != x.size() || A.dim1 != y.size())
     {
	cerr<<"Error | kp_cscmv | diminsion problem"<<endl;
	exit(EXIT_FAILURE);
     }

  kp_gemv_core(NotTrans, &A.dim1, &A.dim2, &alpha, matDescra,
		  A.values, A.rowind, A.colind, &A.nnz,
		  x.getData(),  &beta, y.getData());

  delete[]  NotTrans ; NotTrans=NULL;
  delete[]  matDescra ; matDescra=NULL;
}

template <typename real>
void kp_gemm(char op_A, real alpha, kp_smatrix<real>& A, kp_matrix<real>& B, real beta, kp_matrix<real>& C)
{

char* transA = new char[2];
transA[0] = op_A;
transA[1] = '\0';

char* matDescra = new char[8];
matDescra[0] = 'G';
for (int i=1 ; i<=6 ;i++)
	matDescra[i] = 'F';
matDescra[7] = '\0';

   int B_dim1 = B.getDim1();
   int C_dim1 = C.getDim1();
   int C_dim2 = C.getDim2();


   int opA_dim1, opA_dim2;
   kp_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
   
   if (C.getDim1() != opA_dim1 || opA_dim2 != B_dim1 || C.getDim2() != B.getDim2())
     {
	cerr<<"Error | kp_gemm (sparce) | diminsion problem"<<endl;
	exit(EXIT_FAILURE);
     }


   kp_gemm_core(transA, &A.dim1, &C_dim2,  &A.dim2, &alpha, matDescra, 
	      A.values, A.rowind, A.colind, &A.nnz, 
	      B.getData(), &B_dim1, &beta, C.getData(), & C_dim1); 
   
  delete[]  transA ; transA=NULL;
  delete[]  matDescra ; matDescra=NULL;
}
//                                                                                           
template <typename real>
void kp_check_op_set_dim(int op, const kp_smatrix<real>&M, int& dim1, int& dim2)
{
   if (op == 'N')
     {
	dim1 = M.dim1;
	dim2 = M.dim2;
     }
   else if (op == 'T')
     {
	dim1 = M.dim2;
	dim2 = M.dim1;
     }
   else
     {
	cerr<<"Error | kp_check_op_set_dim | op should ge either N either T"<<endl;
	exit(EXIT_FAILURE);
     }
}
//                                                                                      
template <typename real>
void kp_smatrix<real>::resize2rowMajor()
{
    int i;
    int indice;
    vector<vector<pair<int,int> > > position(dim1, vector<pair<int,int> >());
    int* rowind2 = new int[nnz];
    int* colind2 = new int[nnz];
    real* values2 = new real[nnz];

	//cout<<"dim1="<<this->dim1<<" dim2="<<this->dim2<<" nnz="<<this->nnz<<endl;
	//cout<<"i="<<172248<<" rowind[i]="<<rowind[172248]<<" colind[i]="<<colind[172248]<<endl;

    for (i=0 ; i < nnz ; i++)
    {
	//cout<<"i="<<i<<" rowind[i]="<<rowind[i]-1<<" colind[i]="<<colind[i]-1<<endl;
        position[rowind[i]-1].push_back(make_pair(colind[i],i));

    }

    indice = 0;
    for (i = 0 ; i < dim1 ; i++)
    {
        if(position[i].size() > 0)
        {
            std::sort(position[i].begin(), position[i].end());
            for (unsigned int j = 0 ; j < position[i].size() ; j++)
            {
                rowind2[indice] = rowind[(position[i][j]).second] ;
                colind2[indice] = colind[(position[i][j]).second];
                values2[indice] = values[(position[i][j]).second];
                indice++;
            }
        }
    }
    if (indice != nnz)
    {
        cerr<<"Erreur | kp_smatrix::resize2rowMajor | erreur lors de la conversion."<<endl;
        exit(EXIT_FAILURE);
    }

    delete[] rowind ; rowind = rowind2;
    delete[] colind ; colind = colind2;
    delete[] values ; values = values2;
   
    majorDim= 'R';

}
template <typename real>
void kp_smatrix<real>::resize2colMajor()
{
    int i;
    int indice;
    vector<vector<pair<int,int> > > position(dim2, vector<pair<int,int> >());
    int* rowind2 = new int[nnz];
    int* colind2 = new int[nnz];
    real* values2 = new real[nnz];;

    for (i=0 ; i < nnz ; i++)
    {
        position[colind[i]-1].push_back(make_pair(rowind[i],i));
    }

    indice = 0;
    for (i = 0 ; i < dim2 ; i++)
    {
        if(position[i].size() > 0)
        {
            std::sort(position[i].begin(), position[i].end());
            for (unsigned int j = 0 ; j < position[i].size() ; j++)
            {
                colind2[indice] = colind[(position[i][j]).second] ;
                rowind2[indice] = rowind[(position[i][j]).second];
                values2[indice] = values[(position[i][j]).second];
                indice++;
            }
        }
    }
    if (indice != nnz)
    {
        cerr<<"Erreur | kp_smatrix::resize2colMajor | erreur lors de la conversion."<<endl;
        exit(EXIT_FAILURE);
    }

    delete[] rowind ; rowind = rowind2;
    delete[] colind ; colind = colind2;
    delete[] values ; values = values2;
   
    majorDim= 'C';

}

#endif
