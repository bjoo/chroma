// -*- C++ -*-
// $Id: containers.h,v 1.12 2008-04-02 02:55:58 kostas Exp $

#ifndef _INV_CONTAINERS__H
#define _INV_CONTAINERS__H

#include "chromabase.h"

namespace Chroma
{
  namespace LinAlg
  {

    //--- OPT eigcg space ---//
    class OptEigInfo{
    public:
      int N   ; // vector dimension (large)
      int ncurEvals ;
      int lde ;
      multi1d<Real> evals ;
      multi1d<Complex> evecs ;
      multi1d<Complex> H     ;
      multi1d<Complex> HU    ;
      //OptEigInfo(){}
      void init(int ldh,int lde_, int nn){
	ncurEvals = 0 ;
	N=nn;
	lde = lde_ ;
	evals.resize(ldh);
	evecs.resize(lde*ldh);
	H.resize(ldh*ldh);
	HU.resize(ldh*ldh);
      }
    } ;

    //------------------------------------------------------------------------------
    //! Hold vectors
    /*! \ingroup invert */
    template<class T> class Vectors
    {
    public:
      multi1d<T> vec;
      int N; // number of active vectors size of vec must be larger or equal to N

      Vectors():N(0){} 
      Vectors(const multi1d<T>& v):vec(v),N(v.size()){} 
      Vectors(int size){resize(size); } 
    
      ~Vectors(){}
    
      void AddVector(const T& v,const Subset& s){
	if(N<vec.size()){
	  vec[N][s] = v;
	  N++;
	}
      }

      void NormalizeAndAddVector(const T& v,const Double& inorm, 
				 const Subset& s){
	if(N<vec.size()){// inorm is the inverse of the norm
	  vec[N][s] = v;
	  vec[N][s] *= inorm;
	  N++;
	}
      }
      void AddOrReplaceVector(const T& v,const Subset& s){
	if(N<vec.size()){
	  vec[N][s] = v;
	  N++;
	}
	else{// replace the last vector
	  vec[N-1] = v;
	}
      }

      // This will only add as many vectors as they fit
      void AddVectors(multi1d<T>& v,const Subset& s){
	for(int i(0);i<v.size();i++)
	  AddVector(v[i],s);
      }

    
      void resize(int n) {N=0;vec.resize(n);} 
      int size() const { return vec.size();}
      int Nvecs() const { return N; } 
      T& operator[](int i){ return vec[i];}
    };


    //------------------------------------------------------------------------------
    //! Holds eigenvalues and eigenvectors
    /*! \ingroup invert */
    template<class T> class RitzPairs
    {
    public:
      Vectors<Double> eval;
      Vectors<T>      evec;
      int Neig;

      RitzPairs() {init(0);}
      RitzPairs(int N) {init(N);}

      void init(int N) {
	eval.resize(N);
	evec.resize(N);
	Neig = 0;
      }

      void AddVector(const Double& e, const T& v,const Subset& s){
	eval.AddVector(e,s);
	evec.AddVector(v,s);
	if (eval.N != evec.N)
	{
	  QDPIO::cerr << __func__ << ": length of value and vector arrays are not the same" << endl;
	  QDP_abort(1);
	}
	Neig = evec.N;
      }

      // This will only add as many vectors as they fit
      void AddVectors(const multi1d<Double>& e,const multi1d<T>& v,const Subset& s){
	for(int i(0);i<e.size();i++)
	  eval.AddVector(e[i],s);
	for(int i(0);i<v.size();i++)
	  evec.AddVector(v[i],s);
	if (eval.N != evec.N)
	{
	  QDPIO::cerr << __func__ << ": length of value and vector arrays are not the same" << endl;
	  QDP_abort(1);
	}
	Neig = evec.N;
      }
    };


    //------------------------------------------------------------------------------
    //! This is a square matrix
    /*! \ingroup invert */
    template<class T> class Matrix
    {
    public:
      multi2d<T> mat;
      int N;  // active size 

      Matrix():N(0){} 
      Matrix(const multi2d<T>& v):mat(v),N(v.size1())
	{
	  if(v.size1() != v.size2())
	    QDPIO::cerr<<"WARNING!!! Matrix should be square!! CHECK YOUR CODE!\n";
	}

      Matrix(int size) {resize(size);}
    
      ~Matrix(){}
    
      void resize(int size) {N=0; mat.resize(size,size);}
      int size() const { return mat.size1();}
      int ld() const { return mat.size1();}
      int Nmat() const { return N; } 
      T& operator()(int i,int j){ return mat(i,j);}
    };
  

    //------------------------------------------------------------------------------
    //! Hold vectors
    /*! \ingroup invert */
    template<class T> class VectorArrays
    {
    public:
      multi2d<T> vec;
      int N;   // number of active vectors size of vec must be larger or equal to N
      int Ls;  // size of 5D arrays

      VectorArrays() : N(0) {} 
      VectorArrays(const multi2d<T>& v) : vec(v),N(v.size2()),Ls(v.size1()) {} 
      VectorArrays(int size2,int size1){resize(size2,size1);} 
    
      ~VectorArrays(){}
    
      void AddVector(const multi1d<T>& v,const Subset& s){
	if(N<vec.size2())
	{
	  if (v.size() != Ls)
	  {
	    QDPIO::cerr << "VectorArrays:" << __func__ << ": size of 5D array inconsistent with stored value: Ls=" << Ls 
			<< "  v.size=" << v.size() << endl;
	    QDP_abort(1);
	  }
	  for(int k=0; k < Ls; ++k)
	    vec[N][k][s] = v[k];
	  N++;
	}
      }

      void NormalizeAndAddVector(const multi1d<T>& v,const Double& inorm, 
				 const Subset& s){
	if(N<vec.size2())
	{
	  // inorm is the inverse of the norm
	  if (v.size() != Ls)
	  {
	    QDPIO::cerr << "VectorArrays:" << __func__ << ": size of 5D array inconsistent with stored value: Ls=" << Ls 
			<< "  v.size=" << v.size() << endl;
	    QDP_abort(1);
	  }
	  for(int k=0; k < Ls; ++k)
	  {
	    vec[N][k][s] = v[k];
	    vec[N][k][s] *= inorm;
	  }
	  N++;
	}
      }
      void AddOrReplaceVector(const multi1d<T>& v,const Subset& s){
	if(N<vec.size2())
	{
	  if (v.size() != Ls)
	  {
	    QDPIO::cerr << "VectorArrays:" << __func__ << ": size of 5D array inconsistent with stored value: Ls=" << Ls 
			<< "  v.size=" << v.size() << endl;
	    QDP_abort(1);
	  }
	  for(int k=0; k < Ls; ++k)
	    vec[N][k][s] = v[k];
	  N++;
	}
	else{// replace the last vector
	  for(int k=0; k < Ls; ++k)
	    vec[N-1][k][s] = v[k];
	}
      }

      // This will only add as many vectors as they fit
      void AddVectors(multi2d<T>& v,const Subset& s){
	if (v.size() != Ls)
	{
	  QDPIO::cerr << "VectorArrays:" << __func__ << ": size of 5D array inconsistent with stored value: Ls=" << Ls 
		      << "  v.size=" << v.size() << endl;
	  QDP_abort(1);
	}
	for(int i(0);i<v.size2();i++)
	  AddVector(v[i],s);
      }

    
      void resize(int n2, int n1) {N=Ls=0;vec.resize(n2,n1);} 
      int size() const { return vec.size2();}
      int Nvecs() const { return N; } 
      int N5d() const { return Ls; } 
      multi1d<T> operator[](int i){ return vec[i];}
    };


    //------------------------------------------------------------------------------
    //! Holds eigenvalues and arrays of eigenvectors for use in 5D work
    /*! \ingroup invert */
    template<class T> class RitzPairsArray
    {
    public:
      Vectors<Double>   eval;
      VectorArrays<T>   evec;
      int Neig;

      RitzPairsArray() {init(0,0);}
      RitzPairsArray(int N, int Ls) {init(N,Ls);}

      void init(int N, int Ls) {
	eval.resize(N);
	evec.resize(N,Ls);
	Neig = 0;
      }

      void AddVector(const Double& e, const multi1d<T>& v,const Subset& s){
	eval.AddVector(e,s);
	evec.AddVector(v,s);
	if (eval.N != evec.N)
	{
	  QDPIO::cerr << __func__ << ": length of value and vector arrays are not the same" << endl;
	  QDP_abort(1);
	}
	Neig = evec.N;
      }

      // This will only add as many vectors as they fit
      void AddVectors(const multi1d<Double>& e,const multi2d<T>& v,const Subset& s){
	for(int i(0);i<e.size();i++)
	  eval.AddVector(e[i],s);
	for(int i(0);i<v.size2();i++)
	  evec.AddVector(v[i],s);
	if (eval.N != evec.N)
	{
	  QDPIO::cerr << __func__ << ": length of value and vector arrays are not the same" << endl;
	  QDP_abort(1);
	}
	Neig = evec.N;
      }
    };

  } // namespace LinAlg

} // namespace Chroma

#endif 
