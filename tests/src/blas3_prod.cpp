/* =========================================================================
   Copyright (c) 2010-2011, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

//#define NDEBUG
//#define VIENNACL_DEBUG_BUILD

//
// *** System
//
#include <iostream>

//
// *** Boost
//
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

//
// *** ViennaCL
//
//#define VIENNACL_DEBUG_ALL
//#define VIENNACL_DEBUG_BUILD
#define VIENNACL_HAVE_UBLAS 1
#include "viennacl/scalar.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/direct_solve.hpp"
#include "examples/tutorial/Random.hpp"
//
// -------------------------------------------------------------
//
using namespace boost::numeric;
//
// -------------------------------------------------------------
//
template <typename ScalarType>
ScalarType diff(ScalarType & s1, viennacl::scalar<ScalarType> & s2) 
{
   if (s1 != s2)
      return (s1 - s2) / std::max(fabs(s1), fabs(s2));
   return 0;
}

template <typename ScalarType>
ScalarType diff(ublas::vector<ScalarType> & v1, viennacl::vector<ScalarType> & v2)
{
   ublas::vector<ScalarType> v2_cpu(v2.size());
   viennacl::copy(v2.begin(), v2.end(), v2_cpu.begin());

   for (std::size_t i=0;i<v1.size(); ++i)
   {
      if ( std::max( fabs(v2_cpu[i]), fabs(v1[i]) ) > 0 )
         v2_cpu[i] = fabs(v2_cpu[i] - v1[i]) / std::max( fabs(v2_cpu[i]), fabs(v1[i]) );
      else
         v2_cpu[i] = 0.0;
   }

   return norm_inf(v2_cpu);
}


template <typename ScalarType, typename VCLMatrixType>
ScalarType diff(ublas::matrix<ScalarType> & mat1, VCLMatrixType & mat2)
{
   ublas::matrix<ScalarType> mat2_cpu(mat2.size1(), mat2.size2());
   viennacl::copy(mat2, mat2_cpu);
   double ret = 0;
   double act = 0;

    for (unsigned int i = 0; i < mat2_cpu.size1(); ++i)
    {
      for (unsigned int j = 0; j < mat2_cpu.size2(); ++j)
      {
         act = fabs(mat2_cpu(i,j) - mat1(i,j)) / std::max( fabs(mat2_cpu(i, j)), fabs(mat1(i,j)) );
         if (act > ret)
           ret = act;
      }
    }
   //std::cout << ret << std::endl;
   return ret;
}






//
// Part 1: Matrix-matrix multiplications
//


template< typename NumericT, typename Epsilon,
          typename ReferenceMatrixTypeA, typename ReferenceMatrixTypeB, typename ReferenceMatrixTypeC,
          typename MatrixTypeA, typename MatrixTypeB, typename MatrixTypeC>
int test_prod(Epsilon const& epsilon,
              
              ReferenceMatrixTypeA const & A, ReferenceMatrixTypeA const & A_trans,
              ReferenceMatrixTypeB const & B, ReferenceMatrixTypeB const & B_trans,
              ReferenceMatrixTypeC & C,
              
              MatrixTypeA const & vcl_A, MatrixTypeA const & vcl_A_trans,
              MatrixTypeB const & vcl_B, MatrixTypeB const & vcl_B_trans,
              MatrixTypeC & vcl_C
             )
{
   int retval = EXIT_SUCCESS;
   NumericT act_diff = 0;
   

   // Test: C +-= A * B --------------------------------------------------------------------------       
   C     = viennacl::linalg::prod(A, B);
   vcl_C = viennacl::linalg::prod(vcl_A, vcl_B);
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C = A * B passed!" << std::endl;
   
   
   C     += viennacl::linalg::prod(A, B);
   vcl_C += viennacl::linalg::prod(vcl_A, vcl_B);
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C += A * B passed!" << std::endl;

   C     -= viennacl::linalg::prod(A, B);
   vcl_C -= viennacl::linalg::prod(vcl_A, vcl_B);
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C -= A * B passed!" << std::endl;
   
   
   
   
   
   // Test: C +-= A * trans(B) --------------------------------------------------------------------------       
   C     = boost::numeric::ublas::prod(A, trans(B_trans));
   vcl_C = viennacl::linalg::prod(vcl_A, trans(vcl_B_trans));
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C = A * trans(B) passed!" << std::endl;


   C     += boost::numeric::ublas::prod(A, trans(B_trans));
   vcl_C += viennacl::linalg::prod(vcl_A, trans(vcl_B_trans));
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C += A * trans(B) passed!" << std::endl;

   
   C     -= boost::numeric::ublas::prod(A, trans(B_trans));
   vcl_C -= viennacl::linalg::prod(vcl_A, trans(vcl_B_trans));
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C -= A * trans(B) passed!" << std::endl;
   
   
   
   // Test: C +-= trans(A) * B --------------------------------------------------------------------------       
   C     = boost::numeric::ublas::prod(trans(A_trans), B);
   vcl_C = viennacl::linalg::prod(trans(vcl_A_trans), vcl_B);
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C = trans(A) * B passed!" << std::endl;
   
   
   C     += boost::numeric::ublas::prod(trans(A_trans), B);
   vcl_C += viennacl::linalg::prod(trans(vcl_A_trans), vcl_B);
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C += trans(A) * B passed!" << std::endl;
   
   
   C     -= boost::numeric::ublas::prod(trans(A_trans), B);
   vcl_C -= viennacl::linalg::prod(trans(vcl_A_trans), vcl_B);
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C -= trans(A) * B passed!" << std::endl;
   
   

   
   
   // Test: C +-= trans(A) * trans(B) --------------------------------------------------------------------------       
   C     = boost::numeric::ublas::prod(trans(A_trans), trans(B_trans));
   vcl_C = viennacl::linalg::prod(trans(vcl_A_trans), trans(vcl_B_trans));
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C = trans(A) * trans(B) passed!" << std::endl;
   
   C     += boost::numeric::ublas::prod(trans(A_trans), trans(B_trans));
   vcl_C += viennacl::linalg::prod(trans(vcl_A_trans), trans(vcl_B_trans));
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C += trans(A) * trans(B) passed!" << std::endl;
   
   
   C     -= boost::numeric::ublas::prod(trans(A_trans), trans(B_trans));
   vcl_C -= viennacl::linalg::prod(trans(vcl_A_trans), trans(vcl_B_trans));
   act_diff = fabs(diff(C, vcl_C));
   
   if( act_diff > epsilon )
   {
     std::cout << "# Error at operation: matrix-matrix product" << std::endl;
     std::cout << "  diff: " << act_diff << std::endl;
     retval = EXIT_FAILURE;
   }
   else
     std::cout << "Test C -= trans(A) * trans(B) passed!" << std::endl;
   
   
   
   
   return retval;
}



template< typename NumericT, typename F_A, typename F_B, typename F_C, typename Epsilon >
int test_prod(Epsilon const& epsilon)
{
  int ret;
  
  long matrix_size1 = 131;  //some odd number, not too large
  long matrix_size2 = 67;  //some odd number, not too large
  long matrix_size3 = 73;  //some odd number, not too large
  //long matrix_size1 = 128;  //some odd number, not too large
  //long matrix_size2 = 64;  //some odd number, not too large
  //long matrix_size3 = 128;  //some odd number, not too large
  
  // --------------------------------------------------------------------------       
  
  // ublas reference:
  ublas::matrix<NumericT> A(matrix_size1, matrix_size2);
  ublas::matrix<NumericT> big_A = ublas::scalar_matrix<NumericT>(4*matrix_size1, 4*matrix_size2, 3.1415);
  
  ublas::matrix<NumericT> B(matrix_size2, matrix_size3);
  ublas::matrix<NumericT> big_B = ublas::scalar_matrix<NumericT>(4*matrix_size2, 4*matrix_size3, 42.0);
  
  ublas::matrix<NumericT> C(matrix_size1, matrix_size3);

  //fill A and B:
  for (unsigned int i = 0; i < A.size1(); ++i)
    for (unsigned int j = 0; j < A.size2(); ++j)
        A(i,j) = static_cast<NumericT>(0.1) * random<NumericT>();
  for (unsigned int i = 0; i < B.size1(); ++i)
    for (unsigned int j = 0; j < B.size2(); ++j)
        B(i,j) = static_cast<NumericT>(0.1) * random<NumericT>();

  ublas::matrix<NumericT>     A_trans = trans(A);
  ublas::matrix<NumericT> big_A_trans = trans(big_A);
  
  ublas::matrix<NumericT>     B_trans = trans(B);
  ublas::matrix<NumericT> big_B_trans = trans(big_B);
  
  //
  // ViennaCL objects
  //
  
  // A
  viennacl::range range1_A(matrix_size1, 2*matrix_size1);
  viennacl::range range2_A(matrix_size2, 2*matrix_size2);
  viennacl::slice slice1_A(matrix_size1, 2, matrix_size1);
  viennacl::slice slice2_A(matrix_size2, 3, matrix_size2);
  
  viennacl::matrix<NumericT, F_A>    vcl_A(matrix_size1, matrix_size2);
  viennacl::copy(A, vcl_A);
  
  viennacl::matrix<NumericT, F_A>    vcl_big_range_A(4*matrix_size1, 4*matrix_size2);
  viennacl::matrix_range<viennacl::matrix<NumericT, F_A> > vcl_range_A(vcl_big_range_A, range1_A, range2_A);
  viennacl::copy(A, vcl_range_A);
  
  viennacl::matrix<NumericT, F_A>    vcl_big_slice_A(4*matrix_size1, 4*matrix_size2);
  viennacl::matrix_slice<viennacl::matrix<NumericT, F_A> > vcl_slice_A(vcl_big_slice_A, slice1_A, slice2_A);
  viennacl::copy(A, vcl_slice_A);

  
  // A^T
  viennacl::matrix<NumericT, F_A>    vcl_A_trans(matrix_size2, matrix_size1);
  viennacl::copy(A_trans, vcl_A_trans);

  viennacl::matrix<NumericT, F_A>    vcl_big_range_A_trans(4*matrix_size2, 4*matrix_size1);
  viennacl::matrix_range<viennacl::matrix<NumericT, F_A> > vcl_range_A_trans(vcl_big_range_A_trans, range2_A, range1_A);
  viennacl::copy(A_trans, vcl_range_A_trans);
  
  viennacl::matrix<NumericT, F_A>    vcl_big_slice_A_trans(4*matrix_size2, 4*matrix_size1);
  viennacl::matrix_slice<viennacl::matrix<NumericT, F_A> > vcl_slice_A_trans(vcl_big_slice_A_trans, slice2_A, slice1_A);
  viennacl::copy(A_trans, vcl_slice_A_trans);
  
  
  
  // B
  viennacl::range range1_B(2*matrix_size2, 3*matrix_size2);
  viennacl::range range2_B(2*matrix_size3, 3*matrix_size3);
  viennacl::slice slice1_B(matrix_size2, 3, matrix_size2);
  viennacl::slice slice2_B(matrix_size3, 2, matrix_size3);
  
  viennacl::matrix<NumericT, F_B>    vcl_B(matrix_size2, matrix_size3);
  viennacl::copy(B, vcl_B);
  
  viennacl::matrix<NumericT, F_B>    vcl_big_range_B(4*matrix_size2, 4*matrix_size3);
  viennacl::matrix_range<viennacl::matrix<NumericT, F_B> > vcl_range_B(vcl_big_range_B, range1_B, range2_B);
  viennacl::copy(B, vcl_range_B);
  
  viennacl::matrix<NumericT, F_B>    vcl_big_slice_B(4*matrix_size2, 4*matrix_size3);
  viennacl::matrix_slice<viennacl::matrix<NumericT, F_B> > vcl_slice_B(vcl_big_slice_B, slice1_B, slice2_B);
  viennacl::copy(B, vcl_slice_B);

  
  // B^T
  
  viennacl::matrix<NumericT, F_B>    vcl_B_trans(matrix_size3, matrix_size2);
  viennacl::copy(B_trans, vcl_B_trans);
  
  viennacl::matrix<NumericT, F_B>    vcl_big_range_B_trans(4*matrix_size3, 4*matrix_size2);
  viennacl::matrix_range<viennacl::matrix<NumericT, F_B> > vcl_range_B_trans(vcl_big_range_B_trans, range2_B, range1_B);
  viennacl::copy(B_trans, vcl_range_B_trans);
  
  viennacl::matrix<NumericT, F_B>    vcl_big_slice_B_trans(4*matrix_size3, 4*matrix_size2);
  viennacl::matrix_slice<viennacl::matrix<NumericT, F_B> > vcl_slice_B_trans(vcl_big_slice_B_trans, slice2_B, slice1_B);
  viennacl::copy(B_trans, vcl_slice_B_trans);
  

  // C
  
  viennacl::range range1_C(matrix_size1-1, 2*matrix_size1-1);
  viennacl::range range2_C(matrix_size3-1, 2*matrix_size3-1);
  viennacl::slice slice1_C(matrix_size1-1, 3, matrix_size1);
  viennacl::slice slice2_C(matrix_size3-1, 3, matrix_size3);
  
  viennacl::matrix<NumericT, F_C>    vcl_C(matrix_size1, matrix_size3);
  
  viennacl::matrix<NumericT, F_C>    vcl_big_range_C(4*matrix_size1, 4*matrix_size3);
  viennacl::matrix_range<viennacl::matrix<NumericT, F_C> > vcl_range_C(vcl_big_range_C, range1_C, range2_C);
  
  viennacl::matrix<NumericT, F_C>    vcl_big_slice_C(4*matrix_size1, 4*matrix_size3);
  viennacl::matrix_slice<viennacl::matrix<NumericT, F_C> > vcl_slice_C(vcl_big_slice_C, slice1_C, slice2_C);

  
  std::cout << "--- Part 1: Testing matrix-matrix products ---" << std::endl;
  
  //////
  //////  A: matrix
  //////
  
  //
  //
  std::cout << "Now using A=matrix, B=matrix, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=matrix, B=matrix, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=matrix, B=matrix, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  

  //
  //
  std::cout << "Now using A=matrix, B=range, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=matrix, B=range, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=matrix, B=range, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  
  //
  //
  std::cout << "Now using A=matrix, B=slice, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=matrix, B=slice, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=matrix, B=slice, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_A, vcl_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  
  //////
  //////  A: range
  //////
  
  //
  //
  std::cout << "Now using A=range, B=matrix, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=range, B=matrix, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=range, B=matrix, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  

  //
  //
  std::cout << "Now using A=range, B=range, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=range, B=range, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=range, B=range, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  
  //
  //
  std::cout << "Now using A=range, B=slice, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=range, B=slice, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=range, B=slice, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_range_A, vcl_range_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;
  
  
  
  //////
  //////  A: slice
  //////
  
  //
  //
  std::cout << "Now using A=slice, B=matrix, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=slice, B=matrix, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=slice, B=matrix, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_B, vcl_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  

  //
  //
  std::cout << "Now using A=slice, B=range, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=slice, B=range, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=slice, B=range, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_range_B, vcl_range_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  
  //
  //
  std::cout << "Now using A=slice, B=slice, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_C);
  if (ret != EXIT_SUCCESS)
    return ret;


  //
  //
  std::cout << "Now using A=slice, B=slice, C=range" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_range_C);
  if (ret != EXIT_SUCCESS)
    return ret;

  //
  //
  std::cout << "Now using A=slice, B=slice, C=slice" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            vcl_slice_A, vcl_slice_A_trans,
                            vcl_slice_B, vcl_slice_B_trans,
                            vcl_slice_C);
  if (ret != EXIT_SUCCESS)
    return ret;
  
  
  return ret;
  
}


//
// Control functions
//



template< typename NumericT, typename Epsilon >
int test(Epsilon const& epsilon)
{
  int ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=row, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::row_major, viennacl::row_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;
  
  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=row, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::row_major, viennacl::column_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=col, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::column_major, viennacl::row_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=col, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::column_major, viennacl::column_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=row, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::row_major, viennacl::row_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;
  
  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=row, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::row_major, viennacl::column_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=col, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::column_major, viennacl::row_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=col, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::column_major, viennacl::column_major>(epsilon);
  if (ret != EXIT_SUCCESS)
    return ret;

  
  
  return ret;
}

//
// -------------------------------------------------------------
//
int main()
{
   std::cout << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << "## Test :: BLAS 3 routines" << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << std::endl;

   int retval = EXIT_SUCCESS;

   std::cout << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << std::endl;
   {
      typedef float NumericT;
      NumericT epsilon = NumericT(1.0E-3);
      std::cout << "# Testing setup:" << std::endl;
      std::cout << "  eps:     " << epsilon << std::endl;
      std::cout << "  numeric: float" << std::endl;
      retval = test<NumericT>(epsilon);
      if( retval == EXIT_SUCCESS )
        std::cout << "# Test passed" << std::endl;
      else
        return retval;
   }
   std::cout << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << std::endl;
#ifdef VIENNACL_HAVE_OPENCL   
   if( viennacl::ocl::current_device().double_support() )
#endif
   {
      {
        typedef double NumericT;
        NumericT epsilon = 1.0E-11;
        std::cout << "# Testing setup:" << std::endl;
        std::cout << "  eps:     " << epsilon << std::endl;
        std::cout << "  numeric: double" << std::endl;
        retval = test<NumericT>(epsilon);
        if( retval == EXIT_SUCCESS )
          std::cout << "# Test passed" << std::endl;
        else
          return retval;
      }
      std::cout << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << std::endl;
   }
   
   std::cout << std::endl;
   std::cout << "------- Test completed --------" << std::endl;
   std::cout << std::endl;
   
   
   return retval;
}
