/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

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
#include "viennacl/generator/custom_operation.hpp"

#include "list"
//
// -------------------------------------------------------------
//
using namespace boost::numeric;
//
// -------------------------------------------------------------
//
static const unsigned int min_large_block_size = 32;
static const unsigned int max_large_block_size = 256;
static const unsigned int n_large_blocks = std::log(max_large_block_size/min_large_block_size)/std::log(2)+1;

static const unsigned int min_alignment = 1;
static const unsigned int max_alignment = 8;

static const unsigned int max_small_block_size = max_alignment;

//
// -------------------------------------------------------------

template <typename ScalarType>
ScalarType diff(ScalarType & s1, viennacl::scalar<ScalarType> & s2)
{
   viennacl::backend::finish();
   if (s1 != s2)
      return (s1 - s2) / std::max(fabs(s1), fabs(s2));
   return 0;
}

template <typename ScalarType, typename VCLMatrixType>
ScalarType diff(ublas::matrix<ScalarType> & mat1, VCLMatrixType & mat2)
{
   ublas::matrix<ScalarType> mat2_cpu(mat2.size1(), mat2.size2());
   viennacl::backend::finish();  //workaround for a bug in APP SDK 2.7 on Trinity APUs (with Catalyst 12.8)
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
          typename MatrixTypeA, typename MatrixTypeB, typename MatrixTypeC,
          typename Profiles>
int test_prod(Epsilon const& epsilon,

              ReferenceMatrixTypeA const & A, ReferenceMatrixTypeA const & A_trans,
              ReferenceMatrixTypeB const & B, ReferenceMatrixTypeB const & B_trans,
              ReferenceMatrixTypeC & C,

              MatrixTypeA const & vcl_A, MatrixTypeA const & vcl_A_trans,
              MatrixTypeB const & vcl_B, MatrixTypeB const & vcl_B_trans,
              MatrixTypeC & vcl_C,

              Profiles const & profiles
             )
{
   int retval = EXIT_SUCCESS;
   NumericT act_diff = 0;
   bool failed;



   // Test: C = A * B --------------------------------------------------------------------------
   C     = viennacl::linalg::prod(A, B);

   failed = false;
   for(typename Profiles::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){
       viennacl::generator::custom_operation op;
       if(it!=profiles.begin()) op.operations_manager().override_blas3_model(*it);
       op.add(vcl_C = viennacl::generator::prod(vcl_A,vcl_B));
       op.execute();
       viennacl::ocl::get_queue().finish();
       act_diff = fabs(diff(C, vcl_C.mat()));
       if( act_diff > epsilon )
       {
         std::cout << "# Error at operation: matrix-matrix product for profile " << *it << std::endl;
         std::cout << "  diff: " << act_diff << std::endl;
         retval = EXIT_FAILURE;
         failed = true;
       }
   }
   if(!failed) std::cout << "Test C = A * B passed!" << std::endl;

   // Test: C = trans(A) * B --------------------------------------------------------------------------
   C     = boost::numeric::ublas::prod(trans(A_trans), B);

   failed = false;
   for(typename Profiles::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){
//       std::cout << *it << std::endl;
       viennacl::generator::custom_operation op;
       op.operations_manager().override_blas3_model(*it);
       op.add(vcl_C = viennacl::generator::prod(trans(vcl_A_trans),vcl_B));
       op.execute();
       viennacl::ocl::get_queue().finish();
       act_diff = fabs(diff(C, vcl_C.mat()));
       if( act_diff > epsilon )
       {
         std::cout << "# Error at operation: matrix-matrix product for profile " << *it << std::endl;
         std::cout << "  diff: " << act_diff << std::endl;
         retval = EXIT_FAILURE;
         failed = true;
       }
   }
   if(!failed) std::cout << "Test C = trans(A) * B passed!" << std::endl;


   // Test: C = A * trans(B) --------------------------------------------------------------------------
   C     = boost::numeric::ublas::prod(A, trans(B_trans));

   failed = false;
   for(typename Profiles::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){

       viennacl::generator::custom_operation op;
       op.operations_manager().override_blas3_model(*it);
       op.add(vcl_C = viennacl::generator::prod(vcl_A,trans(vcl_B_trans)));
       op.execute();
       viennacl::ocl::get_queue().finish();
       act_diff = fabs(diff(C, vcl_C.mat()));
       if( act_diff > epsilon )
       {
         std::cout << "# Error at operation: matrix-matrix product for profile " << *it << std::endl;
         std::cout << "  diff: " << act_diff << std::endl;
         retval = EXIT_FAILURE;
         failed = true;
       }
   }
   if(!failed) std::cout << "Test C = A * trans(B) passed!" << std::endl;

   // Test: C = trans(A) * trans(B) --------------------------------------------------------------------------
   C     = boost::numeric::ublas::prod(trans(A_trans), trans(B_trans));

   failed = false;
   for(typename Profiles::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){
       viennacl::generator::custom_operation op;
       op.operations_manager().override_blas3_model(*it);
       op.add(vcl_C = viennacl::generator::prod(trans(vcl_A_trans),trans(vcl_B_trans)));
       op.execute();
       viennacl::ocl::get_queue().finish();
       act_diff = fabs(diff(C, vcl_C.mat()));
       if( act_diff > epsilon )
       {
         std::cout << "# Error at operation: matrix-matrix product for profile " << *it << std::endl;
         std::cout << "  diff: " << act_diff << std::endl;
         retval = EXIT_FAILURE;
         failed = true;
       }
   }
   if(!failed) std::cout << "Test C = trans(A) * trans(B) passed!" << std::endl;

   // Test: C = (A+A) * B --------------------------------------------------------------------------
   ReferenceMatrixTypeA  new_A(A+A);
   C     = viennacl::linalg::prod(new_A, B);

   failed = false;
   for(typename Profiles::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){
       viennacl::generator::custom_operation op;
       op.operations_manager().override_blas3_model(*it);
       op.add(vcl_C = viennacl::generator::prod(vcl_A+vcl_A,vcl_B));
       op.execute();
       viennacl::ocl::get_queue().finish();
       act_diff = fabs(diff(C, vcl_C.mat()));
       if( act_diff > epsilon )
       {
         std::cout << "# Error at operation: matrix-matrix product for profile " << *it << std::endl;
         std::cout << "  diff: " << act_diff << std::endl;
         retval = EXIT_FAILURE;
         failed = true;
       }
   }
   if(!failed) std::cout << "Test C = (A+A) * B passed!" << std::endl;


   // Test: C = A * (B+B) --------------------------------------------------------------------------
   ReferenceMatrixTypeA  new_B(B+B);
   C     = viennacl::linalg::prod(A, new_B);

   failed = false;
   for(typename Profiles::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){
       viennacl::generator::custom_operation op;
       op.operations_manager().override_blas3_model(*it);
       op.add(vcl_C = viennacl::generator::prod(vcl_A,vcl_B + vcl_B));
       op.execute();
       viennacl::ocl::get_queue().finish();
       act_diff = fabs(diff(C, vcl_C.mat()));
       if( act_diff > epsilon )
       {
         std::cout << "# Error at operation: matrix-matrix product for profile " << *it << std::endl;
         std::cout << "  diff: " << act_diff << std::endl;
         retval = EXIT_FAILURE;
         failed = true;
       }
   }
   if(!failed) std::cout << "Test C = A * (B+B) passed!" << std::endl;

   return retval;
}

template< typename NumericT, typename F_A, typename F_B, typename F_C, typename Epsilon, typename Profiles>
int test_prod(Epsilon const& epsilon, Profiles const & profiles)
{
  int ret;

  long matrix_size1 = 5*max_large_block_size;
  long matrix_size2 = 2*max_large_block_size;
  long matrix_size3 = 3*max_large_block_size;

  typedef viennacl::generator::dummy_matrix< viennacl::matrix<NumericT, F_A> > DummyMatrixA;
  typedef viennacl::generator::dummy_matrix< viennacl::matrix<NumericT, F_B> > DummyMatrixB;
  typedef viennacl::generator::dummy_matrix< viennacl::matrix<NumericT, F_C> > DummyMatrixC;

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
  viennacl::matrix<NumericT, F_A>    vcl_A(matrix_size1, matrix_size2);
  viennacl::copy(A, vcl_A);
  DummyMatrixA dummy_vcl_A(vcl_A);

  // A^T
  viennacl::matrix<NumericT, F_A>    vcl_A_trans(matrix_size2, matrix_size1);
  viennacl::copy(A_trans, vcl_A_trans);
  DummyMatrixA dummy_vcl_A_trans(vcl_A_trans);


  // B
  viennacl::matrix<NumericT, F_B>    vcl_B(matrix_size2, matrix_size3);
  viennacl::copy(B, vcl_B);
  DummyMatrixB dummy_vcl_B(vcl_B);

  // B^T
  viennacl::matrix<NumericT, F_B>    vcl_B_trans(matrix_size3, matrix_size2);
  viennacl::copy(B_trans, vcl_B_trans);
  DummyMatrixB dummy_vcl_B_trans(vcl_B_trans);

  // C
  viennacl::matrix<NumericT, F_C>    vcl_C(matrix_size1, matrix_size3);
  DummyMatrixC dummy_vcl_C(vcl_C);

  std::cout << "--- Part 1: Testing matrix-matrix products ---" << std::endl;

  //////
  //////  A: matrix
  //////

  //
  //
  std::cout << "Now using A=matrix, B=matrix, C=matrix" << std::endl;
  ret = test_prod<NumericT>(epsilon,
                            A, A_trans, B, B_trans, C,
                            dummy_vcl_A, dummy_vcl_A_trans,
                            dummy_vcl_B, dummy_vcl_B_trans,
                            dummy_vcl_C,
                            profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  return EXIT_SUCCESS;
}

template< typename NumericT, typename Epsilon >
int test(Epsilon const& epsilon, unsigned int n_profiles)
{
  int ret;

  //Creates profiles
  std::list<viennacl::generator::code_generation::blas3_optimization_profile> profiles;
  while(profiles.size() < n_profiles){
    for(unsigned int alignment = min_alignment ; alignment <= max_alignment ; alignment*=2){
          unsigned int ml = pow(2,rand()%n_large_blocks)*min_large_block_size;
          unsigned int kl = pow(2,rand()%n_large_blocks)*min_large_block_size;
          unsigned int nl = pow(2,rand()%n_large_blocks)*min_large_block_size;
          unsigned int n_small_blocks = std::log(max_small_block_size/alignment)/std::log(2)+1;
          unsigned int ms = pow(2,rand()%n_small_blocks)*alignment;
          unsigned int ks = pow(2,rand()%n_small_blocks)*alignment;
          unsigned int ns = pow(2,rand()%n_small_blocks)*alignment;
          bool use_LHS_shared = rand()%2;
          bool use_RHS_shared = rand()%2;
          viennacl::generator::code_generation::blas3_optimization_profile prof(ml,kl,nl,
                                                                                     ms,ks,ns,
                                                                                     use_LHS_shared, use_RHS_shared,
                                                                                     alignment);

          double lmem_size = 0;
          if(use_LHS_shared) lmem_size += (double)(kl+1)*(ml+1)*sizeof(NumericT)/1024;
          if(use_RHS_shared) lmem_size += (double)(nl+1)*(kl+1)*sizeof(NumericT)/1024;
          if( lmem_size > 32.0) continue;
          if(prof.local_work_size().first*prof.local_work_size().second > viennacl::ocl::current_device().max_workgroup_size()) continue;

          profiles.push_back(prof);

      }
  }

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=row, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::row_major, viennacl::row_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=row, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::row_major, viennacl::column_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=col, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::column_major, viennacl::row_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=row, B=col, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::row_major, viennacl::column_major, viennacl::column_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=row, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::row_major, viennacl::row_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=row, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::row_major, viennacl::column_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=col, C=row ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::column_major, viennacl::row_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "/// Now testing A=col, B=col, C=col ///" << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  ret = test_prod<NumericT, viennacl::column_major, viennacl::column_major, viennacl::column_major>(epsilon,profiles);
  if (ret != EXIT_SUCCESS)
    return ret;



  return ret;
}

int main()
{
   std::cout << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << "## Test :: Generated BLAS 3 routines" << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << std::endl;

   int retval = EXIT_SUCCESS;

   srand(time(NULL));
   unsigned int n_profiles = 10;

   std::cout << std::endl;
   std::cout << "----------------------------------------------" << std::endl;
   std::cout << std::endl;
   {
      typedef float NumericT;
      NumericT epsilon = NumericT(1.0E-3);
      std::cout << "# Testing setup:" << std::endl;

      std::cout << viennacl::ocl::current_device().info() << std::endl;

      std::cout << "  n_profiles: " << n_profiles << std::endl;
      std::cout << "  eps:     " << epsilon << std::endl;
      std::cout << "  numeric: float" << std::endl;
      retval = test<NumericT>(epsilon,n_profiles);
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
        std::cout << "  n_profiles: " << n_profiles << std::endl;
        std::cout << "  eps:     " << epsilon << std::endl;
        std::cout << "  numeric: double" << std::endl;
        retval = test<NumericT>(epsilon,n_profiles);
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

