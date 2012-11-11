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

#define VIENNACL_HAVE_UBLAS
//#define NDEBUG
//#define VIENNACL_BUILD_INFO

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <time.h>
//#include "../benchmarks/benchmark-utils.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"
/*#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/io/matrix_market.hpp"*/
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/vector_proxy.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/io.hpp"


template <typename MatrixType, typename VCLMatrixType>
bool check_for_equality(MatrixType const & ublas_A, VCLMatrixType const & vcl_A, double epsilon)
{
  typedef typename MatrixType::value_type   value_type;
  
  boost::numeric::ublas::matrix<value_type> vcl_A_cpu(vcl_A.size1(), vcl_A.size2());
  viennacl::copy(vcl_A, vcl_A_cpu);
  
  for (std::size_t i=0; i<ublas_A.size1(); ++i)
  {
    for (std::size_t j=0; j<ublas_A.size2(); ++j)
    {
      if (ublas_A(i,j) != vcl_A_cpu(i,j))
      {
        if ( (std::abs(ublas_A(i,j) - vcl_A_cpu(i,j)) / std::max(ublas_A(i,j), vcl_A_cpu(i,j)) > epsilon) || (vcl_A_cpu(i,j) != vcl_A_cpu(i,j)) ) 
        {
          std::cout << "Error at index (" << i << ", " << j << "): " << ublas_A(i,j) << " vs " << vcl_A_cpu(i,j) << std::endl;
          std::cout << std::endl << "TEST failed!" << std::endl;
          return false;
        }
      }
    }
  }
  
  std::cout << "PASSED!" << std::endl;
  return true;
}



           
template <typename UBLASMatrixType,
          typename ViennaCLMatrixType1, typename ViennaCLMatrixType2, typename ViennaCLMatrixType3>
int run_test(double epsilon,
             UBLASMatrixType & ublas_A, UBLASMatrixType & ublas_B, UBLASMatrixType & ublas_C,
             ViennaCLMatrixType1 & vcl_A, ViennaCLMatrixType2 & vcl_B, ViennaCLMatrixType3 vcl_C)
{

  typedef typename viennacl::result_of::cpu_value_type<typename ViennaCLMatrixType1::value_type>::type  cpu_value_type;
  
  cpu_value_type alpha = 3.1415;
  viennacl::scalar<cpu_value_type>   gpu_alpha = alpha;

  cpu_value_type beta = 2.7182;
  viennacl::scalar<cpu_value_type>   gpu_beta = beta;
  
  
  
  std::cout << std::endl;
  //std::cout << "//" << std::endl;
  //std::cout << "////////// Test: Assignments //////////" << std::endl;
  //std::cout << "//" << std::endl;

  std::cout << "Testing matrix assigned to range... ";
  ublas_A = ublas_B;
  vcl_A = vcl_B;
  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;



  //std::cout << std::endl;
  //std::cout << "//" << std::endl;
  //std::cout << "////////// Test 1: Copy to GPU //////////" << std::endl;
  //std::cout << "//" << std::endl;

  ublas_A = ublas_B;
  viennacl::copy(ublas_B, vcl_A);
  std::cout << "Testing upper left copy to GPU... ";
  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;


  ublas_C = ublas_B;
  viennacl::copy(ublas_B, vcl_C);
  std::cout << "Testing lower right copy to GPU... ";
  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;


  //std::cout << std::endl;
  //std::cout << "//" << std::endl;
  //std::cout << "////////// Test 2: Copy from GPU //////////" << std::endl;
  //std::cout << "//" << std::endl;

  std::cout << "Testing upper left copy to A... ";
  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;

  std::cout << "Testing lower right copy to C... ";
  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  
  
  //std::cout << "//" << std::endl;
  //std::cout << "////////// Test 3: Addition //////////" << std::endl;
  //std::cout << "//" << std::endl;
  viennacl::copy(ublas_C, vcl_C);

  std::cout << "Inplace add: ";
  ublas_C += ublas_C;
  vcl_C += vcl_C;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  std::cout << "Scaled inplace add: ";
  ublas_C += beta * ublas_A;
  vcl_C += gpu_beta * vcl_A;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;
  
  std::cout << "Add: ";
  ublas_C = ublas_A + ublas_B;
  vcl_C = vcl_A + vcl_B;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  std::cout << "Scaled add: ";
  ublas_C = alpha * ublas_A + beta * ublas_B;
  vcl_C = alpha * vcl_A + beta * vcl_B;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  std::cout << "Scaled add: ";
  vcl_C = gpu_alpha * vcl_A + gpu_beta * vcl_B;
  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  
  
  //std::cout << "//" << std::endl;
  //std::cout << "////////// Test 4: Subtraction //////////" << std::endl;
  //std::cout << "//" << std::endl;
  viennacl::copy(ublas_C, vcl_C);

  std::cout << "Inplace sub: ";
  ublas_C -= ublas_B;
  vcl_C -= vcl_B;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  std::cout << "Scaled Inplace sub: ";
  ublas_C -= alpha * ublas_B;
  vcl_C -= alpha * vcl_B;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  
  
  
  std::cout << "Sub: ";
  ublas_C = ublas_A - ublas_B;
  vcl_C = vcl_A - vcl_B;

  if (!check_for_equality(ublas_C, vcl_C, epsilon))
    return EXIT_FAILURE;

  std::cout << "Scaled sub: ";
  ublas_B = alpha * ublas_A - beta * ublas_C;
  vcl_B = alpha * vcl_A - vcl_C * beta;

  if (!check_for_equality(ublas_B, vcl_B, epsilon))
    return EXIT_FAILURE;

  std::cout << "Scaled sub: ";
  vcl_B = gpu_alpha * vcl_A - vcl_C * gpu_beta;
  if (!check_for_equality(ublas_B, vcl_B, epsilon))
    return EXIT_FAILURE;

  
  
  
  //std::cout << "//" << std::endl;
  //std::cout << "////////// Test 5: Scaling //////////" << std::endl;
  //std::cout << "//" << std::endl;
  viennacl::copy(ublas_A, vcl_A);

  std::cout << "Multiplication with CPU scalar: ";
  ublas_A *= alpha;
  vcl_A   *= alpha;

  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;

  std::cout << "Multiplication with GPU scalar: ";
  ublas_A *= beta;
  vcl_A *= gpu_beta;

  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;


  std::cout << "Division with CPU scalar: ";
  ublas_A /= alpha;
  vcl_A /= alpha;

  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;

  std::cout << "Division with GPU scalar: ";
  ublas_A /= beta;
  vcl_A /= gpu_beta;

  if (!check_for_equality(ublas_A, vcl_A, epsilon))
    return EXIT_FAILURE;

  std::cout << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << std::endl;


  return EXIT_SUCCESS;
}    



           
template <typename T, typename ScalarType>
int run_test(double epsilon)
{
    //typedef float               ScalarType;
    typedef boost::numeric::ublas::matrix<ScalarType>       MatrixType;
    typedef boost::numeric::ublas::vector<ScalarType>       VectorType;
    
    typedef viennacl::matrix<ScalarType, T>    VCLMatrixType;
    
    //std::size_t dim_rows = 151;
    //std::size_t dim_cols = 37;
    std::size_t dim_rows = 5;
    std::size_t dim_cols = 3;
    
    //setup ublas objects:
    MatrixType ublas_A(dim_rows, dim_cols);
    MatrixType ublas_B(dim_rows, dim_cols);
    MatrixType ublas_C(dim_rows, dim_cols);
    
    for (std::size_t i=0; i<ublas_A.size1(); ++i)
      for (std::size_t j=0; j<ublas_A.size2(); ++j)
      {
        ublas_A(i,j) = ScalarType((i+2) + (j+1)*(i+2));
        ublas_B(i,j) = ScalarType((j+2) + (j+1)*(j+2));
        ublas_C(i,j) = ScalarType((i+1) + (i+1)*(i+2));
      }

    MatrixType ublas_A_large(4 * dim_rows, 4 * dim_cols);
    for (std::size_t i=0; i<ublas_A_large.size1(); ++i)
      for (std::size_t j=0; j<ublas_A_large.size2(); ++j)
        ublas_A_large(i,j) = ScalarType(i * ublas_A_large.size2() + j);

    //Setup ViennaCL objects    
    VCLMatrixType vcl_A_full(4 * dim_rows, 4 * dim_cols);
    VCLMatrixType vcl_B_full(4 * dim_rows, 4 * dim_cols);
    VCLMatrixType vcl_C_full(4 * dim_rows, 4 * dim_cols);
    
    viennacl::copy(ublas_A_large, vcl_A_full);
    viennacl::copy(ublas_A_large, vcl_B_full);
    viennacl::copy(ublas_A_large, vcl_C_full);

    //
    // Create A
    //
    VCLMatrixType vcl_A(dim_rows, dim_cols);
    
    viennacl::range vcl_A_r1(2 * dim_rows, 3 * dim_rows);
    viennacl::range vcl_A_r2(dim_cols, 2 * dim_cols);
    viennacl::matrix_range<VCLMatrixType>   vcl_range_A(vcl_A_full, vcl_A_r1, vcl_A_r2);
    
    viennacl::slice vcl_A_s1(2, 3, dim_rows);
    viennacl::slice vcl_A_s2(2 * dim_cols, 2, dim_cols);
    viennacl::matrix_slice<VCLMatrixType>   vcl_slice_A(vcl_A_full, vcl_A_s1, vcl_A_s2);
    
    
    //
    // Create B
    //
    VCLMatrixType vcl_B(dim_rows, dim_cols);
    
    viennacl::range vcl_B_r1(dim_rows, 2 * dim_rows);
    viennacl::range vcl_B_r2(2 * dim_cols, 3 * dim_cols);
    viennacl::matrix_range<VCLMatrixType>   vcl_range_B(vcl_B_full, vcl_B_r1, vcl_B_r2);
    
    viennacl::slice vcl_B_s1(2 * dim_rows, 2, dim_rows);
    viennacl::slice vcl_B_s2(dim_cols, 3, dim_cols);
    viennacl::matrix_slice<VCLMatrixType>   vcl_slice_B(vcl_B_full, vcl_B_s1, vcl_B_s2);
    
    
    //
    // Create C
    //
    VCLMatrixType vcl_C(dim_rows, dim_cols);
    
    viennacl::range vcl_C_r1(2 * dim_rows, 3 * dim_rows);
    viennacl::range vcl_C_r2(3 * dim_cols, 4 * dim_cols);
    viennacl::matrix_range<VCLMatrixType>   vcl_range_C(vcl_C_full, vcl_C_r1, vcl_C_r2);
    
    viennacl::slice vcl_C_s1(dim_rows, 2, dim_rows);
    viennacl::slice vcl_C_s2(0, 3, dim_cols);
    viennacl::matrix_slice<VCLMatrixType>   vcl_slice_C(vcl_C_full, vcl_C_s1, vcl_C_s2);
    
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_A, vcl_slice_A);

    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_B, vcl_slice_B);

    viennacl::copy(ublas_C, vcl_C);
    viennacl::copy(ublas_C, vcl_range_C);
    viennacl::copy(ublas_C, vcl_slice_C);


    std::cout << std::endl;
    std::cout << "//" << std::endl;
    std::cout << "////////// Test: Copy CTOR //////////" << std::endl;
    std::cout << "//" << std::endl;

    {
      std::cout << "Testing matrix created from range... ";
      VCLMatrixType vcl_temp = vcl_range_A;
      if (check_for_equality(ublas_A, vcl_temp, epsilon))
        std::cout << "PASSED!" << std::endl;
      else
      {
        std::cout << "ublas_A: " << ublas_A << std::endl;
        std::cout << "vcl_temp: " << vcl_temp << std::endl;
        std::cout << "vcl_range_A: " << vcl_range_A << std::endl;
        std::cout << "vcl_A: " << vcl_A << std::endl;
        std::cout << std::endl << "TEST failed!" << std::endl;
        return EXIT_FAILURE;
      }

      std::cout << "Testing matrix created from slice... ";
      VCLMatrixType vcl_temp2 = vcl_range_B;
      if (check_for_equality(ublas_B, vcl_temp2, epsilon))
        std::cout << "PASSED!" << std::endl;
      else
      {
        std::cout << std::endl << "TEST failed!" << std::endl;
        return EXIT_FAILURE;
      }
    }
    
    //
    // run operation tests:
    //
    
    /////// A=matrix:
    std::cout << "Testing A=matrix, B=matrix, C=matrix ..." << std::endl;
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=matrix, B=matrix, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=matrix, B=matrix, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=matrix, B=range, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_range_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=matrix, B=range, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_range_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=matrix, B=range, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_range_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    

    std::cout << "Testing A=matrix, B=slice, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_slice_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=matrix, B=slice, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_slice_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=matrix, B=slice, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_A, vcl_slice_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    


    /////// A=range:
    std::cout << "Testing A=range, B=matrix, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=range, B=matrix, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=range, B=matrix, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    

    std::cout << "Testing A=range, B=range, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_range_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=range, B=range, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_range_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=range, B=range, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_range_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    

    std::cout << "Testing A=range, B=slice, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_slice_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=range, B=slice, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_slice_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=range, B=slice, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_range_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_range_A, vcl_slice_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    
    /////// A=slice:
    std::cout << "Testing A=slice, B=matrix, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=slice, B=matrix, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=slice, B=matrix, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    

    std::cout << "Testing A=slice, B=range, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_range_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=slice, B=range, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_range_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=slice, B=range, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_range_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_range_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    

    
    std::cout << "Testing A=slice, B=slice, C=matrix ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_slice_B, vcl_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    
    std::cout << "Testing A=slice, B=slice, C=range ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_range_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_slice_B, vcl_range_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    std::cout << "Testing A=slice, B=slice, C=slice ..." << std::endl;
    viennacl::copy(ublas_A, vcl_slice_A);
    viennacl::copy(ublas_B, vcl_slice_B);
    viennacl::copy(ublas_C, vcl_slice_C);
    if (run_test(epsilon,
                 ublas_A, ublas_B, ublas_C,
                 vcl_slice_A, vcl_slice_B, vcl_slice_C) != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }
    

    return EXIT_SUCCESS;
}

int main (int, const char **)
{
  std::cout << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "## Test :: Matrix Range" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << std::endl;
   
  double epsilon = 1e-4;
  std::cout << "# Testing setup:" << std::endl;
  std::cout << "  eps:     " << epsilon << std::endl;
  std::cout << "  numeric: float" << std::endl;
  if (run_test<viennacl::row_major, float>(epsilon) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  if (run_test<viennacl::column_major, float>(epsilon) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  
  
#ifdef VIENNACL_HAVE_OPENCL   
   if( viennacl::ocl::current_device().double_support() )
#endif
  {
    double epsilon = 1e-12;
    std::cout << "# Testing setup:" << std::endl;
    std::cout << "  eps:     " << epsilon << std::endl;
    std::cout << "  numeric: double" << std::endl;
    
    if (run_test<viennacl::row_major, double>(epsilon) != EXIT_SUCCESS)
      return EXIT_FAILURE;
    if (run_test<viennacl::column_major, double>(epsilon) != EXIT_SUCCESS)
      return EXIT_FAILURE;
  }

   std::cout << std::endl;
   std::cout << "------- Test completed --------" << std::endl;
   std::cout << std::endl;
   
  
  return EXIT_SUCCESS;
}

