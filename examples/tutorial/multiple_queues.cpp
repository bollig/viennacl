/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

// activate ublas support in ViennaCL
#define VIENNACL_HAVE_UBLAS
#define VIENNACL_DEBUG_ALL

//
// include necessary system headers
//
#include <iostream>
#include <string>

//
// ViennaCL includes
//
#include "viennacl/scalar.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/vector_proxy.hpp"

//
// Boost includes
//
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/io.hpp"

/*
*   Tutorial: Use of matrix ranges
*/

int main (int argc, const char * argv[])
{
  typedef float                                           ScalarType;    //feel free to change this to 'double' if supported by your hardware
  typedef boost::numeric::ublas::vector<ScalarType>       VectorType;

  typedef viennacl::vector<ScalarType>                    VCLVectorType;

  std::size_t dim_large = 7;
  std::size_t dim_small = 3;

  //
  // Setup ublas objects and fill with data:
  //
  VectorType ublas_v1(dim_large);
  VectorType ublas_v2(dim_small);

  for (std::size_t i=0; i<ublas_v1.size(); ++i)
    ublas_v1(i) = i+1;

  for (std::size_t i=0; i<ublas_v2.size(); ++i)
    ublas_v2(i) = dim_large + i;


  //
  // Extract submatrices using the ranges in ublas
  //
  boost::numeric::ublas::range ublas_r1(0, dim_small); //the first 'dim_small' entries
  boost::numeric::ublas::range ublas_r2(dim_small - 1, 2*dim_small - 1); // 'dim_small' entries somewhere from the middle
  boost::numeric::ublas::range ublas_r3(dim_large - dim_small, dim_large); // the last 'dim_small' entries
  boost::numeric::ublas::vector_range<VectorType> ublas_v1_sub1(ublas_v1, ublas_r1); // front part of vector v_1
  boost::numeric::ublas::vector_range<VectorType> ublas_v1_sub2(ublas_v1, ublas_r2); // center part of vector v_1
  boost::numeric::ublas::vector_range<VectorType> ublas_v1_sub3(ublas_v1, ublas_r3); // tail of vector v_1


  //
  // Setup ViennaCL objects
  //
  VCLVectorType vcl_v1(dim_large);
  VCLVectorType vcl_v2(dim_small);

  std::cout << "HERE HERE HERE1\n";

  // Add a new queue for processing (will output debug msg)
  viennacl::ocl::current_context().add_queue(viennacl::ocl::current_device());


  // Need to add Event profiling
  viennacl::ocl::switch_queue(0);
  viennacl::copy(ublas_v1, vcl_v1);
  std::cout << "HERE HERE HERE2\n";
  viennacl::ocl::switch_queue(1);
  viennacl::copy(ublas_v2, vcl_v2);

  std::cout << "HERE HERE HERE3\n";

  //
  // Extract submatrices using the ranges in ViennaCL
  //
  viennacl::range vcl_r1(0, dim_small); //the first 'dim_small' entries
  viennacl::range vcl_r2(dim_small - 1, 2*dim_small - 1);  // 'dim_small' entries somewhere from the middle
  viennacl::range vcl_r3(dim_large - dim_small, dim_large);  // the last 'dim_small' entries
  viennacl::vector_range<VCLVectorType>   vcl_v1_sub1(vcl_v1, vcl_r1); // front part of vector v_1
  viennacl::vector_range<VCLVectorType>   vcl_v1_sub2(vcl_v1, vcl_r2); // center part of vector v_1
  viennacl::vector_range<VCLVectorType>   vcl_v1_sub3(vcl_v1, vcl_r3); // tail of vector v_1


  //
  // Copy from ublas to submatrices and back:
  //

  ublas_v1_sub1 = ublas_v2;
  viennacl::copy(ublas_v2, vcl_v1_sub1);
  viennacl::copy(vcl_v1_sub1, ublas_v2);

  //
  // Addition:
  //

  ublas_v1_sub1 += ublas_v1_sub1;
  vcl_v1_sub1 += vcl_v1_sub1;

  ublas_v1_sub2 += ublas_v1_sub2;
  vcl_v1_sub2 += vcl_v1_sub2;

  ublas_v1_sub3 += ublas_v1_sub3;
  vcl_v1_sub3 += vcl_v1_sub3;

  //
  // print vectors:
  //
  std::cout << "ublas:    " << ublas_v1 << std::endl;
  std::cout << "ViennaCL: " << vcl_v1 << std::endl;

  //
  //  That's it.
  //
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

  return EXIT_SUCCESS;
}

