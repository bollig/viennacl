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

//disable debug mechanisms to have a fair comparison with ublas:
#ifndef NDEBUG
 #define NDEBUG
#endif


//
// include necessary system headers
//
#include <iostream>

//
// ublas includes
//
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_HAVE_UBLAS 1

//
// ViennaCL includes
//
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"

// Some helper functions for this tutorial:
#include "Random.hpp"
#include "vector-io.hpp"

#include "../benchmarks/benchmark-utils.hpp"

/*
*   Tutorial: BLAS level 3 functionality
*   
*/

#define BLAS3_MATRIX_SIZE   700

using namespace boost::numeric;


#ifndef VIENNACL_HAVE_OPENCL
  struct dummy
  {
    std::size_t size() const { return 1; }
  };
#endif
    

int main()
{
  typedef float     ScalarType;

  Timer timer;
  double exec_time;

  //
  // Set up some ublas objects
  //
  ublas::matrix<ScalarType> ublas_A(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  ublas::matrix<ScalarType, ublas::column_major> ublas_B(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  ublas::matrix<ScalarType> ublas_C(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  ublas::matrix<ScalarType> ublas_C1(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  ublas::matrix<ScalarType> ublas_C2(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);

  //
  // One alternative: Put the matrices into a contiguous block of memory (allows to use viennacl::fast_copy(), avoiding temporary memory)
  //
  std::vector<ScalarType> stl_A(BLAS3_MATRIX_SIZE * BLAS3_MATRIX_SIZE);
  std::vector<ScalarType> stl_B(BLAS3_MATRIX_SIZE * BLAS3_MATRIX_SIZE);
  std::vector<ScalarType> stl_C(BLAS3_MATRIX_SIZE * BLAS3_MATRIX_SIZE);

  //
  // Fill the matrix
  //
  for (unsigned int i = 0; i < ublas_A.size1(); ++i)
    for (unsigned int j = 0; j < ublas_A.size2(); ++j)
    {
      ublas_A(i,j) = random<ScalarType>();
      stl_A[i*ublas_A.size2() + j] = ublas_A(i,j);
    }

  for (unsigned int i = 0; i < ublas_B.size1(); ++i)
    for (unsigned int j = 0; j < ublas_B.size2(); ++j)
    {
      ublas_B(i,j) = random<ScalarType>();
      stl_B[i + j*ublas_B.size1()] = ublas_B(i,j);
    }

  //
  // Set up some ViennaCL objects
  //
  //viennacl::ocl::set_context_device_type(0, viennacl::ocl::gpu_tag());  //uncomment this is you wish to use GPUs only
  viennacl::matrix<ScalarType> vcl_A(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  viennacl::matrix<ScalarType, viennacl::column_major> vcl_B(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  viennacl::matrix<ScalarType> vcl_C(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
  
  /////////////////////////////////////////////////
  //////////// Matrix-matrix products /////////////
  /////////////////////////////////////////////////
  
  //
  // Compute reference product using ublas:
  //
  std::cout << "--- Computing matrix-matrix product using ublas ---" << std::endl;
  timer.start();
  ublas_C = ublas::prod(ublas_A, ublas_B);
  exec_time = timer.get();
  std::cout << " - Execution time: " << exec_time << std::endl;
  
  //
  // Now iterate over all OpenCL devices in the context and compute the matrix-matrix product
  //
  std::cout << std::endl << "--- Computing matrix-matrix product on each available compute device using ViennaCL ---" << std::endl;
#ifdef VIENNACL_HAVE_OPENCL
  std::vector<viennacl::ocl::device> devices = viennacl::ocl::current_context().devices();
#else
  dummy devices;
#endif
  
  for (std::size_t i=0; i<devices.size(); ++i)
  {
#ifdef VIENNACL_HAVE_OPENCL
    viennacl::ocl::current_context().switch_device(devices[i]);
    std::cout << " - Device Name: " << viennacl::ocl::current_device().name() << std::endl;
#endif

    //viennacl::copy(ublas_A, vcl_A);
    //viennacl::copy(ublas_B, vcl_B);
    viennacl::fast_copy(&(stl_A[0]),
                        &(stl_A[0]) + stl_A.size(),
                        vcl_A);
    viennacl::fast_copy(&(stl_B[0]),
                        &(stl_B[0]) + stl_B.size(),
                        vcl_B);
    vcl_C = viennacl::linalg::prod(vcl_A, vcl_B);
#ifdef VIENNACL_HAVE_OPENCL
    viennacl::ocl::get_queue().finish();
#endif
    timer.start();
    vcl_C = viennacl::linalg::prod(vcl_A, vcl_B);
#ifdef VIENNACL_HAVE_OPENCL
    viennacl::ocl::get_queue().finish();
#endif
    exec_time = timer.get();
    std::cout << " - Execution time on device (no setup time included): " << exec_time << std::endl;
    
    //
    // Verify the result
    //
    //viennacl::copy(vcl_C, ublas_C1);
    viennacl::fast_copy(vcl_C, &(stl_C[0]));
    for (std::size_t i = 0; i < ublas_C1.size1(); ++i)
      for (std::size_t j = 0; j < ublas_C1.size2(); ++j)
        ublas_C1(i,j) = stl_C[i * ublas_C1.size2() + j];

    std::cout << " - Checking result... ";
    bool check_ok = true;
    for (std::size_t i = 0; i < ublas_A.size1(); ++i)
    {
      for (std::size_t j = 0; j < ublas_A.size2(); ++j)
      {
        if ( fabs(ublas_C1(i,j) - ublas_C(i,j)) / ublas_C(i,j) > 1e-4 )
        {
          check_ok = false;
          break;
        }
      }
      if (!check_ok)
        break;
    }
    if (check_ok)
      std::cout << "[OK]" << std::endl << std::endl;
    else
      std::cout << "[FAILED]" << std::endl << std::endl;
      
  }

  //
  //  That's it. 
  //
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
  return EXIT_SUCCESS;
}

