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
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_HAVE_UBLAS 1


//
// ViennaCL includes
//
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"


// Some helper functions for this tutorial:
#include "Random.hpp"
#include "vector-io.hpp"


/*
*
*   Tutorial:  Modification of sparse matrices
*   
*/
using namespace boost::numeric;


int main()
{
  typedef float       ScalarType;
  
  std::size_t size = 5;
  
  //
  // Set up some ublas objects
  //
  ublas::vector<ScalarType> rhs(size, size);
  ublas::compressed_matrix<ScalarType> ublas_matrix(size, size);
  
  ublas_matrix(0,0) =  2.0; ublas_matrix(0,1) = -1.0;
  ublas_matrix(1,0) = -1.0; ublas_matrix(1,1) =  2.0; ublas_matrix(1,2) = -1.0;
  ublas_matrix(2,1) = -1.0; ublas_matrix(2,2) =  2.0; ublas_matrix(2,3) = -1.0;
  ublas_matrix(3,2) = -1.0; ublas_matrix(3,3) =  2.0; ublas_matrix(3,4) = -1.0;
  ublas_matrix(4,3) = -1.0; ublas_matrix(4,4) =  2.0;

  std::cout << "ublas: " << ublas_matrix << std::endl;
  
  //
  // Set up some ViennaCL objects
  //
  viennacl::vector<ScalarType> vcl_rhs(size); 
  viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(size, size);

  viennacl::copy(ublas_matrix, vcl_compressed_matrix);

  // just get the data directly from the GPU and print it:
  ublas::compressed_matrix<ScalarType> temp(size, size);
  viennacl::copy(vcl_compressed_matrix, temp);
  std::cout << "ViennaCL: " << temp << std::endl;

  // now modify GPU data directly:
  std::cout << "Modifying vcl_compressed_matrix a bit: " << std::endl;
  vcl_compressed_matrix(0, 0) = 3.0;
  vcl_compressed_matrix(2, 3) = -3.0;
  vcl_compressed_matrix(4, 2) = -3.0;  //this is a new nonzero entry
  vcl_compressed_matrix(4, 3) = -3.0;
  
  // and print it again:
  viennacl::copy(vcl_compressed_matrix, temp);
  std::cout << "ViennaCL: " << temp << std::endl;
  
  //
  //  That's it.
  //
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
  
  return 0;
}

