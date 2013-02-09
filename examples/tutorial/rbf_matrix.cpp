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

/*
*
*   Tutorial: Computes an RBF Matrix :
*
*   C(i,j) = norm2( A(i,:) - B(:,j) )
*
*/

//#define VIENNACL_HAVE_OPENCL

// include necessary system headers
#include <iostream>

#include "viennacl/generator_fromscratch/custom_operation.hpp"
#include "viennacl/matrix.hpp"

#include "boost/numeric/ublas/matrix.hpp"

#include "Random.hpp"
#include "vector"

static const unsigned int SIZE = 512;
typedef float NumericT;
using namespace viennacl::generator;

int main(){
    viennacl::matrix<NumericT> A(SIZE,SIZE);
    viennacl::matrix<NumericT> B(SIZE,SIZE);
    viennacl::matrix<NumericT> C(SIZE,SIZE);

    typedef viennacl::generator::dummy_matrix<NumericT, viennacl::row_major> dm_t;

    boost::numeric::ublas::matrix<NumericT> cpu_val(SIZE,SIZE);
    for(unsigned int i=0; i<cpu_val.size1(); ++i)
        for(unsigned int j=0; j<cpu_val.size2(); ++j)
            cpu_val(i,j) = random<NumericT>();

    viennacl::copy(cpu_val, B);
    viennacl::copy(cpu_val, C);

    viennacl::generator::custom_operation op;
    op.add(dm_t(A) = prod_based<add_type>(dm_t(B), dm_t(C),"pow(#1 - #2,2)"));
    op.execute();
    std::cout << op.source_code() << std::endl;




}
