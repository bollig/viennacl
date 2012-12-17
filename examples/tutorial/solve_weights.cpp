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

#define VIENNACL_HAS_UBLAS 1

//
// include necessary system headers
//
#include <iostream>
#include <string>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
//
// ViennaCL includes
//
#include "viennacl/ocl/backend.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/norm_2.hpp"


// Some helper functions for this tutorial:
#include "Random.hpp"


/*
*
*   Tutorial:  Custom compute kernels
*
*/


//
// Custom compute kernels which compute an elementwise product/division of two vectors
// Input: v1 ... vector
//        v2 ... vector
// Output: result ... vector
//
// Algorithm: set result[i] <- v1[i] * v2[i]
//            or  result[i] <- v1[i] / v2[i]
//            (in MATLAB notation this is something like 'result = v1 .* v2' and 'result = v1 ./ v2');
//
const char * my_compute_program =
"__kernel void elementwise_prod(\n"
"          __global const float * vec1,\n"
"          __global const float * vec2, \n"
"          __global float * result,\n"
"          unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    result[i] = vec1[i] * vec2[i];\n"
"};\n\n"
"__kernel void elementwise_div(\n"
"          __global const float * vec1,\n"
"          __global const float * vec2, \n"
"          __global float * result,\n"
"          unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    result[i] = vec1[i] / vec2[i];\n"
"};\n\n"
" \n"
"// Assume the input matrix is dense blocks (array). Then elements[0->sten_size*sten_size] is block 1.\n"
"__kernel void lu_factor_sten_blocks(\n"
"        __global float * _elements,\n"
"        unsigned int _rows,\n"
"        unsigned int _sten_size\n"
"        )\n"
"{\n"
"    float temp;\n"
"    unsigned int b = 0;\n"
"\n"
"    while (b < _rows) {\n"
"        unsigned int offset = b * _sten_size * _sten_size;\n"
"        for (unsigned int i=1; i < _sten_size; ++i)\n"
"        {\n"
"            for (unsigned int k=0; k<i; ++k)\n"
"            {\n"
"                // This logic assumes the elements are sorted by row then column.\n"
"                if (get_global_id(0) == 0) {\n"
"                    _elements[offset + i + k*_sten_size] /= _elements[offset + k + k*_sten_size];\n"
"                }\n"
"                barrier(CLK_GLOBAL_MEM_FENCE);\n"
"                temp = _elements[offset + i + k*_sten_size];\n"
"\n"
"                //parallel subtraction:\n"
"                for (unsigned int j=k+1 + get_global_id(0); j<_sten_size; j += get_global_size(0)) {\n"
" //                   if (j < _sten_size)\n"
"                    _elements[offset + i + j*_sten_size] -= temp * _elements[offset + k + j*_sten_size];\n"
"                }\n"
"            }\n"
"        }\n"
"        barrier(CLK_GLOBAL_MEM_FENCE);\n"
"        b++;\n"
"    }\n"
"};\n";

using namespace boost::numeric;

int main()
{
  typedef float       ScalarType;

  unsigned int n = 3;
  unsigned int N = 2;

  ublas::compressed_matrix<ScalarType> A_host(N*n,N*n);
  ublas::vector<ScalarType> B_host(N*n);

  viennacl::compressed_matrix<ScalarType> A_dev(N*n,N*n);
  viennacl::vector<ScalarType> B_dev(N*n);

  // Dense block (could also be sparse but wont be for our purposes)
  A_host(0, 0) = 3;
  A_host(0, 1) = 2;
  A_host(0, 2) = 1;
  A_host(1, 0) = 2;
  A_host(1, 1) = 3;
  A_host(1, 2) = 2;
  A_host(2, 0) = 1;
  A_host(2, 1) = 2;
  A_host(2, 2) = 3;

  B_host(0) = 1;
  B_host(1) = 1;
  B_host(2) = 1;

#if 1
  A_host(3, 3) = 3;
  A_host(3, 4) = 2;
  A_host(3, 5) = 1;
  A_host(4, 3) = 2;
  A_host(4, 4) = 3;
  A_host(4, 5) = 2;
  A_host(5, 3) = 1;
  A_host(5, 4) = 2;
  A_host(5, 5) = 3;

  B_host(3) = 1;
  B_host(4) = 1;
  B_host(5) = 1;
#endif
  std::cout << "GPU MAT has nnz = " << A_host.nnz() << std::endl;

  viennacl::copy(A_host, A_dev);
  viennacl::copy(B_host, B_dev);

  std::cout << "GPU MAT has nnz = " << A_dev.nnz() << std::endl;
  std::cout << "GPU_MAT = " << A_host << std::endl;


  //
  // Initialize OpenCL vectors:
  //
  unsigned int vector_size = 10;
  viennacl::scalar<ScalarType>  s = 1.0; //dummy
  viennacl::vector<ScalarType>  vec1(vector_size);
  viennacl::vector<ScalarType>  vec2(vector_size);
  viennacl::vector<ScalarType>  result_mul(vector_size);
  viennacl::vector<ScalarType>  result_div(vector_size);

  //
  // fill the operands vec1 and vec2:
  //
  for (unsigned int i=0; i<vector_size; ++i)
  {
    vec1[i] = static_cast<ScalarType>(i);
    vec2[i] = static_cast<ScalarType>(vector_size-i);
  }

  //
  // Set up the OpenCL program given in my_compute_kernel:
  // A program is one compilation unit and can hold many different compute kernels.
  //
  viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(my_compute_program, "my_compute_program");
  my_prog.add_kernel("elementwise_prod");  //register elementwise product kernel
  my_prog.add_kernel("elementwise_div");   //register elementwise division kernel

  my_prog.add_kernel("lu_factor_sten_blocks");       //register Sparse Block solve

  //
  // Now we can get the kernels from the program 'my_program'.
  // (Note that first all kernels need to be registered via add_kernel() before get_kernel() can be called,
  // otherwise existing references might be invalidated)
  //
  viennacl::ocl::kernel & my_kernel_mul = my_prog.get_kernel("elementwise_prod");
  viennacl::ocl::kernel & my_kernel_div = my_prog.get_kernel("elementwise_div");
  viennacl::ocl::kernel & my_kernel_solve = my_prog.get_kernel("lu_factor_sten_blocks");

  //
  // Launch the kernel with 'vector_size' threads in one work group
  // Note that size_t might differ between host and device. Thus, a cast to cl_uint is necessary for the forth argument.
  //
  viennacl::ocl::enqueue(my_kernel_mul(vec1, vec2, result_mul, static_cast<cl_uint>(vec1.size())));
  viennacl::ocl::enqueue(my_kernel_div(vec1, vec2, result_div, static_cast<cl_uint>(vec1.size())));


  // Call custom LU
  viennacl::ocl::enqueue(my_kernel_solve(A_dev.handle(), static_cast<cl_uint>(A_dev.size1()), static_cast<cl_uint>(n)));

  viennacl::copy(A_dev, A_host);

    std::cout << "OUTPUT: " << A_host << std::endl;

  //
  // Print the result:
  //
  std::cout << "        vec1: " << vec1 << std::endl;
  std::cout << "        vec2: " << vec2 << std::endl;
  std::cout << "vec1 .* vec2: " << result_mul << std::endl;
  std::cout << "vec1 /* vec2: " << result_div << std::endl;
  std::cout << "norm_2(vec1 .* vec2): " << viennacl::linalg::norm_2(result_mul) << std::endl;
  std::cout << "norm_2(vec1 /* vec2): " << viennacl::linalg::norm_2(result_div) << std::endl;

  //
  //  That's it.
  //
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

  return 0;
}

