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

#define VIENNACL_USE_SCHEDULER
#define VIENNACL_ENABLE_AUTOTUNE
#define VIENNACL_PROFILING_ENABLE
//#define VIENNACL_DEBUG_ALL

//#define VIENNACL_WITH_OPENCL

// include necessary system headers
#include <iostream>

#include "CL/cl.h"
#include "../benchmarks/benchmark-utils.hpp"

#include "viennacl/distributed/multi_matrix.hpp"
//#include "viennacl/distributed/multi_vector.hpp"

#include "viennacl/distributed/scheduler.hpp"
#include "viennacl/distributed/fission.hpp"

#include "viennacl/linalg/prod.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/event.hpp"

//Boost includes
#include "boost/numeric/ublas/matrix.hpp"

int main(){
    unsigned int size1 = 11264;
    unsigned int size2 = 11264;
    typedef float ScalarType;

    typedef viennacl::distributed::multi_matrix<ScalarType,viennacl::row_major,1> gpu_mat_t;
//    typedef viennacl::distributed::multi_vector<ScalarType,1> gpu_vec_t;

    viennacl::distributed::scheduler::add_all_available_devices(CL_DEVICE_TYPE_GPU);

    gpu_mat_t mat(size1,size2);
    gpu_mat_t mat2(size1,size2);
    gpu_mat_t mat3(size1,size2);

//    gpu_vec_t vec(size2);
//    gpu_vec_t res(size2);

//    for(unsigned int i = 0 ; i < size1 ; ++i){
//        for(unsigned int j=0 ; j < size2 ; ++j){
//            mat(i,j) = 0;
//            mat2(i,j) = i;
//            mat3(i,j) = i;
//        }
//    }


    Timer t;
    t.start();
//    res = viennacl::linalg::prod(mat,vec);

    viennacl::generator::dummy_matrix<gpu_mat_t> dummy_mat(mat);
    viennacl::generator::dummy_matrix<gpu_mat_t> dummy_mat2(mat2);
    viennacl::generator::dummy_matrix<gpu_mat_t> dummy_mat3(mat3);

    mat=viennacl::generator::prod(mat2,mat3);
    viennacl::distributed::scheduler::finish();

    std::cout << "\n\n////////////\n\n" << std::endl;
    std::cout << "Total Time : " << t.get() << std::endl;
}
