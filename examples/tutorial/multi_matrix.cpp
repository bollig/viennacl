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

#define VIENNACL_PROFILING_ENABLE
#define VIENNACL_DEBUG_SCHEDULER
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
    unsigned int size = 15360;
    typedef float ScalarType;

    typedef viennacl::distributed::multi_matrix<ScalarType> gpu_mat_t;

    viennacl::distributed::scheduler::add_all_available_devices(CL_DEVICE_TYPE_GPU);

    gpu_mat_t C(size,size);
    gpu_mat_t A(size,size);
    gpu_mat_t B(size,size);

    gpu_mat_t D(size,size);
    gpu_mat_t E(size,size);


    Timer t;
    t.start();

//    C=viennacl::generator::prod(A,B);
//    viennacl::distributed::scheduler::finish();

    D = A+B;
    viennacl::distributed::scheduler::finish();
    E = A-B;
    viennacl::distributed::scheduler::finish();
    C = viennacl::generator::prod(D,E);
    viennacl::distributed::scheduler::finish();




    std::cout << "\n\n////////////\n\n" << std::endl;
    std::cout << "Total Time : " << t.get() << std::endl;
}
