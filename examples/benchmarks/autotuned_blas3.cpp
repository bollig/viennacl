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
*   Benchmark: BLAS level 3 functionality for dense matrices (blas3.cpp and blas3.cu are identical, the latter being required for compilation using CUDA nvcc)
*   
*/

//#define VIENNACL_DEBUG_ALL

//disable debug mechanisms to have a fair benchmark environment
#ifndef NDEBUG
#define NDEBUG
#endif

//
// include necessary system headers
//
#include <iostream>

//
// ViennaCL includes
//
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/linalg/lu.hpp"

#include "viennacl/generator/custom_operation.hpp"

// Some helper functions for this tutorial:
#include "../tutorial/Random.hpp"


#include "benchmark-utils.hpp"


#define BLAS3_MATRIX_SIZE   2048

template<typename ScalarType, class FA, class FB>
int run_benchmark()
{
    Timer timer;
    double exec_time;

    //
    // One alternative: Put the matrices into a contiguous block of memory (allows to use viennacl::fast_copy(), avoiding temporary memory)
    //
    std::vector<ScalarType> stl_A(BLAS3_MATRIX_SIZE * BLAS3_MATRIX_SIZE);
    std::vector<ScalarType> stl_B(BLAS3_MATRIX_SIZE * BLAS3_MATRIX_SIZE);
    std::vector<ScalarType> stl_C(BLAS3_MATRIX_SIZE * BLAS3_MATRIX_SIZE);

    //
    // Fill the matrix
    //
    for (unsigned int i = 0; i < BLAS3_MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < BLAS3_MATRIX_SIZE; ++j)
            stl_A[i*BLAS3_MATRIX_SIZE + j] = random<ScalarType>();

    for (unsigned int i = 0; i < BLAS3_MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < BLAS3_MATRIX_SIZE; ++j)
            stl_B[i + j*BLAS3_MATRIX_SIZE] = random<ScalarType>();

    typedef viennacl::matrix<ScalarType> MatTypeC;
    typedef viennacl::matrix<ScalarType,FA> MatTypeA;
    typedef viennacl::matrix<ScalarType,FB> MatTypeB;
    //viennacl::ocl::current_context().build_options("-cl-mad-enable -cl-fast-relaxed-math");   //uncomment for additional optimizations
    viennacl::ocl::current_context().build_options("-cl-opt-disable");                        //uncomment to get poor performance
    MatTypeC vcl_C(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
    MatTypeA vcl_A(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);
    MatTypeB vcl_B(BLAS3_MATRIX_SIZE, BLAS3_MATRIX_SIZE);


    /////////////////////////////////////////////////
    //////////// Matrix-matrix products /////////////
    /////////////////////////////////////////////////


    viennacl::fast_copy(&(stl_C[0]),
                        &(stl_C[0]) + stl_C.size(),
                        vcl_C);
    viennacl::fast_copy(&(stl_A[0]),
                        &(stl_A[0]) + stl_A.size(),
                        vcl_A);
    viennacl::fast_copy(&(stl_B[0]),
                        &(stl_B[0]) + stl_B.size(),
                        vcl_B);

    viennacl::generator::dummy_matrix< MatTypeC > dvcl_C(vcl_C);
    viennacl::generator::dummy_matrix< MatTypeA > dvcl_A(vcl_A);
    viennacl::generator::dummy_matrix< MatTypeB  > dvcl_B(vcl_B);



    std::cout << " ----- AA ----- " << std::endl;
    {
    viennacl::generator::custom_operation op( dvcl_C = viennacl::generator::prod( dvcl_A, dvcl_B));
    op.execute();
    viennacl::backend::finish();
    timer.start();
    op.execute();
    viennacl::backend::finish();
    exec_time = timer.get();
    std::cout << " - Execution time on device (no setup time included): " << exec_time << std::endl;
    std::cout << " - GFLOPs (counting multiply&add as one operation): " << (vcl_A.size1() / 1000.0) * (vcl_A.size2() / 1000.0) * (vcl_B.size2() / 1000.0) / exec_time << std::endl;
    std::cout << std::endl;
    }

//    {
//    std::cout << " ----- TA ----- " << std::endl;
//    viennacl::generator::custom_operation op( dvcl_C = viennacl::generator::prod( trans(dvcl_A), dvcl_B) );
//    op.execute();
//    viennacl::backend::finish();
//;
//    timer.start();
//    op.execute();
//    viennacl::backend::finish();
//;
//    exec_time = timer.get();
//    std::cout << " - Execution time on device (no setup time included): " << exec_time << std::endl;
//    std::cout << " - GFLOPs (counting multiply&add as one operation): " << (vcl_A.size1() / 1000.0) * (vcl_A.size2() / 1000.0) * (vcl_B.size2() / 1000.0) / exec_time << std::endl;
//    std::cout << std::endl;
//    }

//    {
//    std::cout << " ----- AT ----- " << std::endl;
//    viennacl::generator::custom_operation op( dvcl_C = viennacl::generator::prod( dvcl_A, trans(dvcl_B)) );
//    op.execute();
//    viennacl::backend::finish();
//;
//    timer.start();
//    op.execute();
//    viennacl::backend::finish();
//;
//    exec_time = timer.get();
//    std::cout << " - Execution time on device (no setup time included): " << exec_time << std::endl;
//    std::cout << " - GFLOPs (counting multiply&add as one operation): " << (vcl_A.size1() / 1000.0) * (vcl_A.size2() / 1000.0) * (vcl_B.size2() / 1000.0) / exec_time << std::endl;
//    std::cout << std::endl;
//    }

//    {
//    std::cout << " ----- TT ----- " << std::endl;
//    viennacl::generator::custom_operation op( dvcl_C = viennacl::generator::prod( trans(dvcl_A), trans(dvcl_B)) );
//    op.execute();
//    viennacl::backend::finish();
//;
//    timer.start();
//    op.execute();
//    viennacl::backend::finish();
//;
//    exec_time = timer.get();
//    std::cout << " - Execution time on device (no setup time included): " << exec_time << std::endl;
//    std::cout << " - GFLOPs (counting multiply&add as one operation): " << (vcl_A.size1() / 1000.0) * (vcl_A.size2() / 1000.0) * (vcl_B.size2() / 1000.0) / exec_time << std::endl;
//    std::cout << std::endl;
//    }

    return EXIT_SUCCESS;
}

int main()
{

    typedef std::vector< viennacl::ocl::platform > platforms_type;
    typedef std::vector<viennacl::ocl::device> devices_type;
    typedef std::vector<cl_device_id> cl_devices_type;

    platforms_type platforms = viennacl::ocl::get_platforms();
    size_t num_platforms = platforms.size();
    for(unsigned int k=0 ; k < num_platforms ; ++k)
    {
        viennacl::ocl::platform pf(k);
        viennacl::ocl::set_context_platform_index(k,k);
        viennacl::ocl::switch_context(k);
        devices_type dev = viennacl::ocl::current_context().devices();
        for(devices_type::iterator it = dev.begin() ; it != dev.end() ; ++it){
            if(it!=dev.begin()) continue;
            viennacl::ocl::switch_device(*it);
            cl_device_id dev_id = it->id();

            if(viennacl::ocl::info<CL_DEVICE_TYPE>(dev_id)==CL_DEVICE_TYPE_GPU){

            std::cout << std::endl;
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "               Device Info" << std::endl;
            std::cout << "----------------------------------------------" << std::endl;

#ifdef VIENNACL_WITH_OPENCL
            std::cout << viennacl::ocl::current_device().info() << std::endl;
#endif  



            std::cout << std::endl;
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "## Benchmark :: Dense Matrix-Matrix product " << std::endl;
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << std::endl;
            std::cout << "   -------------------------------" << std::endl;
            std::cout << "   # benchmarking single-precision" << std::endl;
            std::cout << "   -------------------------------" << std::endl;
            std::cout << "   LHS : Row-Major | RHS : Row-Major"<< std::endl;
            std::cout << "   ---------------------------------"<< std::endl;
            run_benchmark<float,viennacl::row_major,viennacl::row_major>();
            std::cout << "   ---------------------------------"<< std::endl;
            std::cout << "   LHS : Column-Major | RHS : Row-Major"<< std::endl;
            std::cout << "   ---------------------------------"<< std::endl;
            run_benchmark<float,viennacl::column_major,viennacl::row_major>();
//            std::cout << "   ---------------------------------"<< std::endl;
//            std::cout << "   LHS : Row-Major | RHS : Column-Major"<< std::endl;
//            std::cout << "   ---------------------------------"<< std::endl;
//            run_benchmark<float,viennacl::row_major,viennacl::column_major>();
//            std::cout << "   ---------------------------------"<< std::endl;
//            std::cout << "   LHS : Column-Major | RHS : Column-Major"<< std::endl;
//            std::cout << "   ---------------------------------"<< std::endl;
//            run_benchmark<float,viennacl::column_major,viennacl::column_major>();

            //#ifdef VIENNACL_WITH_OPENCL
            //  if( viennacl::ocl::current_device().double_support() )
            //#endif
            //  {
            //    std::cout << std::endl;
            //    std::cout << "   -------------------------------" << std::endl;
            //    std::cout << "   # benchmarking double-precision" << std::endl;
            //    std::cout << "   -------------------------------" << std::endl;
            //    run_benchmark<double,viennacl::row_major,viennacl::row_major>();
            //    std::cout << "   ---------------------------------"<< std::endl;
            //    run_benchmark<double,viennacl::column_major,viennacl::row_major>();
            //    std::cout << "   ---------------------------------"<< std::endl;
            //    run_benchmark<double,viennacl::row_major,viennacl::column_major>();
            //    std::cout << "   ---------------------------------"<< std::endl;
            //    run_benchmark<double,viennacl::column_major,viennacl::column_major>();
            //  }
            }
        }
    }
    return 0;
}
