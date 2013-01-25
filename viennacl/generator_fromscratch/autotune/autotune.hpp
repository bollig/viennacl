#ifndef VIENNACL_GENERATOR_AUTOTUNE_HPP
#define VIENNACL_GENERATOR_AUTOTUNE_HPP

#include "viennacl/generator_fromscratch/forwards.h"
//#include "viennacl/generator_fromscratch/autotune/benchmark-utils.hpp"
#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"

#include "boost/numeric/ublas/matrix.hpp"


namespace viennacl{

namespace generator{

namespace autotune{

typedef std::map<double, viennacl::tools::shared_ptr<viennacl::generator::code_generation::optimization_profile> > timings_map_t;
typedef typename timings_map_t::value_type  timing_t;

timing_t benchmark_blas3(){
    typedef float NumericT;

    static const int BENCHMARK_RUNS = 10;

    timings_map_t timings;
    typedef viennacl::generator::dummy_matrix<NumericT, viennacl::row_major> dm;

    unsigned int size=2048;

    viennacl::matrix<NumericT,viennacl::row_major,16> A(size,size);
    viennacl::matrix<NumericT,viennacl::row_major,16> B(size,size);
    viennacl::matrix<NumericT,viennacl::row_major,16> C(size,size);

    boost::numeric::ublas::matrix<NumericT> cpu_A(size,size);
    boost::numeric::ublas::matrix<NumericT> cpu_B(size,size);
    boost::numeric::ublas::matrix<NumericT> cpu_C(size,size);
    for(unsigned int i=0; i<size; ++i){
        for(unsigned int j=0 ; j<size ; ++j){
            cpu_A(i,j)=0;
            cpu_B(i,j)=i+j;
            cpu_C(i,j)=j;
        }
    }

    viennacl::copy(cpu_A,A);
    viennacl::copy(cpu_B,B);
    viennacl::copy(cpu_C,C);
    viennacl::ocl::get_queue().finish();

    viennacl::generator::custom_operation op;
    op.add(dm(A) = prod(dm(B),dm(C)));
    viennacl::generator::code_generation::blas3_optimization_profile prof(16,256,256,4,4,4,true,false,4);
    op.operations_manager().blas3_model() = prof;
    op.init();

    double exec_time;
    for(unsigned int n=0; n<2+BENCHMARK_RUNS ; ++n){
        op.execute();
        Timer t; t.start();
        viennacl::ocl::get_queue().finish();
        if(n>2) exec_time+=t.get();
    }
    exec_time = exec_time/BENCHMARK_RUNS;

    timings.insert(std::make_pair(exec_time, new viennacl::generator::code_generation::blas3_optimization_profile(prof)));
    return *timings.begin();
}

}

}

}
#endif // AUTOTUNE_HPP
