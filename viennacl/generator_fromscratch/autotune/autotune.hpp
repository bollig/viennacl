#ifndef VIENNACL_GENERATOR_AUTOTUNE_HPP
#define VIENNACL_GENERATOR_AUTOTUNE_HPP

#include "viennacl/generator_fromscratch/autotune/benchmark-utils.hpp"
#include "ctime"
#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"

#include "boost/numeric/ublas/matrix.hpp"


namespace viennacl{

namespace generator{

namespace autotune{

typedef std::map<double, viennacl::tools::shared_ptr<viennacl::generator::code_generation::optimization_profile> > timings_map_t;
typedef typename timings_map_t::value_type  timing_t;

template<class ConfigT>
timing_t benchmark_blas3(ConfigT const & config){
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

    for(unsigned int ml = config.ml_min ; ml <= config.ml_max; ml*=2){
        for(unsigned int kl = config.kl_min ; kl <= config.kl_max; kl*=2){
            for(unsigned int nl = config.nl_min ; nl <= config.nl_max; nl*=2){
                for(unsigned int ms = config.ms_min ; ms <= config.ms_max; ms*=2){
                    for(unsigned int ks = config.ks_min ; ks <= config.ks_max; ks*=2){
                        for(unsigned int ns = config.ns_min ; ns <= config.ns_max; ns*=2){
                            for(unsigned int alignment = 1 ; alignment <= config.alignment_max; alignment *=2){
                                if(alignment>ks || alignment>ns) continue;
                                std::ostringstream oss;
//                                oss << "p" << ml << "_" << kl << "_"  << nl << "_"  << ms << "_"  << ks << "_"  << ns << "_"  << alignment << "_"  << std::endl;
                                std::cout << oss.str() << std::endl;
                                viennacl::generator::custom_operation op;
                                op.add(dm(A) = prod(dm(B),dm(C)));
                                viennacl::generator::code_generation::blas3_optimization_profile prof(ml,kl,nl,ms,ks,ns,true,false,alignment);
                                op.operations_manager().blas3_model() = prof;
                                op.init();
                                double exec_time;
                                for(unsigned int n=0; n<2+BENCHMARK_RUNS ; ++n){
                                    op.execute(oss.str());
                                    Timer t; t.start();
                                    viennacl::ocl::get_queue().finish();
                                    if(n>2) exec_time+=t.get();
                                }
                                exec_time = exec_time/BENCHMARK_RUNS;
                                std::cout << exec_time <<  " " << prof << std::endl;
                                timings.insert(std::make_pair(exec_time, new viennacl::generator::code_generation::blas3_optimization_profile(prof)));
                            }
                        }
                    }
                }
            }
        }
    }

    return *timings.begin();
}

}

}

}
#endif // AUTOTUNE_HPP
