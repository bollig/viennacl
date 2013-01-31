#ifndef VIENNACL_GENERATOR_AUTOTUNE_HPP
#define VIENNACL_GENERATOR_AUTOTUNE_HPP

#include "viennacl/generator_fromscratch/autotune/benchmark-utils.hpp"
#include "ctime"
#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"

namespace viennacl{

namespace generator{

namespace autotune{

typedef std::map<double, viennacl::tools::shared_ptr<viennacl::generator::code_generation::optimization_profile> > timings_t;

template<class ConfigT>
void benchmark_blas3(timings_t & timings, ConfigT const & config){

    viennacl::ocl::device const & dev = viennacl::ocl::current_device();

    unsigned int n_runs =config.n_runs;


    for(unsigned int ml = config.ml_min ; ml <= config.ml_max; ml*=2){
        for(unsigned int kl = config.kl_min ; kl <= config.kl_max; kl*=2){
            for(unsigned int nl = config.nl_min ; nl <= config.nl_max; nl*=2){
                for(unsigned int ms = config.ms_min ; ms <= config.ms_max; ms*=2){
                    for(unsigned int ks = config.ks_min ; ks <= config.ks_max; ks*=2){
                        for(unsigned int ns = config.ns_min ; ns <= config.ns_max; ns*=2){
                            for(unsigned int alignment = config.alignment_min ; alignment <= config.alignment_max; alignment *=2){
                                viennacl::generator::code_generation::blas3_optimization_profile prof(ml,kl,nl,ms,ks,ns,false,true,alignment);

                                if(alignment>ks || alignment>ns) continue;
                                if((double)ml*(kl+1)*4/1024 > 32.0) continue;
                                if(prof.local_work_size(0)*prof.local_work_size(1) > dev.max_workgroup_size()) continue;

                                std::ostringstream oss;
//                                oss << "p" << ml << "_" << kl << "_"  << nl << "_"  << ms << "_"  << ks << "_"  << ns << "_"  << alignment ;
                                viennacl::generator::custom_operation op;
                                config.add_op(op);
                                op.operations_manager().blas3_model() = prof;
                                op.init();
//                                op.program_name(oss.str());
                                viennacl::ocl::program & pgm = op.program();

                                viennacl::ocl::kernel & k = pgm.get_kernel("_k0");

                                //Anticipates kernel failure
                                size_t max_workgroup_size = k.max_workgroup_size(dev.id());
                                if(prof.local_work_size(0)*prof.local_work_size(1) > max_workgroup_size)
                                    continue;

                                //Doesn't execute because it would likelily be a waste of time
                                size_t prefered_workgroup_size_multiple = k.prefered_work_group_size_multiple(dev.id());
                                if( (prof.local_work_size(0)*prof.local_work_size(1)) % prefered_workgroup_size_multiple > 0)
                                    continue;

                                op.execute();
                                viennacl::ocl::get_queue().finish();
                                op.execute();
                                viennacl::ocl::get_queue().finish();

                                double exec_time = 0;
                                for(unsigned int n=0; n<n_runs ; ++n){
                                    op.execute();
                                    Timer t; t.start();
                                    viennacl::ocl::get_queue().flush();
                                    t.start();
                                    viennacl::ocl::get_queue().finish();
                                    exec_time+=t.get();
                                }
                                exec_time = exec_time/n_runs;
                                std::cout << exec_time <<  " " << prof << std::endl;
                                timings.insert(std::make_pair(exec_time, new viennacl::generator::code_generation::blas3_optimization_profile(prof)));
                            }
                        }
                    }
                }
            }
        }
    }
}

}

}

}
#endif // AUTOTUNE_HPP
