#ifndef VIENNACL_GENERATOR_AUTOTUNE_HPP
#define VIENNACL_GENERATOR_AUTOTUNE_HPP

#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/generator_fromscratch/autotune/benchmark-utils.hpp"
#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"

namespace viennacl{

namespace generator{

namespace autotune{

static const int BENCHMARK_RUNS = 10;

class config
{
  public:
    config() {}

    unsigned int min_work_groups() const { return min_work_groups_; }
    void min_work_groups(unsigned int i) { min_work_groups_ = i; }
    unsigned int max_work_groups() const { return max_work_groups_; }
    void max_work_groups(unsigned int i) { max_work_groups_ = i; }


    unsigned int min_local_size() const { return min_local_size_; }
    void min_local_size(unsigned int i) { min_local_size_ = i; }
    unsigned int max_local_size() const { return max_local_size_; }
    void max_local_size(unsigned int i) { max_local_size_ = i; }

    unsigned int min_alignment() const { return min_alignment_; }
    void min_alignment(unsigned int i) { min_alignment_ = i; }
    unsigned int max_alignment() const { return max_alignment_; }
    void max_alignment(unsigned int i) { max_alignment_ = i; }

  private:
    unsigned int min_work_groups_;
    unsigned int max_work_groups_;
    unsigned int min_local_size_;
    unsigned int max_local_size_;
    unsigned int min_alignment_;
    unsigned int max_alignment_;
};

template<class TestConfig>
void benchmark_timings(std::vector<std::list<infos_base*> > kernels, TestConfig & config){
    Timer t;
    for(std::vector<std::list<infos_base*> >::iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
        std::string program_name = "k"+to_string(it-kernels.begin())+"_benchmark";
        std::map<double, viennacl::generator::code_generation::optimization_profile> timings;
        std::ostringstream oss;
        code_generation::utils::kernel_generation_stream kss(oss);
        std::map<std::string, generator::code_generation::kernel_infos_t> kernels_infos;
        for (unsigned int work_groups = config.min_work_groups(); work_groups <= config.max_work_groups(); work_groups *= 2){           //iterate over number of work groups (compute units)
          for (unsigned int local_workers = config.min_local_size(); local_workers <= config.max_local_size(); local_workers *= 2){   //iterate over local thread number
              for(unsigned int alignment = config.min_alignment() ; alignment <= config.max_alignment() ; alignment *= 2){
                  std::string kernel_name(program_name + to_string(work_groups)+to_string(local_workers)+to_string(alignment));
                  code_generation::kernel_infos_t & kernel_infos = kernels_infos[kernel_name];
                  kernel_infos.profile().local_work_size(0,local_workers);
                  kernel_infos.profile().global_work_size(0,local_workers*work_groups);
                  kernel_infos.profile().alignment(alignment);
                  code_generation::kernel_generator kg(*it,kernel_name,kss,kernel_infos);
                  kg.generate();
              }
          }
        }
        viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(oss.str(),"benchmark_test");
        for(std::map<std::string, generator::code_generation::kernel_infos_t>::iterator it = kernels_infos.begin() ; it != kernels_infos.end() ; ++it){
            viennacl::ocl::kernel & k = prog.add_kernel(it->first);
            k.global_work_size(0, it->second.profile().global_work_size(0));
            k.global_work_size(1, it->second.profile().global_work_size(1));
            k.local_work_size(0, it->second.profile().local_work_size(0));
            k.local_work_size(1, it->second.profile().local_work_size(1));
            set_arguments(k,it->second.arguments());
            viennacl::ocl::enqueue(k);
            viennacl::ocl::get_queue().finish();
            t.start();
            for(unsigned int i = 0 ; i < BENCHMARK_RUNS ; ++i){
                viennacl::ocl::enqueue(k);
                viennacl::ocl::get_queue().finish();
            }
            double execution_time = t.get();
            timings[execution_time] = it->second.profile();
        }

        for(std::map<double, viennacl::generator::code_generation::optimization_profile>::iterator it = timings.begin() ; it != timings.end() ; ++it){
            std::cout << it->first << std::endl;
        }
        code_generation::optimization_profile & best_profile = timings.begin()->second;
        std::cout << " Best : ["
                  << "Alignment :" << best_profile.alignment()
                  << ", Local Sork Size 0 :" << best_profile.local_work_size(0)
                  << ", Global Work Size 0 :" << best_profile.global_work_size(0)
                  << "]" << std::endl;
    }


}

template<class TestConfig>
void benchmark_timings(custom_operation& op, TestConfig & config){
    benchmark_timings(op.kernels_list(), config);
}

}

}

}
#endif // AUTOTUNE_HPP
