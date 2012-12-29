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
    config(viennacl::ocl::device const & dev) : device_(dev){
        max_local_size_ = dev.max_work_group_size();


        // GPU specific test setup:
        if (dev.type() == CL_DEVICE_TYPE_GPU)
        {
            min_unroll_ = 1;
            max_unroll_ = 4;
            if(dev.vendor_id()==4318){
                min_work_groups_ = 32;
            }
            else{
                unsigned int units=1;
                do  units *= 2; while (2 * units < dev.compute_units());
                min_work_groups_ = units;
            }
            max_work_groups_ = 512; //reasonable upper limit on current GPUs
            min_local_size_ = 16; //less than 16 threads per work group is unlikely to have any impact

        }
        else if (dev.type() == CL_DEVICE_TYPE_CPU)// CPU specific test setup
        {
            min_unroll_ = 4;
            max_unroll_ = 64;
            min_work_groups_ = 1;
            max_work_groups_ = 2*dev.compute_units(); //reasonable upper limit on current CPUs - more experience needed here!
            min_local_size_ = 1;
        }
        else
        {
            std::cerr << "Unknown device type (neither CPU nor GPU)! Aborting..." << std::endl;
            exit(0);
        }
        cl_uint vector_width_char, vector_width_short, vector_width_int, vector_width_long, vector_width_float, vector_width_double, vector_width_half;
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,sizeof(cl_uint),(void*)&vector_width_char,NULL);
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,sizeof(cl_uint),(void*)&vector_width_short,NULL);
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,sizeof(cl_uint),(void*)&vector_width_int,NULL);
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,sizeof(cl_uint),(void*)&vector_width_long,NULL);
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(cl_uint),(void*)&vector_width_float,NULL);
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),(void*)&vector_width_double,NULL);
        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF,sizeof(cl_uint),(void*)&vector_width_half,NULL);

        min_alignments_["char"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_char/2));
        max_alignments_["char"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_char*2));

        min_alignments_["short"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_short/2));
        max_alignments_["short"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_short*2));

        min_alignments_["int"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_int/4));
        max_alignments_["int"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_int*4));

        min_alignments_["long"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_long/4));
        max_alignments_["long"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_long*4));

        min_alignments_["float"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_float/4));
        max_alignments_["float"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_float*4));

//        min_alignments_["float"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_float/4));
//        max_alignments_["float"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_float*4));
//        std::cout <<  min_alignments_["float"] << std::endl;
//        std::cout <<  max_alignments_["float"] << std::endl;
        min_alignments_["double"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_double/4));
        max_alignments_["double"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_double*4));

        min_alignments_["half"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_half/4));
        max_alignments_["half"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_half*4));
    }

    unsigned int min_work_groups() const { return min_work_groups_; }
    unsigned int max_work_groups() const { return max_work_groups_; }

    unsigned int min_unroll() const { return min_unroll_; }
    unsigned int max_unroll() const { return max_unroll_; }

    unsigned int min_local_size() const { return min_local_size_; }
    unsigned int max_local_size() const { return max_local_size_; }

    unsigned int min_alignment(std::string const & scalartype) const { return min_alignments_.at(scalartype); }
    unsigned int max_alignment(std::string const & scalartype) const { return max_alignments_.at(scalartype); }

    viennacl::ocl::device const & device() const{
        return device_;
    }

  private:
    unsigned int min_work_groups_;
    unsigned int max_work_groups_;
    unsigned int min_local_size_;
    unsigned int max_local_size_;
    unsigned int min_unroll_;
    unsigned int max_unroll_;
    std::map<std::string, unsigned int> min_alignments_;
    std::map<std::string, unsigned int> max_alignments_;
    viennacl::ocl::device const & device_;
};

template<class TestConfig>
std::pair<double, code_generation::optimization_profile> benchmark_timings(std::vector<std::list<infos_base*> > kernels, TestConfig & config){
    for(std::vector<std::list<infos_base*> >::iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
        std::string program_name = "k"+to_string(it-kernels.begin())+"_bench_";
        std::map<double, viennacl::generator::code_generation::optimization_profile> timings;
        std::ostringstream oss;
        code_generation::utils::kernel_generation_stream kss(oss);
        kss << "#if defined(cl_khr_fp64)\n";
        kss <<  "#pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
        kss <<  "#elif defined(cl_amd_fp64)\n";
        kss <<  "#pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
        kss <<  "#endif\n";
        std::map<std::string, generator::code_generation::kernel_infos_t> kernels_infos;
        unsigned int n_kernels=0;
        for (unsigned int local_workers = config.min_local_size(); local_workers <= config.max_local_size(); local_workers *= 2){   //iterate over local thread number
            for(unsigned int alignment = config.min_alignment("float") ; alignment <= config.max_alignment("float") ; alignment *= 2){
                for(unsigned int unroll = config.min_unroll() ; unroll <= config.max_unroll() ; unroll *=2){
                    std::string kernel_name(program_name + "l"+to_string(local_workers)+"a" + to_string(alignment) + "u"+ to_string(unroll));
                    code_generation::kernel_infos_t & kernel_infos = kernels_infos[kernel_name];
                    kernel_infos.profile().local_work_size(0,local_workers);
                    kernel_infos.profile().alignment(alignment);
                    kernel_infos.profile().loop_unroll(unroll);
                    code_generation::kernel_generator kg(*it,kernel_name,kss,kernel_infos);
                    kg.generate();
                    n_kernels++;
                }
            }
        }
        std::cout << "Compiling " << n_kernels << " kernels" << std::endl;
        viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(oss.str(),"benchmark_test");
        std::cout << "Compiled !" << std::endl;
        for(std::map<std::string, generator::code_generation::kernel_infos_t>::iterator it = kernels_infos.begin() ; it != kernels_infos.end() ; ++it){
            viennacl::ocl::kernel & k = prog.add_kernel(it->first);
            code_generation::optimization_profile& profile = it->second.profile();

            //Anticipates kernel failure
            size_t max_workgroup_size = k.max_workgroup_size(config.device().id());
            if(profile.local_work_size(0) > max_workgroup_size || profile.local_work_size(1) > max_workgroup_size)
                continue;

            //Doesn't execute because it would likelily be a waste of time
            size_t prefered_workgroup_size_multiple = k.prefered_work_group_size_multiple(config.device().id());
            if(profile.local_work_size(0) < prefered_workgroup_size_multiple || ((0 < profile.local_work_size(1)) &&  (profile.local_work_size(1)< prefered_workgroup_size_multiple)))
                continue;

            set_arguments(k,it->second.arguments());
            k.local_work_size(0, it->second.profile().local_work_size(0));
            k.local_work_size(1, it->second.profile().local_work_size(1));
            std::cout << "." << std::flush ;
            for (unsigned int work_groups = config.min_work_groups(); work_groups <= config.max_work_groups(); work_groups *= 2){
                k.global_work_size(0, k.local_work_size(0)*work_groups);
                k.global_work_size(1, k.local_work_size(1)*work_groups);
                profile.global_work_size(0,k.global_work_size(0));
                profile.global_work_size(1,k.global_work_size(1));

                Timer t;
                viennacl::ocl::enqueue(k);
                viennacl::ocl::get_queue().finish();
                t.start();
                for(unsigned int i = 0 ; i < BENCHMARK_RUNS ; ++i){
                    viennacl::ocl::enqueue(k);
                    viennacl::ocl::get_queue().finish();
                }
                double execution_time = t.get();
                timings[execution_time] = profile;
            }


        }
//        for(std::map<double, viennacl::generator::code_generation::optimization_profile>::iterator it = timings.begin() ; it!= timings.end() ; ++it){
//            std::cout << it->first << " <== " << it->second << std::endl;
//        }
        std::cout << std::endl;
        return *timings.begin();
    }


}

template<class TestConfig>
std::pair<double, code_generation::optimization_profile> benchmark_timings(custom_operation& op, TestConfig & config){
    return benchmark_timings(op.kernels_list(), config);
}

}

}

}
#endif // AUTOTUNE_HPP
