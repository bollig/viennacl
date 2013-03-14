//~ #define VIENNACL_DEBUG_BUILD
#define VIENNACL_WITH_OPENCL
//#define VIENNACL_DEBUG_ALL

#include <iostream>
#include "CL/cl.hpp"

#include "viennacl/vector.hpp"
#include "viennacl/generator/custom_operation.hpp"
#include "viennacl/generator/dummy_types.hpp"
#include "viennacl/generator/autotune/autotune.hpp"
#include "viennacl/linalg/norm_2.hpp"

#define N_RUNS 5

typedef float ScalarType;
typedef std::vector< viennacl::ocl::platform > platforms_type;
typedef std::vector<viennacl::ocl::device> devices_type;
typedef std::vector<cl_device_id> cl_devices_type;

static const unsigned int size = 5120;


template<class TimingsT, class OpT>
void fill_timings(TimingsT & timings, OpT const & op_template){
    viennacl::ocl::device const & dev = viennacl::ocl::current_device();
    cl_device_id id = dev.id();

    unsigned int min_unroll = 1;
    unsigned int max_unroll = 16;

    unsigned int min_alignment=1;
    unsigned int max_alignment=4;

    unsigned int min_local_size = 64;
    unsigned int max_local_size = viennacl::ocl::info<CL_DEVICE_MAX_WORK_GROUP_SIZE>(id);


    viennacl::generator::autotune::Timer tim;

    for(unsigned int u = max_unroll ; u >= min_unroll ; u/=2){
        for(unsigned int a = min_alignment ; a <= max_alignment ; a*=2){
            for(unsigned int lsize = min_local_size ; lsize <= max_local_size ; lsize *= 2){
                    viennacl::generator::custom_operation op;
                    viennacl::generator::code_generation::blas1_optimization_profile prof(a,u,lsize);
                    op.operations_manager().override_blas1_model(prof);
                    op.add(op_template);
                    op.execute();
                    viennacl::backend::finish();
                    double total_time = 0;
                    for(unsigned int n = 0 ; n < N_RUNS ; ++n){
                        tim.start();
                        op.execute();
                        viennacl::backend::finish();
                        total_time+=tim.get();
                    }
                    total_time /= N_RUNS;
//                    std::cout << total_time << " " << prof << std::endl;
//                    std::cout << "." << std::flush;
                    timings.insert(std::make_pair(total_time,prof));

            }
        }
    }
}

void autotune(){
    std::map<double, viennacl::generator::code_generation::blas1_optimization_profile> timings;

    std::vector<ScalarType> cpu_v1(size), cpu_v2(size), cpu_v3(size), cpu_v4(size);
    for(unsigned int i=0; i<size; ++i){
        cpu_v1[i]=i;
        cpu_v2[i]=2*i;
        cpu_v3[i]=3*i;
        cpu_v4[i]=4*i;
    }

    viennacl::vector<ScalarType> res(size), v1(size), v2(size), v3(size), v4(size);
    viennacl::copy(cpu_v1,v1);
    viennacl::copy(cpu_v2,v2);
    viennacl::copy(cpu_v3,v3);
    viennacl::copy(cpu_v4,v4);


    typedef viennacl::generator::dummy_vector<ScalarType> dv;

    fill_timings(timings,dv(v1) = dv(v2) + dv(v3));

    for(std::map<double,viennacl::generator::code_generation::blas1_optimization_profile>::iterator it = timings.begin(); it!=timings.end(); ++it){
        unsigned int n = std::distance(timings.begin(), it);
        std::cout << n << "th Best : " << it->first << " <=> " << it->second << std::endl;
        if(n == 20) break;
    }
}

int main(){
    platforms_type platforms = viennacl::ocl::get_platforms();
    size_t num_platforms = platforms.size();
	for(unsigned int k=0 ; k < num_platforms ; ++k)
	{
		viennacl::ocl::platform pf(k);
		viennacl::ocl::set_context_platform_index(k,k);
		viennacl::ocl::switch_context(k);
		devices_type dev = viennacl::ocl::current_context().devices();
		for(devices_type::iterator it = dev.begin() ; it != dev.end() ; ++it){
			viennacl::ocl::switch_device(*it);
			std::cout << "-------------------" << std::endl;
            std::cout << "Recording timings for : " << viennacl::ocl::current_device().name() << std::endl;
            autotune();
		}
	}
	
	
}
