//~ #define VIENNACL_DEBUG_BUILD
#define VIENNACL_WITH_OPENCL

#include <iostream>
#include "CL/cl.hpp"

#include "viennacl/vector.hpp"
#include "viennacl/generator_fromscratch/custom_operation.hpp"
#include "viennacl/generator_fromscratch/dummy_types.hpp"
#include "viennacl/generator_fromscratch/autotune/autotune.hpp"
#include "viennacl/linalg/norm_2.hpp"

typedef float ScalarType;
typedef std::vector< viennacl::ocl::platform > platforms_type;
typedef std::vector<viennacl::ocl::device> devices_type;
typedef std::vector<cl_device_id> cl_devices_type;

static const unsigned int size = 10;

std::pair<double, viennacl::generator::code_generation::optimization_profile> autotune(){
std::vector<ScalarType> cpu_v1(size), cpu_v2(size), cpu_v3(size), cpu_v4(size);
for(unsigned int i=0; i<size; ++i){
	cpu_v1[i]=i;
	cpu_v2[i]=2*i;
	cpu_v3[i]=3*i;
	cpu_v4[i]=4*i;
}
	
viennacl::generator::custom_operation op("test");
viennacl::vector<ScalarType,16> res(size), v1(size), v2(size), v3(size), v4(size);
viennacl::copy(cpu_v1,v1);
viennacl::copy(cpu_v2,v2);
viennacl::copy(cpu_v3,v3);
viennacl::copy(cpu_v4,v4);

typedef viennacl::generator::dummy_vector<ScalarType,16> dv; 
op.add(dv(v1) = inner_prod(dv(v2),dv(v3)) );
op.init();
op.execute();
viennacl::ocl::get_queue().finish();
std::cout << v1[0] << std::endl;
std::cout << viennacl::linald::inner_prod(v2,v3) << std::endl;
viennacl::generator::autotune::config conf(viennacl::ocl::current_device());
return viennacl::generator::autotune::benchmark_timings(op,conf);
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
			std::pair<double, viennacl::generator::code_generation::optimization_profile> best = autotune();
			std::cout << "Best : " << best.first << " <=> " << best.second << std::endl;
			
		}
	}
	
	
}
