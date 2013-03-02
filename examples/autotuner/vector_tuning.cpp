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

#define N_RUNS 3

typedef float ScalarType;
typedef std::vector< viennacl::ocl::platform > platforms_type;
typedef std::vector<viennacl::ocl::device> devices_type;
typedef std::vector<cl_device_id> cl_devices_type;

static const unsigned int size = 1024*1024;

//class config
//{
//   public:
//     config() {}
//     config(viennacl::ocl::device const & dev){
//        max_local_size_ = dev.max_work_group_size();

//        min_unroll_ = 1;
//        max_unroll_ = 32;

//        // GPU specific test setup:
//        if (dev.type() == CL_DEVICE_TYPE_GPU)
//        {
//            unsigned int units = 1;
//            do
//              units *= 2;
//            while (2 * units < dev.compute_units());
//            min_work_groups_ = units;
//            max_work_groups_ = 256; //reasonable upper limit on current GPUs
//            min_local_size_ = 32; //less than 32 threads per work group is unlikely to have any impact

//        }
//        else if (dev.type() == CL_DEVICE_TYPE_CPU)// CPU specific test setup
//        {
//            min_work_groups_ = 1;
//            max_work_groups_ = 2*dev.compute_units(); //reasonable upper limit on current CPUs - more experience needed here!
//            min_local_size_ = 1;
//        }
//        else
//        {
//            std::cerr << "Unknown device type (neither CPU nor GPU)! Aborting..." << std::endl;
//            exit(0);
//        }
//        cl_uint vector_width_char, vector_width_short, vector_width_int, vector_width_long, vector_width_float, vector_width_double, vector_width_half;
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,sizeof(cl_uint),(void*)&vector_width_char,NULL);
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,sizeof(cl_uint),(void*)&vector_width_short,NULL);
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,sizeof(cl_uint),(void*)&vector_width_int,NULL);
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,sizeof(cl_uint),(void*)&vector_width_long,NULL);
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(cl_uint),(void*)&vector_width_float,NULL);
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),(void*)&vector_width_double,NULL);
//        clGetDeviceInfo(dev.id(),CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF,sizeof(cl_uint),(void*)&vector_width_half,NULL);

//        min_alignments_["char"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_char/2));
//        max_alignments_["char"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_char*2));

//        min_alignments_["short"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_short/2));
//        max_alignments_["short"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_short*2));

//        min_alignments_["int"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_int/2));
//        max_alignments_["int"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_int*2));

//        min_alignments_["long"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_long/2));
//        max_alignments_["long"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_long*2));

//        min_alignments_["float"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_float/2));
//        max_alignments_["float"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_float*2));

//        min_alignments_["double"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_double/2));
//        max_alignments_["double"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_double*2));

//        min_alignments_["half"] = static_cast<unsigned int>(std::max(cl_uint(1),vector_width_half/2));
//        max_alignments_["half"] = static_cast<unsigned int>(std::min(cl_uint(16),vector_width_half*2));
//    }

//    unsigned int min_work_groups() const { return min_work_groups_; }
//    unsigned int max_work_groups() const { return max_work_groups_; }

//    unsigned int min_unroll() const { return min_unroll_; }
//    unsigned int max_unroll() const { return max_unroll_; }

//    unsigned int min_local_size() const { return min_local_size_; }
//    unsigned int max_local_size() const { return max_local_size_; }

//    unsigned int min_alignment(std::string const & scalartype) const { return min_alignments_.at(scalartype); }
//    unsigned int max_alignment(std::string const & scalartype) const { return max_alignments_.at(scalartype); }

//  private:
//    unsigned int min_work_groups_;
//    unsigned int max_work_groups_;
//    unsigned int min_local_size_;
//    unsigned int max_local_size_;
//    unsigned int min_unroll_;
//    unsigned int max_unroll_;
//    std::map<std::string, unsigned int> min_alignments_;
//    std::map<std::string, unsigned int> max_alignments_;
//};

template<class TimingsT, class OpT>
void fill_timings(TimingsT & timings, OpT const & op_template){
    viennacl::ocl::device const & dev = viennacl::ocl::current_device();
    cl_device_id id = dev.id();

    unsigned int min_unroll = 1;
    unsigned int max_unroll = 64;

    unsigned int min_alignment=1;
    unsigned int max_alignment=16;

    unsigned int min_local_size = 64;
    unsigned int max_local_size = viennacl::ocl::info<CL_DEVICE_MAX_WORK_GROUP_SIZE>(id);

//    unsigned int min_work_groups;
//    unsigned int max_work_groups;

//    // GPU specific test setup:
//    if (dev.type() == CL_DEVICE_TYPE_GPU)
//    {
//        unsigned int units = 1;
//        do
//          units *= 2;
//        while (2 * units < dev.compute_units());
//        min_work_groups = units;
//        max_work_groups = 256; //reasonable upper limit on current GPUs
//        min_local_size = 128; //less than 128 threads per work group is unlikely to have any impact

//    }
//    else if (dev.type() == CL_DEVICE_TYPE_CPU)// CPU specific test setup
//    {
//        min_work_groups = 1;
//        max_work_groups = 2*dev.compute_units(); //reasonable upper limit on current CPUs - more experience needed here!
//        min_local_size = 1;
//    }
//    else
//    {
//        std::cerr << "Unknown device type (neither CPU nor GPU)! Aborting..." << std::endl;
//        exit(0);
//    }


    viennacl::generator::autotune::Timer tim;

    for(unsigned int u = min_unroll ; u <= max_unroll ; u*=2){
        for(unsigned int a = min_alignment ; a <= max_alignment ; a*=2){
            for(unsigned int lsize = min_local_size ; lsize <= max_local_size ; lsize *= 2){
                    viennacl::generator::custom_operation op;
                    viennacl::generator::code_generation::blas1_optimization_profile prof(a,u,lsize);
                    op.operations_manager().override_blas1_model(prof);
                    op.add(op_template);
                    op.execute();
                    viennacl::ocl::get_queue().finish();
                    double total_time = 0;
                    for(unsigned int n = 0 ; n < N_RUNS ; ++n){
                        tim.start();
                        op.execute();
                        viennacl::ocl::get_queue().finish();
                        total_time+=tim.get();
                    }
                    total_time /= N_RUNS;
//                    std::cout << total_time << " " << prof << std::endl;
                    std::cout << "." << std::flush;
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
