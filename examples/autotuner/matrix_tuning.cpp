//#define VIENNACL_DEBUG_BUILD
#define VIENNACL_WITH_OPENCL
//#define VIENNACL_DEBUG_ALL

#define NDEBUG

#include <iostream>
#include "CL/cl.hpp"
#include <sys/time.h>

#include "boost/numeric/ublas/matrix.hpp"

#include "viennacl/matrix.hpp"
#include "viennacl/generator_fromscratch/custom_operation.hpp"
#include "viennacl/generator_fromscratch/dummy_types.hpp"
#include "viennacl/generator_fromscratch/autotune/autotune.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"



#include "../tutorial/Random.hpp"


typedef float  NumericT;
typedef std::vector< viennacl::ocl::platform > platforms_type;
typedef std::vector<viennacl::ocl::device> devices_type;
typedef std::vector<cl_device_id> cl_devices_type;



struct config{

    unsigned int size;
    unsigned int n_runs;

    unsigned int ml_min;
    unsigned int kl_min;
    unsigned int nl_min;
    unsigned int ms_min;
    unsigned int ks_min;
    unsigned int ns_min;

    unsigned int ml_max;
    unsigned int kl_max;
    unsigned int nl_max;
    unsigned int ms_max;
    unsigned int ks_max;
    unsigned int ns_max;

    unsigned int alignment_min;
    unsigned int alignment_max;

    std::vector<bool> LHS_storages;
    std::vector<bool> RHS_storages;
};

template <typename VCLMatrixType1, typename VCLMatrixType2>
double diff(VCLMatrixType1 & mat1, VCLMatrixType2 & mat2)
{
   typedef typename VCLMatrixType1::value_type::value_type ScalarType1;
   typedef typename VCLMatrixType2::value_type::value_type ScalarType2;

   boost::numeric::ublas::matrix<ScalarType1> mat1_cpu(mat1.size1(), mat1.size2());
   boost::numeric::ublas::matrix<ScalarType2> mat2_cpu(mat2.size1(), mat2.size2());

   viennacl::backend::finish();  //workaround for a bug in APP SDK 2.7 on Trinity APUs (with Catalyst 12.8)
   viennacl::copy(mat2, mat2_cpu);
   viennacl::copy(mat1, mat1_cpu);
   viennacl::ocl::get_queue().finish();
   double ret = 0;
   double act = 0;

    for (unsigned int i = 0; i < mat2_cpu.size1(); ++i)
    {
      for (unsigned int j = 0; j < mat2_cpu.size2(); ++j)
      {
         act = fabs(mat2_cpu(i,j) - mat1_cpu(i,j)) / std::max( fabs(mat2_cpu(i, j)), fabs(mat1_cpu(i,j)) );
         if (act > ret)
           ret = act;
      }
    }
   //std::cout << ret << std::endl;
   return ret;
}

template<class RefTypeA,
         class ResTypeA,
         class OpT>
bool test_blas3(RefTypeA const & reference, ResTypeA const & dummyA, OpT const & operation, std::list<viennacl::generator::code_generation::blas3_optimization_profile> const & profiles){
    bool res = true;
    for(std::list<viennacl::generator::code_generation::blas3_optimization_profile>::const_iterator it = profiles.begin(); it!=profiles.end(); ++it){
        viennacl::generator::custom_operation op;
        op.operations_manager().blas3_model() = *it;
        op.add(operation);
        op.execute();
        viennacl::ocl::get_queue().finish();
        double delta = diff(reference,dummyA);
        if(delta>1e-5){
            std::cout << "Failed for " << *it << "| Diff : " << delta << std::endl;
//            std::cout << op.source_code() << std::endl;
            res=false;
        }
    }
    return res;
}


template<class MatTypeA, class MatTypeB, class MatTypeC>
void fill_matrix(MatTypeA & A, MatTypeB & B, MatTypeC & C){
    typedef NumericT ScalarTypeA;
    typedef NumericT ScalarTypeB;
    typedef NumericT ScalarTypeC;

    boost::numeric::ublas::matrix<ScalarTypeA> cpu_A(A.size1(),A.size2());
    boost::numeric::ublas::matrix<ScalarTypeB> cpu_B(B.size1(),B.size2());
    boost::numeric::ublas::matrix<ScalarTypeC> cpu_C(C.size1(),C.size1());

    srand(time(NULL));
    for(unsigned int i=0; i<A.size1(); ++i){
        for(unsigned int j=0 ; j<A.size2() ; ++j){
            cpu_A(i,j)=0;
            cpu_B(i,j) =static_cast<ScalarTypeB>(rand())/static_cast<ScalarTypeB>(RAND_MAX);
            cpu_C(i,j)=static_cast<ScalarTypeB>(rand())/static_cast<ScalarTypeB>(RAND_MAX);
        }
    }

    viennacl::copy(cpu_A,A);
    viennacl::copy(cpu_B,B);
    viennacl::copy(cpu_C,C);
    viennacl::ocl::get_queue().finish();
}


template<class OpT, class MatTypeA, class MatTypeB, class MatTypeC>

void benchmark(OpT const & operation, config conf, MatTypeA & A, MatTypeB & B, MatTypeC & C,
                    std::list<viennacl::generator::code_generation::blas3_optimization_profile> & fastest_firsts){
    viennacl::generator::autotune::timings_t timings;
    unsigned int size;

    std::list<std::pair<unsigned int, unsigned int> > rounds_config;
    rounds_config.push_back(std::make_pair(512 ,200));
    rounds_config.push_back(std::make_pair(2048,20));
    for(std::list<std::pair<unsigned int, unsigned int> >::iterator it = rounds_config.begin() ; it!= rounds_config.end(); ++it){
        unsigned int k = std::distance(rounds_config.begin(),it);
        timings.clear();
        size=it->first;
        unsigned int n_keep=it->second;
        A.resize(size,size,false);
        B.resize(size,size,false);
        C.resize(size,size,false);
        viennacl::ocl::get_queue().finish();
        fill_matrix(A,B,C);
        viennacl::ocl::get_queue().finish();
        if(k==0)
            viennacl::generator::autotune::benchmark_blas3(timings,operation,conf);
        else{
            viennacl::generator::autotune::benchmark_blas3(timings,operation,fastest_firsts);
        }
        fastest_firsts.clear();
        viennacl::ocl::get_queue().finish();
        for(viennacl::generator::autotune::timings_t::iterator itt = timings.begin(); itt!=timings.end() ; ++itt){
            unsigned int n = std::distance(timings.begin(),itt);
            if(n>n_keep) break;
            fastest_firsts.push_back(*static_cast<viennacl::generator::code_generation::blas3_optimization_profile* >(itt->second.get()));
            if(std::distance(rounds_config.begin(),it)==(int)rounds_config.size()-1){
                std::cout << std::distance(timings.begin(),itt) << "th Best : " << itt->first << "s | " << std::pow((NumericT)size/1000,3)/itt->first << " GFlops : " << *itt->second << std::endl;
            }
        }
    }
}

template<class ScalarTypeA, class LayoutA
         , class ScalarTypeB, class LayoutB
         , class ScalarTypeC, class LayoutC>
void run_autotune(){

    using viennacl::generator::prod;




    typedef viennacl::generator::dummy_matrix<ScalarTypeA, LayoutA> dma_t;
    typedef viennacl::generator::dummy_matrix<ScalarTypeB, LayoutB> dmb_t;
    typedef viennacl::generator::dummy_matrix<ScalarTypeC, LayoutC> dmc_t;

    config conf;


    conf.n_runs = 2;
    conf.ml_min = 64; conf.ml_max=64;
    conf.kl_min = 256; conf.kl_max=256;
    conf.nl_min = 64; conf.nl_max=64;
    conf.ms_min = 4; conf.ms_max=4;
    conf.ks_min = 4; conf.ks_max=4;
    conf.ns_min = 8; conf.ns_max=8 ;
    conf.alignment_min = 4 ; conf.alignment_max = 4 ;
//    conf.LHS_storages.push_back(true);
    conf.LHS_storages.push_back(false);
//    conf.RHS_storages.push_back(true);
    conf.RHS_storages.push_back(false);

    viennacl::matrix<ScalarTypeA,LayoutA> A;
    viennacl::matrix<ScalarTypeB,LayoutB> B;
    viennacl::matrix<ScalarTypeC,LayoutC> C;
    std::list<viennacl::generator::code_generation::blas3_optimization_profile> fastest_firsts;

//    std::cout << "------------AA------------" << std::endl;
//    benchmark(dma_t(A) = prod(dmb_t(B),dmc_t(C)),conf,A,B,C,fastest_firsts);
//    std::cout << "Testing " << fastest_firsts.size() << " best configurations" << std::endl;
//    if(!test_blas3(viennacl::matrix<ScalarTypeA,LayoutA>(viennacl::linalg::prod(B,C)),A,dma_t(A) = prod(dmb_t(B),dmc_t(C)),fastest_firsts)){
//        std::cout << "#Fail" << std::endl;
//    }


    std::cout << "------------TA------------" << std::endl;
    benchmark(dma_t(A) = prod(trans(dmb_t(B)),dmc_t(C)),conf,A,B,C,fastest_firsts);
    std::cout << "Testing " << fastest_firsts.size() << " best configurations" << std::endl;

    viennacl::ocl::get_queue().finish();
    if(!test_blas3(viennacl::matrix<ScalarTypeA,LayoutA>(viennacl::linalg::prod(trans(B),C)),A,dma_t(A) = prod(trans(dmb_t(B)),dmc_t(C)),fastest_firsts)){
        std::cout << "#Fail" << std::endl;
    }

//     std::cout << "------------AT------------" << std::endl;
//     benchmark(dma_t(A) = prod(dmb_t(B),trans(dmc_t(C))),conf,A,B,C,fastest_firsts);
//     std::cout << "Testing " << fastest_firsts.size() << " best configurations" << std::endl;

//     viennacl::ocl::get_queue().finish();
//     if(!test_blas3(viennacl::matrix<ScalarTypeA,LayoutA>(viennacl::linalg::prod(B,trans(C))),A,dma_t(A) = prod(dmb_t(B),trans(dmc_t(C))),fastest_firsts)){
//         std::cout << "#Fail" << std::endl;
//     }

//     std::cout << "------------TT------------" << std::endl;
//     benchmark(dma_t(A) = prod(trans(dmb_t(B)),trans(dmc_t(C))),conf,A,B,C,fastest_firsts);
//     std::cout << "Testing " << fastest_firsts.size() << " best configurations" << std::endl;

//     viennacl::ocl::get_queue().finish();
//     if(!test_blas3(viennacl::matrix<ScalarTypeA,LayoutA>(viennacl::linalg::prod(trans(B),trans(C))),A,dma_t(A) = prod(trans(dmb_t(B)),trans(dmc_t(C))),fastest_firsts)){
//         std::cout << "#Fail" << std::endl;
//     }
}

int main(int argc, char* argv[]){
    std::vector<std::string> args(argv,argv+argc);

    if(argc < 2){
        std::cerr << "Usage : PROGRAM_NAME LAYOUT" << std::endl;
        exit(-1);
    }
    int layout = atoi(args[1].c_str());

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
            if(viennacl::ocl::current_device().type()==CL_DEVICE_TYPE_GPU){
                std::cout << "-------------------" << std::endl;
                std::cout << "Recording timings for : " << viennacl::ocl::current_device().name() << std::endl;

                switch(layout){
                case 0 :
                    std::cout << "====== Row-Major = Row-Major * Row-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::row_major
                             ,NumericT,viennacl::row_major
                             ,NumericT,viennacl::row_major >();
                    break;

                case 1:
                    std::cout << "====== Row-Major = Row-Major * Column-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::row_major
                                 ,NumericT,viennacl::row_major
                                 ,NumericT,viennacl::column_major >();
                    break;

                case 2:
                    std::cout << "====== Column-Major = Row-Major * Column-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::column_major
                                 ,NumericT,viennacl::row_major
                                 ,NumericT,viennacl::column_major >();
                    break;

                case 3:
                    std::cout << "====== Column-Major = Row-Major * Row-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::column_major
                                 ,NumericT,viennacl::row_major
                                 ,NumericT,viennacl::row_major >();
                    break;

                case 4:
                    std::cout << "====== Row-Major = Column-Major * Row-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::row_major
                                 ,NumericT,viennacl::column_major
                                 ,NumericT,viennacl::row_major >();
                    break;

                case 5:
                    std::cout << "====== Row-Major = Column-Major * Column-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::row_major
                                 ,NumericT,viennacl::column_major
                                 ,NumericT,viennacl::column_major >();
                    break;

                case 6:
                    std::cout << "====== Column-Major = Column-Major * Row-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::column_major
                                 ,NumericT,viennacl::column_major
                                 ,NumericT,viennacl::row_major >();
                    break;


                case 7:
                    std::cout << "====== Column-Major = Column-Major * Column-Major ======" << std::endl;
                    run_autotune< NumericT,viennacl::column_major
                                 ,NumericT,viennacl::column_major
                                 ,NumericT,viennacl::column_major >();
                    break;
                }












            }

        }
    }


}
