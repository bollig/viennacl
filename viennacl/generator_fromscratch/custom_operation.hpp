#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP

#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"
#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/tools/shared_ptr.hpp"

namespace viennacl
{
  namespace generator
  {



      typedef std::map<viennacl::backend::mem_handle, shared_infos_t> shared_infos_map_t;
      typedef std::map<kernel_argument*,viennacl::backend::mem_handle,deref_less> temporaries_map_t;

      template<class T>
      struct dummy2exptree_impl
      {
      private:
        typedef typename T::Lhs Lhs;
        typedef typename T::Rhs Rhs;
        typedef typename dummy2exptree_impl<Lhs>::result_type LhsResult;
        typedef typename dummy2exptree_impl<Rhs>::result_type RhsResult;
      public:
        typedef typename get_symbolic_type<T,LhsResult,RhsResult>::type result_type;

          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     T const & t){
              return result_type(dummy2exptree_impl<Lhs>::execute(shared_infos, temporaries, t.lhs())
                                 ,dummy2exptree_impl<Rhs>::execute(shared_infos, temporaries, t.rhs()));
          }
      };

      template<class ScalarType, unsigned int Alignment>
      struct dummy2exptree_impl<dummy_vector<ScalarType, Alignment> >{
          typedef symbolic_vector<ScalarType, Alignment> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries_,
                                     dummy_vector<ScalarType,Alignment> const & v){
              return result_type(shared_infos, v.vec());
          }
      };

      template<class LHS, class RHS>
      struct dummy2exptree_impl<inner_prod_wrapper<LHS,RHS> >{
      private:
          typedef typename dummy2exptree_impl<LHS>::result_type LhsResult;
          typedef typename dummy2exptree_impl<RHS>::result_type RhsResult;
      public:
          typedef inprod_infos<LhsResult,RhsResult> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     inner_prod_wrapper<LHS,RHS> const & v){
              return result_type(shared_infos, temporaries,
                                 dummy2exptree_impl<LHS>::execute(shared_infos,temporaries,v.lhs()),
                                 dummy2exptree_impl<RHS>::execute(shared_infos,temporaries,v.rhs()));
          }
      };




      template<class T1, class T2, class T3, class T4, class T5>
      struct dummy2exptree_impl<function_wrapper_impl<T1,T2,T3,T4,T5> >{
      private:
          template<class U>
          static void handle_function_arg(symbolic_function & fun, U const * t, std::string name
                              , shared_infos_map_t & shared_infos
                              , temporaries_map_t & temporaries)
          {

              fun.add_arg(name, dummy2exptree_impl<U>::execute(shared_infos,temporaries,*t));
          }

          static void handle_function_arg(symbolic_function & fun, void const* t, std::string name
                              , shared_infos_map_t & shared_infos
                              , temporaries_map_t & temporaries)
          { }

      public:
          typedef symbolic_function result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     function_wrapper_impl<T1,T2,T3,T4,T5> func){
              result_type res(func.name,func.expr);
              handle_function_arg(res,func.t1,"_1_",shared_infos,temporaries);
              handle_function_arg(res,func.t2,"_2_",shared_infos,temporaries);
              handle_function_arg(res,func.t3,"_3_",shared_infos,temporaries);
              handle_function_arg(res,func.t4,"_4_",shared_infos,temporaries);
              handle_function_arg(res,func.t5,"_5_",shared_infos,temporaries);
              return res;


          }
      };

      template<class T>
      typename dummy2exptree_impl<T>::result_type dummy2exptree(shared_infos_map_t & shared_infos,
                                                                temporaries_map_t & temporaries,
                                                                T const & t){
          return dummy2exptree_impl<T>::execute(shared_infos,temporaries,t);
      }

  /** @brief A class for making a custom operation */
      class custom_operation
      {

      private:
          void compile_program() const{
              assert(!source_code_.empty() && " Custom Operation not initialized ");
              viennacl::ocl::program& program = viennacl::ocl::current_context().add_program(source_code_, operation_name_);
              for(std::map<std::string, generator::code_generation::kernel_infos_t>::const_iterator it = kernels_infos_.begin() ; it !=kernels_infos_.end() ; ++it){
                program.add_kernel(it->first);
              }
          }

      public :

        custom_operation(std::string const & operation_name) : operation_name_(operation_name){ }

          template<class T>
          void add(T const & op){
              operations_manager_.add(dummy2exptree(shared_infos_,temporaries_,op));
          }

          std::vector<std::list<infos_base*> > kernels_list(){
              return operations_manager_.get_kernels_list();
          }

          void init() {
            std::ostringstream oss;
            code_generation::utils::kernel_generation_stream kss(oss);
            kss << "#if defined(cl_khr_fp64)\n";
            kss <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
            kss <<  "#elif defined(cl_amd_fp64)\n";
            kss <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
            kss <<  "#endif\n";

            typedef std::vector<std::list<infos_base*> >  kernels_t;
            kernels_t kernels(operations_manager_.get_kernels_list());
            for(kernels_t::iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
                std::string name = operation_name_+to_string(it-kernels.begin());
                code_generation::kernel_generator kg(*it,name,kss, kernels_infos_[name]);
                kg.generate() ;
            }
            source_code_ = oss.str();
          }

          void execute(){
              compile_program();
              for(std::map<std::string, generator::code_generation::kernel_infos_t>::iterator it = kernels_infos_.begin() ; it != kernels_infos_.end() ; ++it){
                  viennacl::ocl::kernel& k = viennacl::ocl::get_kernel(operation_name_,it->first);
                  set_arguments(k,it->second.arguments());
                  viennacl::ocl::enqueue(k);
              }
          }


          std::string source_code() const{
              return source_code_;
          }

        private:
          std::string operation_name_;
          code_generation::operations_manager operations_manager_;
          shared_infos_map_t shared_infos_;
          temporaries_map_t temporaries_;
          std::map<std::string, generator::code_generation::kernel_infos_t> kernels_infos_;
          std::string source_code_;

     };


//      template <typename TimingType, typename OP, typename TestConfig, typename TestData>
//      void record_full_timings(TimingType & timings,
//                               OP operation,
//                               TestConfig & config)
//      {
//        typedef typename TestData::value_type  ScalarType;

//        double result = 0;
//        functor(data); //startup run (ensures kernel compilation)
//        for (unsigned int work_groups = config.min_work_groups(); work_groups <= config.max_work_groups(); work_groups *= 2)           //iterate over number of work groups (compute units)
//        {
//          for (unsigned int local_workers = config.min_local_size(); local_workers <= config.max_local_size(); local_workers *= 2)   //iterate over local thread number
//          {
//              for(unsigned int vectorize = config.min_vectorize() ; vectorize <= config.max_vectorize() ; vectorize *= 2){
//                  //set parameter:
//                  generate_kernel(work_groups, local_workers, vectorize);

//                  //std::cout << "Benchmarking kernel " << config.kernel_name() << std::endl;
//                  result = execute(functor, data);

//                  //check for valid result: (kernels have an automatic fallback to smaller values included)
//                  if (!validate_result(config.program_name(), config.kernel_name(), work_groups, local_workers))
//                  {
//                  std::cout << "Kernel start failed for kernel " << config.kernel_name() << " [" << work_groups << " groups, " << local_workers << " per group]" << std::endl;
//                    break;
//                  }
//                  else
//                    timings[result] = std::make_pair(work_groups * local_workers, local_workers);
//              }

//          }
//        }
//      }

  }
}
#endif // CUSTOM_OPERATION_HPP
