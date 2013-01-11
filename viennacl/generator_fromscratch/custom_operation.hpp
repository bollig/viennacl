#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP

#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"
#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/tools/shared_ptr.hpp"
#include <bitset>

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

      template<class ScalarType>
      struct dummy2exptree_impl<dummy_vector<ScalarType> >{
          typedef symbolic_vector<ScalarType, 16> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries_,
                                     dummy_vector<ScalarType> const & v){
              return result_type(shared_infos, v.vec());
          }
      };

      template<class ScalarType, class Layout>
      struct dummy2exptree_impl<dummy_matrix<ScalarType, Layout> >{
          typedef symbolic_matrix<ScalarType,Layout, 16> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries_,
                                     dummy_matrix<ScalarType,Layout> const & m){
              return result_type(shared_infos, m.mat());
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

      template<class LHS, class RHS>
      struct dummy2exptree_impl<matmat_prod_wrapper<LHS,RHS> >{
      private:
          typedef typename dummy2exptree_impl<LHS>::result_type LhsResult;
          typedef typename dummy2exptree_impl<RHS>::result_type RhsResult;
      public:
          typedef inprod_infos<LhsResult,RhsResult> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     matmat_prod_wrapper<LHS,RHS> const & v){
              return result_type(shared_infos, temporaries,
                                 dummy2exptree_impl<LHS>::execute(shared_infos,temporaries,v.lhs()),
                                 dummy2exptree_impl<RHS>::execute(shared_infos,temporaries,v.rhs()));
          }
      };



      template<class T1>
      struct dummy2exptree_impl<function_wrapper_impl<T1> >{
      private:
          typedef symbolic_function<typename dummy2exptree_impl<T1>::result_type> result_type;
      public:
          static result_type execute(shared_infos_map_t & shared_infos,temporaries_map_t & temporaries,function_wrapper_impl<T1> func){
              result_type res(func.name,func.expr);
              res.add_arg("_1_", dummy2exptree_impl<T1>::execute(shared_infos,temporaries,*func.t1));
              return res;
          }
      };

      template<class T1, class T2>
      struct dummy2exptree_impl<function_wrapper_impl<T1,T2> >{
      private:
          typedef symbolic_function<typename dummy2exptree_impl<T1>::result_type
                                    ,typename dummy2exptree_impl<T2>::result_type> result_type;
      public:
          static result_type execute(shared_infos_map_t & shared_infos,temporaries_map_t & temporaries,function_wrapper_impl<T1,T2> func){
              result_type res(func.name,func.expr);
              res.add_arg("_1_", dummy2exptree_impl<T1>::execute(shared_infos,temporaries,*func.t1));
              res.add_arg("_2_", dummy2exptree_impl<T2>::execute(shared_infos,temporaries,*func.t2));
              return res;
          }
      };

      template<class T1, class T2, class T3>
      struct dummy2exptree_impl<function_wrapper_impl<T1,T2, T3> >{
      private:
          typedef symbolic_function<typename dummy2exptree_impl<T1>::result_type
                                    ,typename dummy2exptree_impl<T2>::result_type
                                    ,typename dummy2exptree_impl<T3>::result_type> result_type;

      public:
          static result_type execute(shared_infos_map_t & shared_infos,temporaries_map_t & temporaries,function_wrapper_impl<T1,T2,T3> func){
              result_type res(func.name,func.expr);
              res.add_arg("_1_", dummy2exptree_impl<T1>::execute(shared_infos,temporaries,*func.t1));
              res.add_arg("_2_", dummy2exptree_impl<T2>::execute(shared_infos,temporaries,*func.t2));
              res.add_arg("_3_", dummy2exptree_impl<T3>::execute(shared_infos,temporaries,*func.t3));
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
              viennacl::ocl::program& program = viennacl::ocl::current_context().add_program(source_code_, operations_manager_.repr());
              for(std::map<std::string, generator::code_generation::kernel_infos_t>::const_iterator it = kernels_infos_.begin() ; it !=kernels_infos_.end() ; ++it){
                program.add_kernel(it->first);
              }
          }

      public :

        custom_operation() { }

          template<class T>
          void add(T const & op){
//              std::cout << std::bitset<64>(get_operation_id<T>::value(op)).to_string() << std::endl;
//              std::cout << encode_to_kernel_name(get_operation_id<T>::value(op)) << std::endl;
              operations_manager_.add(dummy2exptree(shared_infos_,temporaries_,op));
          }

          std::vector<code_generation::kernel_representation_t> kernels_list(){
              return operations_manager_.get_kernels_list();
          }

          void init() {
              source_code_ = operations_manager_.get_source_code(kernels_infos_);
          }

          void execute(){
              if(!viennacl::ocl::current_context().has_program(operations_manager_.repr()));
                compile_program();
              viennacl::ocl::program & pgm = viennacl::ocl::current_context().get_program(operations_manager_.repr());

              for(std::map<std::string, generator::code_generation::kernel_infos_t>::iterator it = kernels_infos_.begin() ; it != kernels_infos_.end() ; ++it){
                  viennacl::ocl::kernel& k = pgm.get_kernel(it->first);
                  set_arguments(k,it->second.arguments());
                  k.local_work_size(0,it->second.profile().local_work_size(0));
                  k.local_work_size(1,it->second.profile().local_work_size(1));

                  k.global_work_size(0,it->second.profile().global_work_size(0));
                  k.global_work_size(1,it->second.profile().global_work_size(1));

                  viennacl::ocl::enqueue(k);
              }
          }


          std::string source_code() const{
              return source_code_;
          }

        private:
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
