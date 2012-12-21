#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP

#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"
#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/tools/shared_ptr.hpp"

namespace viennacl
{
  namespace generator
  {


  /** @brief A class for making a custom operation */
      class custom_operation
      {
        friend class kernel_argument;

        public :

        custom_operation(std::string const & operation_name) : operation_name_(operation_name){ }

          template<class T>
          void add(T const & op){
              operations_manager_.add(dummy2exptree(shared_infos_,temporaries_,op));
          }

          std::string generate() const{

            std::ostringstream oss;
            code_generation::utils::kernel_generation_stream kss(oss);
            typedef std::vector<std::list<infos_base*> >  kernels_t;
            kernels_t kernels(operations_manager_.get_kernels_list());
            for(kernels_t::const_iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
                std::string name = operation_name_+to_string(it-kernels.begin());
                code_generation::kernel_generator kg(*it,name,kss);
                kg.generate() ;
            }
            return oss.str();
          }

        private:
          std::string operation_name_;
          code_generation::operations_manager operations_manager_;
          shared_infos_map_t shared_infos_;
          temporaries_map_t temporaries_;
     };
  }
}
#endif // CUSTOM_OPERATION_HPP
