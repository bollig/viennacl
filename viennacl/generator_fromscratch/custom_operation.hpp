#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP


#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"
#include "viennacl/tools/shared_ptr.hpp"

namespace viennacl
{
  namespace generator
  {
  /** @brief A class for making a custom operation */
      class custom_operation
      {
        public :

          custom_operation(std::string const & operation_name) : operation_name_(operation_name){ }

          template<class T>
          void add(T const & op){
              operations_manager_.add(op);
          }

          std::string generate() const{
            std::ostringstream oss;
            typedef std::vector<std::list<infos_base*> >  kernels_t;
            kernels_t kernels(operations_manager_.get_kernels_list());
            for(kernels_t::const_iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
                code_generation::kernel_generator kg(*it,operation_name_+to_string(it-kernels.begin()));
                kg.generate(oss);
            }
            return oss.str();
          }

        private:
          std::string operation_name_;
          code_generation::operations_manager operations_manager_;
     };
  }
}
#endif // CUSTOM_OPERATION_HPP
