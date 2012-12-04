#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP


#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/tree_utils.hpp"

namespace viennacl
{
  namespace generator
  {
  /** @brief A class for making a custom operation */
      class custom_operation
      {

          static bool test_pred(infos_base* p){
              return dynamic_cast<vec_infos_base *>(p);
          }

        public :

          /** @brief CTor for 1 expression
          *
          * @param operation_name the code for this expression will be stored in the program provided by this name
          */
          template<class T0>
          custom_operation ( T0 const & , std::string const & operation_name) : program_name_(operation_name)
          {
              infos_base & tree = T0::get();
              std::tree_utils::extract_if(&tree,test_pred);
          }

        private:
          std::string program_name_;
     };
  }
}
#endif // CUSTOM_OPERATION_HPP
