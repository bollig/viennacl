#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP


#include "viennacl/generator_fromscratch/symbolic_types.hpp"

namespace viennacl
{
  namespace generator
  {
  /** @brief A class for making a custom operation */
      class custom_operation
      {

        public :

          /** @brief CTor for 1 expression
          *
          * @param operation_name the code for this expression will be stored in the program provided by this name
          */
          template<class T0>
          custom_operation ( T0 const & , std::string const & operation_name) : program_name_(operation_name)
          {
              infos_base & tree = T0::get();
              std::cout << tree.name() << std::endl;
          }

        private:
          std::string program_name_;
     };
  }
}
#endif // CUSTOM_OPERATION_HPP
