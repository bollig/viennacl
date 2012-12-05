#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP


#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/code_generation/frontend.hpp"

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
          custom_operation ( T0 t0, std::string const & operation_name) : program_name_(operation_name)
          {
              std::list<infos_base*> ops;
              ops.push_back(static_cast<infos_base*>(&t0));
              code_generation::frontend k(ops,operation_name);
              std::cout << k.generate()<< std::endl;
          }

        private:
          std::string program_name_;
     };
  }
}
#endif // CUSTOM_OPERATION_HPP
