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

          /** @brief CTor for 1 expression
          *
          * @param operation_name the code for this expression will be stored in the program provided by this name
          */
          template<class T0>
          custom_operation ( T0 const & t0, std::string const & operation_name) : program_name_(operation_name), raw_ptrs_(1)
          {
              ops_.push_back(viennacl::tools::shared_ptr<infos_base>(static_cast<infos_base*>(new T0(t0))));
              std::transform(ops_.begin(),ops_.end(),raw_ptrs_.begin(),SharedPtr2Raw<infos_base>());
              code_generation::frontend k(raw_ptrs_,operation_name);
              std::cout << k.generate()<< std::endl;
          }

          template<class T0, class T1>
          custom_operation ( T0 const & t0, T1 const & t1,std::string const & operation_name) : program_name_(operation_name), raw_ptrs_(2)
          {
              ops_.push_back(viennacl::tools::shared_ptr<infos_base>(static_cast<infos_base*>(new T0(t0))));
              ops_.push_back(viennacl::tools::shared_ptr<infos_base>(static_cast<infos_base*>(new T1(t1))));
              std::transform(ops_.begin(),ops_.end(),raw_ptrs_.begin(),SharedPtr2Raw<infos_base>());
              code_generation::frontend k(raw_ptrs_,operation_name);
              std::cout << k.generate()<< std::endl;
          }

        private:
          std::string program_name_;
          std::list<viennacl::tools::shared_ptr<infos_base> > ops_;
          std::list<infos_base*> raw_ptrs_;
     };
  }
}
#endif // CUSTOM_OPERATION_HPP
