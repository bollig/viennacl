#ifndef VIENNACL_GENERATOR_KERNEL_GENERATOR_FRONTEND_HPP
#define VIENNACL_GENERATOR_KERNEL_GENERATOR_FRONTEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/tree_utils.hpp"
#include <algorithm>
namespace viennacl{

    namespace generator{


        namespace code_generation{

            class frontend{
            private:
                static bool sort_pred(kernel_argument* const & lhs, kernel_argument* const & rhs){
                    return lhs->id() < rhs->id();
                }

            public:
                frontend(std::list<infos_base*> const & trees, std::string const & op_name) : trees_(trees), op_name_(op_name){ }

                std::string generate_headers() const{
                    std::string res;
                    res+="__kernel void " + op_name_ + "(";
                    std::list<kernel_argument*> args(tree_utils::extract_type<kernel_argument>(trees_));
                    args.sort(sort_pred);
                    for(std::list<kernel_argument*>::iterator it = args.begin() ; it!= args.end() ; ++it){
                        if(it!=args.begin()) res+=",";
                        res+=(*it)->kernel_arguments();
                    }
                    res+=")";
                    return res;
                }

                std::string generate_sources() const{
                    std::string res;
                    return res;
                }

            private:
                std::list<infos_base*> trees_;
                std::string op_name_;
            };

        }


    }

}
#endif // KERNEL_GENERATOR_FRONTEND_HPP
