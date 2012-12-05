#ifndef VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/tree_utils.hpp"
#include "viennacl/generator_fromscratch/code_generation/backend.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"

#include <algorithm>

namespace viennacl{

    namespace generator{


        namespace code_generation{

            class frontend{
            private:
                static bool sort_pred(kernel_argument* const & lhs, kernel_argument* const & rhs){
                    return lhs->id() < rhs->id();
                }

                std::string generate_headers() const{
                    std::string res;
                    res+="__kernel void " + op_name_ + "(";
                    std::list<kernel_argument*> args(utils::cast<kernel_argument>(utils::extract_if(trees_,code_generation::utils::is_type<kernel_argument>)));
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
                    std::list<infos_base *> vec_exprs(code_generation::utils::extract_if (trees_,code_generation::utils::is_type<vector_expression_infos_base>,false));
                    std::list<infos_base *> scal_exprs(code_generation::utils::extract_if(trees_,code_generation::utils::is_type<scalar_expression_infos_base>,false));
                    std::list<infos_base *> mat_exprs(code_generation::utils::extract_if (trees_,code_generation::utils::is_type<matrix_expression_infos_base>,false));
                    code_generation::blas1_generator gen(vec_exprs,scal_exprs);

                    return gen();
                }

            public:
                frontend(std::list<infos_base*> const & trees, std::string const & op_name) : trees_(trees), op_name_(op_name){ }

                std::string generate() const{
                    return generate_headers() + generate_sources();
                }


            private:
                std::list<infos_base*> trees_;
                std::string op_name_;
            };

        }


    }

}
#endif // KERNEL_GENERATOR_FRONTEND_HPP
