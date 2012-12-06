#ifndef VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/code_generation/backend.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"

#include <algorithm>

#include <typeinfo>

namespace viennacl{

    namespace generator{


        namespace code_generation{

            class kernel_generator{
            private:

                std::string generate_headers() const{
                    std::string res;
                    res+="__kernel void " + kernel_name_ + "(";
                    std::list<kernel_argument*> args(utils::cast<kernel_argument>(utils::filter<utils::EXTRACT_IF>(trees_,utils::is_type<kernel_argument>)));
                    utils::remove_unsorted_duplicates(args);
                    args.sort(utils::deref_less());
                    for(std::list<kernel_argument*>::iterator it = args.begin() ; it!= args.end() ; ++it){
                        if(it!=args.begin()) res+=",";
                        res+=(*it)->kernel_arguments();
                    }
                    res+=")";
                    return res;
                }

                std::string generate_sources() const{
                    std::string res;
                    std::list<infos_base *> vec_exprs;
                    std::list<infos_base *> scal_exprs;
                    std::list<infos_base *> mat_exprs;
                    for(std::list<infos_base*>::const_iterator it = trees_.begin(); it!=trees_.end();++it){
                        if(utils::is_type<vector_expression_infos_base>(*it))
                            vec_exprs.push_back(*it);
                        else if(utils::is_type<scalar_expression_infos_base>(*it))
                            scal_exprs.push_back(*it);
                        else
                            mat_exprs.push_back(*it);
                    }
                    code_generation::blas1_generator gen(vec_exprs,scal_exprs);

                    return gen();
                }

            public:
                kernel_generator(std::list<infos_base*> const & trees, std::string const & kernel_name) : trees_(trees), kernel_name_(kernel_name){ }

                void generate(std::ostringstream & oss) const{
                    oss << generate_headers();
                    oss << generate_sources();
                }


            private:
                std::list<infos_base*> trees_;
                std::string kernel_name_;
            };

            class operations_manager{
            private:
                typedef std::list<viennacl::tools::shared_ptr<infos_base> > operations_t;
            public:
                template<class T>
                void add(T const & op){
                    operations_.push_back(viennacl::tools::shared_ptr<infos_base>(new T(op)));
                }

                void flush(){
                    operations_.clear();
                }

                std::vector<std::list<infos_base*> > get_kernels_list() const{
                    typedef std::vector<std::list<infos_base*> > res_t;
                    res_t res(1);
                    for(typename operations_t::const_iterator it = operations_.begin() ; it!=operations_.end() ; ++it){
                        infos_base* p = it->get();
                        std::list<infos_base*> inprods(utils::filter<utils::EXTRACT_IF>(p,utils::is_type<inprod_infos_base>));
                        if(inprods.size()){
                            res.back().merge(inprods);
                            res.push_back(std::list<infos_base*>());
                        }
                        res.back().push_back(p);

                    }
                    return res;
                }

            private:
                operations_t operations_;
            };

        }


    }

}
#endif // KERNEL_GENERATOR_FRONTEND_HPP
