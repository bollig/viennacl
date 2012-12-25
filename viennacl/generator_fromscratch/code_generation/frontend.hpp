#ifndef VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types_base.hpp"
#include "viennacl/generator_fromscratch/code_generation/backend.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"

#include <algorithm>

#include <typeinfo>

namespace viennacl{

    namespace generator{


        template<class ArgumentsT>
        void set_arguments(viennacl::ocl::kernel & k, ArgumentsT & args){
            unsigned int counter=0;
            for(typename ArgumentsT::iterator iit = args.begin(); iit != args.end() ; ++iit){
                (*iit)->enqueue(counter,k);
            }
        }

        namespace code_generation{

            class kernel_infos_t{
            public:
                typedef std::set<kernel_argument*, deref_less> arguments_t;
                arguments_t & arguments(){ return arguments_; }
                code_generation::optimization_profile & profile() { return optimization_profile_; }
            private:
                arguments_t arguments_;
                code_generation::optimization_profile optimization_profile_;
            };

            class kernel_generator{
            private:

                void generate_headers(){
                    kss_ << "__kernel void " + kernel_name_ + "(";
                    for(std::list<infos_base*>::iterator it = trees_.begin() ; it!= trees_.end() ; ++it){
                        extract_as(*it,kernel_infos_.arguments(),utils::is_type<kernel_argument>());
                    }
                    for(std::set<kernel_argument*, deref_less>::iterator it=kernel_infos_.arguments().begin(); it!=kernel_infos_.arguments().end();++it){
                        if(it!=kernel_infos_.arguments().begin()) kss_ << ',';
                        kss_ << (*it)->arguments_string() << std::endl ;
                    }
                    kss_ << ")" << std::endl;
                }

                void generate_sources(){
                    kss_<<"{"<< std::endl;
                    kss_.inc_tab();
                    std::list<infos_base *> vec_exprs;
                    std::list<infos_base *> scal_exprs;
                    std::list<infos_base *> mat_exprs;
                    for(std::list<infos_base*>::const_iterator it = trees_.begin(); it!=trees_.end();++it){
                        if(utils::is_type<vector_expression_infos_base>()(*it))
                            vec_exprs.push_back(*it);
                        else if(utils::is_type<scalar_expression_infos_base>()(*it)
                                ||utils::is_type<inprod_infos_base>()(*it))
                            scal_exprs.push_back(*it);
                        else
                            mat_exprs.push_back(*it);
                    }
                    kernel_infos_.profile().load(viennacl::ocl::current_device());
                    code_generation::blas1_generator gen(vec_exprs,scal_exprs,kernel_infos_.profile());
                    gen(kss_);
                    kss_.dec_tab();
                    kss_<<"}"<< std::endl;
                }

            public:
                kernel_generator(std::list<infos_base*> & trees
                                 , std::string const & kernel_name
                                 , code_generation::utils::kernel_generation_stream & kss
                                 , kernel_infos_t & kernel_infos) : trees_(trees)
                                                                     , kernel_name_(kernel_name)
                                                                     , kss_(kss) ,kernel_infos_(kernel_infos){
                    kernel_infos_.profile().apply(trees);
                }

                void generate(){
                    generate_headers();
                    generate_sources();
                }


            private:
                std::list<infos_base*> trees_;
                std::string kernel_name_;
                utils::kernel_generation_stream & kss_;
                kernel_infos_t & kernel_infos_;
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
                        std::list<infos_base*> inprods(utils::filter<utils::EXTRACT_IF>(p,utils::is_type<inprod_infos_base>()));
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
