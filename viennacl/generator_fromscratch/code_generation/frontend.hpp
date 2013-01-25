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
        void set_arguments(viennacl::ocl::kernel & k, ArgumentsT const & args){
            unsigned int counter=0;
            for(typename ArgumentsT::const_iterator iit = args.begin(); iit != args.end() ; ++iit){
                (*iit)->enqueue(counter,k);
            }
        }



        namespace code_generation{

//            struct kernel_infos_t{
//                enum type_t{
//                    UNDEFINED,
//                    BLAS1_TYPE,
//                    BLAS2_TYPE,
//                    BLAS3_TYPE
//                };

//                std::list<infos_base*> trees;
//                type_t type;
//            };

            class kernel_infos_t{
            public:
                enum type_t{
                    UNDEFINED,
                    BLAS1_TYPE,
                    BLAS2_TYPE,
                    BLAS3_TYPE
                };

            private:
                optimization_profile* make_profile(){
                    if(type_ == kernel_infos_t::BLAS1_TYPE)
                        return new blas1_optimization_profile();
                    if(type_ == kernel_infos_t::BLAS3_TYPE)
                        return new blas3_optimization_profile(16,256,256,4,4,4,true,false,4);
                }

            public:

                kernel_infos_t(std::list<infos_base*> trees, type_t type) : trees_(trees), type_(type), optimization_profile_(make_profile()){

                }

                kernel_infos_t(infos_base* op, type_t type) : type_(type), optimization_profile_(make_profile()){
                    trees_.push_back(op);

                }


                typedef std::set<kernel_argument*, deref_less> arguments_t;

                type_t type(){ return type_; }

                std::list<infos_base*> & trees(){ return trees_; }

                arguments_t & arguments(){ return arguments_; }

                code_generation::optimization_profile* profile() { return optimization_profile_.get(); }

                void config_nd_range(viennacl::ocl::kernel & k) const{
                    if(type_==BLAS3_TYPE){
                        matrix_expression_infos_base * p = static_cast<matrix_expression_infos_base*>(trees_.front());
                        mat_infos_base * mat = static_cast<mat_infos_base*>(&p->lhs());
                        blas3_optimization_profile* prof = static_cast<blas3_optimization_profile*>(optimization_profile_.get());
                        prof->set_global_sizes(mat->real_size1(),mat->real_size2());
                    }
                    optimization_profile_->config_nd_range(k);
                }

            private:
                arguments_t arguments_;
                std::list<infos_base*> trees_;
                type_t type_;
                viennacl::tools::shared_ptr<code_generation::optimization_profile> optimization_profile_;
            };



            class kernel_generator{
            private:

                void generate_headers(){
                    kss_ << "__kernel void " + kernel_name_ + "(";
                    for(std::list<infos_base*>::iterator it = kernel_infos_.trees().begin() ; it!= kernel_infos_.trees().end() ; ++it){
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
                    for(std::list<infos_base*>::const_iterator it = kernel_infos_.trees().begin(); it!=kernel_infos_.trees().end();++it){
                        if(utils::is_type<vector_expression_infos_base>()(*it))
                            vec_exprs.push_back(*it);
                        else if(utils::is_type<scalar_expression_infos_base>()(*it)
                                ||utils::is_type<inprod_infos_base>()(*it))
                            scal_exprs.push_back(*it);
                        else
                            mat_exprs.push_back(*it);
                    }
                    kernel_infos_.profile()->load(viennacl::ocl::current_device());
                    if(kernel_infos_.type()==kernel_infos_t::BLAS1_TYPE){
                        code_generation::blas1_generator gen(vec_exprs,mat_exprs,scal_exprs,kernel_infos_.profile());
                        gen(kss_);
                    }
                    else if(kernel_infos_.type()==kernel_infos_t::BLAS3_TYPE){
                        code_generation::blas3_generator gen(mat_exprs,static_cast<blas3_optimization_profile*>(kernel_infos_.profile()));
                        gen(kss_);
                    }
                    kss_.dec_tab();
                    kss_<<"}"<< std::endl;
                }

            public:
                kernel_generator(kernel_infos_t & kernel_infos
                                 , std::string const & kernel_name
                                 , code_generation::utils::kernel_generation_stream & kss) : kernel_infos_(kernel_infos)
                                                                     , kernel_name_(kernel_name)
                                                                     , kss_(kss){
                    kernel_infos_.profile()->apply(kernel_infos_.trees());
                }

                void generate(){
                    generate_headers();
                    generate_sources();
                }


            private:
                kernel_infos_t & kernel_infos_;
                std::string kernel_name_;
                utils::kernel_generation_stream & kss_;
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

                std::list<kernel_infos_t> get_kernels_list() const{
                    std::list<kernel_infos_t> res;
                    for(typename operations_t::const_iterator it = operations_.begin() ; it!=operations_.end() ; ++it){
                        infos_base* p = it->get();
                        std::list<infos_base*> inprods(utils::filter<utils::EXTRACT_IF>(p,utils::is_type<inprod_infos_base>()));
                        std::list<infos_base*> matmatprods(utils::filter<utils::EXTRACT_IF>(p,utils::is_type<matmat_prod_infos_base>()));
                        kernel_infos_t* last = NULL;
                        if(inprods.size()){
                            assert(matmatprods.size()==0 && "INVALID KERNEL !");
                            if(res.size()){
                                if(res.back().type() != kernel_infos_t::BLAS1_TYPE)
                                    res.push_back(kernel_infos_t(p,kernel_infos_t::BLAS1_TYPE));
                                else
                                    res.back().trees().merge(inprods);
                            }
                            else{
                                res.push_back(kernel_infos_t(p,kernel_infos_t::BLAS1_TYPE));
                            }
                        }
                        else if(matmatprods.size()){
                            if(res.size()){
                                if(res.back().type() != kernel_infos_t::BLAS3_TYPE)
                                    res.push_back(kernel_infos_t(p,kernel_infos_t::BLAS3_TYPE));
                                else
                                    res.back().trees().push_back(p);
                            }
                            else{
                                res.push_back(kernel_infos_t(p,kernel_infos_t::BLAS3_TYPE));
                            }
                        }
                        else{
                            if(res.size()){
                                if(res.back().type() != kernel_infos_t::BLAS1_TYPE)
                                    res.push_back(kernel_infos_t(p,kernel_infos_t::BLAS3_TYPE));
                                else
                                    res.back().trees().push_back(p);
                            }
                        }

                    }
                    return res;
                }

                std::string repr() const{
                    std::string res;
                    std::list<kernel_infos_t> kernels(get_kernels_list());
                    for(std::list<kernel_infos_t>::iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
                        std::string name;
                        for(std::list<infos_base*>::iterator iit = it->trees().begin() ; iit != it->trees().end() ; ++iit){
                            res += (*iit)->repr();
                        }
                    }
                    return res;
                }

                std::string get_source_code( std::map<std::string, generator::code_generation::kernel_infos_t> & kernels_infos) const{
                    std::ostringstream oss;
                    code_generation::utils::kernel_generation_stream kss(oss);
                    kss << "#if defined(cl_khr_fp64)\n";
                    kss <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
                    kss <<  "#elif defined(cl_amd_fp64)\n";
                    kss <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
                    kss <<  "#endif\n";

                    std::list<kernel_infos_t> kernels(get_kernels_list());
                    for(std::list<kernel_infos_t>::iterator it = kernels.begin() ; it !=kernels.end() ; ++it){
                        std::string name("_k"+to_string(std::distance(kernels.begin(),it)));
                        kernel_infos_t & infos = kernels_infos.insert(std::make_pair(name,*it)).first->second;
                        kss <<  "__attribute__((reqd_work_group_size(" << infos.profile()->local_work_size(0) << "," << infos.profile()->local_work_size(1) << ",1)))" << std::endl;
                        code_generation::kernel_generator kg(infos,name,kss);
                        kg.generate() ;
                    }
                    return oss.str();
                }

            private:
                operations_t operations_;
            };

        }


    }

}
#endif // KERNEL_GENERATOR_FRONTEND_HPP
