#ifndef VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types_base.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"
#include <algorithm>

namespace viennacl{

    namespace generator{

        namespace code_generation{

                struct lid_has_same_size : public std::binary_function<local_memory,local_memory,bool>{
                    bool operator()(local_memory const & lhs, local_memory const & rhs) const{ return lhs.size() == rhs.size(); }
                };


                class optimization_profile{
                    typedef unsigned int size_type;
                public:
                    optimization_profile() : alignment_(4){ }

                    void load(viennacl::ocl::device const & d){

                    }

                    void apply(std::list<infos_base*> & expressions){
                        std::set<kernel_argument*,viennacl::generator::deref_less> kernel_arguments;
                        for(std::list<infos_base*>::iterator it = expressions.begin() ; it!=expressions.end() ; ++it){
                            extract_as(*it,kernel_arguments,utils::is_type<kernel_argument>());
                        }
                        for(std::set<kernel_argument*,viennacl::generator::deref_less>::iterator it=kernel_arguments.begin() ; it!= kernel_arguments.end() ; ++it){
                            (*it)->alignment(alignment_);
                        }
                    }

                    void alignment (unsigned int _alignment){
                        alignment_ = _alignment;
                    }

                    unsigned int alignment() const{
                        return alignment_;
                    }


                    void local_work_size(unsigned int index, size_type val){
                        assert(index==0 || index==1);
                        local_work_size_[index]=val;
                    }

                    size_type local_work_size(unsigned int index) const{
                        assert(index==0 || index==1);
                        return local_work_size_[index];
                    }

                    void global_work_size(unsigned int index, size_type val){
                        assert(index==0 || index==1);
                        global_work_size_[index]=val;
                    }

                    size_type global_work_size(unsigned int index) const{
                        assert(index==0 || index==1);
                        return global_work_size_[index];
                    }
                private:
                    size_type local_work_size_[2];
                    size_type global_work_size_[2];
                    unsigned int alignment_;
                };

                void compute_reductions_samesize(utils::kernel_generation_stream& kss, std::list<local_memory> const & lmems){
                   //Same size assumption
                   assert(std::adjacent_find(lmems.begin(), lmems.end(), std::not2(lid_has_same_size()))==lmems.end() && " Calling the wrong function for reducing inner products of different sizes! ");
                   unsigned int size = lmems.front().size();
                   for(unsigned int stride = size/2 ; stride>0 ; stride /=2){
                       kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                       kss << "if(get_local_id(0) < " << to_string(stride) << "){" << std::endl;
                       kss.inc_tab();
                       for(std::list<local_memory>::const_iterator it = lmems.begin(); it != lmems.end() ; ++it){
                           kss <<  it->access("get_local_id(0)") <<  " += " <<  it->access("get_local_id(0)+" + to_string(stride)) << std::endl;
                       }
                       kss.dec_tab();
                       kss << "}" << std::endl;
                   }
                }

                class blas1_generator{
                private:
                    struct is_compute :  public std::unary_function<inprod_infos_base*, bool>{ bool  operator()(inprod_infos_base* p) const { return p->step() == inprod_infos_base::compute; } };
                    struct is_reduce : public std::unary_function<inprod_infos_base*, bool>{ bool  operator()(inprod_infos_base* p) const{ return p->step() == inprod_infos_base::reduce; } };

                public:
                    blas1_generator(std::list<infos_base * > const & vector_expressions
                                    , std::list<infos_base * > const & scalar_expressions
                                    , optimization_profile const & kernel_config): vector_expressions_(vector_expressions), scalar_expressions_(scalar_expressions), optimization_profile_(kernel_config)
                    {
                        std::set<inprod_infos_base*> inprods;
                        for(std::list<infos_base*>::const_iterator it=vector_expressions.begin() ; it!= vector_expressions.end() ; ++it){
                            extract_as(*it,vectors_,utils::is_type<vec_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                        for(std::list<infos_base*>::const_iterator it=scalar_expressions.begin() ; it!= scalar_expressions.end() ; ++it){
                            extract_as(*it,vectors_,utils::is_type<vec_infos_base>());
                            extract_as(*it,gpu_scalars_,utils::is_type<gpu_scal_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                    }


                    void operator()(utils::kernel_generation_stream& kss){

                        unsigned int local_work_size0 = optimization_profile_.local_work_size(0);
                        unsigned int alignment = optimization_profile_.alignment();
                        std::list<vec_infos_base *> assigned_vec;

                        for(std::list<infos_base*>::iterator it=vector_expressions_.begin(); it!= vector_expressions_.end();++it){
                            vector_expression_infos_base* p=static_cast<vector_expression_infos_base*>(*it);
                            if(p->op().is_assignment()==true) assigned_vec.push_back(dynamic_cast<vec_infos_base*>(&p->lhs()));
                        }


                        std::list<gpu_scal_infos_base *> assigned_scal;
                        for(std::list<infos_base*>::iterator it=scalar_expressions_.begin(); it!= scalar_expressions_.end();++it){
                            if(scalar_expression_infos_base* p=dynamic_cast<scalar_expression_infos_base*>(*it)){
                                if(p->op().is_assignment()==true){
                                    assigned_scal.push_back(dynamic_cast<gpu_scal_infos_base*>(&p->lhs()));
                                }
                            }
                        }

                        code_generation::utils::cache_manager<vec_infos_base> vector_cache(vectors_,assigned_vec,kss);
                        code_generation::utils::cache_manager<gpu_scal_infos_base> scalar_cache(gpu_scalars_,assigned_scal,kss);
                        vec_infos_base * first_vector =  NULL;
                        if(vectors_.size())
                            first_vector = *vectors_.begin();
//                        //Assumes same size...

                        if(inner_prods_reduce_.size()){
                            std::list<local_memory> local_mems;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss <<  (*it)->scalartype() <<  " " << (*it)->sum_name()<< " = 0 ;" << std::endl;
                            }
                            kss << "for(unsigned int i = get_global_id(0) ; i <" << first_vector->size() << " ; i += get_global_size(0)){" << std::endl;
                            kss.inc_tab();
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss << (*it)->sum_name() << " += " << (*it)->name() << "[i]" << std::endl;
                            }
                            kss.dec_tab();
                            kss << "}" << std::endl;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(local_work_size0));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " << (*it)->sum_name()<< ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                (*it)->access_name((*it)->make_local_memory(local_work_size0).access("0"));
                            }
                            kss << "barrier(CLK_LOCAL_MEM_FENCE)" << std::endl;
                        }


                        scalar_cache.fetch_entries("0");
                        if(first_vector){
                            for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                                kss << (*it)->scalartype() << " " << (*it)->sum_name() << " = 0;" << std::endl;
                            }
                            kss << "for(unsigned int i = get_global_id(0)" << "; i <" << first_vector->size() << "/" << first_vector->alignment() << " ; i += get_global_size(0)){" << std::endl;
                            kss.inc_tab();
                            vector_cache.fetch_entries("i");
                            for(std::list<infos_base *>::iterator it=vector_expressions_.begin() ; it!=vector_expressions_.end();++it){
                                kss << (*it)->generate() << ";" << std::endl;
                            }
                            for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                                kss << (*it)->generate() << ";" << std::endl;
                            }
                            vector_cache.writeback_entries("i");
                            kss.dec_tab();
                            kss << "}" << std::endl;
                        }
                        scalar_cache.writeback_entries("0");

                        if(inner_prods_compute_.size()){
                            std::list<local_memory> local_mems;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_compute_.begin(); it != inner_prods_compute_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(64));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " <<  (*it)->sum_name() << ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                        }
                        for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                            (*it)->step(inprod_infos_base::reduce);
                        }
                    }

                private:
                    std::list<infos_base * >  vector_expressions_;
                    std::list<infos_base* > scalar_expressions_;
                    std::set<vec_infos_base *, viennacl::generator::deref_less >  vectors_;
                    std::set<gpu_scal_infos_base *, viennacl::generator::deref_less > gpu_scalars_;
                    std::set<inprod_infos_base *, viennacl::generator::deref_less > inner_prods_compute_;
                    std::set<inprod_infos_base *, viennacl::generator::deref_less > inner_prods_reduce_;
                    optimization_profile optimization_profile_;


                };

        }

    }

}
#endif // BACKEND_HPP
