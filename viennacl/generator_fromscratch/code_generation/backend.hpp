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

                void compute_reductions_samesize(utils::kernel_generation_stream& kss, std::list<local_memory> const & lmems){
                   //Same size assumption
                   assert(std::adjacent_find(lmems.begin(), lmems.end(), std::not2(lid_has_same_size()))==lmems.end() && " Calling the wrong function for reducing inner products of different sizes! ");

                   kss << "for(unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2){" << std::endl;
                   kss.inc_tab();
                   kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                   kss << "if(get_local_id(0) < stride){" << std::endl;
                   kss.inc_tab();
                   for(std::list<local_memory>::const_iterator it = lmems.begin(); it != lmems.end() ; ++it){
                   kss <<  it->access("[get_local_id(0)]") <<  " += " <<  it->access("[get_local_id(0)+stride]") << std::endl;
                   }
                   kss.dec_tab();
                   kss << "}" << std::endl;
                   kss.dec_tab();
                   kss << "}" << std::endl;
                }

                class blas1_generator{
                private:
                    struct is_compute :  public std::unary_function<inprod_infos_base*, bool>{ bool  operator()(inprod_infos_base* p) const { return p->step() == inprod_infos_base::compute; } };
                    struct is_reduce : public std::unary_function<inprod_infos_base*, bool>{ bool  operator()(inprod_infos_base* p) const{ return p->step() == inprod_infos_base::reduce; } };

                public:

                    blas1_generator(std::list<infos_base * > const & vector_expressions, std::list<infos_base * > const & scalar_expressions):
                                                                                                        vector_expressions_(vector_expressions),
                                                                                                        scalar_expressions_(scalar_expressions)
                    {
                        std::list<vec_infos_base * > tmp0(code_generation::utils::extract_cast<vec_infos_base>(scalar_expressions_));
                        vectors_ = code_generation::utils::extract_cast<vec_infos_base>(vector_expressions_);
                        vectors_.merge(tmp0);

                        std::list<gpu_scal_infos_base * > tmp1(code_generation::utils::extract_cast<gpu_scal_infos_base>(scalar_expressions_));
                        gpu_scalars_ = code_generation::utils::extract_cast<gpu_scal_infos_base>(vector_expressions_);
                        gpu_scalars_.merge(tmp1);


                        std::list<inprod_infos_base * > tmp2(code_generation::utils::extract_cast<inprod_infos_base>(scalar_expressions_));
                        std::list<inprod_infos_base * >inner_prods(code_generation::utils::extract_cast<inprod_infos_base>(vector_expressions_));
                        inner_prods.merge(tmp2);

                        code_generation::utils::remove_unsorted_duplicates(vectors_);
                        code_generation::utils::remove_unsorted_duplicates(gpu_scalars_);
                        code_generation::utils::remove_unsorted_duplicates(inner_prods);
                        std::remove_copy_if(inner_prods.begin(),inner_prods.end(),std::back_inserter(inner_prods_reduce_),std::not1(is_reduce()));
                        std::remove_copy_if(inner_prods.begin(),inner_prods.end(),std::back_inserter(inner_prods_compute_),std::not1(is_compute()));
                    }


                    void operator()(utils::kernel_generation_stream& kss){

                        std::set<vec_infos_base *> vector_cached_entries;
                        std::set<gpu_scal_infos_base *> scalar_cached_entries;
                        code_generation::utils::cache_manager<vec_infos_base> vector_cache(vectors_,kss,vector_cached_entries);
                        code_generation::utils::cache_manager<gpu_scal_infos_base> scalar_cache(gpu_scalars_,kss,scalar_cached_entries);
                        vec_infos_base * first_vector =  NULL;
                        if(vectors_.size())
                            first_vector = vectors_.front();
//                        //Assumes same size...
                        if(inner_prods_reduce_.size()){
                            std::list<local_memory> local_mems;
                            for( std::list<inprod_infos_base *>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss <<  (*it)->scalartype() << " sum_" << (*it)->name()<< " = 0 ;" << std::endl;
                            }
                             kss << "for(unsigned int i = get_global_id(0) ; i <" << first_vector->size() << " ; i += get_global_size(0){" << std::endl;
                             kss.inc_tab();
                             for( std::list<inprod_infos_base *>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                 kss << "sum_" << (*it)->name() << " = " << (*it)->name() << "[i]" << std::endl;
                             }
                             kss.dec_tab();
                             kss << "}" << std::endl;
                            for( std::list<inprod_infos_base *>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(64));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " << "sum_" << (*it)->name()<< ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                            for( std::list<inprod_infos_base *>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                (*it)->access_name((*it)->make_local_memory(64).access("0"));
                            }
                            kss << "barrier(CLK_LOCAL_MEM_FENCE)" << std::endl;
                        }

                        scalar_cache.fetch_entries("0");
                        if(first_vector){
                            for(std::list<inprod_infos_base *>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                                kss << (*it)->scalartype() << " " << (*it)->name()+"_sum" << " = 0;" << std::endl;
                            }
                            kss << "for(unsigned int i = get_global_id(0) ; i <" << first_vector->size() << " ; i += get_global_size(0){" << std::endl;
                            kss.inc_tab();
                            vector_cache.fetch_entries("i");
                            for(std::list<infos_base *>::iterator it=vector_expressions_.begin() ; it!=vector_expressions_.end();++it){
                                kss << (*it)->generate() << std::endl;
                            }
                            for(std::list<inprod_infos_base *>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                                kss << (*it)->name()+"_sum" << " += " << "(" << (*it)->lhs().generate() << ")" << " * " << "(" << (*it)->rhs().generate() << ")" << std::endl;
                            }
                            vector_cache.writeback_entries("i");
                            kss.dec_tab();
                            kss << "}" << std::endl;
                        }
                        scalar_cache.writeback_entries("0");
                        if(inner_prods_compute_.size()){
                            std::list<local_memory> local_mems;
                            for( std::list<inprod_infos_base *>::const_iterator it = inner_prods_compute_.begin(); it != inner_prods_compute_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(64));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " <<  (*it)->name() + "_sum" << ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                        }
                        for(std::list<inprod_infos_base *>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                            (*it)->step(inprod_infos_base::reduce);
                        }
                    }

                private:
                    std::list<infos_base * >  vector_expressions_;
                    std::list<infos_base* > scalar_expressions_;
                    std::list<vec_infos_base * >  vectors_;
                    std::list<gpu_scal_infos_base * > gpu_scalars_;
                    std::list<inprod_infos_base * > inner_prods_compute_;
                    std::list<inprod_infos_base * > inner_prods_reduce_;


                };

        }

    }

}
#endif // BACKEND_HPP
