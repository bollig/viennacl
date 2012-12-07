#ifndef VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"
#include <algorithm>

namespace viennacl{

    namespace generator{

        namespace code_generation{

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

                    void compute_reductions(utils::kernel_generation_stream& kss, std::list<inprod_infos_base *> const & inprods){
                       for( std::list<inprod_infos_base *>::const_iterator it = inprods.begin(); it != inprods.end() ; ++it){
                           kss << "local_" << (*it)->name() << "[get_local_id(0)]  = " << "sum_" << (*it)->name() << ";" << std::endl;
                       }
                       kss << "for(unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2){" << std::endl;
                       kss.inc_tab();
                       kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                       kss << "if(get_local_id(0) < stride){" << std::endl;
                       kss.inc_tab();
                       for(std::list<inprod_infos_base *>::const_iterator it = inprods.begin(); it != inprods.end() ; ++it){
                       kss <<  (*it)->name() << "[get_local_id(0)]  += " << "local_" << (*it)->name() << "shared_memory_ptr[get_local_id(0)+stride];" << std::endl;
                       }
                       kss.dec_tab();
                       kss << "}" << std::endl;
                       kss.dec_tab();
                       kss << "}" << std::endl;
                       for(std::list<inprod_infos_base *>::const_iterator it = inprods.begin(); it != inprods.end() ; ++it){
                            (*it)->access_name((*it)->name() + "[0]");
                       }
                       kss << "barrier(CLK_LOCAL_MEM_FENCE)" << std::endl;
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
                        if(inner_prods_reduce_.size()>0)
                            compute_reductions(kss,inner_prods_reduce_);
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
                                kss << (*it)->name()+"_sum" << " = " << "(" << (*it)->lhs().generate() << ")" << " * " << "(" << (*it)->rhs().generate() << ")" << std::endl;
                            }
                            vector_cache.writeback_entries("i");
                            kss.dec_tab();
                            kss << "}" << std::endl;
                        }
                        scalar_cache.writeback_entries("0");
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
