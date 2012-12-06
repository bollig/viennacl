#ifndef VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types.hpp"
#include "viennacl/generator_fromscratch/tree_utils.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"

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
                        vectors_ = code_generation::utils::extract_cast<vec_infos_base>(vector_expressions_);

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

                    void compute_reductions(std::ostringstream & oss, std::list<inprod_infos_base *> const & inprods){
                       for( std::list<inprod_infos_base *>::const_iterator it = inprods.begin(); it != inprods.end() ; ++it){
                           oss << " local_" << (*it)->name() << "[get_local_id(0)]  = " << "sum_" << (*it)->name() << ";\n";
                       }
                       oss << "   for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n";
                       oss << "   {\n";
                       oss << "      barrier(CLK_LOCAL_MEM_FENCE);\n";
                       oss << "      if (get_local_id(0) < stride){\n";
                       for(std::list<inprod_infos_base *>::const_iterator it = inprods.begin(); it != inprods.end() ; ++it){
                       oss << "         " << (*it)->name() << "[get_local_id(0)]  += " << "local_" << (*it)->name() << "shared_memory_ptr[get_local_id(0)+stride];\n";
                       }
                       oss << "      }\n";
                       oss << "   }\n";
                       for(std::list<inprod_infos_base *>::const_iterator it = inprods.begin(); it != inprods.end() ; ++it){
                       oss << "    " << (*it)->generate() << " = local_" << (*it)->name() << "[" << 0 << "];\n";
                       }
                    }

                    std::string operator()(){
                        std::ostringstream oss;

                        std::set<vec_infos_base *> vector_cached_entries;
                        std::set<gpu_scal_infos_base *> scalar_cached_entries;
                        code_generation::utils::cache_manager<vec_infos_base> vector_cache(vectors_,oss,vector_cached_entries);
                        code_generation::utils::cache_manager<gpu_scal_infos_base> scalar_cache(gpu_scalars_,oss,scalar_cached_entries);

                        vec_infos_base * first_vector =  NULL;
                        if(vectors_.size())
                            first_vector = vectors_.front();
//                        //Assumes same size...
                        oss << "{\n";

                        if(inner_prods_reduce_.size()>0)
                            compute_reductions(oss,inner_prods_reduce_);
                        scalar_cache.fetch_entries("0");
                        if(first_vector){
                            oss << "for(unsigned int i = get_global_id(0) ; i <" << first_vector->size() << " ; i += get_global_size(0){\n";
                            vector_cache.fetch_entries("i");
                            for(std::list<infos_base *>::iterator it=vector_expressions_.begin() ; it!=vector_expressions_.end();++it){
                                oss << (*it)->generate() << std::endl;
                            }
                            vector_cache.writeback_entries("i");
                            oss << "}\n";
                        }
                        scalar_cache.writeback_entries("0");
                        oss << "}\n";

                        return oss.str();
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
