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

                    optimization_profile() : alignment_(1), loop_unroll_(1){
                        local_work_size_[0] = 128;
                        global_work_size_[0] = 128*128;
                    }

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

                    unsigned int alignment() const{ return alignment_; }
                    unsigned int loop_unroll(){ return loop_unroll_; }

                    friend std::ostream& operator<<(std::ostream& os, optimization_profile const & );

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
                protected:
                    size_type local_work_size_[2];
                    size_type global_work_size_[2];
                    unsigned int alignment_;
                    unsigned int loop_unroll_;
                };

                class blas1_optimization_profile : public optimization_profile{

                };

                class blas3_optimization_profile : public optimization_profile{
                public:

                    blas3_optimization_profile(unsigned int ml, unsigned int kl, unsigned int nl, unsigned int ms, unsigned int ks, unsigned int ns, bool use_LHS_shared, bool use_RHS_shared) : ml_(ml), kl_(kl), nl_(nl), ms_(ms), ks_(ks), ns_(ns), use_LHS_shared_(use_LHS_shared), use_RHS_shared_(use_RHS_shared) {
                        alignment_=4;

                        local_work_size_[0] =ml_/ms_;
                        local_work_size_[1] = nl_/ns_;

                        global_work_size_[0] = 2048/ml_*local_work_size_[0];
                        global_work_size_[1] = 2048/nl_*local_work_size_[1];

                    }

                    unsigned int ml() const{ return ml_ ; }
                    unsigned int kl() const{ return kl_ ; }
                    unsigned int nl() const{ return nl_ ; }
                    unsigned int ms() const{ return ms_ ; }
                    unsigned int ks() const{ return ks_ ; }
                    unsigned int ns() const{ return ns_ ; }
                    bool use_LHS_shared() const{ return use_LHS_shared_; }
                    bool use_RHS_shared() const{ return use_RHS_shared_; }

                private:
                    unsigned int ml_;
                    unsigned int kl_;
                    unsigned int nl_;

                    unsigned int ms_;
                    unsigned int ks_;
                    unsigned int ns_;

                    bool use_LHS_shared_;
                    bool use_RHS_shared_;
                };

                std::ostream& operator<<(std::ostream& os, optimization_profile const & prof){
                    os << " Alignment : " << prof.alignment_ << " | "
                       << " Unroll : " << prof.loop_unroll_ << " | "
                       << " Local Work Size 0 : " << prof.local_work_size_[0] << " | "
                       << " Global Work Size 0 : " << prof.global_work_size_[0] << " | "
                       << " Local Work Size 1 : " << prof.local_work_size_[1] << " | "
                       << " Global Work Size 1 : " << prof.global_work_size_[1] ;
                    return os;
                }

                void compute_reductions_samesize(utils::kernel_generation_stream& kss, std::list<local_memory> const & lmems){
                   //Same size assumption
                   assert(std::adjacent_find(lmems.begin(), lmems.end(), std::not2(lid_has_same_size()))==lmems.end() && " Calling the wrong function for reducing inner products of different sizes! ");
                   unsigned int size = lmems.front().size();
                   for(unsigned int stride = size/2 ; stride>0 ; stride /=2){
                       kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                       kss << "if(get_local_id(0) < " << to_string(stride) << "){" << std::endl;
                       kss.inc_tab();
                       for(std::list<local_memory>::const_iterator it = lmems.begin(); it != lmems.end() ; ++it){
                           kss <<  it->access("get_local_id(0)") <<  " += " <<  it->access("get_local_id(0)+" + to_string(stride)) << ";" << std::endl;
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
                                    , std::list<infos_base * > const & matrix_expressions
                                    , std::list<infos_base * > const & scalar_expressions
                                    , optimization_profile * kernel_config): vector_expressions_(vector_expressions), matrix_expressions_(matrix_expressions), scalar_expressions_(scalar_expressions), optimization_profile_(kernel_config)
                    {
                        std::set<inprod_infos_base*> inprods;
                        for(std::list<infos_base*>::const_iterator it=vector_expressions.begin() ; it!= vector_expressions.end() ; ++it){
                            extract_as(*it,vectors_,utils::is_type<vec_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                        for(std::list<infos_base*>::const_iterator it=matrix_expressions.begin() ; it!= matrix_expressions.end() ; ++it){
                            extract_as(*it,matrices_,utils::is_type<mat_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                        for(std::list<infos_base*>::const_iterator it=scalar_expressions.begin() ; it!= scalar_expressions.end() ; ++it){
                            extract_as(*it,vectors_,utils::is_type<vec_infos_base>());
                            extract_as(*it,matrices_,utils::is_type<mat_infos_base>());
                            extract_as(*it,gpu_scalars_,utils::is_type<gpu_scal_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                    }


                    void operator()(utils::kernel_generation_stream& kss){

                        unsigned int local_work_size0 = optimization_profile_->local_work_size(0);
                        unsigned int alignment = optimization_profile_->alignment();
                        unsigned int n_unroll = optimization_profile_->loop_unroll();
                        std::list<vec_infos_base *> assigned_vec;
                        for(std::list<infos_base*>::iterator it=vector_expressions_.begin(); it!= vector_expressions_.end();++it){
                            vector_expression_infos_base* p=static_cast<vector_expression_infos_base*>(*it);
                            if(p->op().is_assignment()==true) assigned_vec.push_back(dynamic_cast<vec_infos_base*>(&p->lhs()));
                        }

                        std::list<mat_infos_base *> assigned_mat;
                        for(std::list<infos_base*>::iterator it=matrix_expressions_.begin(); it!= matrix_expressions_.end();++it){
                            matrix_expression_infos_base* p=static_cast<matrix_expression_infos_base*>(*it);
                            if(p->op().is_assignment()==true) assigned_mat.push_back(dynamic_cast<mat_infos_base*>(&p->lhs()));
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
                        code_generation::utils::cache_manager<mat_infos_base> matrix_cache(matrices_,assigned_mat,kss);
                        code_generation::utils::cache_manager<gpu_scal_infos_base> scalar_cache(gpu_scalars_,assigned_scal,kss);
                        vec_infos_base * first_vector =  NULL;
                        mat_infos_base * first_matrix = NULL;
                        if(vectors_.size())
                            first_vector = *vectors_.begin();
                        if(matrices_.size())
                            first_matrix = *matrices_.begin();
//                        //Assumes same size...

                        if(inner_prods_reduce_.size()){
                            std::list<local_memory> local_mems;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss <<  (*it)->scalartype() <<  " " << (*it)->sum_name()<< " = 0 ;" << std::endl;
                            }
                            kss << "for(unsigned int i = get_global_id(0) ; i <" << optimization_profile_->global_work_size(0)/optimization_profile_->local_work_size(0) << " ; i += get_global_size(0)){" << std::endl;
                            kss.inc_tab();
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss << (*it)->sum_name() << " += " << (*it)->name() << "[i];" << std::endl;
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
                                (*it)->access_name(0,(*it)->make_local_memory(local_work_size0).access("0"));
                            }
                            kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                        }


                        scalar_cache.fetch_entries(0,"0");
                        if(first_vector){
                            for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                                kss << (*it)->scalartype() << " " << (*it)->sum_name() << " = 0;" << std::endl;
                            }
                            std::list<infos_base*> expressions;
                            std::transform(inner_prods_compute_.begin(), inner_prods_compute_.end(),std::back_inserter(expressions),UnsafeBase2Target<inprod_infos_base,infos_base>());
                            std::copy(vector_expressions_.begin(),vector_expressions_.end(),std::back_inserter(expressions));
                            utils::unroll_gid_loop_contiguous(kss,n_unroll,first_vector->size() + "/" +to_string(alignment) ,expressions,vector_cache);
                        }

                        if(first_matrix){
                            utils::unroll_gid_loop_contiguous(kss,n_unroll,first_matrix->size1()+"*"+first_matrix->size2()+ "/" +to_string(alignment) , matrix_expressions_, matrix_cache);
                        }
                        scalar_cache.writeback_entries(0,"0");

                        if(inner_prods_compute_.size()){
                            std::list<local_memory> local_mems;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_compute_.begin(); it != inner_prods_compute_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(local_work_size0));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " <<  (*it)->sum_name() << ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                        }
                        for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                            (*it)->step(inprod_infos_base::reduce);
                            kss << "if(get_local_id(0)==0) " << (*it)->name() << "[get_group_id(0)]" << "=" << (*it)->make_local_memory(local_work_size0).access("0") << ";" << std::endl;
                        }
                    }

                private:
                    std::list<infos_base * >  vector_expressions_;
                    std::list<infos_base * >  matrix_expressions_;
                    std::list<infos_base* > scalar_expressions_;
                    std::set<vec_infos_base *, viennacl::generator::deref_less >  vectors_;
                    std::set<mat_infos_base *, viennacl::generator::deref_less >  matrices_;
                    std::set<gpu_scal_infos_base *, viennacl::generator::deref_less > gpu_scalars_;
                    std::set<inprod_infos_base *, viennacl::generator::deref_less > inner_prods_compute_;
                    std::set<inprod_infos_base *, viennacl::generator::deref_less > inner_prods_reduce_;
                    optimization_profile * optimization_profile_;


                };

                class blas3_generator{
                    typedef std::set<matmat_prod_infos_base*,viennacl::generator::deref_less> matmat_prods_t;
                public:
                    blas3_generator(std::list<infos_base * > const & blas3_expressions
                                    , blas3_optimization_profile * kernel_config): blas3_expressions_(blas3_expressions), optimization_profile_(kernel_config)
                    {
                        for(std::list<infos_base*>::const_iterator it=blas3_expressions_.begin() ; it!= blas3_expressions_.end() ; ++it){
                            extract_as(*it, matmat_prods_, utils::is_type<matmat_prod_infos_base>());
                            extract_as(*it ,gpu_scalars_,  utils::is_type<gpu_scal_infos_base>());
                            extract_as(*it,matrices_, utils::is_type<mat_infos_base>());
                        }
                    }

                    void operator()(utils::kernel_generation_stream& kss){

                        unsigned int alignment = optimization_profile_->alignment();

                        unsigned int ml = optimization_profile_->ml();
                        unsigned int kl = optimization_profile_->kl()/alignment;
                        unsigned int nl = optimization_profile_->nl()/alignment;

                        unsigned int ms = optimization_profile_->ms();
                        unsigned int ks = optimization_profile_->ks()/alignment;
                        unsigned int ns = optimization_profile_->ns()/alignment;

                        bool use_LHS_shared = optimization_profile_->use_LHS_shared();
                        bool use_RHS_shared = optimization_profile_->use_RHS_shared();

                        std::list<mat_infos_base*> assigned;
                        std::set<mat_infos_base*,viennacl::generator::deref_less> lhss;
                        std::set<mat_infos_base*,viennacl::generator::deref_less> rhss;

                        //Fills assigned matrices set
                        for(std::list<infos_base*>::iterator it = blas3_expressions_.begin() ; it!=blas3_expressions_.end(); ++it){
                            matrix_expression_infos_base* p=static_cast<matrix_expression_infos_base*>(*it);
                            if(p->op().is_assignment()) assigned.push_back(static_cast<mat_infos_base*>(&p->lhs()));
                        }

                        //Fills lhs's
                        for(matmat_prods_t::iterator it = matmat_prods_.begin(); it !=matmat_prods_.end(); ++it){
                            extract_as(&(*it)->lhs(),lhss,utils::is_type<mat_infos_base>());
                            extract_as(&(*it)->rhs(),rhss,utils::is_type<mat_infos_base>());
                        }

                        mat_infos_base* first_lhs = *lhss.begin();
                        mat_infos_base* first_rhs = *rhss.begin();
                        mat_infos_base* first_assigned = assigned.front();

                        std::string local_lhs_name("local_" + first_lhs->name());
                        std::string local_rhs_name("local_" + first_rhs->name());

                        //Declare local memory
                        if(use_LHS_shared){
                                kss << "__local " << first_lhs->scalartype() << " " << local_lhs_name << "[" << kl*alignment * ml << "];" << std::endl;
                        }

                        if(use_RHS_shared){
                                kss << "__local " << first_rhs->aligned_scalartype() << " " << local_rhs_name << "[" << kl*alignment << "]" << "[" << nl << "];" << std::endl;
                        }

                        std::string res_table_name(first_assigned->name() + "_res");
//                        std::string rhs_table_name(first_rhs->name() + "_vals");

                        kss << first_assigned->aligned_scalartype() << " " << res_table_name << "[" << ms << "][" << ns << "];" << std::endl;

                        for(unsigned int m=0; m< ms; ++m){
                            for(unsigned int n=0; n < ns ; ++n){
                             kss << res_table_name << "[" << m << "][" << n << "] = 0 ;" << std::endl;
                            }
                        }
//                        kss << first_rhs->aligned_scalartype() << " " << rhs_table_name << "[" << ns << "];" << std::endl;


                        //Helper variables
                        kss << "unsigned int aligned_size2_lhs = " << first_lhs->internal_size2() << "/" << alignment<< ";" << std::endl;
                        kss << "unsigned int aligned_size2_rhs = " << first_rhs->internal_size2() << "/" << alignment<< ";" << std::endl;
                        kss << "unsigned int offset_m = get_local_id(0)*" << ms << ";" << std::endl;
                        kss << "unsigned int offset_n = get_local_id(1)*" << ns << ";" << std::endl;
                        kss << "unsigned int block_num = ( " << "aligned_size2_lhs" << " + " << kl -1 << ")/" << kl << ";" << std::endl;
                        kss << "for(unsigned int bl=0 ; bl<block_num ; ++bl){" << std::endl;
                        kss.inc_tab();
                        kss << "unsigned int offsetLHS = get_group_id(0)*" << "aligned_size2_lhs" << "*" << ml << "+ bl*" << kl << ";" << std::endl;
                        kss << "unsigned int offsetRHS = bl*" << "aligned_size2_rhs" << "*" << kl*alignment << "+ get_group_id(1)*" << nl << ";" << std::endl;
                        kss << "unsigned int n_subblock = min(" << kl/ks << ", (int)(" << "aligned_size2_lhs" << " - bl*" << kl << ")/ " << ks << ");" << std::endl;

                        if(use_LHS_shared){
                            kss << "for(unsigned int i = get_local_id(0)" << " ; i < " << ml << "; i+= get_local_size(0)){" << std::endl;
                            kss.inc_tab();
                            kss << "for(unsigned int j = get_local_id(1)" << " ; j < " << kl << "; j+= get_local_size(1)){" << std::endl;
                            kss.inc_tab();
                             kss << first_lhs->aligned_scalartype()  << " val_lhs = " << first_lhs->name() <<  "[offsetLHS + j  + aligned_size2_lhs*i];" << std::endl;
                           for(unsigned int a = 0 ; a < alignment ; ++a){
                               kss << local_lhs_name << "[(j*" << alignment << " +" <<a << ")*"  << ml << "+i] =  val_lhs.s" << a << ";" << std::endl;
                            }

                            kss.dec_tab();
                            kss << "}" << std::endl;
                            kss.dec_tab();
                            kss << "}" << std::endl;

                        }

                        if(use_RHS_shared){
                            kss << "for(unsigned int i = get_local_id(0)" << " ; i < " << kl*alignment << "; i+= get_local_size(0)){" << std::endl;
                            kss.inc_tab();
                            kss << "for(unsigned int j = get_local_id(1)" << " ; j < " << nl << "; j+= get_local_size(1)){" << std::endl;
                            kss.inc_tab();
                            kss << local_rhs_name << "[i][j] = " << first_rhs->name() << "[offsetRHS + j + aligned_size2_rhs*i];" << std::endl;
                            kss.dec_tab();
                            kss << "}" << std::endl;
                            kss.dec_tab();
                            kss << "}" << std::endl;

                        }
//                        kss << "#pragma unroll 2" << std::endl;
                        kss << "for(unsigned int bs=0 ; bs < " << kl/ks << " ; ++bs){" << std::endl;
                        kss.inc_tab();
                        if(use_LHS_shared){
                            kss << "__local " << first_lhs->scalartype() << "*lhs_ptr = " << local_lhs_name << "+ bs*" << ks*alignment << "*" << ml << "+ offset_m;" << std::endl;
                        }
                        kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                        if(!use_RHS_shared){
                            kss << "unsigned int smallOffsetRHS = offsetRHS + bs*" << ks*alignment << "*" << "aligned_size2_rhs" << "+ offset_n" << ";" << std::endl;
                        }
                        for(unsigned int k = 0 ; k < ks ; ++k){
                          for(unsigned int a=0; a<alignment; ++a){
                            for(unsigned int n=0 ; n < ns ; ++n){
                                if(use_RHS_shared){
                                   kss << first_rhs->aligned_scalartype() << " val_rhs_" << k*alignment + a << "_" << n << " = " <<  local_rhs_name << "[bs*" << ks*alignment << "+" << k*alignment + a << "][offset_n + " << n << "];" << std::endl;
                               }
                               else{
                                   kss << first_rhs->aligned_scalartype() << " val_rhs_" << k*alignment + a << "_" << n << " = " << first_rhs->name() << "[smallOffsetRHS  + aligned_size2_rhs*" << k*alignment+a<< " + " << n << "];" << std::endl;

                              }
                            }
                          }
                        }
                        for(unsigned int k = 0 ; k < ks ; ++k){
                            for(unsigned int n=0 ; n < ns ; ++n){
                                    for(unsigned int a=0; a<alignment; ++a){
                                      for(unsigned int m=0 ; m < ms ; ++m){
                                            kss << res_table_name<< "["<<m<<"][" << n << "] += " ;
                                            kss << "* lhs_ptr++";
                                            kss << "*";
                                            kss <<" val_rhs_" << k*alignment + a << "_" << n << ";" << std::endl;

                                      }
                                      kss << "lhs_ptr +=" << ml - ms << ";" << std::endl;
                                }
                            }

//       kss << res_table_name<< "["<<m<<"][" << n << "] += " << local_lhs_name << "[bs*" << ks*alignment << "+" << k*alignment + a << "][offset_m + " << m << "]"

    //                            kss << "smallOffsetRHS += " << "aligned_size2_rhs" << ";" << std::endl;
                        }
                        kss.dec_tab();
                        kss << "}" << std::endl;
                        kss.dec_tab();
                        kss << "}" << std::endl;

                        kss << "unsigned int offsetRes = get_group_id(0)*aligned_size2_rhs*" << ml << " + get_group_id(1)*" << nl << " + offset_m*aligned_size2_rhs + offset_n ;"  << std::endl;
                        for(unsigned int m=0 ; m < ms ; ++m){
                            for(unsigned int n=0 ; n < ns ; ++n){
                                kss << first_assigned->name() << "[offsetRes + aligned_size2_rhs*" << m << " + " << n << "]" << "=" << res_table_name << "[" << m << "][" << n << "];" << std::endl;
                            }
                        }
                    }

                private:
                    std::list<infos_base*>  blas3_expressions_;
                    matmat_prods_t matmat_prods_;
                    std::set<mat_infos_base *, viennacl::generator::deref_less >  matrices_;
                    std::set<gpu_scal_infos_base *, viennacl::generator::deref_less > gpu_scalars_;
                    blas3_optimization_profile * optimization_profile_;
                };

        }

    }

}
#endif // BACKEND_HPP
