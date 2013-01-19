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

                    blas3_optimization_profile(unsigned int ml, unsigned int kl, unsigned int nl, unsigned int ms, unsigned int ks, unsigned int ns) : ml_(ml), kl_(kl), nl_(nl), ms_(ms), ks_(ks), ns_(ns) {
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

                private:
                    unsigned int ml_;
                    unsigned int kl_;
                    unsigned int nl_;

                    unsigned int ms_;
                    unsigned int ks_;
                    unsigned int ns_;
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

                        unsigned int ml = optimization_profile_->ml();
                        unsigned int kl = optimization_profile_->kl();
                        unsigned int nl = optimization_profile_->nl();

                        unsigned int ms = optimization_profile_->ms();
                        unsigned int ks = optimization_profile_->ks();
                        unsigned int ns = optimization_profile_->ns();

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

                        //Declare all local memories
//                        if(use_LHS_shared){
//                            for(std::set<mat_infos_base*,viennacl::generator::deref_less>::iterator itm=lhss.begin() ; itm!=lhss.end(); ++itm){
//                                std::string lmem_name("local_"+(*itm)->name());
//                                kss << "__local " << (*itm)->scalartype()
//                                    << " " << lmem_name << "[" << ml << "]" << "[" << kl << "];"
//                                    << std::endl;
//                            }
//                        }
//                        for(std::set<mat_infos_base*,viennacl::generator::deref_less>::iterator itm=rhss.begin() ; itm!=rhss.end(); ++itm){
//                            std::string lmem_name("local_"+(*itm)->name());
//                            kss << "__local_block_ " << (*itm)->scalartype()
//                                << " " << lmem_name << "[" << lsize0 << "]" << "[" << lsize1 << "]"
//                                << std::endl;
//                        }

                        //Declare all assignment registers
//                        for(std::list<mat_infos_base*>::iterator it=assigned.begin(); it!=assigned.end(); ++it){
//                            std::string lmem_name( (*it)->name() + "_block");
//                            kss << "__local " << (*it)->scalartype() << " " << lmem_name << "[" << ml << "]" << "[" << nl << "];" << std::endl;
//                        }

                        std::string res_table_name(first_assigned->name() + "_res");
                        std::string lhs_table_name(first_lhs->name() + "_vals");
                        std::string rhs_table_name(first_rhs->name() + "_vals");

                        kss << first_assigned->scalartype() << " " << res_table_name << "[" << ms << "][" << ns << "] = {0};" << std::endl;
                        kss << first_lhs->scalartype() << " " << lhs_table_name << "[" << ms << "][" << ks << "] = {0};" << std::endl;
                        kss << first_rhs->scalartype() << " " << rhs_table_name << "[" << ks << "][" << ns << "] = {0};" << std::endl;

                        //                        kss << "for(unsigned int m = get_local_id(0)*" << ms << " ; m < (get_local_id(0)+1)*" << ms << "; ++m){" << std::endl;
//                        kss.inc_tab();
//                        kss << "for(unsigned int n = get_local_id(1)*" << ns << " ; n < (get_local_id(1)+1)*" << ns << "; ++n){" << std::endl;
//                        kss.inc_tab();
//                        for(std::list<mat_infos_base*>::iterator it=assigned.begin(); it!=assigned.end(); ++it){
//                            std::string lmem_name( (*it)->name() + "_block");
//                            kss << lmem_name << "[m][n] = 0;" << std::endl;
//                        }
//                        kss.dec_tab();
//                        kss << "}" << std::endl;
//                        kss.dec_tab();
//                        kss << "}" << std::endl;


                        //Helper variables
                        kss << "unsigned int offset_m = get_local_id(0)*" << ms << ";" << std::endl;
                        kss << "unsigned int offset_n = get_local_id(1)*" << ns << ";" << std::endl;
                        kss << "size_t block_num = ( " << first_lhs->size2() << " + " << kl -1 << ")/" << kl << ";" << std::endl;
                        kss << "for(unsigned int bl=0 ; bl<block_num ; ++bl){" << std::endl;
                        kss.inc_tab();
                        kss << "unsigned int offsetLHS = get_group_id(0)*" << first_lhs->internal_size2() << "*" << ml << "+ bl*" << kl << ";" << std::endl;
                        kss << "unsigned int offsetRHS = bl*" << first_rhs->internal_size2() << "*" << kl << "+ get_group_id(1)*" << nl << ";" << std::endl;
                        kss << "unsigned int n_subblock = min(" << kl/ks << ", (int)(" << first_lhs->internal_size2() << " - bl*" << kl << ")/ " << ns << ");" << std::endl;


                        kss << "for(unsigned int bs=0 ; bs < n_subblock ; ++bs){" << std::endl;
                        kss.inc_tab();

                        kss << "unsigned int smallOffsetLHS = offsetLHS + bs*" << ns << "+ offset_m*" << first_lhs->internal_size2() << ";" << std::endl;
                        kss << "unsigned int smallOffsetRHS = offsetRHS + bs*" << ks << "*" << first_rhs->internal_size2()<< "+ offset_n" << ";" << std::endl;


                        for(unsigned int m=0 ; m < ms ; ++m){
                            for(unsigned int k = 0 ; k < ks ; ++k){

                                kss <<  lhs_table_name << "[" << m << "][" << k << "]"
                                     << "="
                                     << first_lhs->name() <<  "[smallOffsetLHS + " << m << "*" << first_lhs->internal_size2() << " + " << k << "];"
                                     << std::endl;
                            }
                        }

                        for(unsigned int k=0 ; k < ks ; ++k){
                            for(unsigned int n = 0 ; n < ns ; ++n){
                                kss <<  rhs_table_name << "[" << k << "][" << n << "]"
                                     << "="
                                     << first_rhs->name() << "[smallOffsetRHS + " << n << " + " << k << "*" << first_rhs->internal_size2() << "];"
                                     << std::endl;
                            }
                        }

                        for(unsigned int m=0 ; m < ms ; ++m){
                            for(unsigned int n=0 ; n < ns ; ++n){
                                for(unsigned int k = 0 ; k < ks ; ++k){
                                         for(std::list<mat_infos_base*>::iterator it=assigned.begin(); it!=assigned.end(); ++it){
                                             kss << res_table_name<< "["<<m<<"][" << n << "] += " << lhs_table_name << "[" << m << "][" << k << "]"
                                                              << "*"
                                                              << rhs_table_name << "[" << k << "]" << "[" << n << "];" << std::endl;
                                         }
                                }
                            }
//                            kss << "smallOffsetLHS +=" << first_lhs->internal_size2() << ";" << std::endl;
                        }
                        kss.dec_tab();
                        kss << "}" << std::endl;
                        kss.dec_tab();
                        kss << "}" << std::endl;

                        kss << "unsigned int offsetRes = get_group_id(0)*" << first_assigned->internal_size2() << "*" << ml << " + get_group_id(1)*" << nl << ";" << std::endl;


                        kss << "for(unsigned int m = 0 ; m < " << ms << "; ++m){" << std::endl;
                        kss.inc_tab();
                        kss << "for(unsigned int n = 0 ; n < " << ns << "; ++n){" << std::endl;
                        kss.inc_tab();
                        for(std::list<mat_infos_base*>::iterator it=assigned.begin(); it!=assigned.end(); ++it){
                            kss << "if (offset_n+n < " << (*it)->internal_size2() << " && offset_m + m < " << (*it)->internal_size1() << " )" << std::endl;
                            kss << (*it)->name() << "[offsetRes + (offset_m+m)*" << first_assigned->internal_size2() << "+ offset_n+n] =" << res_table_name << "[m][n];" << std::endl;

                        }
                        kss.dec_tab();
                        kss << "}" << std::endl;
                        kss.dec_tab();
                        kss << "}" << std::endl;

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
