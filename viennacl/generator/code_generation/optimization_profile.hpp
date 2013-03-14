#ifndef VIENNACL_GENERATOR_CODE_GENERATION_OPTIMIZATION_PROFILE
#define VIENNACL_GENERATOR_CODE_GENERATION_OPTIMIZATION_PROFILE

#include "viennacl/generator/symbolic_types_base.hpp"

namespace viennacl{

namespace generator{

namespace code_generation{

class optimization_profile{
protected:
    typedef unsigned int size_type;


public:

    optimization_profile() : alignment_(1){ }

    void load(viennacl::ocl::device const & d){

    }

    virtual std::string repr() const = 0;

    virtual void config_nd_range(viennacl::ocl::kernel & k, infos_base * p) = 0;

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

    virtual std::pair<size_t,size_t> local_work_size() const = 0;

//    friend std::ostream& operator<<(std::ostream& os, optimization_profile const & );

    virtual ~optimization_profile(){ }

//    void local_work_size(unsigned int index, size_type val){
//        assert(index==0 || index==1);
//        local_work_size_[index]=val;
//    }

//    size_type local_work_size(unsigned int index) const{
//        assert(index==0 || index==1);
//        return local_work_size_[index];
//    }

//    void global_work_size(unsigned int index, size_type val){
//        assert(index==0 || index==1);
//        global_work_size_[index]=val;
//    }

//    size_type global_work_size(unsigned int index) const{
//        assert(index==0 || index==1);
//        return global_work_size_[index];
//    }


protected:
//    size_type local_work_size_[2];
//    size_type global_work_size_[2];
    unsigned int alignment_;
};

class blas1_optimization_profile : public optimization_profile{
public:

    blas1_optimization_profile(){
        alignment_ = 1;
        loop_unroll_ = 1;
        group_size0_ = 128;
    }

    blas1_optimization_profile(unsigned int alignment, unsigned int loop_unroll, size_t group_size0){
        alignment_ = alignment;
        loop_unroll_ = loop_unroll;
        group_size0_ = group_size0;
    }

    unsigned int loop_unroll() const{
        return loop_unroll_;
    }

    std::pair<size_t,size_t> local_work_size() const{
        return std::make_pair(group_size0_,1);
    }

    void config_nd_range(viennacl::ocl::kernel & k, infos_base* p){
        k.local_work_size(0,group_size0_);
        if(vec_infos_base* vec = dynamic_cast<vec_infos_base*>(p)){
            k.global_work_size(0,vec->real_size()/(alignment_*loop_unroll_));
        }
        else if(mat_infos_base * mat = dynamic_cast<mat_infos_base*>(p)){
            k.global_work_size(0,mat->real_size1() * mat->real_size2()/(alignment_*loop_unroll_));
        }
    }


    virtual std::string repr() const{
        std::ostringstream oss;
        oss << " A" << alignment_
           << " U" << loop_unroll_
           << " GROUP" << group_size0_;
        return oss.str();
    }

private:
    unsigned int loop_unroll_;
    unsigned int group_size0_;
};

class blas3_optimization_profile : public optimization_profile{
private:

public:

    blas3_optimization_profile(){
         ml_ = 32;
         kl_ = 32;
         nl_ = 32;

         ms_ = 4;
         ks_ = 4;
         ns_ = 4;

         use_LHS_shared_ = true;
         use_RHS_shared_ = false;

         unroll_ = 1;
    }

    blas3_optimization_profile(unsigned int ml, unsigned int kl, unsigned int nl
                               , unsigned int ms, unsigned int ks, unsigned int ns
                               , bool use_LHS_shared, bool use_RHS_shared
                               , unsigned int alignment
                               , unsigned int unroll) : ml_(ml), kl_(kl), nl_(nl), ms_(ms), ks_(ks), ns_(ns), use_LHS_shared_(use_LHS_shared), use_RHS_shared_(use_RHS_shared){
        ml_= ml; kl_=kl ; nl_=nl;
        ms_ = ms; ks_=ks; ns_=ns;
        use_LHS_shared_ = use_LHS_shared ; use_RHS_shared = use_RHS_shared_;
        alignment_ = alignment;
        unroll_ = unroll;
    }

    virtual std::string repr() const{
        std::ostringstream oss;
        oss << "ML" << ml_
           << "KL" << kl_
           << "NL" << nl_
           << "MS" << ms_
           << "KS" << ks_
           << "NS" << ns_
           << "LHS" << use_LHS_shared_
           << "RHS" << use_RHS_shared_
           << "A" << alignment_
           << "U" << unroll_;

        return oss.str();
    }


    void config_nd_range(viennacl::ocl::kernel & k, infos_base* p){
        mat_infos_base* mat = static_cast<mat_infos_base*>(p);
        k.local_work_size(0, ml_/ms_);
        k.global_work_size(0, mat->real_size1()/ms_);
        k.local_work_size(1, nl_/ns_);
        k.global_work_size(1, mat->real_size2()/ns_);
    }

    std::pair<size_t,size_t> local_work_size() const{
        return std::make_pair(ml_/ms_, nl_/ns_);
    }

    unsigned int ml() const{ return ml_ ; }
    unsigned int kl() const{ return kl_ ; }
    unsigned int nl() const{ return nl_ ; }
    unsigned int ms() const{ return ms_ ; }
    unsigned int ks() const{ return ks_ ; }
    unsigned int ns() const{ return ns_ ; }
    bool use_LHS_shared() const{ return use_LHS_shared_; }
    bool use_RHS_shared() const{ return use_RHS_shared_; }
    unsigned int unroll() const { return unroll_; }
private:
    unsigned int ml_;
    unsigned int kl_;
    unsigned int nl_;

    unsigned int ms_;
    unsigned int ks_;
    unsigned int ns_;

    bool use_LHS_shared_;
    bool use_RHS_shared_;

    unsigned int unroll_;
};

static std::ostream& operator<<(std::ostream& os, optimization_profile const & prof){
    os << prof.repr();
    return os;
}

}

}

}

#endif
