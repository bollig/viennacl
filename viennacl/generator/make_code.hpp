#ifndef VIENNACL_GENERATOR_MAKE_CODE_HPP
#define VIENNACL_GENERATOR_MAKE_CODE_HPP

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file viennacl/generator/make_code.hpp
 *  @brief Definition of code generation policies. Experimental.
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "viennacl/generator/forwards.h"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/symbolic_types.hpp"
#include "viennacl/generator/result_of.hpp"
#include "viennacl/generator/tree_operations.hpp"
#include <list>
#include <set>

namespace viennacl
{
  namespace generator
  {


    /** @brief Inline code for an expression from scalars
        @tparam T Tree specifying the expression
    */
    template <class T>
    struct make_expression_code
    {
      static std::string value(std::string const & loop_accessor)
      {
        return T::name();
      }
    };

    template <long VAL>
    struct make_expression_code<symbolic_constant<VAL> >
    {
      static std::string value(std::string const & )
      {
        return to_string(VAL);
      }
    };

    template <class LHS, class RHS>
    struct make_expression_code<compound_node<LHS,inner_prod_type,RHS > >
    {
      private:
        typedef compound_node<LHS,inner_prod_type,RHS> T;

      public:
        static std::string value(std::string const & )
        {
            return T::name()+"_val";
        }
    };

    template< >
    struct make_expression_code< NullType >
    {
      static std::string value(std::string const & )
      {
        return "";
      }
    };

    template<class T>
    struct make_expression_code< elementwise_modifier<T> >
    {
      static std::string value ( std::string const & loop_accessor)
      {
        return elementwise_modifier<T>::modify(make_expression_code<typename T::PRIOR_TYPE>::value(loop_accessor));
      }
    };

    template<class LHS, class OP, class RHS >
    struct make_expression_code<compound_node<LHS, OP, RHS> >
    {
      static const std::string value(std::string const & loop_accessor)
      {
        return " ( " + make_expression_code<LHS>::value(loop_accessor)
                + OP::expression_string()
                + make_expression_code<RHS>::value(loop_accessor) + " ) ";
      }
    };


    static void replace_all_string(std::string & str, std::string const & from, std::string const & to){
        size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                 str.replace(start_pos, from.length(), to);
                 start_pos += to.length();
        }
    }

    /**
     * @brief The infos_base class
     */
    class infos_base{
    public:
        std::string const & name() const{ return name_; }
        std::string const & scalartype() const{ return scalartype_; }
        void access_name(std::string const & new_name) { access_name_ = new_name; }
        std::string access_name() const { return access_name_; }
        virtual ~infos_base(){ }
        bool is_modified(){ return is_modified_;}
        void is_modified(bool val){ is_modified_ = val; }
    protected:
        infos_base(std::string const & scalartype,
                   std::string const & name): scalartype_(scalartype), name_(name), is_modified_(false){ }
    private:
        std::string scalartype_;
        std::string name_;
        std::string access_name_;
        bool is_modified_;
    };



    class mat_infos_base : public infos_base{
    public:
        std::string const & size1() const{ return size1_; }
        std::string const & size2() const{ return size2_; }
        bool const is_rowmajor() const { return is_rowmajor_; }
        bool const is_transposed() const { return is_transposed_; }
        virtual ~mat_infos_base() { }
    protected:
        mat_infos_base(std::string const & scalartype
                       ,std::string const & name
                       ,std::string const & size1
                       ,std::string const & size2
                       ,bool is_rowmajor
                       ,bool is_transposed) : infos_base(scalartype,name)
                                              ,  size1_(size1)
                                              , size2_(size2)
                                              , is_rowmajor_(is_rowmajor)
                                              , is_transposed_(is_transposed){ }
    private:
        std::string size1_;
        std::string size2_;
        bool is_rowmajor_;
        bool is_transposed_;
    };


    /**
     * @brief The mat_infos class
     */
    template<class T>
    class mat_infos : public mat_infos_base{
    public:
        mat_infos() : mat_infos_base(print_type<typename T::ScalarType,1>::value()
                                     ,T::name()
                                     ,T::size1_name()
                                     ,T::size2_name()
                                     ,result_of::is_row_major<T>::value
                                     ,false){ }
        static infos_base & get(){
            static mat_infos<T> res;
            return res;
        }
    };


    class vec_infos_base : public infos_base{
    public:
        std::string const & size() const{ return size_; }
        virtual ~vec_infos_base(){ }
    protected:
        vec_infos_base(std::string const & scalartype, std::string const & name, std::string const & size) :
                                                          infos_base(scalartype,name)
                                                         ,size_(size){ }
    private:
        std::string size_;
    };



    /**
     * @brief The vec_infos class
     */

    template<class T>
    class vec_infos : public vec_infos_base{
    public:
        vec_infos() : vec_infos_base(print_type<typename T::ScalarType,1>::value(),T::name(),T::size2_name()) { }
        static infos_base & get(){
            static vec_infos<T> res;
            return res;
        }
        virtual ~vec_infos(){ }
    };


    class scal_infos_base : public infos_base{
    protected:
        scal_infos_base(std::string const & scalartype, std::string const & name) : infos_base(scalartype,name){ }
    };

    class cpu_scal_infos_base : public scal_infos_base{
    protected:
        cpu_scal_infos_base(std::string const & scalartype, std::string const & name) : scal_infos_base(scalartype,name){ }
    };

    class gpu_scal_infos_base : public scal_infos_base{
    protected:
        gpu_scal_infos_base(std::string const & scalartype, std::string const & name) : scal_infos_base(scalartype,name){ }
    };

    class constant_scal_infos_base : public scal_infos_base{
    protected:
        constant_scal_infos_base(std::string const & scalartype, std::string const & name) : scal_infos_base(scalartype,name){ }
    };

    /**
     * @brief The scal_infos class
     */
    template<class T>
    class cpu_scal_infos : public cpu_scal_infos_base{
    public:
        cpu_scal_infos() : cpu_scal_infos_base(print_type<typename T::ScalarType,1>::value() ,T::name()) { }

        static infos_base & get(){
            static cpu_scal_infos<T> res;
            return res;
        }
    };

    template<class T>
    class gpu_scal_infos : public gpu_scal_infos_base{
    public:
        gpu_scal_infos() : gpu_scal_infos_base(print_type<typename T::ScalarType,1>::value() ,T::name()) { }

        static infos_base & get(){
            static gpu_scal_infos<T> res;
            return res;
        }
    };

    template<class T>
    class constant_scal_infos : public constant_scal_infos_base{
    public:
        constant_scal_infos() : constant_scal_infos_base(print_type<typename T::ScalarType,1>::value() ,T::name()) { }

        static infos_base & get(){
            static constant_scal_infos<T> res;
            return res;
        }
    };

    template <unsigned int ID, typename SCALARTYPE>
    static infos_base &  get_infos(cpu_symbolic_scalar<ID,SCALARTYPE> const &){
        typedef cpu_symbolic_scalar<ID,SCALARTYPE>  U;
        return cpu_scal_infos<U>::get();
    }

    template<class T, class Enable=void>
    struct inprod_infos;

    template <unsigned int ID, typename SCALARTYPE>
    static infos_base &  get_infos(gpu_symbolic_scalar<ID,SCALARTYPE> const &){
        typedef gpu_symbolic_scalar<ID,SCALARTYPE>  U;
        return gpu_scal_infos<U>::get();
    }

    template<long T>
    static infos_base & get_infos(symbolic_constant<T> const &){
        typedef symbolic_constant<T> U;
        return constant_scal_infos<U>::get();
    }

    template <class T>
    static infos_base &  get_infos(inner_prod_impl_t<T> const &){
        typedef inner_prod_impl_t<T> U;
        return inprod_infos<U>::get();
    }

    template <class LHS, class RHS>
    static infos_base &  get_infos(compound_node<LHS, inner_prod_type, RHS> const &){
        typedef compound_node<LHS, inner_prod_type, RHS> U;
        return inprod_infos<U>::get();
    }


    template <unsigned int ID, typename SCALARTYPE, unsigned int A>
    static infos_base &  get_infos(symbolic_vector<ID,SCALARTYPE,A> const &){
        typedef symbolic_vector<ID,SCALARTYPE,A> U;
        return vec_infos<U>::get();
    }

    template <unsigned int ID, typename SCALARTYPE, class F, unsigned int A>
    static infos_base &  get_infos(symbolic_matrix<ID,SCALARTYPE,F,A> const &){
        typedef symbolic_matrix<ID,SCALARTYPE,F,A> U;
        return mat_infos<U>::get();
    }


    template <class T, class Enable=void>
    struct wrap_expr{
    public:
        static void execute(std::list<infos_base *> & expr){
            expr.push_back(& get_infos(T()));
        }
    };

    template<class T>
    struct wrap_expr<T, typename viennacl::enable_if<result_of::is_arithmetic_compound<T>::value && !result_of::is_assignment_compound<T>::value>::type>
    {
        static void execute(std::list<infos_base *> & expr){
            wrap_expr<typename T::LHS>::execute(expr);
            wrap_expr<typename T::RHS>::execute(expr);
        }
    };

    template<class T>
    struct wrap_expr<T, typename viennacl::enable_if<result_of::is_assignment_compound<T>::value>::type>
    {
        static void execute(std::list<infos_base *> & expr){
            expr.push_back(&get_infos(typename T::LHS()));
            expr.back()->is_modified(true);
            wrap_expr<typename T::RHS>::execute(expr);
        }
    };

    /**
     * The expr_infos_base class
     */
    class expr_infos{
    public:

        typedef std::list<infos_base *> data_t;

        data_t & data() { return data_; }

        std::string generate() {
            std::string res(expression_string_base_);
            for(typename data_t::iterator it = data_.begin() ; it!= data_.end() ; ++it){
                infos_base * p = *it;
                replace_all_string(res,p->name(),p->access_name());
            }
            return res;
        }

        template<class U>
        void find_all(std::list<U* > & res) {
            for(data_t::iterator it = data_.begin() ; it != data_.end() ; ++it ){
                if(U* p = dynamic_cast<U*>(*it)){
                    res.push_back(p);
                }
            }
        }

        virtual ~expr_infos(){ }

    protected:
        expr_infos(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }
        data_t data_;
        std::string expression_string_base_;
    };




    /**
     * @brief The mat_expr_infos class
     */
    class mat_expr_infos : public expr_infos{

    };

    /**
     * @brief The vec_expr_infos_base class
     */
    class vec_expr_infos_base : public expr_infos{
    public:
        typedef vec_infos_base infos_t;
    protected:
        vec_expr_infos_base(std::string const & expression_string_base) : expr_infos(expression_string_base){ }
    };

    /**
     * @brief The vec_expr_infos class
     */
    template<class T>
    class vec_expr_infos : public vec_expr_infos_base{
    public:
        vec_expr_infos() : vec_expr_infos_base(make_expression_code<T>::value("")){
            wrap_expr<T>::execute(data_);
        }
        static vec_expr_infos_base & get(){
            static vec_expr_infos<T> res;
            return res;
        }
    };

    class scal_expr_infos_base : public expr_infos{
    public:
        typedef scal_infos_base infos_t;
    protected:
        scal_expr_infos_base(std::string const & expression_string_base) : expr_infos(expression_string_base){ }
    };



    /**
     * @brief The scal_expr_infos class
     */
    template<class T>
    class scal_expr_infos : public scal_expr_infos_base{
    public:
        scal_expr_infos() : scal_expr_infos_base(make_expression_code<T>::value("")){
            wrap_expr<T>::execute(data_);
        }
        static scal_expr_infos_base & get(){
            static scal_expr_infos<T> res;
            return res;
        }
    };

    class inprod_infos_base : public scal_infos_base{
    public:
        enum step_t{reduce, compute};
        vec_expr_infos_base & lhs(){ return lhs_; }
        vec_expr_infos_base & rhs(){ return rhs_; }
        step_t step(){ return step_; }
    protected:
        inprod_infos_base(std::string const & scalartype, std::string const & name
                               , vec_expr_infos_base & lhs, vec_expr_infos_base & rhs
                               ,step_t step) : scal_infos_base(scalartype,name), lhs_(lhs), rhs_(rhs), step_(step){        }
    private:
        vec_expr_infos_base & lhs_;
        vec_expr_infos_base & rhs_;
        step_t step_;
    };



    template<class T>
    class inprod_infos<T, typename viennacl::enable_if<result_of::is_inner_product_impl<T>::value>::type> : public inprod_infos_base{
    typedef typename T::Type U;
    public:
        inprod_infos() : inprod_infos_base(print_type<typename U::ScalarType,1>::value(), T::name()
                                                   ,vec_expr_infos<typename U::LHS>::get()
                                                   , vec_expr_infos<typename U::RHS>::get()
                                                   ,inprod_infos_base::compute){        }
        static inprod_infos_base & get(){
            static inprod_infos<T> res;
            return res;
        }
    };

    template<class T>
    class inprod_infos<T, typename viennacl::enable_if<result_of::is_inner_product_leaf<T>::value>::type> : public inprod_infos_base{
    public:
        inprod_infos() : inprod_infos_base(print_type<typename T::ScalarType,1>::value("__local"), T::name()
                                                   ,vec_expr_infos<typename T::LHS>::get()
                                                   , vec_expr_infos<typename T::RHS>::get()
                                                   ,inprod_infos_base::reduce){        }
        static inprod_infos_base & get(){
            static inprod_infos<T> res;
            return res;
        }
    };



    template<class T>
    class cache_manager{
    public:
        cache_manager(std::list<T * > const & expressions,  std::ostringstream & oss, std::set<T *>& cached_entries) : expressions_(expressions)
                                                                                                                                   , oss_(oss)
          ,cached_entries_(cached_entries){         }

        void fetch_entries(std::string const & idx){
            for(typename std::list<T * >::iterator it = expressions_.begin() ; it != expressions_.end() ; ++it){
                T * p = *it;
                if(cached_entries_.insert(p).second){
                    p->access_name(p->name()+"_val");
                    oss_ << p->scalartype() << " " << p->access_name() << " = " << p->name() << "[" << idx << "];\n";
                }
            }
        }

        void writeback_entries(std::string const & idx){
            for(typename std::list<T * >::iterator it = expressions_.begin() ; it != expressions_.end() ; ++it){
                T * p = *it;
                if(p->is_modified())
                    oss_<< p->name() << "[" << idx << "]"<< " = "  << p->access_name() << ";\n";
            }
        }

    private:
        std::list<T * > expressions_;
        std::ostringstream & oss_;
        std::set<T *> & cached_entries_;

    };

    struct blas1_generator{
    public:

        blas1_generator(std::list<vec_expr_infos_base * > & vector_expressions, std::list<scal_expr_infos_base * > & scalar_expressions): vector_expressions_(vector_expressions),
                                                                                                                                         scalar_expressions_(scalar_expressions){
            for(std::list<scal_expr_infos_base * >::const_iterator it = scalar_expressions.begin() ; it != scalar_expressions.end() ; ++it){
                for(expr_infos::data_t::iterator iit = (*it)->data().begin() ; iit != (*it)->data().end() ; ++iit){
                    if(inprod_infos_base *p = dynamic_cast<inprod_infos_base *>(*iit)){
                        if(p->step() == inprod_infos_base::compute) inner_prods_compute_.push_back(p);
                        else inner_prods_reduce_.push_back(p);
                        p->lhs().find_all(vectors_);
                        p->lhs().find_all(gpu_scalars_);
                        p->rhs().find_all(vectors_);
                        p->rhs().find_all(gpu_scalars_);
                    }
                }
                (*it)->find_all(gpu_scalars_);
            }
            for(std::list<vec_expr_infos_base *>::const_iterator it = vector_expressions.begin(); it!=vector_expressions.end() ; ++it){
                (*it)->find_all(vectors_);
                (*it)->find_all(gpu_scalars_);
            }
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
           oss << "    " << (*it)->access_name() << " = local_" << (*it)->name() << "[" << 0 << "];\n";
           }
        }

        std::string operator()(){
            std::ostringstream oss;
            std::set<vec_infos_base *> vector_cached_entries;
            std::set<gpu_scal_infos_base *> scalar_cached_entries;
            std::set<inprod_infos_base *> inprod_cached_entries;

            cache_manager<vec_infos_base> vector_cache(vectors_,oss,vector_cached_entries);
            cache_manager<gpu_scal_infos_base> scalar_cache(gpu_scalars_,oss,scalar_cached_entries);
            cache_manager<inprod_infos_base> inprod_cache(inner_prods_reduce_,oss,inprod_cached_entries);

            vec_infos_base * first_vector =  NULL;
            if(vectors_.size())
                first_vector = vectors_.front();
            //Assumes same size...
            oss << "{\n";

            if(inner_prods_reduce_.size()>0)
                compute_reductions(oss,inner_prods_reduce_);
            scalar_cache.fetch_entries("0");
            if(first_vector){
                oss << "for(unsigned int i = get_global_id(0) ; i <" << first_vector->size() << " ; i += get_global_size(0){\n";
                vector_cache.fetch_entries("i");
                vector_cache.writeback_entries("i");
                oss << "}\n";
            }
            scalar_cache.writeback_entries("0");
            oss << "}\n";

            return oss.str();
        }

    private:
        std::list<vec_expr_infos_base * >  vector_expressions_;
        std::list<scal_expr_infos_base * > scalar_expressions_;
        std::list<inprod_infos_base * > inner_prods_compute_;
        std::list<inprod_infos_base * > inner_prods_reduce_;
        std::list<vec_infos_base * > vectors_;
        std::list<gpu_scal_infos_base * > gpu_scalars_;
    };

//    template<class T>
//    class scal_infos<T, typename viennacl::enable_if<result_of::is_inner_product_leaf>::type> : public infos_base{
//        typedef typename tree_utils::extract_if<T,result_of::is_inner_product_leaf >::Result::Head InProd;
//        typedef typename tree_utils::remove_if<T,result_of::is_inner_product_leaf,false >::Result Scalars;
//        typedef typename InProd::LHS InProdLHS;
//        typedef typename InProd::RHS InProdRHS;
//    public:
//        scal_infos() : infos_base(print_type<typename T::ScalarType,1>::value()
//                                  ,T::name()) { }
//        static infos_base & get(){
//            static scal_infos<T> res;
//            return res;
//        }
//    private:
//        vec_expr_infos lhs;
//    };

//    template<class LhsT, class RhsT, class AddExprT>
//    class prod_infos_base : public infos_base{
//    public:
//        prod_infos_base(std::string const & name
//                                  ,infos_base const & assigned
//                                  ,std::string const & op_str
//                                  ,LhsT const & lhs
//                                  ,RhsT const & rhs
//                                  ,AddExprT const & additional_expression
//                                  ,std::string const & expression_string_base) : infos_base(assigned.scalartype(),name)
//                                                                              , assigned_(assigned)
//                                                                              , lhs_(lhs)
//                                                                              , op_str_(op_str)
//                                                                              , rhs_(rhs)
//                                                                              , additional_expression_(additional_expression)
//                                                                              , expression_string_base_(expression_string_base){ }
//        LhsT const & lhs() const { return lhs_; }
//        RhsT const & rhs() const { return rhs_; }
//        AddExprT const & additional_expression() const { return additional_expression_; }
//        infos_base const & assigned() const { return assigned_; }
//        std::string generate(){
//            std::string tmp(expression_string_base_);
//            replace_all_string(tmp,name_,access_name_);
//            return assigned_.access_name() + op_str_ + tmp + "+" + additional_expression_.generate();
//        }

//    private:
//        std::string op_str_;
//        infos_base const & assigned_;
//        LhsT lhs_;
//        RhsT rhs_;
//        AddExprT additional_expression_;
//        mutable std::string expression_string_base_;
//    };

//    /**
//     * @brief The matvec_prod_infos struct
//     */
//    class matvec_prod_infos : public prod_infos_base<vec_infos, mat_expr_infos, vec_expr_infos, vec_expr_infos>{
//        typedef prod_infos_base<vec_infos, mat_expr_infos, vec_expr_infos, vec_expr_infos> BaseT;
//    public:
//        matvec_prod_infos(std::string const & name
//                          ,vec_infos const & assigned
//                          ,std::string const & op_str
//                          ,mat_expr_infos const & lhs
//                          ,vec_expr_infos const & rhs
//                          ,vec_expr_infos const & additional_expression
//                          ,std::string const & expression_string_base): BaseT(name,assigned,op_str,lhs,rhs,additional_expression,expression_string_base){ }
//    };

//    /**
//     * @brief The inprod_infos struct
//     */




//    template<class U>
//    static vec_infos wrap_vec(){
//        return vec_infos(print_type<typename U::ScalarType,1>::value(),
//                                 U::name(),
//                                 U::internal_size2_name());
//    }

//    template<class U>
//    static mat_infos wrap_mat(){
//       return mat_infos(print_type<typename U::ScalarType,1>::value(),
//                                      U::name(),
//                                      U::size1_name(),
//                                      U::size2_name(),
//                                      result_of::is_row_major<U>::value,
//                                      false);
//    }

//    template<class U>
//    static scal_infos wrap_scal(){
//       return scal_infos(print_type<typename U::ScalarType,1>::value(),
//                                      U::name());
//    }

//    template<class U>
//    struct foreach_functor{
//        static void execute(mat_expr_infos & res){ res.add_infos(mat_infos<U>::get()); }
//        static void execute(vec_expr_infos & res){ res.add_infos(vec_infos<U>::get()); }
//        static void execute(scal_expr_infos & res){ res.add_infos(scal_infos<U>::get()); }
//    };

//    /**
//     * @brief wrap_vec_expr class
//     */
//    template<class U>
//    static vec_expr_infos wrap_vecexpr(){
//            typedef typename tree_utils::extract_if<U,result_of::is_symbolic_vector>::Result List;
//            vec_expr_infos res(make_expression_code<U>::value(""));
//            typelist_utils::ForEach<List,foreach_functor>::execute(res);
//            return res;
//    }

//    template<class U>
//    static mat_expr_infos wrap_matexpr(){
//            typedef typename tree_utils::extract_if<U,result_of::is_symbolic_matrix>::Result List;
//            mat_expr_infos res(make_expression_code<U>::value(""));
//            typelist_utils::ForEach<List,foreach_functor>::execute(res);
//            return res;
//    }

//    template<class U>
//    static scal_expr_infos wrap_scalexpr(){
//            typedef typename tree_utils::extract_if<U,result_of::is_symbolic_scalar>::Result List;
//            scal_expr_infos res(make_expression_code<U>::value(""));
//            typelist_utils::ForEach<List,foreach_functor>::execute(res);
//            return res;
//    }

//    template<unsigned int>
//    class inprod_infos;

//    template<>
//    class inprod_infos<1>: public prod_infos_base<infos_base,vec_expr_infos,vec_expr_infos, scal_expr_infos>{
//        typedef prod_infos_base<infos_base,vec_expr_infos,vec_expr_infos, scal_expr_infos> BaseT;
//    public:
//        inprod_infos(std::string const & name
//                          ,vec_expr_infos const & lhs
//                     ,vec_expr_infos const & rhs): BaseT(name,infos_base("",""),"",lhs,rhs,wrap_scalexpr<NullType>(),""){ }
//    };


//    template<>
//    class inprod_infos<2>: public prod_infos_base<scal_infos,vec_expr_infos,vec_expr_infos, scal_expr_infos>{
//        typedef prod_infos_base<scal_infos,vec_expr_infos,vec_expr_infos, scal_expr_infos> BaseT;
//    public:
//        inprod_infos(std::string const & name
//                          ,scal_infos const & assigned
//                          ,std::string const & op_str
//                          ,vec_expr_infos const & lhs
//                          ,vec_expr_infos const & rhs
//                          ,scal_expr_infos const & additional_expression
//                     ,std::string const & expression_string_base): BaseT(name,assigned,op_str,lhs,rhs,additional_expression,expression_string_base){ }
//    };

//    template<class T>
//    static matvec_prod_infos wrap_matvec_prod(){
//        typedef typename tree_utils::extract_if<T,result_of::is_product_leaf>::Result::Head Product;
//        typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vectors;
//        typedef typename Product::LHS ProdLHS;
//        typedef typename Product::RHS ProdRHS;
//        typedef typename result_of::expression_type<T>::Result    IntermediateType;
//        static const unsigned int Alignment = IntermediateType::Alignment;
//            VIENNACL_STATIC_ASSERT(Alignment==1,AlignmentNotSupported);
//            return matvec_prod_infos(Product::name()
//                                     ,wrap_vec<typename T::LHS>()
//                                     ,T::OP::expression_string()
//                                     ,wrap_matexpr<ProdLHS>()
//                                     ,wrap_vecexpr<ProdRHS>()
//                                     ,wrap_vecexpr<Vectors>()
//                                     ,make_expression_code<Product>::value(""));
//    }/*

//    template<class T>
//    struct wrap_inprod{

//            template<class U>
//            struct functor_compute{
//                typedef typename U::Type InProd;
//                typedef typename InProd::LHS InProdLHS;
//                typedef typename InProd::RHS InProdRHS;

//                static void execute(std::vector<inprod_infos<1> > & arg){
//                    arg.push_back(inprod_infos<1>(InProd::name()
//                                                  ,wrap_vecexpr<InProdLHS>()
//                                                  ,wrap_vecexpr<InProdRHS>()));
//                }
//            };

//            template<class U>
//            struct functor_reduce{
//                typedef typename tree_utils::extract_if<T,result_of::is_inner_product_leaf >::Result::Head InProd;
//                typedef typename tree_utils::remove_if<T,result_of::is_inner_product_leaf,false >::Result Scalars;
//                typedef typename InProd::LHS InProdLHS;
//                typedef typename InProd::RHS InProdRHS;
//                typedef typename result_of::expression_type<T>::Result    IntermediateType;


//                static void execute(std::vector<inprod_infos<2> > & arg){
//                    arg.push_back(inprod_infos<2>(InProd::name()
//                                                  ,wrap_scal<typename T::LHS>()
//                                                  ,T::OP::expression_string()
//                                                  ,wrap_vecexpr<InProdLHS>()
//                                                  ,wrap_vecexpr<InProdRHS>()
//                                                  ,wrap_scalexpr<Scalars>()
//                                                  ,make_expression_code<InProd>::value("")));
//                }
//            };

//            static void execute(std::vector<inprod_infos<1> > & arg_compute, std::vector<inprod_infos<2> > & arg_reduce){
//                typelist_utils::ForEach< typename tree_utils::extract_if<T,result_of::is_inner_product_impl>::Result, functor_compute>::execute(arg_compute);
//                typelist_utils::ForEach< typename tree_utils::extract_if<T,result_of::is_inner_product_leaf>::Result, functor_reduce>::execute(arg_reduce);
//            }
//    };

//    static std::string generate_matvec_prod_backend(std::vector<matvec_prod_infos> const & infos){
//        std::string res;
//        mat_infos const & mat = infos.front().lhs().data().front();
//        vec_infos const & vec = infos.front().rhs().data().front();
//        res += "   __local shared_memory_ptr[64];\n";
//        res += "   unsigned int row_gid = get_global_id(0)/get_local_size(0);\n";
//        res += "   unsigned int col_gid = get_global_id(0)%get_local_size(0);\n";
//        res += "   unsigned int lid = get_local_id(0);\n";
//        res += "   for(unsigned int row = row_gid ; row < " + mat.size1() + " ; row+=get_num_groups(0)){\n";
//        res += "       " + mat.scalartype() + " sum = 0;\n";
//        res += "       for(unsigned int col = col_gid ; col < " + mat.size2() + " ; col+=get_local_size(0)){\n";
//        mat.access_name(mat.name()+"row");
//        vec.access_name("col");
//        res += "            sum +=  " + infos.front().lhs().generate() + "*" +  infos.front().rhs().generate() + ";\n";
//        res += "       }\n";
//        res += "       shared_memory_ptr[lid]=sum;\n";
//        res += "       for(unsigned int stride=get_local_size(0)/2 ; stride>0 ; stride>>=1){\n";
//        res += "           barrier(CLK_LOCAL_MEM_FENCE);\n";
//        res += "           if(lid < stride) shared_memory_ptr[lid]+=shared_memory_ptr[lid+stride];\n";
//        res += "       }\n";
//        res += "       if(lid==0) shared_memory_ptr[0] ;\n";
//        res+=";\n";
//        res += "   }\n";
//        res += "}\n";
//        return res;
//    }


//    struct blas1_generator{

//        std::string operator()(std::vector<inprod_infos<2> > & inprods_to_reduce,
//                                std::vector<vec_expr_infos> & vec_exprs,
//                                std::vector<inprod_infos<1> > & inprods_to_compute){
//            bool has_inprods_to_reduce = inprods_to_reduce.size();
//            bool has_vec_exprs = vec_exprs.size();
//            bool has_inprods_to_compute = inprods_to_compute.size();
//            std::string size;
//            if(has_vec_exprs)
//                size = vec_exprs.front().data().front().size() ;
//            else if(has_inprods_to_compute)
//                size = inprods_to_compute.front().lhs().data().front().size();
//            else
//                size = inprods_to_reduce.front().lhs().data().front().size();


//            std::ostringstream oss;
//            oss << "{\n";

//            //Reduces inner products previously computed...
//            if(has_inprods_to_reduce){
//                for( std::vector<inprod_infos<2> >::const_iterator it = inprods_to_reduce.begin(); it != inprods_to_reduce.end() ; ++it){
//                    it->access_name(it->name()+"_val");
//                    oss << "   float " << it->access_name() << " = 0;\n " ;
//                }
//                oss << "   {\n";
//                for( std::vector<inprod_infos<2> >::const_iterator it = inprods_to_reduce.begin(); it != inprods_to_reduce.end() ; ++it){
//                    oss << "   float sum_" << it->name() << "  = 0;\n";
//                    oss << "   __local local_" << it->name() << "[64];\n";
//                }
//                oss << "   for (unsigned int i = get_local_id(0) ; i<get_num_groups(0) ; i+=get_local_size(0))\n";
//                oss << "   {\n";
//                for( std::vector<inprod_infos<2> >::const_iterator it = inprods_to_reduce.begin(); it != inprods_to_reduce.end() ; ++it){
//                    oss << " sum_" << it->name() << "  += " << it->name() << "[i];\n";
//                }
//                oss << "   }\n";
//                for( std::vector<inprod_infos<2> >::const_iterator it = inprods_to_reduce.begin(); it != inprods_to_reduce.end() ; ++it){
//                    oss << " local_" << it->name() << "[get_local_id(0)]  = " << "sum_" << it->name() << ";\n";
//                }
//                oss << "   for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n";
//                oss << "   {\n";
//                oss << "      barrier(CLK_LOCAL_MEM_FENCE);\n";
//                oss << "      if (get_local_id(0) < stride){\n";
//                for( std::vector<inprod_infos<2> >::const_iterator it = inprods_to_reduce.begin(); it != inprods_to_reduce.end() ; ++it){
//                oss << "         local_" << it->name() << "[get_local_id(0)]  += " << "local_" << it->name() << "shared_memory_ptr[get_local_id(0)+stride];\n";
//                }
//                oss << "      }\n";
//                oss << "   }\n";
//                for( std::vector<inprod_infos<2> >::const_iterator it = inprods_to_reduce.begin(); it != inprods_to_reduce.end() ; ++it){
//                oss << "    " << it->access_name() << " = local_" << it->name() << "[" << 0 << "];\n";
//                }
//            }

//            if(has_inprods_to_compute || has_vec_exprs){
//                oss << "for(unsigned int i = get_global_id(0) ; i < " << size << " i+= get_global_size(0){\n";
//                for( std::vector<inprod_infos<1> >::const_iterator it = inprods_to_compute.begin(); it != inprods_to_compute.end() ; ++it){
//                    it->access_name(it->name()+"_val");
//                    oss << it->scalartype() << " " << it->access_name() << " = " << it->name() << "[i];\n";
//                }
//                oss << "}\n";
//            }

//            oss << "}\n";
//            return oss.str();
//        }


//    };

    };

    template <class T>
    struct make_code<InProdToken<T, 0> >
    {
      template<class U>
      struct generate_code
      {
        private:
          typedef typename U::LHS LHS;
          typedef typename U::RHS RHS;

        public:

          static void execute(std::string & generated_code)
          {
              generated_code+=
                      "{\n"
                      "   __local shared_memory_ptr[64];\n"
                      "   float sum = 0;\n"
                      "   for (unsigned int i = get_local_id(0) ; i<get_num_groups(0) ; i+=get_local_size(0))\n"
                      "   {\n"
                      "      sum+= " +U::name() +"[i];\n"
                      "   };\n"
                      "   shared_memory_ptr[get_local_id(0)]=sum;\n"
                      "   for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
                      "   {\n"
                      "      barrier(CLK_LOCAL_MEM_FENCE);\n"
                      "      if (get_local_id(0) < stride)\n"
                      "      shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
                      "   }\n"
                      "   if(get_local_id(0)==0)\n"
                      "       "+U::local_value() + " = shared_memory_ptr[0];\n"
                      "   barrier(CLK_LOCAL_MEM_FENCE);\n"
                      "}\n";
          }
      };

      static std::string value()
      {
          std::string res;
          typelist_utils::ForEach<T,generate_code>::execute(res);
          return res;
      }
    };

    void replace_all_string(std::string & str, std::string const & from, std::string const & to){
        size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                 str.replace(start_pos, from.length(), to);
                 start_pos += to.length();
        }
    }

    template<class T, class OP, class Assigned>
    struct make_code<MatVecToken<T,OP,Assigned> >
    {
      private:

        typedef typename tree_utils::extract_if<T,result_of::is_product_leaf>::Result::Head Products;
        typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vectors;

        typedef typename result_of::expression_type<T>::Result    IntermediateType;
        static const unsigned int Alignment = IntermediateType::Alignment;

        std::string const & size1() const{ return size1_; }
        std::string const & size2() const{ return size2_; }
        std::string const & name() const{ return name_; }
        std::string const & scalartype() const{ return scalartype_; }
        bool const is_rowmajor() const { return is_rowmajor_; }
        bool const is_transposed() const { return is_transposed_; }

    private:
        std::string scalartype_;
        std::string name_;
        std::string size1_;
        std::string size2_;
        bool is_rowmajor_;
        bool is_transposed_;
        mutable std::string expression_string_base_;
    };

    class vec_infos{
    public:
        vec_infos(std::string const & scalartype,
                          std::string const & name,
                          std::string const & size,
                    std::string const & expression_string_base) : scalartype_(scalartype), name_(name), size_(size), expression_string_base_(expression_string_base) { }




    template<class T>
    struct make_code<MatVecToken<T> >
    {

      template<class U>
      struct fill_infos{
          static void execute(std::vector<matvec_prod_infos> & res){
              res.push_back(wrap_matvec_prod<U>::create());
          }
      };

//    /** @brief Functor to make code from a token
//        @tparam TOKEN The token to make the code for
//    */
//    template <class TOKEN>
//    struct make_code;

//    template<>
//    struct make_code<NullType>
//    {
//      static std::string value() { return ""; }

//      static std::string sum() { return ""; }

//      static std::string reduction() { return ""; }
//    };

//    template <class EXPR>
//    struct make_code<ArithmeticToken<EXPR> >
//    {
//      static std::string value()
//      {
//        std::string res;
//        res+="\n//Arithmetic Token\n";
//        res+=make_expression_code<EXPR>::value("gid") + ";\n";
//        return res;
//      }
//    };

//    template <class T>
//    struct make_code<InProdToken<T, 1>  >
//    {
//      template<class U>
//      struct generate_code_sum
//      {
//        private:
//          typedef typename U::Type ARG;
//          typedef typename ARG::LHS LHS;
//          typedef typename ARG::RHS RHS;
//        public:
//          static void execute(std::string & generated_code)
//          {
////              generated_code+= U::private_value() + " += " + dot_product<LHS,RHS>::value("gid","gid") + ";\n";
//          }
//      };


//      template<class U>
//      struct generate_code_reduction
//      {
//        private:
//          typedef typename U::Type ARG;
//          typedef typename ARG::LHS LHS;
//          typedef typename ARG::RHS RHS;
//        public :

//          static void execute(std::string & generated_code)
//          {
//              generated_code+=
//                      "shared_memory_ptr[get_local_id(0)] = " + U::private_value() + ";\n"
//                      "for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
//                      "{\n"
//                      "  barrier(CLK_LOCAL_MEM_FENCE);\n"
//                      "  if (get_local_id(0) < stride)\n"
//                      "  shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
//                      "}\n"
//                      +ARG::name() + "[get_group_id(0)] = shared_memory_ptr[0];";
//          }
//      };

//      static std::string sum()
//      {
//          std::string res;
//          typelist_utils::ForEach<T,generate_code_sum>::execute(res);
//          return res;
//      }

//      static std::string reduction()
//      {
//          std::string res;
//          typelist_utils::ForEach<T,generate_code_reduction>::execute(res);
//          return res;
//      }

//    };

//    template <class T>
//    struct make_code<InProdToken<T, 0> >
//    {
//      template<class U>
//      struct generate_code
//      {
//        private:
//          typedef typename U::LHS LHS;
//          typedef typename U::RHS RHS;

//        public:

//          static void execute(std::string & generated_code)
//          {
//              generated_code+=
//                      "{\n"
//                      "   __local shared_memory_ptr[64];\n"
//                      "   float sum = 0;\n"
//                      "   for (unsigned int i = get_local_id(0) ; i<get_num_groups(0) ; i+=get_local_size(0))\n"
//                      "   {\n"
//                      "      sum+= " +U::name() +"[i];\n"
//                      "   };\n"
//                      "   shared_memory_ptr[get_local_id(0)]=sum;\n"
//                      "   for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
//                      "   {\n"
//                      "      barrier(CLK_LOCAL_MEM_FENCE);\n"
//                      "      if (get_local_id(0) < stride)\n"
//                      "      shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
//                      "   }\n"
//                      "   if(get_local_id(0)==0)\n"
//                      "       "+U::local_value() + " = shared_memory_ptr[0];\n"
//                      "   barrier(CLK_LOCAL_MEM_FENCE);\n"
//                      "}\n";
//          }
//      };

//      static std::string value()
//      {
//          std::string res;
//          typelist_utils::ForEach<T,generate_code>::execute(res);
//          return res;
//      }
//    };





//    template<class T>
//    struct make_code<MatVecToken<T> >
//    {

//      template<class U>
//      struct fill_infos{
//          static void execute(std::vector<matvec_prod_infos> & res){
//              res.push_back(wrap_matvec_prod<U>());
//          }
//      };

//      public:
//        static std::string value()
//        {
//          std::vector<matvec_prod_infos> infos;
//          typelist_utils::ForEach<T,fill_infos>::execute(infos);
//          return generate_matvec_prod_backend(infos);
//        }

//    };


  }

}

#endif

