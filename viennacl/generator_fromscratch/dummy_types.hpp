#ifndef VIENNACL_GENERATOR_DUMMY_TYPES_HPP
#define VIENNACL_GENERATOR_DUMMY_TYPES_HPP

#include "viennacl/meta/enable_if.hpp"
#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/vector.hpp"
#include <set>

namespace viennacl{

namespace generator{

template<class LHS, class OP, class RHS>
class compile_time_beast{
public:
    typedef LHS Lhs;
    typedef OP Op;
    typedef RHS Rhs;
    LHS const & lhs() const{return lhs_;}
    RHS const & rhs() const{return rhs_;}
    OP const & op() const{ return op_; }
protected:
    compile_time_beast(LHS const & lhs, RHS const & rhs) : lhs_(lhs), rhs_(rhs){ }
private:
    LHS const & lhs_;
    RHS const & rhs_;
    OP  op_;
};

struct inprod_type{ };


template<class LHS, class OP, class RHS>
class vector_expression_wrapper : public compile_time_beast<LHS,OP,RHS>{
public: vector_expression_wrapper(LHS const & lhs, RHS const & rhs) : compile_time_beast<LHS,OP,RHS>(lhs,rhs){ }
};
template<class LHS, class OP, class RHS>
class scalar_expression_wrapper : public compile_time_beast<LHS,OP,RHS>{
public: scalar_expression_wrapper(LHS const & lhs, RHS const & rhs) : compile_time_beast<LHS,OP,RHS>(lhs,rhs){ }
};
template<class LHS, class OP, class RHS>
class matrix_expression_wrapper : public compile_time_beast<LHS,OP,RHS>{
public: matrix_expression_wrapper(LHS const & lhs, RHS const & rhs) : compile_time_beast<LHS,OP,RHS>(lhs,rhs){ }
};
template<class LHS, class RHS>
class inner_prod_wrapper : public compile_time_beast<LHS,inprod_type,RHS>{
public: inner_prod_wrapper(LHS const & lhs, RHS const & rhs) : compile_time_beast<LHS,inprod_type, RHS>(lhs,rhs){ }
};


template<class T1, class T2=void, class T3=void>
struct function_wrapper_impl{
    function_wrapper_impl(unsigned int _func_id, std::string const & _name, std::string const & _expr) : func_id(_func_id), name(_name), expr(_expr), t1(NULL), t2(NULL), t3(NULL){ }
    unsigned int func_id;
    std::string name;
    std::string expr;
    T1 const * t1;
    T2 const * t2;
    T3 const * t3;
};




class function_wrapper{
private:
    static unsigned int create_func_id(){
        static unsigned long i = 0;
        return i++;
    }
public:
    function_wrapper(std::string const & name
                     ,std::string const & expr) : name_(name), func_id_(create_func_id()), expr_(expr){
        n_args_ = 0;
        std::cout << func_id_ << std::endl;
        bool keep_going = true;
        while(keep_going){
            std::string current_arg = "_"+to_string(n_args_+1)+"_";
            if(expr_.find(current_arg)!=std::string::npos)
                ++n_args_;
            else
                keep_going=false;
        }
        assert(n_args_>0 && "\nNo argument specified for the function\n"
                            "\nRecall : 1st arg : _1_\n"
                            "\n         2nd arg : _2_\n"
                                      "...");
    }

    template<class T1>
    function_wrapper_impl<T1> operator()(T1 const & t1){
        assert(n_args_==1);
        function_wrapper_impl<T1> res(func_id_,name_,expr_);
        res.t1 = &t1;
        return res;
    }

    template<class T1, class T2>
    function_wrapper_impl<T1,T2> operator()(T1 const & t1, T2 const & t2){
        assert(n_args_==2);
        function_wrapper_impl<T1, T2> res(func_id_,name_,expr_);
        res.t1 = &t1; res.t2 = &t2;
        return res;
    }

    template<class T1, class T2, class T3>
    function_wrapper_impl<T1,T2, T3> operator()(T1 const & t1, T2 const & t2, T3 const & t3){
        assert(n_args_==3);
        function_wrapper_impl<T1, T2,T3> res(func_id_,name_,expr_);
        res.t1 = &t1; res.t2 = &t2; res.t3 = &t3;
        return res;
    }

private:
    std::string name_;
    unsigned int func_id_;
    std::string expr_;
    unsigned int n_args_;
};

template<typename SCALARTYPE, unsigned int Alignment=1>
class dummy_vector{
    typedef dummy_vector<SCALARTYPE> self_type;
    typedef viennacl::vector<SCALARTYPE,Alignment> vcl_vec_t;
    vcl_vec_t const & vec_;
public:

    dummy_vector(vcl_vec_t const & vec): vec_(vec){ }

    vcl_vec_t const & vec() const{ return vec_; }

    template<typename RHS_TYPE>
    vector_expression_wrapper<self_type, assign_type, RHS_TYPE >
    operator= ( RHS_TYPE const & rhs ){
      return vector_expression_wrapper<self_type,assign_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    vector_expression_wrapper<self_type, inplace_scal_mul_type, RHS_TYPE >
    operator*= ( RHS_TYPE const & rhs ){
      return vector_expression_wrapper<self_type,inplace_scal_mul_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    vector_expression_wrapper<self_type, inplace_scal_div_type, RHS_TYPE >
    operator/= ( RHS_TYPE const & rhs ){
      return vector_expression_wrapper<self_type,inplace_scal_div_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    vector_expression_wrapper<self_type, inplace_add_type, RHS_TYPE >
    operator+= ( RHS_TYPE const & rhs ){
      return vector_expression_wrapper<self_type,inplace_add_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    vector_expression_wrapper<self_type, inplace_sub_type, RHS_TYPE >
    operator-= ( RHS_TYPE const & rhs ){
      return vector_expression_wrapper<self_type,inplace_sub_type,RHS_TYPE >(*this,rhs);
    }
};

template<class ScalarType>
class dummy_scalar{
    typedef dummy_scalar<ScalarType> self_type;
public:
    template<typename RHS_TYPE>
    scalar_expression_wrapper<self_type, assign_type, RHS_TYPE >
    operator= ( RHS_TYPE const & rhs ){
      return scalar_expression_wrapper<self_type,assign_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    scalar_expression_wrapper<self_type, inplace_scal_mul_type, RHS_TYPE >
    operator*= ( RHS_TYPE const & rhs ){
      return scalar_expression_wrapper<self_type,inplace_scal_mul_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    scalar_expression_wrapper<self_type, inplace_scal_div_type, RHS_TYPE >
    operator/= ( RHS_TYPE const & rhs ){
      return scalar_expression_wrapper<self_type,inplace_scal_div_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    scalar_expression_wrapper<self_type, inplace_add_type, RHS_TYPE >
    operator+= ( RHS_TYPE const & rhs ){
      return scalar_expression_wrapper<self_type,inplace_add_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    scalar_expression_wrapper<self_type, inplace_sub_type, RHS_TYPE >
    operator-= ( RHS_TYPE const & rhs ){
      return scalar_expression_wrapper<self_type,inplace_sub_type,RHS_TYPE >(*this,rhs);
    }
};


template<class ScalarType>
class dummy_matrix{
    typedef dummy_matrix<ScalarType> self_type;
public:
    template<typename RHS_TYPE>
    matrix_expression_wrapper<self_type, assign_type, RHS_TYPE >
    operator= ( RHS_TYPE const & rhs ){
      return matrix_expression_wrapper<self_type,assign_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    matrix_expression_wrapper<self_type, inplace_scal_mul_type, RHS_TYPE >
    operator*= ( RHS_TYPE const & rhs ){
      return matrix_expression_wrapper<self_type,inplace_scal_mul_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    matrix_expression_wrapper<self_type, inplace_scal_div_type, RHS_TYPE >
    operator/= ( RHS_TYPE const & rhs ){
      return matrix_expression_wrapper<self_type,inplace_scal_div_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    matrix_expression_wrapper<self_type, inplace_add_type, RHS_TYPE >
    operator+= ( RHS_TYPE const & rhs ){
      return matrix_expression_wrapper<self_type,inplace_add_type,RHS_TYPE >(*this,rhs);
    }

    template<typename RHS_TYPE>
    matrix_expression_wrapper<self_type, inplace_sub_type, RHS_TYPE >
    operator-= ( RHS_TYPE const & rhs ){
      return matrix_expression_wrapper<self_type,inplace_sub_type,RHS_TYPE >(*this,rhs);
    }
};


template<class T>
struct is_vector_expression_t{ enum { value = 0 }; };
template<typename ScalarType, unsigned int Alignment>
struct is_vector_expression_t<dummy_vector<ScalarType,Alignment> >{ enum { value = 1}; };
template<class LHS, class OP, class RHS>
struct is_vector_expression_t<vector_expression_wrapper<LHS,OP,RHS> >{ enum { value = 1}; };

template<class T>
struct is_scalar_expression_t{ enum { value = 0 }; };
template<class ScalarType>
struct is_scalar_expression_t<dummy_scalar<ScalarType> >{ enum { value = 1}; };
template<class LHS, class OP, class RHS>
struct is_scalar_expression_t<scalar_expression_wrapper<LHS,OP,RHS> >{ enum { value = 1}; };
template<class LHS, class RHS>
struct is_scalar_expression_t<inner_prod_wrapper<LHS,RHS> >{ enum { value = 1}; };
template<class T1, class T2, class T3>
struct is_scalar_expression_t<function_wrapper_impl<T1,T2,T3> >{ enum { value = 1}; };


template<class T>
struct is_matrix_expression_t{ enum { value = 0 }; };
template<class ScalarType>
struct is_matrix_expression_t<dummy_matrix<ScalarType> >{ enum { value = 1}; };
template<class LHS, class OP, class RHS>
struct is_matrix_expression_t<matrix_expression_wrapper<LHS,OP,RHS> >{ enum { value = 1}; };

template<class LHS, class OP, class RHS, bool create_vector, bool create_scalar, bool create_matrix>
struct convert_to_expr;
template<class LHS, class OP, class RHS>
struct convert_to_expr<LHS,OP,RHS,true,false,false>{ typedef vector_expression_wrapper<LHS,OP,RHS> type; };
template<class LHS, class OP, class RHS>
struct convert_to_expr<LHS,OP,RHS,false,true,false>{ typedef scalar_expression_wrapper<LHS,OP,RHS> type; };
template<class LHS, class OP, class RHS>
struct convert_to_expr<LHS,OP,RHS,false,false,true>{ typedef matrix_expression_wrapper<LHS,OP,RHS> type; };

template<class T>
struct is_operator{ enum{ value = 0}; };
template<> struct is_operator<assign_type>{ enum { value = 1}; };
template<> struct is_operator<add_type>{ enum { value = 1}; };
template<> struct is_operator<inplace_add_type>{ enum { value = 1}; };
template<> struct is_operator<sub_type>{ enum { value = 1}; };
template<> struct is_operator<inplace_sub_type>{ enum { value = 1}; };
template<> struct is_operator<scal_mul_type>{ enum { value = 1}; };
template<> struct is_operator<inplace_scal_mul_type>{ enum { value = 1}; };
template<> struct is_operator<scal_div_type>{ enum { value = 1}; };
template<> struct is_operator<inplace_scal_div_type>{ enum { value = 1}; };

template<class T>
struct is_leaf{ enum{ value = 0}; };
template<class ScalarType> struct is_leaf<dummy_vector<ScalarType> >{ enum { value = 1 }; };
template<class ScalarType> struct is_leaf<dummy_scalar<ScalarType> >{ enum { value = 1 }; };
template<class ScalarType> struct is_leaf<dummy_matrix<ScalarType> >{ enum { value = 1 }; };

//template<class T>
//unary_minus<T> operator -(T const &)
//{
//  return unary_minus<T>();
//}


template<class LHS, class RHS> struct create_vector{
    enum{  value= (is_vector_expression_t<LHS>::value && is_scalar_expression_t<RHS>::value)
         || (is_scalar_expression_t<LHS>::value && is_vector_expression_t<RHS>::value)
         || (is_vector_expression_t<LHS>::value && is_vector_expression_t<RHS>::value) };
};

template<class LHS, class RHS> struct create_scalar{
    enum{  value= (is_scalar_expression_t<LHS>::value && is_scalar_expression_t<RHS>::value) };
};


template<class LHS, class RHS> struct create_matrix{
    enum{  value= (is_matrix_expression_t<LHS>::value && is_scalar_expression_t<RHS>::value)
         || (is_scalar_expression_t<LHS>::value && is_matrix_expression_t<RHS>::value)
         || (is_matrix_expression_t<LHS>::value && is_matrix_expression_t<RHS>::value) };
};


template<class LHS, class RHS>
typename viennacl::enable_if<is_vector_expression_t<LHS>::value && is_vector_expression_t<RHS>::value
                            ,inner_prod_wrapper<LHS,RHS> >::type
inner_prod(LHS const & lhs, RHS const & rhs)
{
    return inner_prod_wrapper<LHS,RHS>(lhs,rhs);
}

template<class LHS, class RHS>
typename viennacl::enable_if< (is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value)
                             ||(is_vector_expression_t<LHS>::value && is_vector_expression_t<RHS>::value)
                             ||(is_matrix_expression_t<LHS>::value && is_matrix_expression_t<RHS>::value)
                            ,typename convert_to_expr<LHS,add_type,RHS
                                                    ,create_vector<LHS,RHS>::value
                                                    ,create_scalar<LHS,RHS>::value
                                                    ,create_matrix<LHS,RHS>::value>::type>::type
operator+(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,add_type,RHS
            ,create_vector<LHS,RHS>::value
            ,create_scalar<LHS,RHS>::value
            ,create_matrix<LHS,RHS>::value>::type(lhs,rhs);
}

template<class LHS, class RHS>
typename viennacl::enable_if< (is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value)
                             ||(is_vector_expression_t<LHS>::value && is_vector_expression_t<RHS>::value)
                             ||(is_matrix_expression_t<LHS>::value && is_matrix_expression_t<RHS>::value)
                            ,typename convert_to_expr<LHS,sub_type,RHS
                                                    ,create_vector<LHS,RHS>::value
                                                    ,create_scalar<LHS,RHS>::value
                                                    ,create_matrix<LHS,RHS>::value>::type>::type
operator-(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,sub_type,RHS
            ,create_vector<LHS,RHS>::value
            ,create_scalar<LHS,RHS>::value
            ,create_matrix<LHS,RHS>::value>::type(lhs,rhs);
}
template<class LHS, class RHS>
typename viennacl::enable_if< is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
                            ,typename convert_to_expr<LHS,scal_mul_type,RHS
                            ,create_vector<LHS,RHS>::value
                            ,create_scalar<LHS,RHS>::value
                            ,create_matrix<LHS,RHS>::value>::type>::type
operator*(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,scal_mul_type,RHS
                                    ,create_vector<LHS,RHS>::value
                                    ,create_scalar<LHS,RHS>::value
                                    ,create_matrix<LHS,RHS>::value>::type(lhs,rhs);
}

template<class LHS, class RHS>
typename viennacl::enable_if< is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
                            ,typename convert_to_expr<LHS,scal_div_type,RHS
                                                    ,is_vector_expression_t<LHS>::value || is_vector_expression_t<RHS>::value
                                                    ,is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
                                                    ,is_matrix_expression_t<LHS>::value || is_matrix_expression_t<RHS>::value>::type>::type
operator/(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,scal_div_type,RHS
            ,is_vector_expression_t<LHS>::value || is_vector_expression_t<RHS>::value
            ,is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
            ,is_matrix_expression_t<LHS>::value || is_matrix_expression_t<RHS>::value>::type(lhs,rhs);
}

//static const unsigned int N_BITS_LEAF = 6;
//static const unsigned int N_BITS_OP = 4;

//template<class T, class Enable=void> struct n_bits;
//template<class T> struct n_bits<T, typename enable_if<is_operator<T>::value>::type> { enum{ value = N_BITS_OP}; };
//template<class T> struct n_bits<T, typename enable_if<is_leaf<T>::value>::type> { enum{ value = N_BITS_LEAF}; };
//template<class T1> struct n_bits<function_wrapper_impl<T1> >{ enum { value = N_BITS_LEAF*2}; };
//template<class T1, class T2> struct n_bits<function_wrapper_impl<T1, T2> >{ enum { value = N_BITS_LEAF*3}; };
//template<class T1, class T2, class T3> struct n_bits<function_wrapper_impl<T1, T2, T3> >{ enum { value = N_BITS_LEAF*4}; };

//template<class TREE_T>
//struct get_operation_id{
//    typedef typename TREE_T::Lhs Lhs;
//    typedef typename TREE_T::Rhs Rhs;
//    typedef typename TREE_T::Op Op;
//    static unsigned long value(TREE_T const & tree){
//        return get_operation_id<Lhs>::value(tree.lhs())
//                | get_operation_id<Op>::value(tree.op()) << (n_bits<Lhs>::value)
//                | get_operation_id<Rhs>::value(tree.rhs())<< (n_bits<Lhs>::value + n_bits<Op>::value);
//    }
//};

//#define DEFINE_ID_FOR(TYPE, VALUE) \
//    struct get_operation_id< TYPE >{\
//        static unsigned long value(TYPE const & t){ return VALUE; }\
//    }

////Defines Operators id
//template<> DEFINE_ID_FOR(assign_type, 0);
//template<> DEFINE_ID_FOR(add_type, 1);
//template<> DEFINE_ID_FOR(inplace_add_type,2);


//template<> DEFINE_ID_FOR(dummy_vector<float>, 0);
//template<> DEFINE_ID_FOR(dummy_vector<double>, 1);
//template<class T1>
//struct get_operation_id<function_wrapper_impl<T1> >{
//    static unsigned long value(function_wrapper_impl<T1> const & t){
//        return 2 | get_operation_id<T1>::value(*t.t1) << N_BITS_LEAF;
//    }
//};
//template<class T1, class T2>
//struct get_operation_id<function_wrapper_impl<T1, T2> >{
//    static unsigned long value(function_wrapper_impl<T1,T2> const & t){
//        return 3
//                | get_operation_id<T1>::value(*t.t1) << N_BITS_LEAF
//                | get_operation_id<T2>::value(*t.t2) << N_BITS_LEAF*2;
//    }
//};
//template<class T1, class T2, class T3>
//struct get_operation_id<function_wrapper_impl<T1, T2, T3> >{
//    static unsigned long value(function_wrapper_impl<T1,T2, T3> const & t){
//        return 4
//                | get_operation_id<T1>::value(*t.t1) << N_BITS_LEAF
//                | get_operation_id<T2>::value(*t.t2) << N_BITS_LEAF*2
//                | get_operation_id<T2>::value(*t.t3) << N_BITS_LEAF*3;
//    }
//};

std::string encode_to_kernel_name(unsigned long id){
    const char *digs = "01234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabzdefghijklmnopqrstuvwxyz_";
    unsigned int n = strlen(digs);
    char buffer[16];
    int i = 0;
    do  buffer[i++] = digs[id%n];  while((id /= n)>0);
    buffer[i] = 0;
    return std::string(buffer);
}

/*
template<> static unsigned long get_operation_id<assign_type>(assign_type const &){ return 0; }
template<> static unsigned long get_operation_id<add_type>(add_type const &){ return 1; }


template<> static unsigned long get_operation_id<dummy_vector<float> >(dummy_vector<float> const &){ return 0; }
template<> static unsigned long get_operation_id<dummy_matrix<float> >(dummy_matrix<float> const &){ return 1; }
template<> static unsigned long get_operation_id<dummy_scalar<float> >(dummy_scalar<float> const &){ return 2; }
template<class T1> static unsigned long get_operation_id<function_wrapper_impl<T1> >(function_wrapper_impl<T1> const &){ return 3; }
*/


#define MAKE_BUILTIN_FUNCTION1(name) static function_wrapper name = function_wrapper(#name,#name "(_1_)")
#define MAKE_BUILTIN_FUNCTION2(name) static function_wrapper name = function_wrapper(#name,#name "(_1_,_2_)")
#define MAKE_BUILTIN_FUNCTION3(name) static function_wrapper name = function_wrapper(#name,#name "(_1_,_2_,_3_)")

MAKE_BUILTIN_FUNCTION1(acos);
MAKE_BUILTIN_FUNCTION1(acosh);
MAKE_BUILTIN_FUNCTION1(acospi);
MAKE_BUILTIN_FUNCTION1(asin);
MAKE_BUILTIN_FUNCTION1(asinh);
MAKE_BUILTIN_FUNCTION1(asinpi);
MAKE_BUILTIN_FUNCTION1(atan);
MAKE_BUILTIN_FUNCTION2(atan2);
MAKE_BUILTIN_FUNCTION1(atanh);
MAKE_BUILTIN_FUNCTION1(atanpi);
MAKE_BUILTIN_FUNCTION2(atan2pi);
MAKE_BUILTIN_FUNCTION1(cbrt);
MAKE_BUILTIN_FUNCTION1(ceil);
MAKE_BUILTIN_FUNCTION2(copysign);
MAKE_BUILTIN_FUNCTION1(cos);
MAKE_BUILTIN_FUNCTION1(cosh);
MAKE_BUILTIN_FUNCTION1(cospi);
MAKE_BUILTIN_FUNCTION1(erfc);
MAKE_BUILTIN_FUNCTION1(erf);
MAKE_BUILTIN_FUNCTION1(exp);
MAKE_BUILTIN_FUNCTION1(exp2);
MAKE_BUILTIN_FUNCTION1(exp10);
MAKE_BUILTIN_FUNCTION1(expm1);
MAKE_BUILTIN_FUNCTION1(fabs);
MAKE_BUILTIN_FUNCTION2(fdim);
MAKE_BUILTIN_FUNCTION1(floor);
MAKE_BUILTIN_FUNCTION3(fma);
MAKE_BUILTIN_FUNCTION2(fmax);
MAKE_BUILTIN_FUNCTION2(fmin);
MAKE_BUILTIN_FUNCTION2(fmod);
//    MAKE_BUILTIN_FUNCTION1(fract);
//    MAKE_BUILTIN_FUNCTION1(frexp);
MAKE_BUILTIN_FUNCTION2(hypot);
MAKE_BUILTIN_FUNCTION1(ilogb);
MAKE_BUILTIN_FUNCTION2(ldexp);
MAKE_BUILTIN_FUNCTION1(lgamma);
//    MAKE_BUILTIN_FUNCTION1(lgamma_r);
MAKE_BUILTIN_FUNCTION1(log);
MAKE_BUILTIN_FUNCTION1(log2);
MAKE_BUILTIN_FUNCTION1(log10);
MAKE_BUILTIN_FUNCTION1(log1p);
MAKE_BUILTIN_FUNCTION1(logb);
MAKE_BUILTIN_FUNCTION3(mad);
//    MAKE_BUILTIN_FUNCTION1(modf);
MAKE_BUILTIN_FUNCTION1(nan);
MAKE_BUILTIN_FUNCTION2(nextafter);
MAKE_BUILTIN_FUNCTION2(pow);
MAKE_BUILTIN_FUNCTION2(pown);
MAKE_BUILTIN_FUNCTION2(powr);
MAKE_BUILTIN_FUNCTION2(remainder);
//    MAKE_BUILTIN_FUNCTION1(remquo);
MAKE_BUILTIN_FUNCTION1(rint);
MAKE_BUILTIN_FUNCTION1(rootn);
MAKE_BUILTIN_FUNCTION1(round);
MAKE_BUILTIN_FUNCTION1(rsqrt);
MAKE_BUILTIN_FUNCTION1(sin);
//    MAKE_BUILTIN_FUNCTION1(sincos);
MAKE_BUILTIN_FUNCTION1(sinh);
MAKE_BUILTIN_FUNCTION1(sinpi);
MAKE_BUILTIN_FUNCTION1(sqrt);
MAKE_BUILTIN_FUNCTION1(tan);
MAKE_BUILTIN_FUNCTION1(tanh);
MAKE_BUILTIN_FUNCTION1(tanpi);
MAKE_BUILTIN_FUNCTION1(tgamma);
MAKE_BUILTIN_FUNCTION1(trunc);



}

}


#endif // DUMMY_TYPES_HPP
