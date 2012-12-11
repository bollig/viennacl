#ifndef VIENNACL_GENERATOR_DUMMY_TYPES_HPP
#define VIENNACL_GENERATOR_DUMMY_TYPES_HPP

#include "viennacl/meta/enable_if.hpp"
#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/vector.hpp"

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

class inline_function_wrapper{
public:
    inline_function_wrapper(std::string const & name
                    ,std::string const & expr) : name_(name), expr_(expr){ }
private:
    std::string name_;
    std::string expr_;
};

template<typename ScalarType, unsigned int Alignment=1>
class dummy_vector{
    typedef dummy_vector self_type;
    typedef viennacl::vector<ScalarType,Alignment> vcl_vec_t;
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

class dummy_scalar{
    typedef dummy_scalar self_type;
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

class dummy_matrix{
    typedef dummy_matrix self_type;
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
template<>
struct is_scalar_expression_t<dummy_scalar>{ enum { value = 1}; };
template<class LHS, class OP, class RHS>
struct is_scalar_expression_t<scalar_expression_wrapper<LHS,OP,RHS> >{ enum { value = 1}; };

template<class T>
struct is_matrix_expression_t{ enum { value = 0 }; };
template<>
struct is_matrix_expression_t<dummy_matrix>{ enum { value = 1}; };
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


//template<class T>
//unary_minus<T> operator -(T const &)
//{
//  return unary_minus<T>();
//}

template<class LHS, class RHS>
typename viennacl::enable_if< (is_vector_expression_t<LHS>::value && is_vector_expression_t<RHS>::value)
                             ||(is_scalar_expression_t<LHS>::value && is_scalar_expression_t<RHS>::value)
                             ||(is_matrix_expression_t<LHS>::value && is_matrix_expression_t<RHS>::value)
                            ,typename convert_to_expr<LHS,add_type,RHS
                                                    ,is_vector_expression_t<LHS>::value
                                                    ,is_scalar_expression_t<LHS>::value
                                                    ,is_matrix_expression_t<LHS>::value>::type>::type
operator+(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,add_type,RHS, is_vector_expression_t<LHS>::value, is_scalar_expression_t<LHS>::value, is_matrix_expression_t<LHS>::value>::type(lhs,rhs);
}

template<class LHS, class RHS>
typename viennacl::enable_if< (is_vector_expression_t<LHS>::value && is_vector_expression_t<RHS>::value)
                             ||(is_scalar_expression_t<LHS>::value && is_scalar_expression_t<RHS>::value)
                             ||(is_matrix_expression_t<LHS>::value && is_matrix_expression_t<RHS>::value)
                            ,typename convert_to_expr<LHS,sub_type,RHS
                                                    ,is_vector_expression_t<LHS>::value
                                                    ,is_scalar_expression_t<LHS>::value
                                                    ,is_matrix_expression_t<LHS>::value>::type>::type
operator-(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,sub_type,RHS, is_vector_expression_t<LHS>::value, is_scalar_expression_t<LHS>::value, is_matrix_expression_t<LHS>::value>::type(lhs,rhs);
}

template<class LHS, class RHS>
typename viennacl::enable_if< is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
                            ,typename convert_to_expr<LHS,scal_mul_type,RHS
                                                    ,is_vector_expression_t<LHS>::value || is_vector_expression_t<RHS>::value
                                                    ,is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
                                                    ,is_matrix_expression_t<LHS>::value || is_matrix_expression_t<RHS>::value>::type>::type
operator*(LHS const & lhs, RHS const & rhs){
    return typename convert_to_expr<LHS,scal_mul_type,RHS
            ,is_vector_expression_t<LHS>::value || is_vector_expression_t<RHS>::value
            ,is_scalar_expression_t<LHS>::value || is_scalar_expression_t<RHS>::value
            ,is_matrix_expression_t<LHS>::value || is_matrix_expression_t<RHS>::value>::type(lhs,rhs);
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



#define MAKE_BUILTIN_FUNCTION1(name) static inline_function_wrapper name = inline_function_wrapper(#name,#name "(_1_)")
#define MAKE_BUILTIN_FUNCTION2(name) static inline_function_wrapper name = inline_function_wrapper(#name,#name "(_1_,_2_)")
#define MAKE_BUILTIN_FUNCTION3(name) static inline_function_wrapper name = inline_function_wrapper(#name,#name "(_1_,_2_,_3_)")

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
