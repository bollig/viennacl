#ifndef VIENNACL_GENERATOR_FORWARDS_H
#define VIENNACL_GENERATOR_FORWARDS_H

namespace viennacl{

namespace generator{

class custom_operation;
class infos_base;
class assign_type;
class add_type;
class inplace_add_type;
class sub_type;
class inplace_sub_type;
class scal_mul_type;
class inplace_scal_mul_type;
class scal_div_type;
class inplace_scal_div_type;
class elementwise_prod_type;
class elementwise_div_type;

template<class LHS, class OP, class RHS>
class compile_time_beast;

template<typename ScalarType>
class dummy_vector;

template<typename ScalarType>
class dummy_scalar;

template<typename ScalarType, class Layout>
class dummy_matrix;

}

}
#endif // FORWARDS_H
