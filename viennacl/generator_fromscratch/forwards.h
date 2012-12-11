#ifndef VIENNACL_GENERATOR_FORWARDS_H
#define VIENNACL_GENERATOR_FORWARDS_H

namespace viennacl{

namespace generator{

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
template<typename ScalarType, unsigned int Alignment>
class dummy_vector;
class dummy_scalar;
class dummy_matrix;

}

}
#endif // FORWARDS_H
