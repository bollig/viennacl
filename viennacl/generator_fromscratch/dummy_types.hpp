#ifndef VIENNACL_GENERATOR_DUMMY_TYPES_HPP
#define VIENNACL_GENERATOR_DUMMY_TYPES_HPP

namespace viennacl{

namespace generator{

class dummy_vector{

};

class dummy_scalar{

};

class dummy_matrix{

};

template<class LHS, class OP, class RHS>
class compile_time_beast{
    compile_time_beast(LHS const & lhs, RHS const & rhs) : lhs_(lhs), rhs_(rhs){ }
    LHS const & lhs(){return lhs_;}
    RHS const & rhs(){return rhs_;}
private:
    LHS const & lhs_;
    RHS const & rhs_;
};


}

}


#endif // DUMMY_TYPES_HPP
