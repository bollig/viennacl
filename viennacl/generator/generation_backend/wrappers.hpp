#ifndef VIENNACL_GENERATOR_BACKEND_WRAPPERS_HPP
#define VIENNACL_GENERATOR_BACKEND_HPP

#include <list>

namespace viennacl{

    namespace generator{

            static void replace_all_string(std::string & str, std::string const & from, std::string const & to){
                size_t start_pos = 0;
                while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                         str.replace(start_pos, from.length(), to);
                         start_pos += to.length();
                }
            }

            /**
             * @brief The mat_infos class
             */
            class mat_infos{
            public:
                mat_infos(std::string const & scalartype,
                          std::string const & name,
                          std::string const & size1,
                          std::string const & size2,
                          bool const & is_rowmajor,
                          bool const & is_transposed) : scalartype_(scalartype), name_(name)
                                                                      ,  size1_(size1), size2_(size2)
                                                                      , is_rowmajor_(is_rowmajor), is_transposed_(is_transposed){ }
                std::string const & size1() const{ return size1_; }
                std::string const & size2() const{ return size2_; }
                std::string const & name() const{ return name_; }
                std::string const & scalartype() const{ return scalartype_; }
                bool const is_rowmajor() const { return is_rowmajor_; }
                bool const is_transposed() const { return is_transposed_; }
                void set_access_name(std::string const & new_name) const{ access_name_ = new_name; }
                std::string access_name() const { return access_name_; }
            private:
                mutable std::string access_name_;
                std::string scalartype_;
                std::string name_;
                std::string size1_;
                std::string size2_;
                bool is_rowmajor_;
                bool is_transposed_;
            };

            /**
             * @brief The mat_expr_infos class
             */
            class mat_expr_infos{
                mat_expr_infos(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }

                void add_infos(std::string const & scalartype,
                               std::string const & name,
                               std::string const & size1,
                               std::string const & size2,
                               bool const & is_rowmajor,
                               bool const & is_transposed){
                    mats_.push_back(mat_infos(scalartype,name,size1,size2,is_rowmajor,is_transposed));
                }

                std::string generate(){
                    std::string res(expression_string_base_);
                    for(std::list<mat_infos>::iterator it = mats_.begin() ; it!= mats_.end() ; ++it){
                        replace_all_string(res,it->name(),it->access_name());
                    }
                }

            private:
                std::list<mat_infos> mats_;
                mutable std::string expression_string_base_;
            };

            /**
             * @brief The vec_infos class
             */
            class vec_infos{
            public:
                vec_infos(std::string const & scalartype,
                                  std::string const & name,
                                  std::string const & size) : scalartype_(scalartype), name_(name), size_(size){ }
                std::string const & scalartype(){ return scalartype_; }
                std::string const & name(){ return name_; }
                std::string const & size(){ return size_; }
                void set_access_name(std::string const & new_name) const{ access_name_ = new_name; }
                std::string access_name() const { return access_name_; }

            private:
                mutable std::string access_name_;
                std::string scalartype_;
                std::string name_;
                std::string size_;
            };

            /**
             * @brief The vec_expr_infos class
             */
            class vec_expr_infos{

                vec_expr_infos(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }

                void add_infos(std::string const & scalartype,
                               std::string const & name,
                               std::string const & size){
                    vecs_.push_back(vec_infos(scalartype,name,size));
                }

                std::string generate(){
                    std::string res(expression_string_base_);
                    for(std::list<vec_infos>::iterator it = vecs_.begin() ; it!= vecs_.end() ; ++it){
                        replace_all_string(res,it->name(),it->access_name());
                    }
                }

            private:
                std::list<vec_infos> vecs_;
                mutable std::string expression_string_base_;
            };

            /**
             * @brief The matvec_prod_infos struct
             */
            struct matvec_prod_infos{
                matvec_prod_infos(mat_infos const & lhs, vec_infos const & rhs, vec_infos const & additional_expression): lhs_(lhs)
                                                                                                                    , rhs_(rhs)
                                                                                                                    , additional_expression_(additional_expression){ }
                mat_expr_infos const & lhs() const { return lhs_; }
                vec_expr_infos const & rhs() const { return rhs_; }
                vec_expr_infos const & additional_expression() const { return additional_expression_; }

            private:
                mat_expr_infos lhs_;
                vec_expr_infos rhs_;
                vec_expr_infos additional_expression_;
            };



            template<class T>
            struct wrap_vec_expr{
                typedef typename tree_utils::extract_if<T,result_of::is_symbolic_vector,false>::Result VectorsList;

                template<class U>
                struct functor{
                    static void execute(vec_expr_infos & res){
                        res.add_infos(print_type<typename U::ScalarType,1>::value(),
                                      U::name(),
                                      U::internal_size2_name());
                    }
                };

                vec_expr_infos create_vector_expression_wrapper(){
                    vec_expr_infos res(make_expression_code<T>::value(""));

                }
            };



            template<class T>
            matvec_prod_infos create_matvec_prod_wrapper(){
                typedef typename tree_utils::remove_if<T,result_of::is_symbolic_vector,false >::Result Product;
                typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vector;
                typedef typename Product::LHS ProdLHS;
                typedef typename Product::RHS ProdRHS;
                typedef typename result_of::expression_type<T>::Result    IntermediateType;
                static const unsigned int Alignment = IntermediateType::Alignment;
                VIENNACL_STATIC_ASSERT(Alignment==1,AlignmentNotSupported);

                mat_expr_infos lhs_infos(make_expression_code<ProdLHS>::value(""));
                vec_expr_infos rhs_infos(make_expression_code<ProdRHS>::value(""));
                vec_expr_infos additional_expression(make_expression_code<Vector>::value(""));

                lhs_infos.add_infos(print_type<typename ProdLHS::ScalarType,1>::value(),
                                    ProdLHS::name(),
                                    ProdLHS::internal_size1_name(),
                                    ProdLHS::internal_size2_name(),
                                    result_of::is_row_major<ProdLHS>::value,
                                    false,
                                    );

                rhs_infos.add_infos(print_type<typename ProdRHS::ScalarType,1>::value(),
                                    ProdRHS::name(),
                                    ProdRHS::internal_size2_name(),
                                   );

                return matvec_prod_infos(lhs_infos,rhs_infos,make_expression_code<T>::value("__OFFSET__"));
            }
    }

}




#endif // WRAPPERS_HPP
