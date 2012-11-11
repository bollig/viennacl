// Evan Bollig 2012
//
// Performs GMRES in parallel using an assumed additive schwarz decomposition.
// The comm_unit provides details of the MPI environment (rank, size, communicator, etc.)
// The domain provides the details of additive schwarz decomposition (implicit restriction operator)
//      - Implicit because we work directly on solution nodes, not on the differentiation matrix.
//      - if R is a restriction operator (eye only where stencils are part of
//      subdomain), then we have
//          \sum_{p=1}^{nproc}(R_p' A R_p)u = \sum_{p=1}^{nproc}(R_p'R_p)F = Au = F.
//      - For now we assume that R can be constructed with domain->Q; in the
//      future we might generalize this so the code is not specific to my
//      decomposition class

#ifndef VIENNACL_PARALLEL_GMRES_HPP_
#define VIENNACL_PARALLEL_GMRES_HPP_

// License related to ViennaCL content:
/* =========================================================================
   Copyright (c) 2010-2011, Institute for Microelectronics,
   Institute for Analysis and Scientific Computing,
   TU Wien.

   -----------------
   ViennaCL - The Vienna Computing Library
   -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file parallel_gmres.hpp
  @brief Implementations of the generalized minimum residual method are in this file.
  */

#include <vector>
#include <cmath>
#include <limits>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/traits/clear.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/meta/result_of.hpp"

// We'll extend the original gmres class:
#include "viennacl/linalg/gmres.hpp"

// But we need our own parallel norms
#include "linalg/parallel_norm_1.hpp"
#include "linalg/parallel_norm_2.hpp"
#include "linalg/parallel_norm_inf.hpp"
#include "linalg/parallel_inner_prod.hpp"

#include "utils/comm/communicator.h"
#include "grids/domain.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/filesystem.hpp>

#include "timer_eb.h"

namespace viennacl
{
    namespace linalg
    {

        /** @brief A tag for the solver GMRES. Used for supplying solver parameters and for dispatching the solve() function
        */
        class parallel_gmres_tag : public gmres_tag      //generalized minimum residual
        {
            public:
                /** @brief The constructor
                 *
                 * @param tol            Relative tolerance for the residual (solver quits if ||r|| < tol * ||r_initial||)
                 * @param max_iterations The maximum number of iterations (including restarts
                 * @param R     The maximum dimension of the Krylov space before restart (number of restarts is found by max_iterations / R)
                 */
                parallel_gmres_tag(Communicator& comm_unit, Domain& decomposition, double tol = 1e-10, unsigned int max_iterations = 300, unsigned int R = 20, unsigned int solution_dim_per_node = 1)
                    : gmres_tag(tol, max_iterations, R),
                    comm_ref(comm_unit), grid_ref(decomposition),
                    sol_dim(solution_dim_per_node)
            {
                std::cout << "GMRES for " << sol_dim << " components\n";
                tlist["allocate"] = new EB::Timer("Allocate Comm Buffers");
                tlist["gmres_alltoall"] = new EB::Timer("MPI_Alltoallv");
                tlist["setR"] = new EB::Timer("Sync set R");
                tlist["setO"] = new EB::Timer("Sync set O");
                allocateCommBuffers();
            };

                ~parallel_gmres_tag() {
                    delete [] sbuf;
                    delete [] sendcounts;
                    delete [] sdispls;
                    delete [] rdispls;
                    delete [] recvcounts;
                    delete [] rbuf;

                    tlist.printAllNonStatic();
                    tlist.clear();
                }

                Communicator const& comm() const { return this->comm_ref; }

                void allocateCommBuffers() {
                    tlist["allocate"]->start();
#if 0
                    double* sbuf;
                    int* sendcounts;
                    int* sdispls;
                    int* rdispls;
                    int* recvcounts;
                    double* rbuf;
#endif
                    this->sdispls = new int[grid_ref.O_by_rank.size()];
                    this->sendcounts = new int[grid_ref.O_by_rank.size()];
                    sdispls[0] = 0;
                    sendcounts[0] = sol_dim*grid_ref.O_by_rank[0].size();
                    unsigned int O_tot = sendcounts[0];
                    for (size_t i = 1; i < grid_ref.O_by_rank.size(); i++) {
                        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
                        sendcounts[i] = sol_dim*grid_ref.O_by_rank[i].size();
                        O_tot += sendcounts[i];
                    }

                    this->rdispls = new int[grid_ref.R_by_rank.size()];
                    this->recvcounts = new int[grid_ref.R_by_rank.size()];
                    rdispls[0] = 0;
                    recvcounts[0] = sol_dim*grid_ref.R_by_rank[0].size();
                    unsigned int R_tot = recvcounts[0];
                    for (size_t i = 1; i < grid_ref.R_by_rank.size(); i++) {
                        recvcounts[i] = sol_dim*grid_ref.R_by_rank[i].size();
                        rdispls[i] = rdispls[i-1] + recvcounts[i-1];
                        R_tot += recvcounts[i];
                    }

                    std::cout << "O_tot = " << O_tot << std::endl;
                    std::cout << "R_tot = " << R_tot << std::endl;

                    // Not sure if we need to use malloc to guarantee contiguous?
                    this->sbuf = new double[O_tot];
                    this->rbuf = new double[R_tot];
                    tlist["allocate"]->stop();
                }

                // TODO: (bug) fix the copy set R and set O so when a value is
                // sent to multiple CPUs it is placed in the correct vector for
                // alltoallv.
                template <typename GPUVectorType>
                    void syncSetR(GPUVectorType& gpu_vec) const {
                        tlist["setR"]->start();
                        //unsigned int nb_nodes = grid_ref.getNodeListSize();
                        //unsigned int set_G_size = grid_ref.G.size();
                        unsigned int set_Q_size = grid_ref.Q.size();
                        //unsigned int set_O_size = grid_ref.O.size();
                        unsigned int set_R_size = grid_ref.R.size();
                        unsigned int nb_bnd = grid_ref.getBoundaryIndicesSize();

                        // OUR SOLUTION IS ARRANGED IN THIS FASHION:
                        //  { Q\B D O R } where B = union(D, O) and Q = union(Q\B D O)

                        // TODO: fix this. We have to maintain an additional index
                        // map to convert from local node indices to the linear
                        // system indices (i.e. filter off the dirichlet boundary
                        // node indices
                        unsigned int offset_to_interior = nb_bnd;
                        unsigned int offset_to_set_R = set_Q_size;

                        viennacl::vector_range< GPUVectorType > setR(gpu_vec, viennacl::range((offset_to_set_R-offset_to_interior) * sol_dim, ((offset_to_set_R-offset_to_interior)+set_R_size) * sol_dim));

                        double* vec = new double[set_R_size * sol_dim];

                        unsigned int k = 0;
                        for (size_t i = 0; i < grid_ref.R_by_rank.size(); i++) {
                            k = this->rdispls[i];
                            for (size_t j = 0; j < grid_ref.R_by_rank[i].size(); j++) {
                                unsigned int r_indx = grid_ref.g2l(grid_ref.R_by_rank[i][j]);
                                // Offset the r_indx to 0 and we'll copy the subset
                                // into the proper range
                                r_indx -= offset_to_set_R;

                                r_indx *= sol_dim;

                                // TODO: need to translate to local
                                // indexing properly. This hack assumes all
                                // boundary are dirichlet and appear first
                                // in the list
                                for (unsigned int d = 0; d < sol_dim; d++) {
                                    vec[r_indx+d] = this->rbuf[k];
                                k++;
                                }
                            }
                        }

                        viennacl::copy(vec, setR, set_R_size * sol_dim);
                        delete [] vec;
                        tlist["setR"]->stop();
                    }

                template <typename GPUVectorType>
                    void syncSetO(GPUVectorType& gpu_vec) const {
                        tlist["setO"]->start();
                        unsigned int set_Q_size = grid_ref.Q.size();
                        unsigned int set_O_size = grid_ref.O.size();
                        unsigned int nb_bnd = grid_ref.getBoundaryIndicesSize();

                        //std::cout << "set_Q_size = " << set_Q_size << ", set_O_size = " << set_O_size << ", nb_bnd = " << nb_bnd << std::endl;

                        // OUR SOLUTION IS ARRANGED IN THIS FASHION:
                        //  { Q\B D O R } where B = union(D, O) and Q = union(Q\B D O)
                        //  Minus 1 because we start indexing at 0

                        // TODO: fix this. We have to maintain an additional index
                        // map to convert from local node indices to the linear
                        // system indices (i.e. filter off the dirichlet boundary
                        // node indices
                        unsigned int offset_to_interior = nb_bnd;
                        unsigned int offset_to_set_O = (set_Q_size - set_O_size);

                        //                    std::cout << "set_Q_size = " << set_Q_size << ", set_O_size = " << set_O_size << ", nb_bnd = " << nb_bnd << std::endl;

                        viennacl::vector_range< GPUVectorType > setO(gpu_vec, viennacl::range((offset_to_set_O - offset_to_interior) * sol_dim, ((offset_to_set_O-offset_to_interior)+set_O_size) * sol_dim));

                        double* vec = new double[set_O_size * sol_dim];

                        viennacl::copy(setO, vec, set_O_size * sol_dim);

                        // Copy elements to sbuf individually in case multiple procs receive the values
                        unsigned int k = 0;
                        for (size_t i = 0; i < grid_ref.O_by_rank.size(); i++) {
                            k = this->sdispls[i];
                            for (size_t j = 0; j < grid_ref.O_by_rank[i].size(); j++) {
                                unsigned int s_indx = grid_ref.g2l(grid_ref.O_by_rank[i][j]);

                                s_indx -= offset_to_set_O;

                                //std::cout << "offset_to_set_O = " << offset_to_set_O << ", s_indx = " << s_indx << std::endl;
                                s_indx *= sol_dim;

                                for (unsigned int d = 0; d < sol_dim; d++) {
                                    //std::cout << "k = " << k << ", s_ind+d = " << s_indx+d << std::endl;
                                    this->sbuf[k] = vec[s_indx+d];
                                k++;
                                }
                            }
                        }

                        delete [] vec;
                        tlist["setO"]->stop();
                    }

                // Generic for GPU (transfer GPU subset to cpu buffer, alltoallv
                // and return to the gpu)
                template <typename VectorType>
                    void
                    alltoall_subset(VectorType& vec, typename VectorType::value_type& dummy) const
                    {

                        // Share data in vector with all other processors. Will only transfer
                        // data associated with nodes in overlap between domains.
                        // Uses MPI_Alltoallv and MPI_Barrier.
                        // Copies data from vec to transfer, then writes received data into vec
                        // before returning.
                        if (comm_ref.getSize() > 1) {
                            tlist["gmres_alltoall"]->start();

                            // GPU array is arranged as {QmB BmO O R}
                            // TODO: however the set O and set R might have values that span multiple CPUs.

                            //std::cout << "VCL TRANSFER " << grid_ref.O_by_rank.size() << "\n";
                            syncSetO(vec);

                            MPI_Alltoallv(this->sbuf, this->sendcounts, this->sdispls, MPI_DOUBLE, this->rbuf, this->recvcounts, this->rdispls, MPI_DOUBLE, comm_ref.getComm());
                            comm_ref.barrier();

                            syncSetR(vec);
                            tlist["gmres_alltoall"]->stop();
                        }
                    }


                // For CPU (copy subset to cpu buffer, alltoallv
                // and return to the gpu)
                template <typename VectorType>
                    void
                    alltoall_subset(VectorType& vec, double& dummy) const
                    {

#if 1
                        // Share data in vector with all other processors. Will only transfer
                        // data associated with nodes in overlap between domains.
                        // Uses MPI_Alltoallv and MPI_Barrier.
                        // Copies data from vec to transfer, then writes received data into vec
                        // before returning.
                        if (comm_ref.getSize() > 1) {

                            //std::cout << "vec size = " << vec.size() << std::endl;
                            tlist["gmres_alltoall"]->start();
                            // Copy elements of set to sbuf
                            unsigned int k = 0;
                            for (size_t i = 0; i < grid_ref.O_by_rank.size(); i++) {
                                k = this->sdispls[i];
                                for (size_t j = 0; j < grid_ref.O_by_rank[i].size(); j++) {
                                    unsigned int s_indx = grid_ref.g2l(grid_ref.O_by_rank[i][j]) - grid_ref.getBoundaryIndicesSize();
                                    s_indx *= sol_dim;
                                    //                                    std::cout << "Sending " << s_indx << "\n";
                                    //std::cout << "s_indx = " << s_indx << ", k = " << k << std::endl;
                                    for (unsigned int d=0; d < sol_dim; d++) {
                                        this->sbuf[k] = vec[s_indx+d];
                                        k++;
                                    }
                                }
                            }

                            MPI_Alltoallv(this->sbuf, this->sendcounts, this->sdispls, MPI_DOUBLE, this->rbuf, this->recvcounts, this->rdispls, MPI_DOUBLE, comm_ref.getComm());

                            comm_ref.barrier();
                            k = 0;
                            for (size_t i = 0; i < grid_ref.R_by_rank.size(); i++) {
                                k = this->rdispls[i];
                                for (size_t j = 0; j < grid_ref.R_by_rank[i].size(); j++) {
                                    unsigned int r_indx = grid_ref.g2l(grid_ref.R_by_rank[i][j]) - grid_ref.getBoundaryIndicesSize();
                                    r_indx *= sol_dim;
                                    //std::cout << "r_indx = " << r_indx << ", k = " << k << std::endl;
                                    //                                    std::cout << "Receiving " << r_indx << "\n";
                                    // TODO: need to translate to local
                                    // indexing properly. This hack assumes all
                                    // boundary are dirichlet and appear first
                                    // in the list
                                    for (unsigned int d=0; d < sol_dim; d++) {
                                        vec[r_indx+d] = this->rbuf[k];
                                        k++;
                                    }
                                }
                            }
                            tlist["gmres_alltoall"]->stop();
                        }
#endif
                    }

                template <typename VectorType>
                    void
                    alltoall_subset(VectorType& v) const
                    {
                        typename VectorType::value_type dummy;
                        alltoall_subset(v, dummy);
                    }




            protected:
                Communicator& comm_ref;
                Domain& grid_ref;
                unsigned int sol_dim;
                mutable double* sbuf;
                mutable int* sendcounts;
                mutable int* sdispls;
                mutable int* rdispls;
                mutable int* recvcounts;
                mutable double* rbuf;
            protected:
                mutable EB::TimerList tlist;
        };


        namespace ublas = boost::numeric::ublas;
        namespace vcl = viennacl;

#if 1
        // ApplyPlaneRotation and Givens rotation based GMRES based on GMRES
        // implementation in CUSP-v0.3
        // Adapted to use ViennaCL primitives and my own communicators

// License from to CUSP implementation

/*
 *  Copyright 2011 The Regents of the University of California
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

        template <typename ValueType>
            void ApplyPlaneRotation(ValueType& dx,
                    ValueType& dy,
                    ValueType& cs,
                    ValueType& sn)
            {
                ValueType temp = cs * dx + sn *dy;
                dy = -sn*dx+cs*dy;
                dx = temp;
            }

        template <typename ValueType>
            void GeneratePlaneRotation(ValueType& dx,
                    ValueType& dy,
                    ValueType& cs,
                    ValueType& sn)
            {
                if(dy == ValueType(0.0)){
                    cs = 1.0;
                    sn = 0.0;
                }else if (abs(dy) > abs(dx)) {
                    ValueType tmp = dx / dy;
                    sn = ValueType(1.0) / sqrt(ValueType(1.0) + tmp*tmp);
                    cs = tmp*sn;
                }else {
                    ValueType tmp = dy / dx;
                    cs = ValueType(1.0) / sqrt(ValueType(1.0) + tmp*tmp);
                    sn = tmp*cs;
                }
            }


        template <typename MatrixType, typename VectorType>
            void PlaneRotation(MatrixType& H,
                    VectorType& cs,
                    VectorType& sn,
                    VectorType& s,
                    int i)
            {
                // Apply previous rotations
                for (int k = 0; k < i; k++){
                    ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
                }
                // Generate new rotation
                GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
                ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
                ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
            }
#endif

        /** @brief Implementation of the GMRES solver.
         *
         * Following the algorithm 2.1 proposed by Walker in "A Simpler GMRES" (1988)
         * (Evan Bollig): I changed variable names to be consistent with the paper and other literature
         *
         * @param matrix     The system matrix
         * @param rhs        The load vector
         * @param tag        Solver configuration tag
         * @param precond    A preconditioner. Precondition operation is done via member function apply()
         * @return The result vector
         */
        template <typename MatrixType, typename VectorType, typename PreconditionerType>
                          VectorType
                          solve(const MatrixType & A, VectorType & b_full, parallel_gmres_tag const & tag, PreconditionerType const & precond)
        {
#if 1
#define GMRES_DEBUG 1
#endif
            std::cout << "INSIDE VCL PARALLEL\n";
            EB::TimerList tlist;
            tlist["inner"] = new EB::Timer("GMRES Inner Iteration");
            tlist["outer"] = new EB::Timer("GMRES Outer Iteration");

            //typedef vcl::vector<double>                                             VectorType;
            typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
            typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

            unsigned int NN = A.size1(); //viennacl::traits::size1(matrix);
            unsigned int MM = A.size2(); //viennacl::traits::size2(matrix);
            int R  = (int)tag.krylov_dim();

            // Solution
            VectorType x_full(MM);
            vcl::vector_range<VectorType> x(x_full, vcl::range(0,NN));
            viennacl::traits::clear(x_full);
            tag.alltoall_subset(x_full);

            vcl::vector_range<VectorType> b(b_full, vcl::range(0,NN));

            // Workspace
            VectorType w_full(MM);
            vcl::vector_range<VectorType> w(w_full, vcl::range(0,NN));

            // Arnoldi Matrix
            std::vector< VectorType > v_full(R+1);
            std::vector< vcl::vector_range<VectorType> * > v(R+1);

            VectorType v0_full(MM);
            vcl::vector_range< VectorType > v0(v0_full, vcl::range(0,NN));

            // Givens rotations (NOTE: the R+1 x R least squares problem is solved on the CPU)
            ublas::vector<double> s(R+1);

            // Hessenberg matrix (if we do the givens rotations properly this ends as upper triangular)
            //std::vector< std::vector<CPU_ScalarType> > H(R+1);
            ublas::matrix<double> H(R+1,R, 0);

            // Rotations (cs = cosine; sn = sine)
            ublas::vector<double> cs(R);
            ublas::vector<double> sn(R);

            //representing the scalar '-1' on the GPU. Prevents blocking write operations
            const CPU_ScalarType gpu_scalar_minus_1 = static_cast<CPU_ScalarType>(-1);
            //representing the scalar '1' on the GPU. Prevents blocking write operations
            const CPU_ScalarType gpu_scalar_1 = static_cast<CPU_ScalarType>(1);
#if 0
            //representing the scalar '2' on the GPU. Prevents blocking write operations
            const CPU_ScalarType gpu_scalar_2 = static_cast<CPU_ScalarType>(2);
#endif

            double beta = 0;
            double rel_resid0 = 0;

            for (int k = 0; k < R+1; ++k)
            {
                //H[k].resize(tag.krylov_dim());
                v_full[k].resize(MM);
                v[k] = new vcl::vector_range<VectorType>(v_full[k], vcl::range(0,NN));
            }

            MPI_Barrier(MPI_COMM_WORLD);

#if GMRES_DEBUG
            if (tag.comm().isMaster())
                std::cout << "Starting Parallel GMRES..." << std::endl;
#endif
            tag.iters(0);

            // Save very first residual norm so we know when to stop
            CPU_ScalarType b_norm = viennacl::linalg::norm_2(b, tag.comm());
            v0 = b_full;
            precond.apply(v0_full);
            CPU_ScalarType resid0 = viennacl::linalg::norm_2(v0, tag.comm()) / b_norm;
            //std::cout << "B_norm = " << b_norm << ", Resid0 = " << resid0 << std::endl;

            do{
                tlist["outer"]->start();
                // compute initial residual and its norm //
                w = viennacl::linalg::prod(A, x_full) - b;                  // V(0) = A*x        //
                tag.alltoall_subset(w_full);
                precond.apply(w_full);                                  // V(0) = M*V(0)     //
                beta = viennacl::linalg::norm_2(w, tag.comm());
                w /= gpu_scalar_minus_1 * beta;                                         // V(0) = -V(0)/beta //

                *(v[0]) = w;

                // First givens rotation
                for (int ii = 0; ii < R+1; ii++) {
                    s[ii] = 0.;
                }
                s[0] = beta;
                int i = -1;

                if (beta / b_norm < tag.tolerance() || (b_norm == CPU_ScalarType(0.0)) )
                {
#if GMRES_DEBUG
                    if (tag.comm().isMaster())
                        std::cout << "Allowed Error reached at begin of loop" << std::endl;
#endif
                    tag.error(beta / b_norm);
                    tlist["outer"]->stop();
                    tlist.printAllNonStatic();
                    tlist.clear();
                    return x_full;
                }

                do{
                    tlist["inner"]->start();
                    ++i;
                    tag.iters(tag.iters() + 1);

                    tag.alltoall_subset(w_full);

                    //apply preconditioner
                    //can't pass in ref to column in V so need to use copy (w)
                    v0 = viennacl::linalg::prod(A,w_full);
                    tag.alltoall_subset(v0_full);
                    //V(i+1) = A*w = M*A*V(i)    //
                    precond.apply(v0_full);
                    w = v0;

                    for (int k = 0; k <= i; k++){
                        //  H(k,i) = <V(i+1),V(k)>    //
                        H(k, i) = viennacl::linalg::inner_prod(w, *(v[k]), tag.comm());
                        // V(i+1) -= H(k, i) * V(k)  //
                        //std::cout << v[k]->size() << "," << w.size() << H(k,i) * *(v[k])  << std::endl;
                        w -= H(k,i) * *(v[k]);
                    }

                    H(i+1,i) = viennacl::linalg::norm_2(w, tag.comm());


#ifdef VIENNACL_GMRES_DEBUG
                    std::cout << "H[" << i << "] = ";
                    for (int j = 0; j < R+1; j++) {
                        for (int k = 0; k < R; k++) {
                            std::cout << H(j,k) << ", ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                    //exit(-1);
#endif

                    // V(i+1) = V(i+1) / H(i+1, i) //
                    w /= gpu_scalar_1*H(i+1,i);
                    v_full[i+1] = w_full;

                    // Rotation takes place on host
                    PlaneRotation(H,cs,sn,s,i);

                    rel_resid0 = fabs(s[i+1]) / resid0;
#if GMRES_DEBUG
                    std::cout << "\t" << tag.iters() << "\t" << rel_resid0 << std::endl;
#endif
                    tag.error(rel_resid0);
                    tlist["inner"]->stop();
                    // We could add absolute tolerance here as well:
                    if (rel_resid0 < b_norm * tag.tolerance() ) {
                        break;
                    }

                }while (i+1 < R && tag.iters()+1 <= tag.max_iterations());
                // -------------------- SOLVE PROCESS ----------------------------------


                // After the Givens rotations, we have H is an upper triangular matrix
                for (int j = i; j >= 0; j--) {
                    s[j] /= H(j,j);
                    for (int k = j-1; k >= 0; k--) {
                        s[k] -= H(k,j) * s[j];
                    }
                }

#ifdef VIENNACL_GMRES_DEBUG
                std::cout << "H[" << i << "] = ";
                for (int j = 0; j < R+1; j++) {
                    for (int k = 0; k < R; k++) {
                        std::cout << H(j,k) << ", ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;

                std::cout << "s[" << i << "] = ";
                for (int j =0 ; j < R+1; j++) {
                    std::cout << s[j] << ", ";
                }
                std::cout << std::endl;
#endif

                // Update our solution
                for (int j = 0 ; j < i; j++) {
                    x +=  s[j] * *(v[j]);
                }
                tag.alltoall_subset(x_full);
                tlist["outer"]->stop();

#ifdef VIENNACL_GMRES_DEBUG
                if (tag.iters()+1 > 2*R) {
                    std::cout << "EXITING GMRES\n";
                    exit(-1);
                }
#endif

            } while (rel_resid0 >= b_norm*tag.tolerance() && tag.iters()+1 <= tag.max_iterations());

            tlist.printAllNonStatic();
            tlist.clear();
            return x_full;
        }


        /** @brief Implementation of the GMRES solver.
         *
         * Following the algorithm 2.1 proposed by Walker in "A Simpler GMRES" (1988)
         * (Evan Bollig): I changed variable names to be consistent with the paper and other literature
         *
         * @param matrix     The system matrix
         * @param rhs        The load vector
         * @param tag        Solver configuration tag
         * @param precond    A preconditioner. Precondition operation is done via member function apply()
         * @return The result vector
         */
        template <typename MatrixType, typename PreconditionerType>
                          ublas::vector<double>
                          solve(const MatrixType & A, ublas::vector<double> & b_full, parallel_gmres_tag const & tag, PreconditionerType const & precond)
        {
            std::cout << "INSIDE UBLAS PARALLEL\n";

            EB::TimerList tlist;
            tlist["inner"] = new EB::Timer("GMRES Inner Iteration");
            tlist["outer"] = new EB::Timer("GMRES Outer Iteration");

            typedef ublas::vector<double>                                             VectorType;
            typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
            typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

            unsigned int NN = A.size1(); //viennacl::traits::size1(matrix);
            unsigned int MM = A.size2(); //viennacl::traits::size2(matrix);
            int R  = (int)tag.krylov_dim();

            // Solution
            VectorType x_full(MM);
            ublas::vector_range<VectorType> x(x_full, ublas::range(0,NN));
            viennacl::traits::clear(x_full);
            tag.alltoall_subset(x_full);

            ublas::vector_range<VectorType> b(b_full, ublas::range(0,NN));

            // Workspace
            VectorType w_full(MM);
            ublas::vector_range<VectorType> w(w_full, ublas::range(0,NN));

            // Arnoldi Matrix
            std::vector< VectorType > v_full(R+1);
            std::vector< ublas::vector_range<VectorType> * > v(R+1);

            VectorType v0_full(MM);
            ublas::vector_range< VectorType > v0(v0_full, ublas::range(0,NN));

            // Givens rotations
            VectorType s(R+1);

            // Hessenberg matrix (if we do the givens rotations properly this ends as upper triangular)
            //std::vector< std::vector<CPU_ScalarType> > H(R+1);
            ublas::matrix<double> H(R+1,R,0);

            // Rotations (cs = cosine; sn = sine)
            VectorType cs(R);
            VectorType sn(R);

#if 0
            //representing the scalar '-1' on the GPU. Prevents blocking write operations
            const CPU_ScalarType gpu_scalar_minus_1 = static_cast<CPU_ScalarType>(-1);
            //representing the scalar '1' on the GPU. Prevents blocking write operations
            const CPU_ScalarType gpu_scalar_1 = static_cast<CPU_ScalarType>(1);
            //representing the scalar '2' on the GPU. Prevents blocking write operations
            const CPU_ScalarType gpu_scalar_2 = static_cast<CPU_ScalarType>(2);
#endif

            double beta = 0;
            double rel_resid0 = 0;

            for (int k = 0; k < R+1; ++k)
            {
                //H[k].resize(tag.krylov_dim());
                v_full[k].resize(MM);
                v[k] = new ublas::vector_range<VectorType>(v_full[k], ublas::range(0,NN));
            }

            MPI_Barrier(MPI_COMM_WORLD);

#if GMRES_DEBUG
            if (tag.comm().isMaster())
                std::cout << "Starting Parallel GMRES..." << std::endl;
#endif
            tag.iters(0);

            // Save very first residual norm so we know when to stop
            double b_norm = viennacl::linalg::norm_2(b, tag.comm());
            v0 = b_full;
            precond.apply(v0);
            double resid0 = viennacl::linalg::norm_2(v0, tag.comm()) / b_norm;
            //std::cout << "B_norm = " << b_norm << ", Resid0 = " << resid0 << std::endl;


            do{
                tlist["outer"]->start();
                // compute initial residual and its norm //
                w = b - viennacl::linalg::prod(A, x_full);                  // V(0) = A*x        //
                tag.alltoall_subset(w_full);
                precond.apply(w_full);                                  // V(0) = M*V(0)     //
                beta = viennacl::linalg::norm_2(w, tag.comm());
                w /= beta;                                         // V(0) = -V(0)/beta //

                *(v[0]) = w;

                // First givens rotation
                for (int ii = 0; ii < R+1; ii++) {
                    s[ii] = 0.;
                }
                s[0] = beta;
                int i = -1;

                if (beta / b_norm < tag.tolerance() || (b_norm == CPU_ScalarType(0.0)) )
                {
#if GMRES_DEBUG
                    if (tag.comm().isMaster())
                        std::cout << "Allowed Error reached at begin of loop" << std::endl;
#endif
                    tag.error(beta / b_norm);
                    tlist["outer"]->stop();
                    tlist.printAllNonStatic();
                    tlist.clear();
                    return x_full;
                }

                do{
                    tlist["inner"]->start();
                    ++i;
                    tag.iters(tag.iters() + 1);

                    tag.alltoall_subset(w_full);

                    //apply preconditioner
                    //can't pass in ref to column in V so need to use copy (w)
                    v0 = viennacl::linalg::prod(A,w_full);
                    tag.alltoall_subset(v0_full);
                    //V(i+1) = A*w = M*A*V(i)    //
                    precond.apply(v0_full);
                    w = v0;

                    for (int k = 0; k <= i; k++){
                        //  H(k,i) = <V(i+1),V(k)>    //
                        H(k, i) = viennacl::linalg::inner_prod(w, *(v[k]), tag.comm());
                        // V(i+1) -= H(k, i) * V(k)  //
                        w -= H(k,i) * *(v[k]);
                    }

                    H(i+1,i) = viennacl::linalg::norm_2(w, tag.comm());

#ifdef VIENNACL_GMRES_DEBUG
                    std::cout << "H[" << i << "] = ";
                    for (int j = 0; j < R+1; j++) {
                        for (int k = 0; k < R; k++) {
                            std::cout << H(j,k) << ", ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                    //exit(-1);
#endif
                    // V(i+1) = V(i+1) / H(i+1, i) //
                    w *= 1.0/H(i+1,i);
                    v_full[i+1] = w;

                    PlaneRotation(H,cs,sn,s,i);

                    rel_resid0 = fabs(s[i+1]) / resid0;
#if GMRES_DEBUG
                    std::cout << "\t" << tag.iters() << "\t" << rel_resid0 << std::endl;
#endif
                    tag.error(rel_resid0);
                    tlist["inner"]->stop();
                    // We could add absolute tolerance here as well:
                    if (rel_resid0 < b_norm * tag.tolerance() ) {
                        break;
                    }

                }while (i+1 < R && tag.iters()+1 <= tag.max_iterations());

                // -------------------- SOLVE PROCESS ----------------------------------


                // After the Givens rotations, we have H is an upper triangular matrix
                for (int j = i; j >= 0; j--) {
                    s[j] /= H(j,j);
                    for (int k = j-1; k >= 0; k--) {
                        s[k] -= H(k,j) * s[j];
                    }
                }

#ifdef VIENNACL_GMRES_DEBUG
                std::cout << "H[" << i << "] = ";
                for (int j = 0; j < R+1; j++) {
                    for (int k = 0; k < R; k++) {
                        std::cout << H(j,k) << ", ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;

                std::cout << "s[" << i << "] = ";
                for (int j =0 ; j < R+1; j++) {
                    std::cout << s[j] << ", ";
                }
                std::cout << std::endl;
#endif

                // Update our solution
                for (int j = 0 ; j < i; j++) {
                    x += *(v[j]) * s[j];
                }
                tag.alltoall_subset(x_full);
                tlist["outer"]->stop();

#ifdef VIENNACL_GMRES_DEBUG
                if (tag.iters()+1 > 2*R) {
                    std::cout << "EXITING GMRES\n";
                    exit(-1);
                }
#endif
            } while (rel_resid0 >= b_norm*tag.tolerance() && tag.iters()+1 <= tag.max_iterations());

            tlist.printAllNonStatic();
            tlist.clear();
            return x_full;
        }

        /** @brief Convenience overload of the solve() function using GMRES. Per default, no preconditioner is used
        */
        template <typename MatrixType, typename VectorType>
            VectorType solve(const MatrixType & matrix, VectorType & rhs, parallel_gmres_tag const & tag)
            {
                std::cout << "CALLING SOLVER\n";
                return solve(matrix, rhs, tag, no_precond());
            }


    }
}

#endif
