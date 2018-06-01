/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


// external
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

// acados
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_gnsf.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "interfaces/acados_c/sim_interface.h"

// crane dae model
#include "examples/c/crane_dae_model/crane_dae_model.h"


extern "C"
{

}

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

sim_solver_t hashitsim_dae(std::string const& inString)
{
    if (inString == "ERK") return ERK;
    if (inString == "IRK") return IRK;
    if (inString == "LIFTED_IRK") return LIFTED_IRK;
    if (inString == "GNSF") return GNSF;
    if (inString == "NEW_LIFTED_IRK") return NEW_LIFTED_IRK;

    return (sim_solver_t) -1;
}

double sim_solver_tolerance_dae(std::string const& inString)
{
    // if (inString == "ERK") return 1e-3;
    if (inString == "IRK") return 1e-4;
    // if (inString == "LIFTED_IRK") return 1e-3;
    if (inString == "GNSF") return 1e-4;
    // if (inString == "NEW_LIFTED_IRK") return 1e-3;

    return -1;
}




TEST_CASE("crane_dae_example", "[integrators]")
{
    vector<std::string> solvers = {"IRK", "GNSF"};
      // {"ERK", "IRK", "LIFTED_IRK", "GNSF", "NEW_LIFTED_IRK"};
    // initialize dimensions
    int ii, jj;

    const int nx = 9;
    const int nu = 2;
    const int nz = 2;
    const int nx1 = 8;  // gnsf split
    const int nx2 = 1;
    const int n_out = 3;
    const int ny = 5;
    const int nuhat = 1;

    // generate x0, u_sim
    double x0[nx];
    double u_sim[nu];

    for (int ii = 0; ii < nx; ii++) {
        x0[ii] = 0.0;
    }
    x0[2] = 0.8;

    u_sim[0] = 40.108149413030752;
    u_sim[1] = -50.446662212534974;

    int NF = nx + nu;  // columns of forward seed

    int nsim0 = 1;  // nsim;

    double T = 0.01;  // simulation time // former values .1, .05
    // reduced for faster test

    double x_sim[nx*(nsim0+2)];
    double x_ref_sol[nx];
    double S_forw_ref_sol[nx*NF];
    double S_adj_ref_sol[NF];
    double z_ref_sol[nz];
    double S_alg_ref_sol[NF*nz];

    double error[nx];
    double error_z[nz];
    double error_S_forw[nx*NF];
    double error_S_adj[NF];
    double error_S_alg[NF*nz];

    double max_error, max_error_forw, max_error_adj, max_error_z, max_error_sens_alg;

    for (ii=0; ii < nx; ii++)
        x_sim[ii] = x0[ii];

/************************************************
* external functions
************************************************/

    // impl_ode_fun
    external_function_casadi impl_ode_fun;
    impl_ode_fun.casadi_fun = &crane_dae_impl_ode_fun;
    impl_ode_fun.casadi_work = &crane_dae_impl_ode_fun_work;
    impl_ode_fun.casadi_sparsity_in = &crane_dae_impl_ode_fun_sparsity_in;
    impl_ode_fun.casadi_sparsity_out = &crane_dae_impl_ode_fun_sparsity_out;
    impl_ode_fun.casadi_n_in = &crane_dae_impl_ode_fun_n_in;
    impl_ode_fun.casadi_n_out = &crane_dae_impl_ode_fun_n_out;
    external_function_casadi_create(&impl_ode_fun);

    // impl_ode_fun_jac_x_xdot
    external_function_casadi impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_fun = &crane_dae_impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_work = &crane_dae_impl_ode_fun_jac_x_xdot_work;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &crane_dae_impl_ode_fun_jac_x_xdot_sparsity_in;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &crane_dae_impl_ode_fun_jac_x_xdot_sparsity_out;
    impl_ode_fun_jac_x_xdot.casadi_n_in = &crane_dae_impl_ode_fun_jac_x_xdot_n_in;
    impl_ode_fun_jac_x_xdot.casadi_n_out = &crane_dae_impl_ode_fun_jac_x_xdot_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_fun = &crane_dae_impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_work = &crane_dae_impl_ode_jac_x_xdot_u_work;
    impl_ode_jac_x_xdot_u.casadi_sparsity_in = &crane_dae_impl_ode_jac_x_xdot_u_sparsity_in;
    impl_ode_jac_x_xdot_u.casadi_sparsity_out = &crane_dae_impl_ode_jac_x_xdot_u_sparsity_out;
    impl_ode_jac_x_xdot_u.casadi_n_in = &crane_dae_impl_ode_jac_x_xdot_u_n_in;
    impl_ode_jac_x_xdot_u.casadi_n_out = &crane_dae_impl_ode_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_jac_x_xdot_u);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_fun = &crane_dae_impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_work = &crane_dae_impl_ode_fun_jac_x_xdot_u_work;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_in =
                            &crane_dae_impl_ode_fun_jac_x_xdot_u_sparsity_in;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_out =
                            &crane_dae_impl_ode_fun_jac_x_xdot_u_sparsity_out;
    impl_ode_fun_jac_x_xdot_u.casadi_n_in = &crane_dae_impl_ode_fun_jac_x_xdot_u_n_in;
    impl_ode_fun_jac_x_xdot_u.casadi_n_out = &crane_dae_impl_ode_fun_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot_u);

    /************************************************
    * external functions (Generalized Nonlinear Static Feedback (GNSF) model)
    ************************************************/
    // phi_fun
    external_function_casadi phi_fun;
    phi_fun.casadi_fun            = &crane_dae_phi_fun;
    phi_fun.casadi_work           = &crane_dae_phi_fun_work;
    phi_fun.casadi_sparsity_in    = &crane_dae_phi_fun_sparsity_in;
    phi_fun.casadi_sparsity_out   = &crane_dae_phi_fun_sparsity_out;
    phi_fun.casadi_n_in           = &crane_dae_phi_fun_n_in;
    phi_fun.casadi_n_out          = &crane_dae_phi_fun_n_out;
    external_function_casadi_create(&phi_fun);

    // phi_fun_jac_y
    external_function_casadi phi_fun_jac_y;
    phi_fun_jac_y.casadi_fun            = &crane_dae_phi_fun_jac_y;
    phi_fun_jac_y.casadi_work           = &crane_dae_phi_fun_jac_y_work;
    phi_fun_jac_y.casadi_sparsity_in    = &crane_dae_phi_fun_jac_y_sparsity_in;
    phi_fun_jac_y.casadi_sparsity_out   = &crane_dae_phi_fun_jac_y_sparsity_out;
    phi_fun_jac_y.casadi_n_in           = &crane_dae_phi_fun_jac_y_n_in;
    phi_fun_jac_y.casadi_n_out          = &crane_dae_phi_fun_jac_y_n_out;
    external_function_casadi_create(&phi_fun_jac_y);

    // phi_jac_y_uhat
    external_function_casadi phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_fun                = &crane_dae_phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_work               = &crane_dae_phi_jac_y_uhat_work;
    phi_jac_y_uhat.casadi_sparsity_in        = &crane_dae_phi_jac_y_uhat_sparsity_in;
    phi_jac_y_uhat.casadi_sparsity_out       = &crane_dae_phi_jac_y_uhat_sparsity_out;
    phi_jac_y_uhat.casadi_n_in               = &crane_dae_phi_jac_y_uhat_n_in;
    phi_jac_y_uhat.casadi_n_out              = &crane_dae_phi_jac_y_uhat_n_out;
    external_function_casadi_create(&phi_jac_y_uhat);

    // f_lo_fun_jac_x1k1uz
    external_function_casadi f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_fun            = &crane_dae_f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_work           = &crane_dae_f_lo_fun_jac_x1k1uz_work;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_in    = &crane_dae_f_lo_fun_jac_x1k1uz_sparsity_in;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_out   = &crane_dae_f_lo_fun_jac_x1k1uz_sparsity_out;
    f_lo_fun_jac_x1k1uz.casadi_n_in           = &crane_dae_f_lo_fun_jac_x1k1uz_n_in;
    f_lo_fun_jac_x1k1uz.casadi_n_out          = &crane_dae_f_lo_fun_jac_x1k1uz_n_out;
    external_function_casadi_create(&f_lo_fun_jac_x1k1uz);

    // get_matrices_fun
    external_function_casadi get_matrices_fun;
    get_matrices_fun.casadi_fun            = &crane_dae_get_matrices_fun;
    get_matrices_fun.casadi_work           = &crane_dae_get_matrices_fun_work;
    get_matrices_fun.casadi_sparsity_in    = &crane_dae_get_matrices_fun_sparsity_in;
    get_matrices_fun.casadi_sparsity_out   = &crane_dae_get_matrices_fun_sparsity_out;
    get_matrices_fun.casadi_n_in           = &crane_dae_get_matrices_fun_n_in;
    get_matrices_fun.casadi_n_out          = &crane_dae_get_matrices_fun_n_out;
    external_function_casadi_create(&get_matrices_fun);


/************************************************
* Create Reference Solution
************************************************/

    sim_solver_plan plan;
    plan.sim_solver = GNSF;  // IRK; -- works but super slow

    sim_solver_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);

    /* set dimensions */
    config->set_nx(dims, nx);
    config->set_nu(dims, nu);
    config->set_nz(dims, nz);

    // GNSF -- set additional dimensions
    sim_gnsf_dims *gnsf_dim;
    if (plan.sim_solver == GNSF)
    {
        gnsf_dim = (sim_gnsf_dims *) dims;
        gnsf_dim->nx1 = nx1;
        gnsf_dim->nx2 = nx2;
        gnsf_dim->ny = ny;
        gnsf_dim->nuhat = nuhat;
        gnsf_dim->n_out = n_out;
    }

    // set opts
    void *opts_ = sim_opts_create(config, dims);
    sim_rk_opts *opts = (sim_rk_opts *) opts_;
    config->opts_initialize_default(config, dims, opts);

    // opts reference solution
    opts->sens_forw = true;
    opts->sens_adj = true;
    opts->sens_algebraic = true;
    opts->output_z = true;
    opts->jac_reuse = false;  // jacobian reuse
    opts->newton_iter = 5;  // number of newton iterations per integration step
    opts->num_steps = 500;  // number of steps
    opts->ns = 8;  // number of stages in rk integrator

    sim_in *in = sim_in_create(config, dims);
    sim_out *out = sim_out_create(config, dims);

    in->T = T;

    // import model matrices
    external_function_generic *get_model_matrices =
            (external_function_generic *) &get_matrices_fun;
    gnsf_model *model = (gnsf_model *) in->model;

    // set model
    switch (plan.sim_solver)
    {
        case IRK:  // IRK
        {
            sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
            sim_set_model(config, in, "impl_ode_fun_jac_x_xdot",
                    &impl_ode_fun_jac_x_xdot);
            sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
            break;
        }
        case GNSF:  // GNSF
        {
            // set model funtions
            sim_set_model(config, in, "phi_fun", &phi_fun);
            sim_set_model(config, in, "phi_fun_jac_y", &phi_fun_jac_y);
            sim_set_model(config, in, "phi_jac_y_uhat", &phi_jac_y_uhat);
            sim_set_model(config, in, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);

            // import model matrices
            external_function_generic *get_model_matrices =
                    (external_function_generic *) &get_matrices_fun;
            gnsf_model *model = (gnsf_model *) in->model;
            sim_gnsf_import_matrices(gnsf_dim, model, get_model_matrices);
            break;
        }
        default :
        {
            printf("\nnot plan.sim_solver not supported!\n");
            exit(1);
        }
    }

    // seeds forw
    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    // seeds adj
    for (ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;
    for (ii = nx; ii < nx + nu; ii++)
        in->S_adj[ii] = 0.0;

    /************************************************
    * sim solver
    ************************************************/

    sim_solver *sim_solver = sim_create(config, dims, opts);


    if (plan.sim_solver == GNSF){  // for gnsf: perform precomputation
        gnsf_model *model = (gnsf_model *) in->model;
        sim_gnsf_precompute(config, gnsf_dim, model, opts,
                    sim_solver->mem, sim_solver->work, in->T);
    }

    int acados_return;

    for (ii=0; ii < nsim0; ii++)
    {
        // x
        for (jj = 0; jj < nx; jj++)
            in->x[jj] = x_sim[ii*nx+jj];

        // u
        for (jj = 0; jj < nu; jj++)
            in->u[jj] = u_sim[ii*nu+jj];

        acados_return = sim_solve(sim_solver, in, out);
        REQUIRE(acados_return == 0);

        for (jj = 0; jj < nx; jj++)
            x_sim[(ii+1)*nx+jj] = out->xn[jj];
    }

    for (jj = 0; jj < nx; jj++)
        x_ref_sol[jj] = out->xn[jj];

    for (jj = 0; jj < nx*NF; jj++)
        S_forw_ref_sol[jj] = out->S_forw[jj];

    for (jj = 0; jj < NF; jj++)
        S_adj_ref_sol[jj] = out->S_adj[jj];

    for (jj = 0; jj < nz; jj++)
        z_ref_sol[jj] = out->zn[jj];

    for (jj = 0; jj < nz*NF; jj++)
        S_alg_ref_sol[jj] = out->S_algebraic[jj];

    // printf("Reference xn \n");
    // d_print_e_mat(1, nx, &x_ref_sol[0], 1);

    // printf("Reference zn \n");
    // d_print_e_mat(1, nz, &z_ref_sol[0], 1);

    // printf("Reference forward sensitivities \n");
    // d_print_e_mat(nx, NF, &S_forw_ref_sol[0], nx);

    // printf("tested adjoint sensitivities \n");
    // d_print_e_mat(1, NF, &S_adj_ref_sol[0], 1);

    /* free */
    free(config);
    free(dims);
    free(opts);

    free(in);
    free(out);
    free(sim_solver);

/************************************************
* test solver loop
************************************************/



    for (int sens_forw = 1; sens_forw < 2; sens_forw++)
    {
    SECTION("sens_forw = " + std::to_string((bool)sens_forw))
    {
        for (int sens_adj = 1; sens_adj < 2; sens_adj++)
        {
        SECTION("sens_adj = " + std::to_string((bool)sens_adj))
        {
            for (int output_z = 1; output_z < 2; output_z++)
            {
            SECTION("output_z = " + std::to_string((bool)output_z))
            {
            for (int sens_alg = 1; sens_alg < 2; sens_alg++)
            {
            SECTION("sens_alg = " + std::to_string((bool)sens_alg))
            {
            for (int num_stages = 3; num_stages < 5; num_stages++)
            {
            SECTION("num_stages = " + std::to_string(num_stages))
            {
            for (int num_steps = 2; num_steps < 7; num_steps += 3)
            {
            SECTION("num_steps = " + std::to_string(num_steps))
            {

            for (std::string solver : solvers)
            {
            SECTION(solver)
            {


                double tol = sim_solver_tolerance_dae(solver);

                plan.sim_solver = hashitsim_dae(solver);

                // create correct config based on plan
                sim_solver_config *config = sim_config_create(plan);

            /* sim dims */
                void *dims = sim_dims_create(config);
                config->set_nx(dims, nx);
                config->set_nu(dims, nu);
                config->set_nz(dims, nz);
                // GNSF -- set additional dimensions
                sim_gnsf_dims *gnsf_dim;
                if (plan.sim_solver == GNSF)
                {
                    gnsf_dim = (sim_gnsf_dims *) dims;
                    gnsf_dim->nx1 = nx1;
                    gnsf_dim->nx2 = nx2;
                    gnsf_dim->ny = ny;
                    gnsf_dim->nuhat = nuhat;
                    gnsf_dim->n_out = n_out;
                }

            /* sim options */

                void *opts_ = sim_opts_create(config, dims);
                sim_rk_opts *opts = (sim_rk_opts *) opts_;
                config->opts_initialize_default(config, dims, opts);

                opts->jac_reuse = false;        // jacobian reuse
                opts->newton_iter = 2;          // number of newton iterations per integration step

                opts->ns                = num_stages;          // number of stages in rk integrator
                opts->num_steps         = num_steps;    // number of steps
                opts->sens_forw         = (bool) sens_forw;
                opts->sens_adj          = (bool) sens_adj;
                opts->output_z          = (bool) output_z;
                opts->sens_algebraic    = (bool) sens_alg;


            /* sim in / out */

                sim_in *in = sim_in_create(config, dims);
                sim_out *out = sim_out_create(config, dims);

                in->T = T;

                // external functions -- model
                switch (plan.sim_solver)
                {
                    case IRK:  // IRK
                    {
                        sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
                        sim_set_model(config, in, "impl_ode_fun_jac_x_xdot",
                                &impl_ode_fun_jac_x_xdot);
                        sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
                        break;
                    }
                    case GNSF:  // GNSF
                    {
                        // set model funtions
                        sim_set_model(config, in, "phi_fun", &phi_fun);
                        sim_set_model(config, in, "phi_fun_jac_y", &phi_fun_jac_y);
                        sim_set_model(config, in, "phi_jac_y_uhat", &phi_jac_y_uhat);
                        sim_set_model(config, in, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);

                        // import model matrices
                        external_function_generic *get_model_matrices =
                                (external_function_generic *) &get_matrices_fun;
                        gnsf_model *model = (gnsf_model *) in->model;
                        sim_gnsf_import_matrices(gnsf_dim, model, get_model_matrices);
                        break;
                    }
                    // case NEW_LIFTED_IRK:  // new_lifted_irk
                    // {
                    //     sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
                    //     sim_set_model(config, in, "impl_ode_fun_jac_x_xdot_u",
                    //              &impl_ode_fun_jac_x_xdot_u);
                    //     break;
                    // }
                    default :
                    {
                        printf("\nnot enough sim solvers implemented!\n");
                        exit(1);
                    }
                }

            /* seeds */
                for (ii = 0; ii < nx * NF; ii++)
                    in->S_forw[ii] = 0.0;
                for (ii = 0; ii < nx; ii++)
                    in->S_forw[ii * (nx + 1)] = 1.0;

                // seeds adj
                for (ii = 0; ii < nx; ii++)
                    in->S_adj[ii] = 1.0;
                for (ii = nx; ii < nx + nu; ii++)
                    in->S_adj[ii] = 0.0;

            /** sim solver  */
                sim_solver = sim_create(config, dims, opts);
                int acados_return;

                if (plan.sim_solver == GNSF){  // for gnsf: perform precomputation
                    gnsf_model *model = (gnsf_model *) in->model;
                    sim_gnsf_precompute(config, gnsf_dim, model, opts,
                             sim_solver->mem, sim_solver->work, in->T);
                }
        
            /* print */
            std::cout << "\n---> testing integrator " << solver;
            std::cout << " OPTS: num_steps = " << opts->num_steps;
            std::cout << ", num_stages = " << opts->ns;
            std::cout << ", jac_reuse = " << opts->jac_reuse;
            std::cout << ", newton_iter = " << opts->newton_iter << ")\n";

                for (ii=0; ii < nsim0; ii++)
                {
                    // x
                    for (jj = 0; jj < nx; jj++)
                        in->x[jj] = x_sim[ii*nx+jj];

                    // u
                    for (jj = 0; jj < nu; jj++)
                        in->u[jj] = u_sim[ii*nu+jj];

                    acados_return = sim_solve(sim_solver, in, out);
                    REQUIRE(acados_return == 0);

                    for (jj = 0; jj < nx; jj++){
                        x_sim[(ii+1)*nx+jj] = out->xn[jj];
                    }

                }

            /************************************************
            * compute error w.r.t. reference solution
            ************************************************/

                // error sim
                for (jj = 0; jj < nx; jj++){
                    error[jj] = fabs(out->xn[jj] - x_ref_sol[jj]);
                    REQUIRE(std::isnan(out->xn[jj]) == false);
                }
                // max_error
                max_error = 0.0;
                for (int ii = 0; ii < nx; ii++)
                    max_error = (error[ii] >= max_error) ? error[ii] : max_error;

                if ( opts->sens_forw ){     // error_S_forw
                    max_error_forw = 0.0;
                    for (jj = 0; jj < nx*NF; jj++){
                        REQUIRE(std::isnan(out->S_forw[jj]) == 0);
                        error_S_forw[jj] = fabs(S_forw_ref_sol[jj] - out->S_forw[jj]);
                        max_error_forw = (error_S_forw[jj] >= max_error_forw)
                                ? error_S_forw[jj] : max_error_forw;
                    }
                }

                if ( opts->sens_adj ){               // error_S_adj
                    for (jj = 0; jj < nx + nu; jj++){
                        REQUIRE(std::isnan(out->S_adj[jj]) == 0);
                        error_S_adj[jj] = S_adj_ref_sol[jj] - out->S_adj[jj];
                    }
                    max_error_adj = 0.0;
                    for (jj = 0; jj < nx + nu; jj++)
                        max_error_adj = (error_S_adj[jj] >= max_error_adj)
                                ? error_S_adj[jj] : max_error_adj;
                }

                if ( opts->output_z ){      // error_z
                    max_error_z = 0.0;
                    for (jj = 0; jj < nz; jj++){
                        error_z[jj] = fabs(out->zn[jj] - z_ref_sol[jj]);
                        REQUIRE(std::isnan(out->zn[jj]) == 0);
                        max_error_z = (error_z[ii] >= max_error_z) ? error_z[ii] : max_error_z;
                    }
                }

                if ( opts->sens_algebraic ){        // error_S_alg
                    max_error_sens_alg = 0.0;
                    for (jj = 0; jj < nz*NF; jj++){
                        REQUIRE(std::isnan(out->S_algebraic[jj]) == 0);
                        error_S_alg[jj] = fabs(out->S_algebraic[jj] - S_alg_ref_sol[jj]);
                    }
                    for (int ii = 0; ii < nz * NF; ii++)
                        max_error_sens_alg = (error_S_alg[ii] >= max_error_sens_alg) ?
                                            error_S_alg[ii] : max_error_sens_alg;
                }




            /************************************************
            * printing
            ************************************************/

                std::cout  << "error_sim        = " << max_error << "\n";
                if ( opts->sens_forw )
                std::cout  << "error_forw       = " << max_error_forw << "\n";
                if ( opts->sens_adj )
                std::cout  << "error_adj        = " << max_error_adj  << "\n";
                if ( opts->sens_algebraic )
                std::cout  << "error_algeb_sens = " << max_error_sens_alg  << "\n";
                if ( opts->output_z )
                std::cout  << "error_z          = " << max_error_z << "\n";

                // printf("tested algebraic sensitivities \n");
                // d_print_e_mat(nz, NF, &out->S_algebraic[0], nz);

            /************************************************
            * asserts on erors
            ************************************************/
                REQUIRE(max_error <= tol);

                if ( opts->sens_forw )
                    REQUIRE(max_error_forw <= tol);

                if ( opts->sens_adj )
                    REQUIRE(max_error_adj <= tol);

                if ( opts->output_z )
                    REQUIRE(max_error_z <= 1e1*tol);

                if ( opts->sens_algebraic )
                    REQUIRE(max_error_sens_alg <= 1e3*tol);

            /************************************************
            * free tested solver
            ************************************************/
                free(config);
                free(dims);
                free(opts);

                free(in);
                free(out);
                free(sim_solver);
            }  // end SECTION
            }  // end for
            }  // end SECTION
            }  // end for
            }  // end SECTION
            }  // end for
            }  // end SECTION
            }  // end for
            }  // end SECTION
            }  // end for num_steps
            }  // end SECTION
            }  // end for num_stages
        }  // end section solver
    }  // END FOR SOLVERS

    // implicit model
    external_function_casadi_free(&impl_ode_fun);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot_u);
    external_function_casadi_free(&impl_ode_jac_x_xdot_u);
    // gnsf functions:
    external_function_casadi_free(&phi_fun);
    external_function_casadi_free(&phi_fun_jac_y);
    external_function_casadi_free(&phi_jac_y_uhat);
    external_function_casadi_free(&f_lo_fun_jac_x1k1uz);
    external_function_casadi_free(&get_matrices_fun);

}  // END_TEST_CASE
