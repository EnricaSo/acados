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

#include "acados/ocp_nlp/ocp_nlp_constraints.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"
#include "acados/ocp_qp/ocp_qp_common.h"



/************************************************
* config
************************************************/

int ocp_nlp_constraints_config_calculate_size()
{

	int size = 0;

	size += sizeof(ocp_nlp_constraints_config);

	return size;

}



ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory)
{

	char *c_ptr = raw_memory;

	ocp_nlp_constraints_config *config = (ocp_nlp_constraints_config *) c_ptr;
	c_ptr += sizeof(ocp_nlp_constraints_config);

	return config;

}



/************************************************
* dims
************************************************/

int ocp_nlp_constraints_dims_calculate_size(void *config_)
{
    int size = sizeof(ocp_nlp_constraints_dims);

    return size;
}



void *ocp_nlp_constraints_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_dims *dims = (ocp_nlp_constraints_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_dims);

    assert((char *) raw_memory + ocp_nlp_constraints_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void  ocp_nlp_constraints_dims_initialize(void *config_, void *dims_, int nx, int nu, int nbx, int nbu, int ng, int nh, int nq, int ns)
{
	ocp_nlp_constraints_dims *dims = dims_;

	dims->nx = nx;
	dims->nu = nu;
	dims->nbx = nbx;
	dims->nbu = nbu;
	dims->nb = nbx+nbu;
	dims->ng = ng;
	dims->nh = nh;
	dims->nq = nq;
	dims->ns = ns;

	return;
}



/************************************************
* linear constraints
************************************************/

/* model */

int ocp_nlp_constraints_model_calculate_size(void *config, void *dims_)
{
	ocp_nlp_constraints_dims *dims = dims_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;

	int size = 0;

    size += sizeof(ocp_nlp_constraints_model);

	size += sizeof(int)*nb;  // idxb
	size += blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // d
	size += blasfeo_memsize_dmat(nu+nx, ng); // DCt

	size += 64; // blasfeo_mem align

	return size;
}



void *ocp_nlp_constraints_model_assign(void *config, void *dims_, void *raw_memory)
{
	ocp_nlp_constraints_dims *dims = dims_;

	char *c_ptr = (char *) raw_memory;

	// extract sizes
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;

	// struct
    ocp_nlp_constraints_model *model = (ocp_nlp_constraints_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_model);

	// dims
//	model->dims = dims;

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// blasfeo_dmat
	// DCt
	assign_and_advance_blasfeo_dmat_mem(nu+nx, ng, &model->DCt, &c_ptr);

	// blasfeo_dvec
	// d
	assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, &model->d, &c_ptr);

    // idxb
    assign_and_advance_int(dims->nbx+dims->nbu, &model->idxb, &c_ptr);

	// h
//	model->h = NULL;

	// assert
    assert((char *) raw_memory + ocp_nlp_constraints_model_calculate_size(config, dims) >= c_ptr);

	return model;
}



/* options */

int ocp_nlp_constraints_opts_calculate_size(void *config_, void *dims_)
{
    int size = 0;

    size += sizeof(ocp_nlp_constraints_opts);

    return size;
}



void *ocp_nlp_constraints_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_opts *opts = (ocp_nlp_constraints_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_opts);

    assert((char*)raw_memory + ocp_nlp_constraints_opts_calculate_size(config_, dims_) >= c_ptr);

    return opts;
}



void ocp_nlp_constraints_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
//	ocp_nlp_constraints_opts *opts = opts_;

	return;

}



void ocp_nlp_constraints_opts_update(void *config_, void *dims_, void *opts_)
{
//	ocp_nlp_constraints_opts *opts = opts_;

	return;

}



/* memory */

int ocp_nlp_constraints_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
	ocp_nlp_constraints_dims *dims = dims_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;

	int size = 0;

    size += sizeof(ocp_nlp_constraints_memory);

	size += 1*blasfeo_memsize_dvec(2*nb+2*ng+2*nh);  // fun
	size += 1*blasfeo_memsize_dvec(nu+nx);  // adj

	size += 1*64;  // blasfeo_mem align

	return size;
}



void *ocp_nlp_constraints_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
	ocp_nlp_constraints_dims *dims = dims_;

	char *c_ptr = (char *) raw_memory;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;

	// struct
    ocp_nlp_constraints_memory *memory = (ocp_nlp_constraints_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_memory);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// fun
	assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, &memory->fun, &c_ptr);
	// adj
	assign_and_advance_blasfeo_dvec_mem(nu+nx, &memory->adj, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_constraints_memory_calculate_size(config_, dims, opts_) >= c_ptr);

	return memory;
}



struct blasfeo_dvec *ocp_nlp_constraints_memory_get_fun_ptr(void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	return &memory->fun;
}



struct blasfeo_dvec *ocp_nlp_constraints_memory_get_adj_ptr(void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	return &memory->adj;
}



void ocp_nlp_constraints_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->ux = ux;
}



void ocp_nlp_constraints_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->lam = lam;
}



void ocp_nlp_constraints_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->DCt = DCt;
}



void ocp_nlp_constraints_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->RSQrq = RSQrq;
}



void ocp_nlp_constraints_memory_set_idxb_ptr(int *idxb, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->idxb = idxb;
}



/* workspace */

int ocp_nlp_constraints_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
	ocp_nlp_constraints_dims *dims = dims_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;
	int nq = dims->nq;

	int size = 0;

    size += sizeof(ocp_nlp_constraints_workspace);

	size += 1*blasfeo_memsize_dvec(nb+ng+nh);  // tmp_nbgh
	if (nq > 0) {
//		size += nq*(nx+nu)*sizeof(double);
//		size += 1*blasfeo_memsize_dmat(nx+nu, nq);
	}

	size += 2*64;  // blasfeo_mem align

	return size;

}



static void ocp_nlp_constraints_cast_workspace(void *config_, void *dims_, void *opts_, void *work_)
{
	ocp_nlp_constraints_dims *dims = dims_;
	ocp_nlp_constraints_workspace *work = work_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;
	int nq = dims->nq;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_constraints_workspace);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// tmp_nbgh
	assign_and_advance_blasfeo_dvec_mem(nb+ng+nh, &work->tmp_nbgh, &c_ptr);
	if (nq > 0) {
//		c_ptr += nq*(nx+nu)*sizeof(double);
//		align_char_to(64, &c_ptr);
//		assign_and_advance_blasfeo_dmat_mem(nx+nu, nq, &work->jacobian_quadratic, &c_ptr);	
	}

    assert((char *)work + ocp_nlp_constraints_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

	return;
}



/* functions */

void ocp_nlp_constraints_initialize(void *config_, void *dims_, void *model_, void *opts, void *memory_, void *work_)
{

	ocp_nlp_constraints_dims *dims = dims_;
	ocp_nlp_constraints_model *model = model_;
	ocp_nlp_constraints_memory *memory = memory_;

	// loop index
	int j;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;

	// initialize idxb
	for (j=0; j<nb; j++)
	{
		memory->idxb[j] = model->idxb[j];
	}

	// initialize general constraints matrix
	blasfeo_dgecp(nu+nx, ng, &model->DCt, 0, 0, memory->DCt, 0, 0);

	return;

}


void ocp_nlp_constraints_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{

	ocp_nlp_constraints_dims *dims = dims_;
	ocp_nlp_constraints_model *model = model_;
	ocp_nlp_constraints_memory *memory = memory_;
	ocp_nlp_constraints_workspace *work = work_;

	ocp_nlp_constraints_cast_workspace(config_, dims, opts_, work_);

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;
	int nq = dims->nq;

	ext_fun_arg_t ext_fun_type_in[2];
	void *ext_fun_in[2]; // XXX large enough ?
	ext_fun_arg_t ext_fun_type_out[2];
	void *ext_fun_out[2]; // XXX large enough ?

	blasfeo_dvecex_sp(nb, 1.0, model->idxb, memory->ux, 0, &work->tmp_nbgh, 0);
	blasfeo_dgemv_t(nu+nx, ng, 1.0, memory->DCt, 0, 0, memory->ux, 0, 0.0, &work->tmp_nbgh, nb, &work->tmp_nbgh, nb);

	if (nh>0)
	{
		// XXX DANGER !!! DO NOT TRY THIS AT HOME !!!
		// Cast of BLASFEO internals, only if you know what is happening !!!
		work->tmp_h = work->tmp_nbgh;
		work->tmp_h.pa = &(BLASFEO_DVECEL(&work->tmp_nbgh, nb+ng));
		work->tmp_Jht = memory->DCt[0];
		work->tmp_Jht.pA = &(BLASFEO_DMATEL(memory->DCt, 0, ng)); // only works for col offsets !!!

		ext_fun_type_in[0] = BLASFEO_VEC;
		ext_fun_in[0] = memory->ux; // ux: nu+nx

		ext_fun_type_out[0] = BLASFEO_VEC;
		ext_fun_out[0] = &work->tmp_h; // fun: nh
		ext_fun_type_out[1] = BLASFEO_MAT;
		ext_fun_out[1] = &work->tmp_Jht; // jac': (nx+nu) * nh

		model->h->evaluate(model->h, ext_fun_type_in, ext_fun_in, ext_fun_type_out, ext_fun_out);

	}

	if (nq>0)
	{
//		if (nh != 1) {
			printf("Not implemented");
			exit(1);
//		}
//		double lam = blasfeo_dvecex1(memory->lam, 2*(nb+ng)+nh) - blasfeo_dvecex1(memory->lam, nb+ng);
//		blasfeo_pack_tran_dmat(nq, nx+nu, work->nl_constraint_output+nh*(1+nx+nu), nq, &work->jacobian_quadratic, 0, 0);
		// NOTE(giaf) here the Hessian in overwritten, not updated. Is this correct?
//		blasfeo_dsyrk_ln(nx+nu, nq, 2*lam, &work->jacobian_quadratic, 0, 0, &work->jacobian_quadratic, 0, 0,
//			0.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);
	}

	blasfeo_daxpy(nb+ng+nh, -1.0, &work->tmp_nbgh, 0, &model->d, 0, &memory->fun, 0);
	blasfeo_daxpy(nb+ng+nh, -1.0, &model->d, nb+ng+nh, &work->tmp_nbgh, 0, &memory->fun, nb+ng+nh);

	// nlp_mem: ineq_adj
	blasfeo_dvecse(nu+nx, 0.0, &memory->adj, 0);
	blasfeo_daxpy(nb+ng+nh, -1.0, memory->lam, nb+ng+nh, memory->lam, 0, &work->tmp_nbgh, 0);
	blasfeo_dvecad_sp(nb, 1.0, &work->tmp_nbgh, 0, model->idxb, &memory->adj, 0);
	blasfeo_dgemv_n(nu+nx, ng+nh, 1.0, memory->DCt, 0, 0, &work->tmp_nbgh, nb, 1.0, &memory->adj, 0, &memory->adj, 0);

	return;

}



void ocp_nlp_constraints_config_initialize_default(void *config_)
{
	ocp_nlp_constraints_config *config = config_;

	config->dims_calculate_size = &ocp_nlp_constraints_dims_calculate_size;
	config->dims_assign = &ocp_nlp_constraints_dims_assign;
	config->dims_initialize = &ocp_nlp_constraints_dims_initialize;
	config->model_calculate_size = &ocp_nlp_constraints_model_calculate_size;
	config->model_assign = &ocp_nlp_constraints_model_assign;
	config->opts_calculate_size = &ocp_nlp_constraints_opts_calculate_size;
	config->opts_assign = &ocp_nlp_constraints_opts_assign;
	config->opts_initialize_default = &ocp_nlp_constraints_opts_initialize_default;
	config->opts_update = &ocp_nlp_constraints_opts_update;
	config->memory_calculate_size = &ocp_nlp_constraints_memory_calculate_size;
	config->memory_assign = &ocp_nlp_constraints_memory_assign;
	config->memory_get_fun_ptr = &ocp_nlp_constraints_memory_get_fun_ptr;
	config->memory_get_adj_ptr = &ocp_nlp_constraints_memory_get_adj_ptr;
	config->memory_set_ux_ptr = &ocp_nlp_constraints_memory_set_ux_ptr;
	config->memory_set_lam_ptr = &ocp_nlp_constraints_memory_set_lam_ptr;
	config->memory_set_DCt_ptr = &ocp_nlp_constraints_memory_set_DCt_ptr;
	config->memory_set_RSQrq_ptr = &ocp_nlp_constraints_memory_set_RSQrq_ptr;
	config->memory_set_idxb_ptr = &ocp_nlp_constraints_memory_set_idxb_ptr;
	config->workspace_calculate_size = &ocp_nlp_constraints_workspace_calculate_size;
	config->initialize = &ocp_nlp_constraints_initialize;
	config->update_qp_matrices = &ocp_nlp_constraints_update_qp_matrices;
	config->config_initialize_default = &ocp_nlp_constraints_config_initialize_default;

	return;

}




