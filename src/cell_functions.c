#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cell_functions.h"
#include "mesh.h"

/* Function pointers for pressure cell operations */
const cell_ops_t cell_press_ops = {
	.cell_compute_beta 	   	= &press_compute_beta,
	.cell_compute_A        	= &press_compute_A,
	.cell_update_interior  	= &press_update_interior,
	.cell_update_boundary  	= &press_update_boundary,
	.cell_update_corner    	= &press_update_corner,
	.cell_compute_robin		= &press_compute_robin
};

/* Function pointers for diffusion cell operations */
const cell_ops_t cell_diff_ops = {
	.cell_compute_beta 	   	= &diff_compute_beta,
	.cell_compute_A        	= &diff_compute_A,
	.cell_update_interior  	= &diff_update_interior,
	.cell_update_boundary  	= &diff_update_boundary,
	.cell_update_corner    	= &diff_update_corner,
	.cell_compute_robin		= &diff_compute_robin
};

/* Retrieves the cell adjacent to the current cell */
/* 0 - up, 1 - right,  2 - down, 3 - left */
int get_adjacent_index(mesh_t *mesh, int direction, int cur_y, int cur_x)
{
    switch (direction) {
        case 3:
            return MESH_INDEX(cur_y, (cur_x - 1));
        case 1:
            return MESH_INDEX(cur_y, (cur_x + 1));
        case 0:
            return MESH_INDEX((cur_y - 1), cur_x);
        case 2:
            return MESH_INDEX((cur_y + 1), cur_x);
    }
    return 0;
}

/* Retrieves the cell diagonal to the current cell */
/* 0 - up left, 1 - up right, 2 - down right, 3 - down left */
int get_diagonal_index(mesh_t *mesh, int direction, int cur_y, int cur_x)
{
	switch (direction) {
		case 0:
			return MESH_INDEX((cur_y - 1), (cur_x - 1));
		case 1:
			return MESH_INDEX((cur_y - 1), (cur_x + 1));
		case 2:
			return MESH_INDEX((cur_y + 1), (cur_x + 1));
		case 3:
			return MESH_INDEX((cur_y + 1), (cur_x - 1));
	}
	return 0;
}

/* Relative permeability of the oil phase, K_ro */
static double rel_perm_o(cell_t *cell, global_mesh_params_t *global)
{
    return pow(1 - (cell->saturation/(1-global->sat_rel_o)), 2);
}

/* Relative permeability of the water phase, K_rw */
static double rel_perm_w(cell_t *cell, global_mesh_params_t *global)
{
    double num, denom;

    num = pow(cell->saturation - global->sat_rel_w, 2);
    denom = pow(1 - global->sat_rel_w, 2);

    return pow(num / denom, 2);
}

/* Total mobility, lambda */
double total_mobility(cell_t *cell, global_mesh_params_t *global)
{
    return (rel_perm_o(cell, global) / global->visc_o +
                rel_perm_w(cell, global) / global->visc_w);
}

/* Phase mobility of the oil phase, lambda_o */
double phase_mobility_o(cell_t *cell, global_mesh_params_t *global)
{
    double total_mob = total_mobility(cell, global);
    double rel_perm = rel_perm_o(cell, global);

    return rel_perm / (global->visc_o * total_mob);
}

/* Phase mobility of the water phase, lambda_w */
double phase_mobility_w(cell_t *cell, global_mesh_params_t *global)
{
    double total_mob = total_mobility(cell, global);
    double rel_perm = rel_perm_w(cell, global);

    return rel_perm / (global->visc_w * total_mob);
}

/* Derivative of the capillary pressure, P_c' */
double cap_pressure_deriv(cell_t *cell, global_mesh_params_t *global)
{
    double z, s_o, pc_deriv;

    s_o = 1 - global->sat_rel_o;
    z = pow(global->sat_rel_o, 2) * pow(s_o - global->sat_rel_w, -2);

    pc_deriv = pow(cell->saturation - global->sat_rel_w, -3);
    pc_deriv += z * pow(1 - cell->saturation, -3);
    pc_deriv *= -2 * global->eta;

    return pc_deriv;
}

/* Computes diffusion at the current cell */
void diff_compute_diffusion(mesh_t *mesh, int cur_y, int cur_x)
{
    double total_mob, w_mob, o_mob, pc_deriv;
    cell_t *cur_cell;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

    total_mob = total_mobility(cur_cell, &mesh->global);
    w_mob = phase_mobility_w(cur_cell, &mesh->global);
    o_mob = phase_mobility_o(cur_cell, &mesh->global);
    pc_deriv = cap_pressure_deriv(cur_cell, &mesh->global);

    cur_cell->diffusion = cur_cell->perm * total_mob * w_mob * o_mob * pc_deriv;
}

/* Combutes beta at the current cell for diffusion problem */
void diff_compute_beta(mesh_t *mesh, int cur_y, int cur_x, double beta_coef)
{
	cell_t *cur_cell, *adj_cell;
	double diff_eff;

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	for (int k = 0; k < 4; k++) {
		adj_cell = &mesh->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		if (adj_cell->diffusion == 0) {
			cur_cell->beta[k] = 0;
		} else {
			diff_eff = 2 * adj_cell->diffusion * cur_cell->diffusion;
			diff_eff /= (adj_cell->diffusion + cur_cell->diffusion);
			cur_cell->beta[k] = beta_coef * mesh->dim.h / diff_eff;
		}
	}
}

/* Computes beta at the current cell for pressure problem */
void press_compute_beta(mesh_t *mesh, int cur_y, int cur_x, double beta_coef)
{
    double perm_eff;
    cell_t *cur_cell, *adj_cell;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

    for (int k = 0; k < 4; k++) {
        adj_cell = &mesh->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
        if (adj_cell->perm == 0) {
            cur_cell->beta[k] = 0;
        } else {
            perm_eff = 2 * adj_cell->perm * cur_cell->perm;
            perm_eff /= (adj_cell->perm + cur_cell->perm);
            cur_cell->beta[k] = beta_coef * mesh->dim.h / perm_eff;
        }

        /* temp fix for weird inf bug */
        if ((cur_cell->beta[k] == INFINITY) || (cur_cell->beta[k] == NAN)) {
            cur_cell->beta[k] = 0;
        }
    }
}

/* Computes robin conditions for the pressure problem */
void press_compute_robin(mesh_t *mesh, int cur_y, int cur_x)
{
	cell_t *cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	for (int k = 0; k < 4; k++) {
		cur_cell->robin[k] = cur_cell->beta[k] * cur_cell->flux_p[k] + cur_cell->l_p[k];
	}
}

/* Computes robin conditions for the diffusion problem */
void diff_compute_robin(mesh_t *mesh, int cur_y, int cur_x)
{
	cell_t *cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	for (int k = 0; k < 4; k++) {
		cur_cell->robin[k] = cur_cell->beta[k] * cur_cell->flux_d[k] + cur_cell->l_d[k];
	}
}

/* Computes A_alpha = xi/(1+beta_alpha*xi), xi = 2k/h */
void press_compute_A(mesh_t *mesh, int cur_y, int cur_x)
{
    double xi;
    cell_t *cur_cell;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
    xi = 2 * cur_cell->perm / mesh->dim.h;

    for (int k = 0; k < 4; k++) {
        cur_cell->A_p[k] = xi / (1 + cur_cell->beta[k] * xi);
    }
}

/* Computes A_alpha = xi/(1+beta_alpha*xi), xi = 2k/h */
void diff_compute_A(mesh_t *mesh, int cur_y, int cur_x)
{
    double xi;
    cell_t *cur_cell;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
    xi = 2 * cur_cell->diffusion / mesh->dim.h;

    for (int k = 0; k < 4; k++) {
        cur_cell->A_d[k] = xi / (1 + cur_cell->beta[k] * xi);
    }
}

/* Updates the pressure, flux, and robin conditions for a cell */
void press_update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x)
{
    int k;
    cell_t *cur_cell, *cur_cell_old, *adj_cell;
    double sum_A, sum_A_R, A;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
    cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

    sum_A = 0;
    sum_A_R = 0;

    for (k = 0; k < 4; k++) {
        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
        sum_A += cur_cell_old->A_p[k];
        sum_A_R += cur_cell_old->A_p[k] * adj_cell->robin[(k + 2) % 4];
    }

    /* Updates the pressure at the current cell on the new mesh */
    cur_cell->pressure = (cur_cell_old->source * mesh->dim.h + sum_A_R) / sum_A;

    /* Updates the flux at the current cell on the new mesh */
    for (k = 0; k < 4; k++) {
        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
        A = cur_cell->A_p[k];
        cur_cell->flux_p[k] = A * (cur_cell->pressure - adj_cell->robin[(k + 2) % 4]);
    }

    /* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
        cur_cell->l_p[k] = cur_cell->beta[k] * cur_cell->flux_p[k] + adj_cell->robin[(k + 2) % 4];
    }
}

/* Updates the saturation and flux for the diffusion problem on interior cells */
void diff_update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x)
{
	int k;
	cell_t *cur_cell, *cur_cell_old, *adj_cell;
	double A, sum_A, sum_A_R, num, denom, phi_h_dt;

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

	sum_A = 0;
	sum_A_R = 0;

	/* Updates saturation for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		sum_A += cur_cell_old->A_d[k];
		sum_A_R += cur_cell_old->A_d[k] * adj_cell->robin[(k + 2) % 4];
	}

	phi_h_dt = (mesh->global.porosity * mesh->dim.h / mesh->dim.dt);

	num = cur_cell_old->source_d * mesh->dim.h + sum_A_R + phi_h_dt * cur_cell_old->saturation_prev;
	denom = phi_h_dt + sum_A;

	cur_cell->saturation = num / denom;

	/* Updates flux for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
        A = cur_cell->A_d[k];
        cur_cell->flux_d[k] = A * (cur_cell->saturation - adj_cell->robin[(k + 2) % 4]);
	}

	/* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
        cur_cell->l_d[k] = cur_cell->beta[k] * cur_cell->flux_d[k] + adj_cell->robin[(k + 2) % 4];
    }
}

void press_update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
        int boundary_side)
{
    int k;
    cell_t *cur_cell, *cur_cell_old, *adj_cell;
    double sum_A, sum_A_R, A;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
    cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

    sum_A = 0;
    sum_A_R = 0;

    for (k = 0; k < 4; k++) {
        if (k != boundary_side) {
            adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
            sum_A += cur_cell_old->A_p[k];
            sum_A_R += cur_cell_old->A_p[k] * adj_cell->robin[(k + 2) % 4];
        }
    }

    /* Updates the pressure at the current cell on the new mesh */
    cur_cell->pressure = (cur_cell_old->source * mesh->dim.h + sum_A_R) / sum_A;

    /* Updates the flux at the current cell on the new mesh */
    for (k = 0; k < 4; k++) {
        if (k != boundary_side) {
            adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
            A = cur_cell->A_p[k];
            cur_cell->flux_p[k] = A * (cur_cell->pressure - adj_cell->robin[(k + 2) % 4]);
        }
    }

    /* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        if (k != boundary_side) {
            adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
            cur_cell->l_p[k] = cur_cell->beta[k] * cur_cell->flux_p[k] +
                                adj_cell->robin[(k + 2) % 4];
        }
    }
}

void diff_update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
        int boundary_side)
{
	int k;
	cell_t *cur_cell, *cur_cell_old, *adj_cell;
	double A, sum_A, sum_A_R, num, denom, phi_h_dt;

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

	sum_A = 0;
	sum_A_R = 0;

	/* Updates saturation for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if (k != boundary_side) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
			sum_A += cur_cell_old->A_d[k];
			sum_A_R += cur_cell_old->A_d[k] * adj_cell->robin[(k + 2) % 4];
		}
	}

	phi_h_dt = (mesh->global.porosity * mesh->dim.h / mesh->dim.dt);

	num = cur_cell_old->source_d * mesh->dim.h + sum_A_R + phi_h_dt * cur_cell_old->saturation_prev;
	denom = phi_h_dt + sum_A;

	cur_cell->saturation = num / denom;

	/* Updates flux for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if (k != boundary_side) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
	        A = cur_cell->A_d[k];
	        cur_cell->flux_d[k] = A * (cur_cell->saturation - adj_cell->robin[(k + 2) % 4]);
		}
	}

	/* Updates the saturation at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
		if (k != boundary_side) {
	        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		    cur_cell->l_d[k] = cur_cell->beta[k] * cur_cell->flux_d[k] + adj_cell->robin[(k + 2) % 4];
		}
	}
}

void press_update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2)
{
    int k;
    cell_t *cur_cell, *cur_cell_old, *adj_cell;
    double sum_A, sum_A_R, A;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
    cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

    sum_A = 0;
    sum_A_R = 0;

    for (k = 0; k < 4; k++) {
        if ((k != boundary_side1) && (k != boundary_side2)) {
            adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
            sum_A += cur_cell_old->A_p[k];
            sum_A_R += cur_cell_old->A_p[k] * adj_cell->robin[(k + 2) % 4];
        }
    }

    /* Updates the pressure at the current cell on the new mesh */
    cur_cell->pressure = (cur_cell_old->source * mesh->dim.h + sum_A_R) / sum_A;

    /* Updates the flux at the current cell on the new mesh */
    for (k = 0; k < 4; k++) {
        if ((k != boundary_side1) && (k != boundary_side2)) {
            adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
            A = cur_cell->A_p[k];
            cur_cell->flux_p[k] = A * (cur_cell->pressure - adj_cell->robin[(k + 2) % 4]);
        }
    }

    /* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        if ((k != boundary_side1) && (k != boundary_side2)) {
            adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
            cur_cell->l_p[k] = cur_cell->beta[k] * cur_cell->flux_p[k] +
                                adj_cell->robin[(k + 2) % 4];
        }
    }
}

void diff_update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2)
{
	int k;
	cell_t *cur_cell, *cur_cell_old, *adj_cell;
	double A, sum_A, sum_A_R, num, denom, phi_h_dt;

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

	sum_A = 0;
	sum_A_R = 0;

	/* Updates saturation for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if ((k != boundary_side1) && (k != boundary_side2)) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
			sum_A += cur_cell_old->A_d[k];
			sum_A_R += cur_cell_old->A_d[k] * adj_cell->robin[(k + 2) % 4];
		}
	}

	phi_h_dt = (mesh->global.porosity * mesh->dim.h / mesh->dim.dt);

	num = cur_cell_old->source_d * mesh->dim.h + sum_A_R + phi_h_dt * cur_cell_old->saturation_prev;
	denom = phi_h_dt + sum_A;

	cur_cell->saturation = num / denom;

	/* Updates flux for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if ((k != boundary_side1) && (k != boundary_side2)) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
	        A = cur_cell->A_d[k];
	        cur_cell->flux_d[k] = A * (cur_cell->saturation - adj_cell->robin[(k + 2) % 4]);
		}
	}

	/* Updates the saturation at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
		if ((k != boundary_side1) && (k != boundary_side2)) {
	        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		    cur_cell->l_d[k] = cur_cell->beta[k] * cur_cell->flux_d[k] + adj_cell->robin[(k + 2) % 4];
		}
	}
}

void diff_update_boundary_dirichlet(mesh_t *mesh, mesh_t *mesh_old, int cur_y,
	int cur_x, int boundary_side)
{
	int k;
	cell_t *cur_cell, *cur_cell_old, *adj_cell;
	double A, sum_A, sum_A_R, num, denom, phi_h_dt, xi;

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

	sum_A = 0;
	sum_A_R = 0;

	xi = 2 * cur_cell->diffusion / mesh->dim.h;

	/* Updates saturation for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if (k != boundary_side) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
			sum_A += cur_cell_old->A_d[k];
			sum_A_R += cur_cell_old->A_d[k] * adj_cell->robin[(k + 2) % 4];
		}
	}

	phi_h_dt = (mesh->global.porosity * mesh->dim.h / mesh->dim.dt);

	num = cur_cell_old->source_d * mesh->dim.h + sum_A_R + phi_h_dt * cur_cell_old->saturation_prev;
	denom = phi_h_dt + sum_A + xi;

	cur_cell->saturation = num / denom;

	/* Updates flux for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if (k != boundary_side) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
	        A = cur_cell->A_d[k];
	        cur_cell->flux_d[k] = A * (cur_cell->saturation - adj_cell->robin[(k + 2) % 4]);
		}
	}

	/* Updates the saturation at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
		if (k != boundary_side) {
	        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		    cur_cell->l_d[k] = cur_cell->beta[k] * cur_cell->flux_d[k] + adj_cell->robin[(k + 2) % 4];
		}
	}
}

void diff_update_corner_dirichlet(mesh_t *mesh, mesh_t *mesh_old, int cur_y,
	int cur_x, int boundary_side1, int boundary_side2)
{
	int k;
	cell_t *cur_cell, *cur_cell_old, *adj_cell;
	double A, sum_A, sum_A_R, num, denom, phi_h_dt, xi;

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh_old->cell[MESH_INDEX(cur_y, cur_x)];

	sum_A = 0;
	sum_A_R = 0;

	xi = 2 * cur_cell->diffusion / mesh->dim.h;

	/* Updates saturation for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if ((k != boundary_side1) && (k != boundary_side2)) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
			sum_A += cur_cell_old->A_d[k];
			sum_A_R += cur_cell_old->A_d[k] * adj_cell->robin[(k + 2) % 4];
		}
	}

	phi_h_dt = (mesh->global.porosity * mesh->dim.h / mesh->dim.dt);

	num = cur_cell_old->source_d * mesh->dim.h + sum_A_R + phi_h_dt * cur_cell_old->saturation_prev;
	denom = phi_h_dt + sum_A + xi;

	cur_cell->saturation = num / denom;

	/* Updates flux for the current cell on the new mesh */
	for (k = 0; k < 4; k++) {
		if ((k != boundary_side1) && (k != boundary_side2)) {
			adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
	        A = cur_cell->A_d[k];
	        cur_cell->flux_d[k] = A * (cur_cell->saturation - adj_cell->robin[(k + 2) % 4]);
		}
	}

	/* Updates the saturation at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
		if ((k != boundary_side1) && (k != boundary_side2)) {
	        adj_cell = &mesh_old->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		    cur_cell->l_d[k] = cur_cell->beta[k] * cur_cell->flux_d[k] + adj_cell->robin[(k + 2) % 4];
		}
	}
}
