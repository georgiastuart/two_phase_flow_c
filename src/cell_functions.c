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

const cell_ops_t cell_trans_ops = {
	.cell_compute_beta		= NULL,
	.cell_compute_A			= NULL,
	.cell_update_interior	= &trans_update_interior,
	.cell_update_boundary	= &trans_update_boundary,
	.cell_update_corner		= &trans_update_corner,
	.cell_compute_robin		= NULL
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

/* Returns the position for method of characteristics calculation during transport */
/* Direction = 0 for x, 1 for y */
double trans_get_old_position(mesh_t *mesh, int cur_y, int cur_x, int direction)
{
	cell_t *cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	double pos;
	pos = -mesh->global.porosity * cur_cell->pm_w_deriv;
	pos *= mesh->dim.dt_transport;

	if (direction == 0)
		return pos * cur_cell->velocity_x;
	else if (direction == 1)
		return pos * cur_cell->velocity_y;
	else
		return 0;
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

	double total_mob, adj_total_mob;

    cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	total_mob = total_mobility(cur_cell, &mesh->global);

    for (int k = 0; k < 4; k++) {
        adj_cell = &mesh->cell[get_adjacent_index(mesh, k, cur_y, cur_x)];
		adj_total_mob = total_mobility(adj_cell, &mesh->global);
		printf("cur y %d, cur x %d, adj total mob %e\n", cur_y, cur_x, adj_total_mob);
        if (adj_cell->perm == 0) {
            cur_cell->beta[k] = 0;
        } else {
            perm_eff = 2 * adj_cell->perm * adj_total_mob * cur_cell->perm * total_mob;
            perm_eff /= (adj_cell->perm * adj_total_mob + cur_cell->perm * total_mob);
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
    xi = 2 * cur_cell->perm * total_mobility(cur_cell, &mesh->global)/ mesh->dim.h;
	// printf("lambda: %e\n", total_mobility(cur_cell, &mesh->global));
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

/* Returns the average saturation at t_old for the transport method of characteristics */
double trans_get_average_sat(cell_t *cur_cell, cell_t *adj_cell_hor, cell_t *adj_cell_vert,
							 cell_t *adj_cell_diag, double y_comp, double x_comp)
{
	double saturation;

	double x_1 = (1.0 - x_comp);
	double y_1 = (1.0 - y_comp);

	/* Finds average saturation weighted to areas of the new cell in each */
	/* adjacent cell */
	saturation = cur_cell->saturation * x_1 * y_1;
	saturation += adj_cell_hor->saturation * x_comp * y_1;
	saturation += adj_cell_vert->saturation * x_1 * y_comp;
	saturation += adj_cell_diag->saturation * x_comp * y_comp;

	return saturation;
}

static int get_quadrant(double y_comp, double x_comp)
{
	double sign_y = copysign(1.0, y_comp), sign_x = copysign(1.0, x_comp);

	/* Selects the quadrant */
	if ((sign_y == 1.0) && (sign_x == -1.0)) {
		return 2;
	} else if ((sign_y == 1.0) && (sign_x == 1.0)) {
		return 1;
	} else if ((sign_y == -1.0) && (sign_x == 1.0)) {
		return 4;
	} else {
		return 3;
	}
}

/* Cell duplication functions for transport step */
static void dup_cell_hor(mesh_t *mesh, cell_t *cur_cell, cell_t **adj_cell_hor,
						 cell_t **adj_cell_vert, cell_t **adj_cell_diag,
						 int cur_y, int cur_x, int side_vert)
{
	*adj_cell_vert = &mesh->cell[get_adjacent_index(mesh, side_vert, cur_y, cur_x)];
	*adj_cell_hor = cur_cell;
	*adj_cell_diag = *adj_cell_vert;
}

static void dup_cell_vert(mesh_t *mesh, cell_t *cur_cell, cell_t **adj_cell_hor,
						 cell_t **adj_cell_vert, cell_t **adj_cell_diag,
						 int cur_y, int cur_x, int side_hor)
{
	*adj_cell_vert = cur_cell;
	*adj_cell_hor = &mesh->cell[get_adjacent_index(mesh, side_hor, cur_y, cur_x)];
	*adj_cell_diag = *adj_cell_hor;
}

static void dup_cell_corner(mesh_t *mesh, cell_t *cur_cell, cell_t **adj_cell_hor,
						 cell_t **adj_cell_vert, cell_t **adj_cell_diag)
{
	*adj_cell_vert = cur_cell;
	*adj_cell_hor = cur_cell;
	*adj_cell_diag = cur_cell;
}

/* Cell assignments for normal quadrants */
static void assign_cells(mesh_t *mesh, cell_t **adj_cell_hor, cell_t **adj_cell_vert,
						 cell_t **adj_cell_diag, int cur_y, int cur_x, int quad)
{
	switch (quad) {
		case 1:
			*adj_cell_vert = &mesh->cell[get_adjacent_index(mesh, 0, cur_y, cur_x)];
			*adj_cell_hor = &mesh->cell[get_adjacent_index(mesh, 1, cur_y, cur_x)];
			*adj_cell_diag = &mesh->cell[get_diagonal_index(mesh, 1, cur_y, cur_x)];
			break;
		case 2:
			*adj_cell_vert = &mesh->cell[get_adjacent_index(mesh, 0, cur_y, cur_x)];
			*adj_cell_hor = &mesh->cell[get_adjacent_index(mesh, 3, cur_y, cur_x)];
			*adj_cell_diag = &mesh->cell[get_diagonal_index(mesh, 0, cur_y, cur_x)];
			break;
		case 3:
			*adj_cell_vert = &mesh->cell[get_adjacent_index(mesh, 2, cur_y, cur_x)];
			*adj_cell_hor = &mesh->cell[get_adjacent_index(mesh, 3, cur_y, cur_x)];
			*adj_cell_diag = &mesh->cell[get_diagonal_index(mesh, 3, cur_y, cur_x)];
			break;
		default:
			*adj_cell_vert = &mesh->cell[get_adjacent_index(mesh, 2, cur_y, cur_x)];
			*adj_cell_hor = &mesh->cell[get_adjacent_index(mesh, 1, cur_y, cur_x)];
			*adj_cell_diag = &mesh->cell[get_diagonal_index(mesh, 2, cur_y, cur_x)];
	}
}

/* Interior update for tranport step */
void trans_update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x)
{
	cell_t *cur_cell, *cur_cell_old, *adj_cell_vert, *adj_cell_hor, *adj_cell_diag;
	double y_comp, x_comp;

	/* Gets the x and y components of the backwards in time vector for method */
	/* of characteristics */
	y_comp = trans_get_old_position(mesh_old, cur_y, cur_x, 1);
	x_comp = trans_get_old_position(mesh_old, cur_y, cur_x, 0);

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	/* Gets the quadrant the foot is in */
	int quad = get_quadrant(y_comp, x_comp);

	/* Sets components to proportions of full cell */
	y_comp = fabs(y_comp / mesh->dim.h);
	x_comp = fabs(x_comp / mesh->dim.h);

	/* Selects the configuration of cells */
	assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert, &adj_cell_diag, cur_y, cur_x, quad);

	cur_cell->saturation = trans_get_average_sat(cur_cell_old, adj_cell_hor,
							adj_cell_vert, adj_cell_diag, y_comp, x_comp);
}

/* Boundary update for tranport step */
void trans_update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
						   int boundary_side)
{
	cell_t *cur_cell, *cur_cell_old, *adj_cell_vert, *adj_cell_hor, *adj_cell_diag;
	double y_comp, x_comp;

	/* Gets the x and y components of the backwards in time vector for method */
	/* of characteristics */
	y_comp = trans_get_old_position(mesh_old, cur_y, cur_x, 1);
	x_comp = trans_get_old_position(mesh_old, cur_y, cur_x, 0);

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	/* Gets the quadrant the foot is in */
	int quad = get_quadrant(y_comp, x_comp);

	/* Sets components to proportions of full cell */
	y_comp = fabs(y_comp / mesh->dim.h);
	x_comp = fabs(x_comp / mesh->dim.h);

	/* Selects the correct configuration of adjacent cells */

	switch (boundary_side) {
		case 0:
			switch (quad) {
				case 1:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
								  &adj_cell_diag, cur_y, cur_x, 1);
					break;
				case 2:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							  &adj_cell_diag, cur_y, cur_x, 3);
					break;
				case 3:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
								 &adj_cell_diag, cur_y, cur_x, quad);
					break;
				default:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
								 &adj_cell_diag, cur_y, cur_x, quad);
					break;
			}
			break;
		case 1:
			switch (quad) {
				case 1:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							  &adj_cell_diag, cur_y, cur_x, 0);
					break;
				case 2:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
							 &adj_cell_diag, cur_y, cur_x, quad);
					break;
				case 3:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
							 &adj_cell_diag, cur_y, cur_x, quad);
					break;
				default:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
						  &adj_cell_diag, cur_y, cur_x, 2);
					break;
			}
			break;
		case 2:
			switch (quad) {
				case 1:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
						 &adj_cell_diag, cur_y, cur_x, quad);
					break;
				case 2:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
						 &adj_cell_diag, cur_y, cur_x, quad);
					break;
				case 3:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							  &adj_cell_diag, cur_y, cur_x, 3);
					break;
				default:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							  &adj_cell_diag, cur_y, cur_x, 1);
					break;
			}
			break;
		default:
			switch (quad) {
				case 1:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
						 &adj_cell_diag, cur_y, cur_x, quad);
					break;
				case 2:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
						  &adj_cell_diag, cur_y, cur_x, 0);
					break;
				case 3:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
						  &adj_cell_diag, cur_y, cur_x, 2);
					break;
				default:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
						 &adj_cell_diag, cur_y, cur_x, quad);
					break;
			}
			break;
	}

	cur_cell->saturation = trans_get_average_sat(cur_cell_old, adj_cell_hor,
							adj_cell_vert, adj_cell_diag, y_comp, x_comp);
}

/* Corner configuration for transport corner calculation */
/* Corners are labeled like quadrants:
	1 - Upper right
	2 - Upper left
	3 - Lower left
	4 - Lower right
*/
static int get_corner_config(int boundary_side1, int boundary_side2)
{
	if ((boundary_side1 == 0) || (boundary_side2 == 0)) {
		if ((boundary_side1 == 3) || (boundary_side2 == 3)) {
			return 2;
		} else {
			return 1;
		}
	} else {
		if ((boundary_side1 == 3) || (boundary_side2 == 3)) {
			return 3;
		} else {
			return 4;
		}
	}
}

void trans_update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2)
{
	cell_t *cur_cell, *cur_cell_old, *adj_cell_vert, *adj_cell_hor, *adj_cell_diag;
	double y_comp, x_comp;

	/* Gets the x and y components of the backwards in time vector for method */
	/* of characteristics */
	y_comp = trans_get_old_position(mesh_old, cur_y, cur_x, 1);
	x_comp = trans_get_old_position(mesh_old, cur_y, cur_x, 0);

	cur_cell = &mesh->cell[MESH_INDEX(cur_y, cur_x)];
	cur_cell_old = &mesh->cell[MESH_INDEX(cur_y, cur_x)];

	/* Gets the quadrant the foot is in */
	int quad = get_quadrant(y_comp, x_comp);
	int corner_type = get_corner_config(boundary_side1, boundary_side2);

	/* Sets components to proportions of full cell */
	y_comp = fabs(y_comp / mesh->dim.h);
	x_comp = fabs(x_comp / mesh->dim.h);
	// printf("ycomp %e, xcomp %e\n", y_comp, x_comp);

	/* Selects the configuration of cells */
	switch (corner_type) {
		case 1:
			switch (quad) {
				case 1:
					dup_cell_corner(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
								&adj_cell_diag);
					break;
				case 2:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							  	&adj_cell_diag, cur_y, cur_x, 3);
					break;
				case 3:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
						 		&adj_cell_diag, cur_y, cur_x, quad);
					break;
				default:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
					  			&adj_cell_diag, cur_y, cur_x, 2);
					break;
			}
			break;
		case 2:
			switch (quad) {
				case 1:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
						  		&adj_cell_diag, cur_y, cur_x, 1);
					break;
				case 2:
					dup_cell_corner(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
								&adj_cell_diag);
					break;
				case 3:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
				  				&adj_cell_diag, cur_y, cur_x, 2);
					break;
				default:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
								&adj_cell_diag, cur_y, cur_x, quad);
					break;
			}
			break;
		case 3:
			switch (quad) {
				case 1:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag, cur_y, cur_x, quad);
					break;
				case 2:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag, cur_y, cur_x, 0);
					break;
				case 3:
					dup_cell_corner(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag);
					break;
				default:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag, cur_y, cur_x, 1);
					break;
			}
			break;
		default:
			switch (quad) {
				case 1:
					dup_cell_hor(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag, cur_y, cur_x, 0);
					break;
				case 2:
					assign_cells(mesh_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag, cur_y, cur_x, quad);
					break;
				case 3:
					dup_cell_vert(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag, cur_y, cur_x, 3);
					break;
				default:
					dup_cell_corner(mesh, cur_cell_old, &adj_cell_hor, &adj_cell_vert,
							&adj_cell_diag);
					break;
			}
	}

	cur_cell->saturation = trans_get_average_sat(cur_cell_old, adj_cell_hor,
							adj_cell_vert, adj_cell_diag, y_comp, x_comp);
	// printf("sat: %e\n", cur_cell->saturation);
}

/* Boundary update for the diffusion test problem */
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

/* Corner update for the diffusion test problem */
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
