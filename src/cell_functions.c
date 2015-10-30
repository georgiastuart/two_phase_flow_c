#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cell_functions.h"

/* Retrieves the cell adjacent to the current cell */
/* 0 - up, 1 - right,  2 - down, 3 - left */
int get_adjacent_index(int direction, int cur_y, int cur_x)
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

/* Computes beta at each gridpoint on the computational domain */
void compute_beta(cell_t *mesh, int cur_y, int cur_x, double beta_coef)
{
    int k;
    double perm_eff;
    cell_t *cur_cell, *adj_cell;

    cur_cell = &mesh[MESH_INDEX(cur_y, cur_x)];

    for (k = 0; k < 4; k++) {
        adj_cell = &mesh[get_adjacent_index(k, cur_y, cur_x)];
        if (adj_cell->perm == 0) {
            cur_cell->beta[k] = 0;
        }
        else {
            perm_eff = 2 * adj_cell->perm * cur_cell->perm;
            perm_eff /= (adj_cell->perm + cur_cell->perm);
            cur_cell->beta[k] = beta_coef * dim.h / perm_eff;
        }
    }
}

/* Computes A_alpha = xi/(1+beta_alpha*xi), xi = 2k/h */
void compute_A(cell_t *mesh, int cur_y, int cur_x)
{
    double xi;
    int k;
    cell_t *cur_cell;

    cur_cell = &mesh[MESH_INDEX(cur_y, cur_x)];
    xi = 2 * cur_cell->perm / dim.h;

    for (k = 0; k < 4; k++) {
        cur_cell->A[k] = xi / (1 + cur_cell->beta[k] * xi);
    }
}

/* Updates the pressure, flux, and robin conditions for a cell */
void update_interior(cell_t *mesh, cell_t *mesh_old, int cur_y, int cur_x)
{
    int k;
    cell_t *cur_cell, *cur_cell_old, *adj_cell;
    double sum_A, sum_A_R, A;

    cur_cell = &mesh[MESH_INDEX(cur_y, cur_x)];
    cur_cell_old = &mesh_old[MESH_INDEX(cur_y, cur_x)];

    sum_A = 0;
    sum_A_R = 0;

    for (k = 0; k < 4; k++) {
        adj_cell = &cur_cell_old[get_adjacent_index(k, cur_y, cur_x)];
        sum_A += cur_cell_old->A[k];
        sum_A_R += cur_cell_old->A[k] * adj_cell->robin[(k + 2) % 4];
    }

    /* Updates the pressure at the current cell on the new mesh */
    cur_cell->pressure = (cur_cell_old->source * dim.h + sum_A_R) / sum_A;

    /* Updates the flux at the current cell on the new mesh */
    for (k = 0; k < 4; k++) {
        adj_cell = &cur_cell_old[get_adjacent_index(k, cur_y, cur_x)];
        A = cur_cell->A[k];
        cur_cell->flux[k] = A * (cur_cell->pressure - adj_cell->robin[(k + 2) % 4]);
    }

    /* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        adj_cell = &mesh_old[get_adjacent_index(k, cur_y, cur_x)];
        cur_cell->l[k] = cur_cell->beta[k] * cur_cell->flux[k] + adj_cell->robin[(k + 2) % 4];
    }
}

void update_boundary(cell_t *mesh, cell_t *mesh_old, int cur_y, int cur_x, int boundary_side)
{
    int k;
    cell_t *cur_cell, *cur_cell_old, *adj_cell;
    double sum_A, sum_A_R, A;

    cur_cell = &mesh[MESH_INDEX(cur_y, cur_x)];
    cur_cell_old = &mesh_old[MESH_INDEX(cur_y, cur_x)];

    sum_A = 0;
    sum_A_R = 0;

    for (k = 0; k < 4; k++) {
        if (k != boundary_side) {
            adj_cell = &cur_cell_old[get_adjacent_index(k, cur_y, cur_x)];
            sum_A += cur_cell_old->A[k];
            sum_A_R += cur_cell_old->A[k] * adj_cell->robin[(k + 2) % 4];
        }
    }

    /* Updates the pressure at the current cell on the new mesh */
    cur_cell->pressure = (cur_cell_old->source * dim.h + sum_A_R) / sum_A;

    /* Updates the flux at the current cell on the new mesh */
    for (k = 0; k < 4; k++) {
        if (k != boundary_side) {
            adj_cell = &cur_cell_old[get_adjacent_index(k, cur_y, cur_x)];
            A = cur_cell->A[k];
            cur_cell->flux[k] = A * (cur_cell->pressure - adj_cell->robin[(k + 2) % 4]);
        }
    }

    /* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        if (k != boundary_side) {
            adj_cell = &mesh_old[get_adjacent_index(k, cur_y, cur_x)];
            cur_cell->l[k] = cur_cell->beta[k] * cur_cell->flux[k] +
                                adj_cell->robin[(k + 2) % 4];
        }
    }
}

void update_corner(cell_t *mesh, cell_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2)
{
    int k;
    cell_t *cur_cell, *cur_cell_old, *adj_cell;
    double sum_A, sum_A_R, A;

    cur_cell = &mesh[MESH_INDEX(cur_y, cur_x)];
    cur_cell_old = &mesh_old[MESH_INDEX(cur_y, cur_x)];

    sum_A = 0;
    sum_A_R = 0;

    for (k = 0; k < 4; k++) {
        if ((k != boundary_side1) && (k != boundary_side2)) {
            adj_cell = &cur_cell_old[get_adjacent_index(k, cur_y, cur_x)];
            sum_A += cur_cell_old->A[k];
            sum_A_R += cur_cell_old->A[k] * adj_cell->robin[(k + 2) % 4];
        }
    }

    /* Updates the pressure at the current cell on the new mesh */
    cur_cell->pressure = (cur_cell_old->source * dim.h + sum_A_R) / sum_A;

    /* Updates the flux at the current cell on the new mesh */
    for (k = 0; k < 4; k++) {
        if ((k != boundary_side1) && (k != boundary_side2)) {
            adj_cell = &cur_cell_old[get_adjacent_index(k, cur_y, cur_x)];
            A = cur_cell->A[k];
            cur_cell->flux[k] = A * (cur_cell->pressure - adj_cell->robin[(k + 2) % 4]);
        }
    }

    /* Updates the pressure at the edges of the current cell in the new mesh */
    for (k = 0; k < 4; k ++) {
        if ((k != boundary_side1) && (k != boundary_side2)) {
            adj_cell = &mesh_old[get_adjacent_index(k, cur_y, cur_x)];
            cur_cell->l[k] = cur_cell->beta[k] * cur_cell->flux[k] +
                                adj_cell->robin[(k + 2) % 4];
        }
    }
}
