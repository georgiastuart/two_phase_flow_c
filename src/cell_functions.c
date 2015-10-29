#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cell_functions.h"

/* Retrieves the cell adjacent to the current cell
    0 - left
    1 - right
    2 - up
    3 - down */
int get_adjacent_index(int direction, int cur_y, int cur_x)
{
    switch (direction) {
        case 0:
            return MESH_INDEX(cur_y, (cur_x - 1));
        case 1:
            return MESH_INDEX(cur_y, (cur_x + 1));
        case 2:
            return MESH_INDEX((cur_y - 1), cur_x);
        case 3:
            return MESH_INDEX((cur_y + 1), cur_x);
    }
    return 0;
}

/* Computes beta at each gridpoint on the computational domain */
void compute_beta(cell_t *mesh, double beta_coef)
{
    int i, j, k;
    double perm_eff, perm, beta;
    cell_t *cur_cell, *adj_cell;

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            cur_cell = &mesh[MESH_INDEX(i, j)];

            for (k = 0; k < 4; k++) {
                adj_cell = &mesh[get_adjacent_index(k, i, j)];
            }
        }
    }
}
