#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_functions.h"
#include "cell_functions.h"
#include "mesh.h"

double rel_perm_o(cell_t *cell, global_mesh_params_t *global)
{
    return pow(1 - (cell->saturation/(1-global->sat_rel_o)), 2);
}

double rel_perm_w(cell_t *cell, global_mesh_params_t *global)
{
    double num, denom;

    num = pow(cell->saturation - global->sat_rel_w, 2);
    denom = pow(1 - global->sat_rel_w, 2);

    return pow(num / denom, 2);
}

double total_mobility(cell_t *cell, global_mesh_params_t *global)
{
    return (rel_perm_o(cell, global) / global->visc_o +
                rel_perm_w(cell, global) / global->visc_w);
}

double phase_mobility_o(cell_t *cell, global_mesh_params_t *global)
{
    double total_mob = total_mobility(cell, global);
    double rel_perm = rel_perm_o(cell, global);

    return rel_perm / (global->visc_o * total_mob);
}

double phase_mobility_w(cell_t *cell, global_mesh_params_t *global)
{
    double total_mob = total_mobility(cell, global);
    double rel_perm = rel_perm_w(cell, global);

    return rel_perm / (global->visc_w * total_mob);
}
