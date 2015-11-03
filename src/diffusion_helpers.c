#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_helpers.h"
#include "cell_functions.h"
#include "mesh.h"

static double rel_perm_o(cell_t *cell, global_mesh_params_t *global)
{
    return pow(1 - (cell->saturation/(1-global->sat_rel_o)), 2);
}

static double rel_perm_w(cell_t *cell, global_mesh_params_t *global)
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
