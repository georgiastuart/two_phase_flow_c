#include <stdlib.h>
#include <math.h>
#include "parameters.h"
#include "cell_functions.h"
#include "mesh.h"

/* Relative permeability of the oil phase, K_ro */
static double rel_perm_o(cell_t *cell, global_mesh_params_t *global)
{
    return pow(1 - (cell->saturation / (1 - global->sat_rel_o)), 2);
}

/* Derivative of the relative permeability of oil, K_ro' */
static double rel_perm_o_deriv(cell_t *cell, global_mesh_params_t *global)
{
    double denom = (1 - global->sat_rel_o);
    return -2.0 / denom * (1 - cell->saturation / denom);
}

/* Relative permeability of the water phase, K_rw */
static double rel_perm_w(cell_t *cell, global_mesh_params_t *global)
{
    double num = pow(cell->saturation - global->sat_rel_w, 2);
    double denom = pow(1 - global->sat_rel_w, 2);

    return pow(num / denom, 2);
}

/* Derivative of the relative permeability of water, K_rw' */
static double rel_perm_w_deriv(cell_t *cell, global_mesh_params_t *global)
{
    double num = 4 * pow(cell->saturation - global->sat_rel_w, 3);
    double denom = pow(1 - global->sat_rel_w, 4);

    return num / denom;
}

/* Total mobility, lambda */
double total_mobility(cell_t *cell, global_mesh_params_t *global)
{
    return 1;
    //return (rel_perm_o(cell, global) / global->visc_o +
    //        rel_perm_w(cell, global) / global->visc_w);
}

/* Derivative of the total mobility, lambda' */
static double total_mobility_deriv(cell_t *cell, global_mesh_params_t *global)
{
    return (rel_perm_o_deriv(cell, global) / global->visc_o +
            rel_perm_w_deriv(cell, global) / global->visc_w);
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

/* Derivative of phase mobility of water, lambda_w' */
double phase_mobility_w_deriv(cell_t *cell, global_mesh_params_t *global)
{
    double visc = global->visc_w;
    double total_mob = total_mobility(cell, global);
    double num = visc * total_mob * rel_perm_w_deriv(cell, global);
    num -= rel_perm_w(cell, global) * visc * total_mobility_deriv(cell, global);
    double denom = pow(visc * total_mob, 2);

    return num / denom;
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
