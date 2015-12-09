#include <stdlib.h>
#include <math.h>
#include "parameters.h"
#include "cell_functions.h"
#include "mesh.h"

/* Relative permeability of the oil phase, K_ro */
static double rel_perm_o_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    return pow(1.0 - (cell->saturation / (1.0 - global->sat_rel_o)), 2);
}

/* Relative permeability of the oil phase for linear transport */
static double rel_perm_o_linear(cell_t *cell, global_mesh_params_t *global)
{
    return 1.0 - cell->saturation;
}

/* Derivative of the relative permeability of oil, K_ro' */
static double rel_perm_o_deriv_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double denom = (1.0 - global->sat_rel_o);
    return -2.0 / denom * (1.0 - cell->saturation / denom);
}

/* Derivative of the relative permeability of oil for linear transport */
static double rel_perm_o_deriv_linear(cell_t *cell, global_mesh_params_t *global)
{
    return -1.0;
}

/* Relative permeability of the water phase, K_rw */
static double rel_perm_w_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double num = pow(cell->saturation - global->sat_rel_w, 2);
    double denom = pow(1.0 - global->sat_rel_w, 2);

    return pow(num / denom, 2);
}

/* Relative permeability of the water phase for linear transport */
static double rel_perm_w_linear(cell_t *cell, global_mesh_params_t *global)
{
    return cell->saturation;
}

/* Derivative of the relative permeability of water, K_rw' */
static double rel_perm_w_deriv_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double num = 4.0 * pow(cell->saturation - global->sat_rel_w, 3);
    double denom = pow(1.0 - global->sat_rel_w, 4);

    return num / denom;
}

/* Derivative of the relative permeability of water for linear transport */
static double rel_perm_w_deriv_linear(cell_t *cell, global_mesh_params_t *global)
{
    return 1.0;
}

/* Pointers for linear and nonlinear parameters */
const linearity_ops_t param_nonlinear_ops = {
    .rel_perm_o         = &rel_perm_o_nonlinear,
    .rel_perm_o_deriv   = &rel_perm_o_deriv_nonlinear,
    .rel_perm_w         = &rel_perm_w_nonlinear,
    .rel_perm_w_deriv   = &rel_perm_w_deriv_nonlinear
};

const linearity_ops_t param_linear_ops = {
    .rel_perm_o         = &rel_perm_o_linear,
    .rel_perm_o_deriv   = &rel_perm_o_deriv_linear,
    .rel_perm_w         = &rel_perm_w_linear,
    .rel_perm_w_deriv   = &rel_perm_w_deriv_linear
};

static const linearity_ops_t *lin_ops;

/* Linearity is 1 if linear, 0 if nonlinear */
void set_linearity(int linearity)
{
    if (linearity)
        lin_ops = &param_linear_ops;
    else
        lin_ops = &param_nonlinear_ops;
}

/* Total mobility, lambda */
double total_mobility(cell_t *cell, global_mesh_params_t *global)
{
    return (lin_ops->rel_perm_o(cell, global) / global->visc_o +
           lin_ops->rel_perm_w(cell, global) / global->visc_w);
}

/* Derivative of the total mobility, lambda' */
static double total_mobility_deriv(cell_t *cell, global_mesh_params_t *global)
{
    return (lin_ops->rel_perm_o_deriv(cell, global) / global->visc_o +
            lin_ops->rel_perm_w_deriv(cell, global) / global->visc_w);
}

/* Phase mobility of the oil phase, lambda_o */
double phase_mobility_o(cell_t *cell, global_mesh_params_t *global)
{
    double total_mob = total_mobility(cell, global);
    double rel_perm = lin_ops->rel_perm_o(cell, global);

    return rel_perm / (global->visc_o * total_mob);
}

/* Phase mobility of the water phase, lambda_w */
double phase_mobility_w(cell_t *cell, global_mesh_params_t *global)
{
    double total_mob = total_mobility(cell, global);
    double rel_perm = lin_ops->rel_perm_w(cell, global);

    return rel_perm / (global->visc_w * total_mob);
}

/* Derivative of phase mobility of water, lambda_w' */
double phase_mobility_w_deriv(cell_t *cell, global_mesh_params_t *global)
{
    double visc = global->visc_w;
    double total_mob = total_mobility(cell, global);
    double num = visc * total_mob * lin_ops->rel_perm_w_deriv(cell, global);
    num -= lin_ops->rel_perm_w(cell, global) * visc * total_mobility_deriv(cell, global);
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
