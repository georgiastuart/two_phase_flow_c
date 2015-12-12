#include <stdlib.h>
#include <math.h>
#include "parameters.h"
#include "cell_functions.h"
#include "mesh.h"

/* Selects the saturation level to guard against saturations that are too high */
static double select_sat(cell_t *cell, global_mesh_params_t *global)
{
    // if (cell->saturation < (global->sat_rel_w + 0.001)) {
    //     return global->sat_rel_w + 0.001;
    // } else if (cell->saturation > (1.0 - (global->sat_rel_o + 0.001))) {
    //     return 1.0 - (global->sat_rel_o + 0.001);
    // } else {
    //     return cell->saturation;
    // }

    return cell->saturation;
}

/* Relative permeability of the oil phase, K_ro */
static double rel_perm_o_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double sat = select_sat(cell, global);
    return pow(1.0 - (sat / (1.0 - global->sat_rel_o)), 2);
}

/* Relative permeability of the oil phase for linear transport */
static double rel_perm_o_linear(cell_t *cell, global_mesh_params_t *global)
{
    double sat = select_sat(cell, global);
    return 1.0 - sat;
}

/* Derivative of the relative permeability of oil, K_ro' */
static double rel_perm_o_deriv_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double sat = select_sat(cell, global);
    double denom = (1.0 - global->sat_rel_o);
    return -2.0 / denom * (1.0 - sat / denom);
}

/* Derivative of the relative permeability of oil for linear transport */
static double rel_perm_o_deriv_linear(cell_t *cell, global_mesh_params_t *global)
{
    return -1.0;
}

/* Relative permeability of the water phase, K_rw */
static double rel_perm_w_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double sat = select_sat(cell, global);
    double num = pow(sat - global->sat_rel_w, 2);
    double denom = pow(1.0 - global->sat_rel_w, 2);

    return num / denom;
}

/* Relative permeability of the water phase for linear transport */
static double rel_perm_w_linear(cell_t *cell, global_mesh_params_t *global)
{
    return select_sat(cell, global);
}

/* Derivative of the relative permeability of water, K_rw' */
static double rel_perm_w_deriv_nonlinear(cell_t *cell, global_mesh_params_t *global)
{
    double sat = select_sat(cell, global);
    double num = 2.0 * (sat - global->sat_rel_w);
    double denom = pow(1.0 - global->sat_rel_w, 2);
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
    // printf("rel_perm_o: %e\n", lin_ops->rel_perm_o(cell, global));
    // printf("rel_perm_w: %e\n", lin_ops->rel_perm_w(cell, global));
    // printf("visc o %e\n", global->visc_o);
    // printf("visc w %e\n", global->visc_w);
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
    // printf("rel perm: %e\n", rel_perm);
    // printf("total mob: %e\n", total_mob);
    // printf("phase mob inside %e\n", rel_perm / (global->visc_w * total_mob));

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
    double sat = select_sat(cell, global);

    s_o = 1.0 - global->sat_rel_o;
    z = pow(global->sat_rel_o, 2) * pow(s_o - global->sat_rel_w, -2);

    pc_deriv = pow(sat - global->sat_rel_w, -3);
    pc_deriv += z * pow(1.0 - sat, -3);
    pc_deriv *= -2.0 * global->eta;

    return pc_deriv;
}
