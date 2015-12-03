#ifndef H_PARAMETERS
#define H_PARAMETERS

struct cell;
typedef struct cell cell_t;

struct global_mesh_params;
typedef struct global_mesh_params global_mesh_params_t;

typedef struct linearity_ops
{
    double (*rel_perm_o)(cell_t*, global_mesh_params_t*);
    double (*rel_perm_w)(cell_t*, global_mesh_params_t*);
    double (*rel_perm_o_deriv)(cell_t*, global_mesh_params_t*);
    double (*rel_perm_w_deriv)(cell_t*, global_mesh_params_t*);
} linearity_ops_t;

/* Parameter calculations. Note: relative permeability functions are static */
void set_linearity(int linearity);
double total_mobility(cell_t *cell, global_mesh_params_t *global);
double phase_mobility_o(cell_t *cell, global_mesh_params_t *global);
double phase_mobility_w(cell_t *cell, global_mesh_params_t *global);
double phase_mobility_w_deriv(cell_t *cell, global_mesh_params_t *global);
double cap_pressure_deriv(cell_t *cell, global_mesh_params_t *global);

#endif /* H_PARAMETERS */
