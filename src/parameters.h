#ifndef H_PARAMETERS
#define H_PARAMETERS

struct cell;
typedef struct cell cell_t;

struct global_mesh_params;
typedef struct global_mesh_params global_mesh_params_t;

/* Parameter calculations. Note: relative permeability functions are static */
double total_mobility(cell_t *cell, global_mesh_params_t *global);
double phase_mobility_o(cell_t *cell, global_mesh_params_t *global);
double phase_mobility_w(cell_t *cell, global_mesh_params_t *global);
double phase_mobility_w_deriv(cell_t *cell, global_mesh_params_t *global);
double cap_pressure_deriv(cell_t *cell, global_mesh_params_t *global);

#endif /* H_PARAMETERS */
