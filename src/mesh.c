#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "cell_functions.h"

/* Deep copies mesh 1 to mesh 2. Dimensions must be the same */
void mesh_copy(mesh_t *mesh_1, mesh_t *mesh_2)
{
    size_t size = sizeof(cell_t) * (mesh_1->dim.xdim + 2) * (mesh_1->dim.ydim + 2);
    memcpy(mesh_2->cell, mesh_1->cell, size);
}

/* Computes beta and A at all mesh points */
void mesh_compute_beta_A(mesh_t *mesh, const cell_ops_t *cell_ops)
{
    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cell_ops->cell_compute_beta(mesh, i, j, mesh->global.beta_coef);
            cell_ops->cell_compute_A(mesh, i, j);
        }
    }
}

/* Sets up the mesh with cell structs at each gridpoint */
mesh_t* mesh_init_mesh(dim_t dim, double *perm, double *source, double *sat, config_t *config)
{
    int i, j;
    mesh_t *mesh;
    cell_t *cur_cell;

    /* Allocates memory for the mesh */
    mesh = malloc(sizeof(mesh_t));
    mesh->cell = malloc((dim.ydim + 2) * (dim.xdim + 2) * sizeof(cell_t));

    /* Sets the mesh dimensions */
    mesh->dim = dim;

    /* Sets global parameters */
    mesh->global.porosity = config->porosity;
    mesh->global.sat_rel_o = config->sat_rel_o;
    mesh->global.sat_rel_w = config->sat_rel_w;
    mesh->global.visc_o = config->visc_o;
    mesh->global.visc_w = config->visc_w;
    mesh->global.eta = config->eta;
    mesh->global.beta_coef = config->beta_coef;

    /* Sets cell permeability and sources */
    for (i = 0; i < (mesh->dim.ydim + 2); i++) {
        for (j = 0; j < (mesh->dim.xdim + 2); j++) {
            cur_cell = &mesh->cell[MESH_INDEX_INC_PAD(i, j)];
            cur_cell->saturation = sat[MESH_INDEX_INC_PAD(i, j)];
            cur_cell->perm = config->perm_scale * exp(config->perm_strength
                    * perm[MESH_INDEX_INC_PAD(i, j)]);
            cur_cell->source = source[MESH_INDEX_INC_PAD(i, j)];
        }
    }

    /* Computes beta and A at all mesh points */
    mesh_compute_beta_A(mesh, &cell_press_ops);

    return mesh;
}

/* One iteration of the algorithm over all mesh points */
/* 0 - up, 1 - right,  2 - down, 3 - left */
void update_9(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Updates the corners */
    cell_ops->cell_update_corner(mesh, mesh_old, 0, 0, 0, 3);
    cell_ops->cell_update_corner(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    /* Updates the boundaries */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);
    }

    for (i = 1; i < (mesh->dim.xdim - 1); i++) {
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    /* Updates the interior cells */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < (mesh->dim.xdim - 1); j++) {
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Iteration for a type 0 block */
void update_0(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, 0, 0, 0, 3);

    for (i = 1; i < mesh->dim.xdim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 1 block */
void update_1(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.xdim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 2 block */
void update_2(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);

    for (i = 0; i < (mesh->dim.xdim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 3 block */
void update_3(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 4 block */
void update_4(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Iteration for a type 5 block */
void update_5(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 6 block */
void update_6(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);

    for (i = 1; i < mesh->dim.xdim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);

    for (i = 0; i < (mesh->dim.ydim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 7 block */
void update_7(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.xdim; i++) {
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 0; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 8 block */
void update_8(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    for (i = 0; i < (mesh->dim.xdim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);

    for (i = 0; i < (mesh->dim.ydim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

void mesh_update(mesh_t *mesh, mesh_t *mesh_old, int block_type, const cell_ops_t *cell_ops)
{
    switch (block_type) {
        case 9:
            update_9(mesh, mesh_old, cell_ops);
            break;
        case 0:
            update_0(mesh, mesh_old, cell_ops);
            break;
        case 1:
            update_1(mesh, mesh_old, cell_ops);
            break;
        case 2:
            update_2(mesh, mesh_old, cell_ops);
            break;
        case 3:
            update_3(mesh, mesh_old, cell_ops);
            break;
        case 4:
            update_4(mesh, mesh_old, cell_ops);
            break;
        case 5:
            update_5(mesh, mesh_old, cell_ops);
            break;
        case 6:
            update_6(mesh, mesh_old, cell_ops);
            break;
        case 7:
            update_7(mesh, mesh_old, cell_ops);
            break;
        case 8:
            update_8(mesh, mesh_old, cell_ops);
            break;
    }
}

/* Checks for convergence at a specified cutoff. Returns 1 if relative error */
/* is less than the convergence cutoff, 0 otherwise */
int mesh_press_convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank)
{
    int i, j;
    double num, denom, p_new, p_old, rel_error, global_num, global_denom;

    num = 0;
    denom = 0;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            p_new = mesh->cell[MESH_INDEX(i, j)].pressure;
            p_old = mesh_old->cell[MESH_INDEX(i, j)].pressure;
            num += pow(p_new - p_old, 2);
            denom += pow(p_new, 2);
        }
    }

    MPI_Reduce(&num, &global_num, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&denom, &global_denom, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        rel_error = sqrt(global_num / global_denom);
    }

    MPI_Bcast(&rel_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // printf("rel error %e\n", rel_error);

    if (rel_error < conv_cutoff) {
        return 1;
    }
    else {
        return 0;
    }
}

/* Checks for convergence at a specified cutoff. Returns 1 if relative error */
/* is less than the convergence cutoff, 0 otherwise */
int mesh_diff_convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank)
{
    double num, denom, s_new, s_old, rel_error, global_num, global_denom;

    num = 0;
    denom = 0;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            s_new = mesh->cell[MESH_INDEX(i, j)].saturation;
            s_old = mesh_old->cell[MESH_INDEX(i, j)].saturation;
            num += pow(s_new - s_old, 2);
            denom += pow(s_new, 2);
        }
    }

    MPI_Reduce(&num, &global_num, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&denom, &global_denom, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        rel_error = sqrt(global_num / global_denom);
    }

    MPI_Bcast(&rel_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rel_error < conv_cutoff) {
        return 1;
    }
    else {
        return 0;
    }
}

/* Ensures the average pressure is 0 */
void mesh_press_impose_0_average(mesh_t *mesh, int rank)
{
    int N, i, j, k;
    double sum, global_sum, avg;

    N = mesh->dim.x_full_dim * mesh->dim.y_full_dim;
    sum = 0;
    global_sum = 0;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            sum += mesh->cell[MESH_INDEX(i, j)].pressure;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        avg = global_sum / N;
    }

    MPI_Bcast(&avg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            mesh->cell[MESH_INDEX(i, j)].pressure -= avg;

            for (k = 0; k < 4; k++) {
                mesh->cell[MESH_INDEX(i, j)].l_p[k] -= avg;
            }
        }
    }
}

/* Updates the robin conditions along the boundaries */
void mesh_update_robin(mesh_t *mesh, const cell_ops_t *cell_ops)
{
    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cell_ops->cell_compute_robin(mesh, i, j);
        }
    }
}

/* Computes the velocity information on the mesh */
void mesh_compute_velocity(mesh_t *mesh)
{
    cell_t *cur_cell;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cur_cell = &mesh->cell[MESH_INDEX(i, j)];

            cur_cell->velocity_y = (cur_cell->flux_p[0] - cur_cell->flux_p[2]) / 2;
            cur_cell->velocity_x = (cur_cell->flux_p[1] - cur_cell->flux_p[3]) / 2;
        }
    }
}

/* Sets the velocity in output mesh equal to the velocity in input mesh */
void mesh_set_velocity(mesh_t *mesh, mesh_t *output_mesh)
{
    cell_t *in_cell, *out_cell;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            in_cell = &mesh->cell[MESH_INDEX(i, j)];
            out_cell = &output_mesh->cell[MESH_INDEX(i, j)];
            out_cell->velocity_x = in_cell->velocity_x;
            out_cell->velocity_y = in_cell->velocity_y;
        }
    }
}

/* For the pressure step in the flow problem */
int mesh_pressure_iteration(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff,
    int block_type, int rank, send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    int itr = 0;
    mesh_t *temp;

    mesh_compute_beta_A(mesh, &cell_press_ops);
    mesh_update_robin(mesh, &cell_press_ops);

    mesh_compute_beta_A(mesh_old, &cell_press_ops);
    mesh_update_robin(mesh_old, &cell_press_ops);

    mpi_comm(mesh_old, send_vec, rec_vec, block_type, rank, ROBIN_MODE);
    mpi_comm(mesh, send_vec, rec_vec, block_type, rank, ROBIN_MODE);

    for (;;) {
        itr++;

        mesh_update(mesh, mesh_old, block_type, &cell_press_ops);

        if (mesh_press_convergence_check(mesh, mesh_old, conv_cutoff, rank)) {
            break;
        }

        // if (rank == 0) {
        //     print_attribute(mesh, "saturation");
        // }

        // break;
        mesh_press_impose_0_average(mesh, rank);
        mesh_update_robin(mesh, &cell_press_ops);
        mpi_comm(mesh, send_vec, rec_vec, block_type, rank, ROBIN_MODE);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;
    }


    mesh_compute_velocity(mesh);
    mesh_set_velocity(mesh, mesh_old);

    mesh_copy(mesh, mesh_old);

    return itr;
}

/* Computes diffusion D at all mesh points */
void mesh_compute_diffusion_coef_source(mesh_t *mesh)
{
    cell_t *cur_cell;
    for (int i = 0; i < (mesh->dim.ydim + 2); i++) {
        for (int j = 0; j < (mesh->dim.xdim + 2); j++) {
            cur_cell = &mesh->cell[MESH_INDEX_INC_PAD(i, j)];
            diff_compute_diffusion(mesh, cur_cell);
            diff_compute_source(mesh, cur_cell);
        }
    }
}

/* Only for the diffusion test */
void diffusion_test_update(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    /* Updates the corners */
    diff_update_corner_dirichlet(mesh, mesh_old, 0, 0, 0, 3);
    diff_update_corner_dirichlet(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);
    diff_update_corner_dirichlet(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);
    diff_update_corner_dirichlet(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    /* Updates the boundaries */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        diff_update_boundary_dirichlet(mesh, mesh_old, i, 0, 3);
        diff_update_boundary_dirichlet(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);
    }

    for (i = 1; i < (mesh->dim.xdim - 1); i++) {
        diff_update_boundary(mesh, mesh_old, 0, i, 0);
        diff_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    /* Updates the interior cells */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < (mesh->dim.xdim - 1); j++) {
            diff_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* For diffusion step in the flow problem */
int mesh_diffusion_iteration(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff,
    int block_type, int rank, send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    int itr = 0;
    mesh_t *temp;

    mesh_compute_diffusion_coef_source(mesh);
    mesh_compute_diffusion_coef_source(mesh_old);

    mesh_compute_beta_A(mesh, &cell_diff_ops);
    mesh_update_robin(mesh, &cell_diff_ops);

    mesh_compute_beta_A(mesh_old, &cell_diff_ops);
    mesh_update_robin(mesh_old, &cell_diff_ops);

    mpi_comm(mesh, send_vec, rec_vec, block_type, rank, ROBIN_MODE);
    mpi_comm(mesh_old, send_vec, rec_vec, block_type, rank, ROBIN_MODE);

    for (;;) {
        itr++;

        // diffusion_test_update(mesh, mesh_old);
        mesh_update(mesh, mesh_old, block_type, &cell_diff_ops);

        // break;
        if (mesh_diff_convergence_check(mesh, mesh_old, conv_cutoff, rank)) {
            break;
        }

        if (!(itr % 100)) {
            printf("Iteration: %d\n", itr);
        }

        mesh_update_robin(mesh, &cell_diff_ops);
        mpi_comm(mesh, send_vec, rec_vec, block_type, rank, ROBIN_MODE);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;
    }

    mesh_copy(mesh, mesh_old);

    return itr;
}

/* Sets saturation_prev to saturation */
void mesh_update_saturation_time(mesh_t *mesh)
{
    cell_t *cur_cell;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cur_cell = &mesh->cell[MESH_INDEX(i, j)];
            cur_cell->saturation_prev = cur_cell->saturation;
        }
    }
}

/* Finds the maxiumum time step for transport */
void mesh_max_time_step(mesh_t *mesh, mesh_t *mesh_old)
{
    double max = 0, mult = 0, vel_norm = 0, global_max = 0;
    cell_t *cur_cell;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cur_cell = &mesh_old->cell[MESH_INDEX(i, j)];
            vel_norm = sqrt(pow(cur_cell->velocity_x, 2) + pow(cur_cell->velocity_y, 2));
            mult = cur_cell->pm_w_deriv / mesh->global.porosity * vel_norm;
            if (mult > max)
                max = mult;
        }
    }

    MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    mesh->dim.dt_transport = (0.95 * mesh->dim.h / 2) / global_max;
    mesh_old->dim.dt_transport = (0.95 * mesh->dim.h / 2) / global_max;
}

/* Calculates phase mobility derivative of the meshes using information from mesh old */
static void set_phase_mob_deriv(mesh_t *mesh, mesh_t *mesh_old)
{
    double lambda_w_deriv;
    cell_t *cur_cell;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cur_cell = &mesh_old->cell[MESH_INDEX(i, j)];
            lambda_w_deriv = phase_mobility_w_deriv(cur_cell, &mesh_old->global);
            cur_cell->pm_w_deriv = lambda_w_deriv;
            mesh->cell[MESH_INDEX(i, j)].pm_w_deriv = lambda_w_deriv;
        }
    }
}

/* Time step for the transport problem */
int mesh_transport_iteration(mesh_t *mesh, mesh_t *mesh_old, int block_type, int rank,
            send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    mesh_t *temp;

    /* Computes phase mobility of water derivative*/
    set_phase_mob_deriv(mesh, mesh_old);

    /* Computes dt_transport for both meshes */
    mesh_max_time_step(mesh, mesh_old);
    double dtt = mesh->dim.dt_transport;
    int num_ts = (int) (mesh->dim.dt / dtt) + 1;
    double remainder_ts = mesh->dim.dt - ((double) (num_ts - 1) * dtt);

    // for (int i = 0; i < (num_ts - 1); i++) {
    for (int i = 0; i < num_ts; i++) {
        mesh_update(mesh, mesh_old, block_type, &cell_trans_ops);
        mpi_comm(mesh, send_vec, rec_vec, block_type, rank, SAT_MODE);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;
    }

    mesh->dim.dt_transport = remainder_ts;
    mesh_old->dim.dt_transport = remainder_ts;

    mesh_update(mesh, mesh_old, block_type, &cell_trans_ops);

    /* sets current s to s_prev */
    mesh_update_saturation_time(mesh);

    mesh_copy(mesh, mesh_old);

    return num_ts;
}

void setup_diffusion_test(mesh_t *mesh)
{
    cell_t *cur_cell;
    double delta_x = mesh->dim.xlen / mesh->dim.xdim;
    double x;

    mesh->global.porosity = 1;

    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            x = delta_x * j;
            cur_cell = &mesh->cell[MESH_INDEX(i, j)];
            cur_cell->saturation = 0;
            cur_cell->saturation_prev = sin(M_PI * x) + 0.5 * sin(3.0 * M_PI * x);
            cur_cell->diffusion = 1;
        }
    }
}

void setup_transport_test(mesh_t *mesh)
{
    cell_t *cur_cell;

    cur_cell = &mesh->cell[MESH_INDEX((mesh->dim.ydim - 1), 0)];
    cur_cell->saturation = 0.84;
}
