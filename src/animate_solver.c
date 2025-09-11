#define SETUP
#define UPDATE
#define RENDER
#define ALMOG_DRAW_LIBRARY_IMPLEMENTATION
#include "./includes/Almog_Draw_Library.h"
#include "./includes/display.c"
#define MATRIX2D_IMPLEMENTATION
#include "./includes/Matrix2D.h"
#include "./includes/mesher.h"
#include "./includes/solver.h"


Figure figure1;
double *x_2Dmat, *y_2Dmat;
double *rho_2Dmat, *u_2Dmat, *v_2Dmat, *e_2Dmat, *M_2Dmat;
double *J_vals_mat, *first_Q, *current_Q, *next_Q, *S, *W, *dxi_dx_mat, *dxi_dy_mat, *deta_dx_mat, *deta_dy_mat, *s2, *rspec, *qv, *dd, *U_mat, *V_mat, *A, *B, *C, *D, *drr, *drp, max_S_norm = 0, current_S_norm, first_S_norm;

int i_LE, i_TEL, i_TEU, j_TEL, j_LE, j_TEU;
double angle_of_attack_rad, epsi;
int max_ni_nj;



bool stop_solving;
int tot_iter;
char iter_counter[256];

void setup(game_state_t *game_state)
{
    game_state->const_fps = 30;
    game_state->to_limit_fps = 0;

    Input_param input_param = game_state->input_param;

    /* allocating the matrices */
    int i_index, j_index, k_index;
    int ni = input_param.ni;
    int nj = input_param.nj;

    rho_2Dmat = (double *)malloc(sizeof(double) * input_param.ni * input_param.nj);
    for (i_index = 0; i_index < input_param.ni; i_index++) {
        for (j_index = 0; j_index < input_param.nj; j_index++) {
            rho_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = 0;
        }
    }
    u_2Dmat = (double *)malloc(sizeof(double) * input_param.ni * input_param.nj);
    for (i_index = 0; i_index < input_param.ni; i_index++) {
        for (j_index = 0; j_index < input_param.nj; j_index++) {
            u_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = 0;
        }
    }
    v_2Dmat = (double *)malloc(sizeof(double) * input_param.ni * input_param.nj);
    for (i_index = 0; i_index < input_param.ni; i_index++) {
        for (j_index = 0; j_index < input_param.nj; j_index++) {
            v_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = 0;
        }
    }
    e_2Dmat = (double *)malloc(sizeof(double) * input_param.ni * input_param.nj);
    for (i_index = 0; i_index < input_param.ni; i_index++) {
        for (j_index = 0; j_index < input_param.nj; j_index++) {
            e_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = 0;
        }
    }
    M_2Dmat = (double *)malloc(sizeof(double) * input_param.ni * input_param.nj);
    for (i_index = 0; i_index < input_param.ni; i_index++) {
        for (j_index = 0; j_index < input_param.nj; j_index++) {
            M_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = 0;
        }
    }

    /* creating mesh */
    printf("[INFO] meshing\n");

    int mesh_rc = create_mesh(&x_2Dmat, &y_2Dmat, input_param.NACA, input_param.ni, input_param.nj, input_param.num_points_on_airfoil, input_param.delta_y, input_param.XSF, input_param.YSF, input_param.r, input_param.omega, input_param.output_dir);
    if (mesh_rc != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating mesh\n", __FILE__, __LINE__);
        exit(1);
    }


    i_LE  = (ni-1) / 2;
    i_TEL = i_LE - input_param.num_points_on_airfoil / 2;
    i_TEU = i_LE + input_param.num_points_on_airfoil / 2;
    j_TEL = 0;
    j_LE  = 0;
    j_TEU = 0;
    angle_of_attack_rad = 2 * PI / 180 * input_param.angle_of_attack_deg;
    epsi   = input_param.epse * 2;
    max_ni_nj = (int)fmax(ni, nj);

    /* allocating the matrices */
    J_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            J_vals_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    dxi_dx_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            dxi_dx_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    dxi_dy_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            dxi_dy_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    deta_dx_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            deta_dx_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    deta_dy_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            deta_dy_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    U_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            U_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    V_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            V_mat[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    s2 = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        s2[i_index] = 0;
    }
    rspec = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        rspec[i_index] = 0;
    }
    qv = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        qv[i_index] = 0;
    }
    dd = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        dd[i_index] = 0;
    }
    drr = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        drr[i_index] = 0;
    }
    drp = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        drp[i_index] = 0;
    }
    W = (double *)malloc(sizeof(double) * max_ni_nj * 4);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            W[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    D = (double *)malloc(sizeof(double) * max_ni_nj * 4);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            D[offset2d_solver(i_index, j_index, ni, nj)] = 0;
        }
    }
    A = (double *)malloc(sizeof(double) * 4 * 4 * max_ni_nj);
    for (i_index = 0; i_index < 4; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            for (k_index = 0; k_index < max_ni_nj; k_index++) {
                A[offset2d_solver(i_index, j_index, ni, nj)] = 0;
            }
        }
    }
    B = (double *)malloc(sizeof(double) * 4 * 4 * max_ni_nj);
    for (i_index = 0; i_index < 4; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            for (k_index = 0; k_index < max_ni_nj; k_index++) {
                B[offset2d_solver(i_index, j_index, ni, nj)] = 0;
            }
        }
    }
    C = (double *)malloc(sizeof(double) * 4 * 4 * max_ni_nj);
    for (i_index = 0; i_index < 4; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            for (k_index = 0; k_index < max_ni_nj; k_index++) {
                C[offset2d_solver(i_index, j_index, ni, nj)] = 0;
            }
        }
    }
    first_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                first_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    current_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                current_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    next_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                next_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    S = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                S[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }

    /*------------------------------------------------------------*/

    initialize_solver(current_Q, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat, deta_dy_mat, x_2Dmat, y_2Dmat, ni, nj, input_param.Mach_inf, angle_of_attack_rad, input_param.environment_pressure, input_param.density, input_param.Gamma);
    copy_3Dmat_to_3Dmat(first_Q, current_Q, ni, nj);

    stop_solving = false;
    tot_iter = 0;

    
    figure1 = adl_alloc_figure(1200, 1200, (Point){300, 50, 0, 0});
    figure1.background_color = 0xFFFFFFFF;
    figure1.to_draw_axis = true;
    figure1.to_draw_max_min_values = true;

}

void update(game_state_t *game_state)
{
    figure1.offset_zoom_param = game_state->offset_zoom_param;
    Input_param input_param = game_state->input_param;
    int ni = input_param.ni;
    int nj = input_param.nj;
    
    if (!stop_solving) {
        // for (int iteration = 0; iteration < input_param.max_iteration; iteration++) {
        for (int iteration = 0; iteration < 200; iteration++) {
            tot_iter++;
            apply_BC(current_Q, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat, deta_dy_mat, ni, nj, i_TEL, i_LE, i_TEU, input_param.Gamma);

            double delta_t_fixed = input_param.delta_t;
            // double delta_t_fixed = 5 * fix_delta_t(current_Q, x_vals_mat, y_vals_mat, ni, nj);
            // dprintD(delta_t_fixed);
            // return 0;

            current_S_norm = step_solver(A, B, C, D, current_Q, S, W, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat, deta_dy_mat, s2, drr, drp, rspec, qv, dd, ni, nj, max_ni_nj, input_param.Mach_inf, delta_t_fixed, input_param.Gamma, input_param.epse, epsi);
            if (max_S_norm < fabs(current_S_norm)) {
                max_S_norm = fabs(current_S_norm);
            }
            if (iteration == 0) {
                first_S_norm = current_S_norm;
            }
            advance_Q(next_Q, current_Q, S, J_vals_mat, ni, nj);
            copy_3Dmat_to_3Dmat(current_Q, next_Q, ni, nj);
            
            if (!(tot_iter % 100)) {
                printf("\r%5d: %13.10f", tot_iter, current_S_norm);
                fflush(stdout);
                sprintf(iter_counter, "%5d: %13.10f", tot_iter, current_S_norm);
            }

            if (fabs(current_S_norm) / first_S_norm < 1e-5 || current_S_norm == 0 || isnan(current_S_norm)) {
                printf("\r%5d: %13.10f\n", tot_iter, current_S_norm);
                fflush(stdout);
                stop_solving = true;
            }
        }
        /* preparing matrixes for export */
        for (int i = 0; i < ni; i++) {
            for (int j = 0; j < nj; j++) {
                rho_2Dmat[offset2d_solver(i, j, ni, nj)]  = current_Q[offset3d(i, j, 0, ni, nj)];
                u_2Dmat[offset2d_solver(i, j, ni, nj)]  = current_Q[offset3d(i, j, 1, ni, nj)] / current_Q[offset3d(i, j, 0, ni, nj)];
                v_2Dmat[offset2d_solver(i, j, ni, nj)] = current_Q[offset3d(i, j, 2, ni, nj)] / current_Q[offset3d(i, j, 0, ni, nj)];
                e_2Dmat[offset2d_solver(i, j, ni, nj)] = current_Q[offset3d(i, j, 3, ni, nj)];
            }
        }

        for (int i_index = 0; i_index < input_param.ni; i_index++) {
            for (int j_index = 0; j_index < input_param.nj; j_index++) {
                double e = e_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];
                double rho = rho_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];
                double u = u_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];
                double v = v_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];

                double p = (input_param.Gamma - 1) * (e - 0.5 * rho * (u * u + v * v));
                double a = sqrt(input_param.Gamma * p / rho);
                M_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = sqrt(u * u + v * v) / a;
            }
        }
    }


}

void render(game_state_t *game_state)
{
    adl_interp_scalar_2D_on_figure(figure1, x_2Dmat, y_2Dmat, M_2Dmat, game_state->input_param.ni, game_state->input_param.nj, "b-r");
    adl_copy_figure_to_screen(game_state->window_pixels_mat, figure1);

    adl_draw_sentence(game_state->window_pixels_mat, game_state->input_param.NACA, strlen(game_state->input_param.NACA), 10, 10, 40, 0xFFFFFFFF, ADL_DEFAULT_OFFSET_ZOOM);
    adl_draw_sentence(game_state->window_pixels_mat, iter_counter, strlen(iter_counter), 1550, 10, 40, 0xFFFFFFFF, ADL_DEFAULT_OFFSET_ZOOM);
}

