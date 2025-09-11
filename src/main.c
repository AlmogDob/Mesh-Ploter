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
void setup(game_state_t *game_state)
{
    game_state->const_fps = 30;
    game_state->to_limit_fps = 0;

    Input_param input_param = game_state->input_param;

    /* allocating the matrices */
    int i_index, j_index;
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

    /* solving flow */
    printf("[INFO] solving flow field\n");
    double final_S_norm;

    int solver_rc = solver(input_param.output_dir, &rho_2Dmat, &u_2Dmat, &v_2Dmat, &e_2Dmat, x_2Dmat, y_2Dmat, &final_S_norm, input_param.ni, input_param.nj, input_param.num_points_on_airfoil, input_param.Mach_inf, input_param.angle_of_attack_deg, input_param.density, input_param.environment_pressure, input_param.delta_t, input_param.Gamma, input_param.epse, input_param.max_iteration);
    if (solver_rc != 0 || isnan(final_S_norm)) {
        fprintf(stderr, "%s:%d: [ERROR] unable to solve the flow\n", __FILE__, __LINE__);
        exit(1);
    }
    
    for (i_index = 0; i_index < input_param.ni; i_index++) {
        for (j_index = 0; j_index < input_param.nj; j_index++) {
            double e = e_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];
            double rho = rho_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];
            double u = u_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];
            double v = v_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)];

            double p = (input_param.Gamma - 1) * (e - 0.5 * rho * (u * u + v * v));
            double a = sqrt(input_param.Gamma * p / rho);
            M_2Dmat[adl_offset2d(i_index, j_index, input_param.ni)] = sqrt(u * u + v * v) / a;
        }
    }

    figure1 = adl_alloc_figure(800, 800, (Point){300, 50, 0, 0});
    figure1.background_color = 0xFFFFFFFF;
    figure1.to_draw_axis = true;
    figure1.to_draw_max_min_values = true;

}

void update(game_state_t *game_state)
{
    figure1.offset_zoom_param = game_state->offset_zoom_param;
}

void render(game_state_t *game_state)
{
    // adl_plot_curves_on_figure(figure1);
    adl_interp_scalar_2D_on_figure(figure1, x_2Dmat, y_2Dmat, M_2Dmat, game_state->input_param.ni, game_state->input_param.nj, "g-r");

    adl_copy_figure_to_screen(game_state->window_pixels_mat, figure1);
    adl_draw_sentence(game_state->window_pixels_mat, game_state->input_param.NACA, strlen(game_state->input_param.NACA), 10, 10, 40, 0xFFFFFFFF, ADL_DEFAULT_OFFSET_ZOOM);
}

