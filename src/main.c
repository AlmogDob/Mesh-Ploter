#define SETUP
#define UPDATE
#define RENDER
#define ALMOG_DRAW_LIBRARY_IMPLEMENTATION
#include "./includes/Almog_Draw_Library.h"
#include "./includes/display.c"
#define MATRIX2D_IMPLEMENTATION
#include "./includes/Matrix2D.h"
#include "./includes/mesher.h"


Figure figure1;
void setup(game_state_t *game_state)
{
    game_state->const_fps = 30;
    // game_state->to_limit_fps = 0;

    /* creating mesh */
    printf("[INFO] meshing\n");
    double *x_2Dmat, *y_2Dmat;

    int mesh_rc = create_mesh(&x_2Dmat, &y_2Dmat, game_state->input_param.NACA, game_state->input_param.ni, game_state->input_param.nj, game_state->input_param.num_points_on_airfoil, game_state->input_param.delta_y, game_state->input_param.XSF, game_state->input_param.YSF, game_state->input_param.r, game_state->input_param.omega, game_state->input_param.output_dir);
    if (mesh_rc != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating mesh\n", __FILE__, __LINE__);
        exit(1);
    }

    figure1 = adl_alloc_figure(800, 800, (Point){300, 50, 0, 0});
    figure1.background_color = 0xFFFFFFFF;
    figure1.to_draw_axis = true;

    Curve points;
    ada_init_array(Point, points);

    Input_param input_param = game_state->input_param;
    for (int i = 0; i < input_param.ni; i++) {
        points.length = 0;
        for (int j = 0; j < input_param.nj; j++) {
            Point current_point = {0};

            current_point.x = x_2Dmat[offset2d_mesher(i, j, input_param.ni)];
            current_point.y = y_2Dmat[offset2d_mesher(i, j, input_param.ni)];

            ada_appand(Point, points, current_point);
        }
        adl_add_curve_to_figure(&figure1, points.elements, points.length, 0x0);
    }
    for (int j = 0; j < input_param.nj; j++) {
        points.length = 0;
        for (int i = 0; i < input_param.ni; i++) {
            Point current_point = {0};

            current_point.x = x_2Dmat[offset2d_mesher(i, j, input_param.ni)];
            current_point.y = y_2Dmat[offset2d_mesher(i, j, input_param.ni)];

            ada_appand(Point, points, current_point);
        }
        adl_add_curve_to_figure(&figure1, points.elements, points.length, 0x0);
    }



}

void update(game_state_t *game_state)
{
    figure1.offset_zoom_param = game_state->offset_zoom_param;
}

void render(game_state_t *game_state)
{
    adl_plot_curves_on_figure(figure1);
    adl_copy_figure_to_screen(game_state->window_pixels_mat, figure1);
    adl_draw_sentence(game_state->window_pixels_mat, game_state->input_param.NACA, strlen(game_state->input_param.NACA), 10, 10, 40, 0xFFFFFFFF, game_state->offset_zoom_param);
}

