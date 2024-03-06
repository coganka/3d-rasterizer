#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include "upng.h"
#include "array.h"
#include "display.h"
#include "clipping.h"
#include "vector.h"
#include "mesh.h"
#include "matrix.h"
#include "light.h"
#include "camera.h"
#include "triangle.h"
#include "texture.h"


#define MAX_TRIANGLES_PER_MESH 10000
triangle_t triangles_to_render[MAX_TRIANGLES_PER_MESH];
int num_triangles_to_render = 0;

bool is_running = false;
int previous_frame_time = 0;
float delta_time = 0;

mat4_t world_matrix;
mat4_t proj_matrix;
mat4_t view_matrix;

void setup(void) {
    // Initialize render mode and triangle culling method
    set_render_method(RENDER_WIRE);
    set_cull_method(RENDER_WIRE);

    // init light direction
    init_light(vec3_new(0, 0, 1));

    // initialize the perspective projection matrix
    float aspecty = (float)get_window_height() / (float)get_window_width();
    float aspectx = (float)get_window_width() / (float)get_window_height();
    float fovy = M_PI / 3.0; // 60 degrees -> 1.047197
    float fovx = atan(tan(fovy / 2) * aspectx) * 2.0;
    float znear = 0.1;
    float zfar = 100.0;
    proj_matrix = mat4_make_perspective(fovy, aspecty, znear, zfar);

    // initialize frustum planes with a point and a normal
    init_frustum_planes(fovx, fovy, znear, zfar);

    load_mesh("./assets/runway.obj", "./assets/runway.png", vec3_new(1, 1, 1), vec3_new(0, -1.5, +23), vec3_new(0, 0, 0));
    load_mesh("./assets/f22.obj", "./assets/f22.png", vec3_new(1, 1, 1), vec3_new(0, -1.3, +5), vec3_new(0, -M_PI/2, 0));
    load_mesh("./assets/efa.obj", "./assets/efa.png", vec3_new(1, 1, 1), vec3_new(-2, -1.3, +9), vec3_new(0, -M_PI/2, 0));
    load_mesh("./assets/f117.obj", "./assets/f117.png", vec3_new(1, 1, 1), vec3_new(+2, -1.3, +9), vec3_new(0, -M_PI/2, 0));
}


void process_input(void) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
            case SDL_QUIT:
                is_running = false;
                break;
            case SDL_KEYDOWN:
                if (event.key.keysym.sym == SDLK_ESCAPE) {
                    is_running = false;
                    break;
                }
                if (event.key.keysym.sym == SDLK_1) {
                    set_render_method(RENDER_WIRE_VERTEX);
                    break;
                }
                if (event.key.keysym.sym == SDLK_2) {
                    set_render_method(RENDER_WIRE);
                    break;
                }
                if (event.key.keysym.sym == SDLK_3) {
                    set_render_method(RENDER_FILL_TRIANGLE);
                    break;
                }
                if (event.key.keysym.sym == SDLK_4) {
                    set_render_method(RENDER_FILL_TRIANGLE_WIRE);
                    break;
                }
                if (event.key.keysym.sym == SDLK_5) {
                    set_render_method(RENDER_TEXTURED);
                    break;
                }
                if (event.key.keysym.sym == SDLK_6) {
                    set_render_method(RENDER_TEXTURED_WIRE);
                    break;
                }
                if (event.key.keysym.sym == SDLK_c) {
                    set_cull_method(CULL_BACKFACE);
                    break;
                }
                if (event.key.keysym.sym == SDLK_x) {
                    set_cull_method(CULL_NONE);
                    break;
                }
                // camera movement, rotation
                if (event.key.keysym.sym == SDLK_UP) {
                    camera.position.y += 3.0 * delta_time;
                    break;
                }
                if (event.key.keysym.sym == SDLK_DOWN) {
                    camera.position.y -= 3.0 * delta_time;
                    break;
                }
                if (event.key.keysym.sym == SDLK_RIGHT) {
                    camera.position.x += 2.0 * delta_time;
                    break;
                }
                if (event.key.keysym.sym == SDLK_LEFT) {
                    camera.position.x -= 2.0 * delta_time;
                    break;
                }
                if (event.key.keysym.sym == SDLK_a) {
                    camera.yaw -= 5.0 * delta_time;
                    break;
                }
                if (event.key.keysym.sym == SDLK_d) {
                    camera.yaw += 5.0 * delta_time;
                    break;
                }
                if (event.key.keysym.sym == SDLK_w) {
                    camera.forward_velocity = vec3_mul(camera.direction, 5.0 * delta_time);
                    camera.position = vec3_add(camera.position, camera.forward_velocity);
                    break;
                }
                if (event.key.keysym.sym == SDLK_s) {
                    camera.forward_velocity = vec3_mul(camera.direction, 5.0 * delta_time);
                    camera.position = vec3_sub(camera.position, camera.forward_velocity);
                    break;
                }
                break;
        }
    }
}


void process_graphics_pipeline_stages(mesh_t* mesh) {
    // create scale, rotation, translation matrices
    mat4_t scale_matrix = mat4_make_scale(mesh->scale.x, mesh->scale.y, mesh->scale.z);
    mat4_t translation_matrix = mat4_make_translation(mesh->translation.x, mesh->translation.y, mesh->translation.z);
    mat4_t rotation_matrix_x = mat4_make_rotation_x(mesh->rotation.x);
    mat4_t rotation_matrix_y = mat4_make_rotation_y(mesh->rotation.y);
    mat4_t rotation_matrix_z = mat4_make_rotation_z(mesh->rotation.z);

    // Loop all triangle faces of our mesh       
    int num_faces = array_length(mesh->faces);
    for (int i = 0; i < num_faces; i++) {
        face_t mesh_face = mesh->faces[i];

        vec3_t face_vertices[3];
        face_vertices[0] = mesh->vertices[mesh_face.a];
        face_vertices[1] = mesh->vertices[mesh_face.b];
        face_vertices[2] = mesh->vertices[mesh_face.c];

        vec4_t transformed_vertices[3];

        // Loop all three vertices of this current face and apply transformations
        for (int j = 0; j < 3; j++) {
            vec4_t transformed_vertex = vec4_from_vec3(face_vertices[j]);

            //create a world matrix combining scale, rotation, translation matrices
            world_matrix = mat4_identity();

            // multiply all matrices
            world_matrix = mat4_mul_mat4(scale_matrix, world_matrix);
            world_matrix = mat4_mul_mat4(rotation_matrix_z, world_matrix);
            world_matrix = mat4_mul_mat4(rotation_matrix_y, world_matrix);
            world_matrix = mat4_mul_mat4(rotation_matrix_x, world_matrix);
            world_matrix = mat4_mul_mat4(translation_matrix, world_matrix);

            // multiply world matrix by the original vector
            transformed_vertex = mat4_mul_vec4(world_matrix, transformed_vertex);

            // multiply the view matrix by the vector to transform the scene to camera space
            transformed_vertex = mat4_mul_vec4(view_matrix, transformed_vertex);
            

            // Save transformed vertex in the array of transformed vertices
            transformed_vertices[j] = transformed_vertex;
        }

        // calculate the triangle face normal
        vec3_t face_normal = get_triangle_normal(transformed_vertices);

        if (is_cull_backface()) {
            // Find the vector between vertex A in the triangle and the camera origin
            vec3_t origin = { 0, 0, 0 };
            vec3_t camera_ray = vec3_sub(origin, vec3_from_vec4(transformed_vertices[0]));

            // Calculate how aligned the camera ray is with the face normal (using dot product)
            float dot_normal_camera = vec3_dot(face_normal, camera_ray);

            // Bypass the triangles that are looking away from the camera
            if (dot_normal_camera < 0) {
                continue;
            }
        }

        // Clipping
        polygon_t polygon = create_polygon_from_triangle(
            vec3_from_vec4(transformed_vertices[0]),
            vec3_from_vec4(transformed_vertices[1]),
            vec3_from_vec4(transformed_vertices[2]),
            mesh_face.a_uv,
            mesh_face.b_uv,
            mesh_face.c_uv
        );

        // clips the polygon and returns a new polygon with potential new vertices
        clip_polygon(&polygon);

        // after clipping break the polygon into triangle
        triangle_t triangles_after_clipping[MAX_NUM_POLY_TRIANGLES];
        int num_triangles_after_clpping = 0;

        triangles_from_polygon(&polygon, triangles_after_clipping, &num_triangles_after_clpping);

        for (int t = 0; t < num_triangles_after_clpping; t++) {
            triangle_t triangle_after_clipping = triangles_after_clipping[t];

            vec4_t projected_points[3];

            // Loop all three vertices to perform projection
            for (int j = 0; j < 3; j++) {
                // Project the current vertex
                projected_points[j] = mat4_mul_vec4_project(proj_matrix, triangle_after_clipping.points[j]);

                // invert the y values to account for flipped screen y coords
                projected_points[j].y *= -1;

                // scale into the view
                projected_points[j].x *= (get_window_width() / 2.0); 
                projected_points[j].y *= (get_window_height() / 2.0); 

                // translate the projected points to the middle of the screen
                projected_points[j].x += (get_window_width() / 2.0);
                projected_points[j].y += (get_window_height() / 2.0);
            }

            //calculate the shade intensity based on the dot product between face normal and light ray
            float light_intensity_factor = -vec3_dot(face_normal, get_light_dir());

            // calculate the triangle color based on the light angle
            uint32_t triangle_color = light_apply_intensity(mesh_face.color, light_intensity_factor);

            triangle_t triangle_to_render = {
                .points = {
                    { projected_points[0].x, projected_points[0].y, projected_points[0].z, projected_points[0].w },
                    { projected_points[1].x, projected_points[1].y, projected_points[1].z, projected_points[1].w },
                    { projected_points[2].x, projected_points[2].y, projected_points[2].z, projected_points[2].w }
                },
                .texcoords = {
                    { triangle_after_clipping.texcoords[0].u, triangle_after_clipping.texcoords[0].v },
                    { triangle_after_clipping.texcoords[1].u, triangle_after_clipping.texcoords[1].v },
                    { triangle_after_clipping.texcoords[2].u, triangle_after_clipping.texcoords[2].v }
                },
                .color = triangle_color,
                .texture = mesh->texture
            };

            // Save the projected triangle in the array of triangles to render
            if (num_triangles_to_render < MAX_TRIANGLES_PER_MESH) {
                triangles_to_render[num_triangles_to_render++] = triangle_to_render;
            }
        }
    }
}

void update(void) {
    // Wait some time until the reach the target frame time in milliseconds
    int time_to_wait = FRAME_TARGET_TIME - (SDL_GetTicks() - previous_frame_time);
    // Only delay execution if we are running too fast
    if (time_to_wait > 0 && time_to_wait <= FRAME_TARGET_TIME) {
        SDL_Delay(time_to_wait);
    }
    // get a delta time factor converted to seconds to be used to update our game objects
    delta_time = (SDL_GetTicks() - previous_frame_time) / 1000.0;
    previous_frame_time = SDL_GetTicks();

    // Initialize the counter of triangles to render for the current face
    num_triangles_to_render = 0;

    for (int mesh_index = 0; mesh_index < get_num_meshes(); mesh_index++) {
        mesh_t* mesh = get_mesh(mesh_index);

        //mesh.rotation.x += 0.5 * delta_time;
        //mesh.rotation.y += 0.8 * delta_time;
        //mesh.translation.z = 5.0;

        // change the cam position per animation frame
        //camera.position.x += 0.5 * delta_time; 
        //camera.position.y += 0.3 * delta_time; 

        vec3_t up_dir = { 0, 1, 0 };

        vec3_t target = { 0, 0, 1 };
        mat4_t camera_yaw_rotation = mat4_make_rotation_y(camera.yaw);
        camera.direction = vec3_from_vec4(mat4_mul_vec4(camera_yaw_rotation, vec4_from_vec3(target)));

        // offset cam position in the direction where the cam is pointing at
        target = vec3_add(camera.position, camera.direction);
        // create the view matrix
        view_matrix = mat4_look_at(camera.position, target, up_dir);

        process_graphics_pipeline_stages(mesh);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Render function to draw objects on the display
///////////////////////////////////////////////////////////////////////////////
void render(void) {
    clear_color_buffer(0xFF000000);
    clear_z_buffer();

    draw_grid();

    // Loop all projected triangles and render them
    for (int i = 0; i < num_triangles_to_render; i++) {
        triangle_t triangle = triangles_to_render[i];

        // Draw filled triangle
        if (should_render_filled_triangles()) {
            draw_filled_triangle(
                triangle.points[0].x, triangle.points[0].y, triangle.points[0].z, triangle.points[0].w, // vertex A
                triangle.points[1].x, triangle.points[1].y, triangle.points[1].z, triangle.points[1].w,// vertex B
                triangle.points[2].x, triangle.points[2].y, triangle.points[2].z, triangle.points[2].w, // vertex C
                triangle.color
            );
        }

        if (should_render_textured_triangles()) {
            draw_textured_triangle(
                triangle.points[0].x, triangle.points[0].y, triangle.points[0].z, triangle.points[0].w, triangle.texcoords[0].u, triangle.texcoords[0].v, // vertex A
                triangle.points[1].x, triangle.points[1].y, triangle.points[1].z, triangle.points[1].w, triangle.texcoords[1].u, triangle.texcoords[1].v, // vertex B
                triangle.points[2].x, triangle.points[2].y, triangle.points[2].z, triangle.points[2].w, triangle.texcoords[2].u, triangle.texcoords[2].v, // vertex C
                triangle.texture
            );
        }

        // Draw triangle wireframe
        if (should_render_wireframe()) {
            draw_triangle(
                triangle.points[0].x, triangle.points[0].y, // vertex A
                triangle.points[1].x, triangle.points[1].y, // vertex B
                triangle.points[2].x, triangle.points[2].y, // vertex C
                0xFFFFFFFF
            );
        }

        // Draw triangle vertex points
        if (should_render_wire_vertex()) {
            draw_rect(triangle.points[0].x - 3, triangle.points[0].y - 3, 6, 6, 0xFFFF0000); // vertex A
            draw_rect(triangle.points[1].x - 3, triangle.points[1].y - 3, 6, 6, 0xFFFF0000); // vertex B
            draw_rect(triangle.points[2].x - 3, triangle.points[2].y - 3, 6, 6, 0xFFFF0000); // vertex C
        }
    }

    render_color_buffer();
}


void free_resources(void) {
    free_meshes();
    destroy_window();
}


int main(void) {
    is_running = initialize_window();

    setup();

    while (is_running) {
        process_input();
        update();
        render();
    }

    free_resources();

    return 0;
}
