// #define SR_DEBUG
#define DEBUGUI
#define SDL_MAIN_HANDLED

#ifdef DEBUGUI
#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_STANDARD_BOOL 
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_INCLUDE_VERTEX_BUFFER_OUTPUT
#define NK_INCLUDE_FONT_BAKING
#define NK_INCLUDE_DEFAULT_FONT
#define NK_IMPLEMENTATION
#define NK_SDL_RENDERER_IMPLEMENTATION
#include "external/nuklear/nuklear.h"
#include "external/nuklear/nuklear_sdl_renderer.h"
#endif // DEBUGUI

#include <SDL.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "simplerenderer.h"
#include "objloader.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define WIDTH 1066
#define HEIGHT 600
#define ROTATION_SPEED 15.0f // degrees per second
#define VIEWPORT_WIDTH (WIDTH /(float) HEIGHT)
#define VIEWPORT_HEIGHT 1.0f
#define VIEWPORT_DISTANCE 1.0f
#define PIXEL_DEPTH 4
#define PITCH (PIXEL_DEPTH * WIDTH)


typedef struct game_state_t {
    // Loop control
    bool         running;
    double       elapsedTime;
    uint64_t     lastTime;

    // Rendering
    SDL_Event     event;
    SDL_Window*   window;
    SDL_Renderer* renderer;
    SDL_Texture*  texture;
    canvas_t      canvas;
    color_t       backgroundColor;
    bool          drawLights;
    bool          draw3DObjects;
    bool          draw2DObjects;
    bool          drawWire;
    bool          drawFilled;
    bool          diffuseLighting;
    bool          specularLighting;
    bool          shaded;
    int           shadingMode;
    bool          backfaceCulling;
    bool          bilinearFiltering;

    // Game objects
    int              numMeshes;
    mesh_t*          meshes;
    int              numObjects;
    object3D_t*      objects;
    int              numAmbientLights;
    ambient_light_t* ambientLights;
    int              numDirectionalLights;
    dir_light_t*     directionalLights;
    int              numPointLights;
    point_light_t*   pointLights;
    object3D_t*      pointLightObjects;

    // Camera
    camera_t     camera;
    
    // Animation
    float        rotationSpeed;

    // Input
    const uint8_t* keys;

    // GUI
    #ifdef DEBUGUI
    struct nk_context* nuklearContext;
    bool               showGUI;
    bool               toggleGUIKeyPressed;
    #endif // DEBUG
} game_state_t;


// TODO: Code here is a bit repeated between directional and point lights. Maybe refactor?
// TODO: Make it independent of game and move it to simplerenderer.h
float shadeVertex(vec3_t v, vec3_t normal, float invMagnitudeNormal, float specularExponent, game_state_t* game) {
    float diffuseIntensity  = 0.0;
    float specularIntensity = 0.0;
    float ambientIntensity  = 0.0;
    
    // Directional lights
    for (int i = 0; i < game->numDirectionalLights; i++) {
        vec3_t lightDirection = game->directionalLights[i].direction;
        float magnitudeLightDirection = magnitude(lightDirection);
        float invMagnitudeLightDirection = 1.0f / magnitudeLightDirection;
        if (game->diffuseLighting) {
            float cos_alpha = -dot(lightDirection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            diffuseIntensity += MAX(cos_alpha, 0.0f) * game->directionalLights[i].intensity;
        }

        if (game->specularLighting) {
            vec3_t reflection = sub(mulScalarV3(2 * -dot(lightDirection, normal), normal), lightDirection);
            float cos_beta = -dot(reflection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            specularIntensity += pow(MAX(cos_beta, 0.0f), specularExponent) * game->directionalLights[i].intensity;
        }
    }

    // Point lights
    for (int i = 0; i < game->numPointLights; i++) {
        vec3_t lightDirection = sub(game->pointLights[i].position, v);
        float invMagnitudeLightDirection = 1.0f / magnitude(lightDirection);
        if (game->diffuseLighting) {
            float cos_alpha = dot(lightDirection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            diffuseIntensity += MAX(cos_alpha, 0) * game->pointLights[i].intensity;
        }

        if (game->specularLighting) {
            vec3_t reflection = sub(mulScalarV3(2 * dot(lightDirection, normal), normal), lightDirection);
            float cos_beta = dot(reflection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            specularIntensity += pow(MAX(cos_beta, 0), specularExponent) * game->pointLights[i].intensity;
        }
    }

    // Ambient light
    for (int i = 0; i < game->numAmbientLights; i++) {
        ambientIntensity += game->ambientLights[i].intensity;
    }

    return (diffuseIntensity + specularIntensity + ambientIntensity);
}

// TODO: Move to simplerenderer.h
int edgeCross(int ax, int ay, int bx, int by, int px, int py) {
  point_t ab = { bx - ax, by - ay };
  point_t ap = { px - ax, py - ay };
  return ab.x * ap.y - ab.y * ap.x;
}

// TODO: Move to simplerenderer.h
void drawTriangleWireframe(int x0, int x1, int x2,
                           int y0, int y1, int y2,
                           color_t color, canvas_t canvas) {
    drawLine(x0, x1, y0, y1, color, canvas);
    drawLine(x1, x2, y1, y2, color, canvas);
    drawLine(x2, x0, y2, y0, color, canvas);
}

// TODO: Remove dependency with game_state_t and move to simplerenderer.h
void drawTriangleFilled(int x0, int x1, int x2,
                        int y0, int y1, int y2,
                        float invz0, float invz1, float invz2,
                        float i0, float i1, float i2,
                        vec3_t n0, vec3_t n1, vec3_t n2,
                        vec3_t t0, vec3_t t1, vec3_t t2,
                        color_t c0, color_t c1, color_t c2,
                        float specularExponent,
                        uint32_t* texture, int textureWidth, int textureHeight,
                        mat4x4_t invCameraTransform,
                        int area,
                        game_state_t* game) {
    canvas_t canvas = game->canvas;
    camera_t camera = game->camera;

    int x_min = MAX(MIN(MIN(x0, x1), x2), 0);
    int x_max = MIN(MAX(MAX(x0, x1), x2), WIDTH - 1); 
    int y_min = MAX(MIN(MIN(y0, y1), y2), 0);
    int y_max = MIN(MAX(MAX(y0, y1), y2), HEIGHT - 1);

    float invArea = 1.0f / area;

    // Compute the constant delta_s that will be used for the horizontal and vertical steps
    int delta_w0_col = (y1 - y2);
    int delta_w1_col = (y2 - y0);
    int delta_w2_col = (y0 - y1);
    int delta_w0_row = (x2 - x1);
    int delta_w1_row = (x0 - x2);
    int delta_w2_row = (x1 - x0);

    // Compute the edge functions for the fist (top-left) point
    int w0_row = edgeCross(x1, y1, x2, y2, x_min, y_min);
    int w1_row = edgeCross(x2, y2, x0, y0, x_min, y_min);
    int w2_row = edgeCross(x0, y0, x1, y1, x_min, y_min);

    // Illuminate each vertex
    if (game->shaded && !(game->shadingMode == 2)) {
        c0 = mulScalarColor(i0, c0);
        c1 = mulScalarColor(i1, c1);
        c2 = mulScalarColor(i2, c2);
    }
    
    for (int y = y_min; y <= y_max; y++) {
        bool was_inside = false;
        int w0 = w0_row;
        int w1 = w1_row;
        int w2 = w2_row;
        for (int x = x_min; x <= x_max; x++) {
            bool is_inside = (w0 | w1 | w2) >= 0;
            if (is_inside) {
                was_inside = true;
            
                float alpha = w0 * invArea;
                float beta  = w1 * invArea;
                float gamma = w2 * invArea;
                float invz = alpha * invz0 + beta * invz1 + gamma * invz2;
                if (invz > canvas.depthBuffer[y * WIDTH + x]) {
                    uint32_t color = colorToUint32(c0); // Fallback in case of no texture and no shading
                    float light = 1;

                    if (game->shadingMode == 2) { // phong shading
                        point_t p = {x, y, invz};
                        vec3_t v = mulMV3(invCameraTransform, unprojectPoint(p, canvas, camera));
                        vec3_t normal = add(add(mulScalarV3(alpha, n0), mulScalarV3(beta, n1)), mulScalarV3(gamma, n2));
                        light = shadeVertex(v, normal , 1/magnitude(normal), specularExponent, game);
                    } else if (game->shaded) {
                        light = alpha * i0 + beta * i1 + gamma * i2;
                    }
                    
                    if (textureWidth != 0 && textureHeight != 0) {
                        // Interpolate u/z and v/z to get perspective correct texture coordinates
                        float u_over_z = alpha * (t0.x * invz0) + beta * (t1.x * invz1) + gamma * (t2.x * invz2);
                        float v_over_z = alpha * (t0.y * invz0) + beta * (t1.y * invz1) + gamma * (t2.y * invz2);
                        color_t color_typed;
                        // TODO: Fix crash when we have overflow here
                        if (game->bilinearFiltering) {
                            float tex_u = u_over_z/invz;
                            if (tex_u < 0) {
                                tex_u = 1 + tex_u;
                            }
                            tex_u = MIN(tex_u * textureWidth, textureWidth - 1);

                            float tex_v = v_over_z/invz;
                            if (tex_v < 0) {
                                tex_v = 1 + tex_v;
                            }
                            tex_v = MIN(tex_v * textureHeight, textureHeight - 1);

                            int floor_u = floor(tex_u);
                            int floor_v = floor(tex_v);
                            int next_u = MIN(floor_u + 1, textureWidth - 1);
                            int next_v = MIN(floor_v + 1, textureHeight - 1);
                            float frac_u = tex_u - floor_u;
                            float frac_v = tex_v - floor_v;
                            color_t color_tl = colorFromUint32(texture[floor_v * textureWidth + floor_u]);
                            color_t color_tr = colorFromUint32(texture[floor_v * textureWidth + next_u]);
                            color_t color_bl = colorFromUint32(texture[next_v * textureWidth + floor_u]);
                            color_t color_br = colorFromUint32(texture[next_v * textureWidth + next_u]);
                            color_t color_b = sumColors(mulScalarColor(1 - frac_u, color_bl), mulScalarColor(frac_u, color_br));
                            color_t color_tp = sumColors(mulScalarColor(1 - frac_u, color_tl), mulScalarColor(frac_u, color_tr));
                            color_typed = sumColors(mulScalarColor(1 - frac_v, color_b), mulScalarColor(frac_v, color_tp));
                        } else {
                            int tex_x = MIN(abs((int)((u_over_z/invz) * textureWidth)), textureWidth - 1);
                            int tex_y = MIN(abs((int)((v_over_z/invz) * textureHeight)), textureHeight - 1);
                            color_typed = colorFromUint32(texture[tex_y * textureWidth + tex_x]);
                        }
                        
                        color_t color_shaded = mulScalarColor(light, color_typed);
                        color = colorToUint32(color_shaded);
                    } else if (game->shaded) {
                        color_t color_typed = {
                            (uint8_t) (c0.r * alpha + c1.r * beta + c2.r * gamma),
                            (uint8_t) (c0.g * alpha + c1.g * beta + c2.g * gamma),
                            (uint8_t) (c0.b * alpha + c1.b * beta + c2.b * gamma)
                        };
                        color_t color_shaded = mulScalarColor(light, color_typed);
                        color = colorToUint32(color_shaded);
                    }
                    drawPixelDepthBuffer(x, y, invz, color, canvas);
                }
            }

            // Go to next row if we jumped outside the triangle
            if (!is_inside && was_inside) {
                break;
            }

            w0 += delta_w0_col;
            w1 += delta_w1_col;
            w2 += delta_w2_col;
        }
        w0_row += delta_w0_row;
        w1_row += delta_w1_row;
        w2_row += delta_w2_row;
    }
}

// TODO: Maybe move to simplerenderer.h
void drawObject(object3D_t* object, game_state_t* game) {
    // TODO: This function is very dumb, as it is allocating transformed
    //       vertices just to allocate the projections right after.

    mesh_t* mesh = object->mesh;
    camera_t camera = game->camera;
    material_t* materials = mesh->materials;
    mat4x4_t invCameraTransform = inverseM4(camera.transform);
    
    // Transform the bounding sphere and don't draw when object is fully out of the camera volume
    mat4x4_t transform = mulMM4(camera.transform, object->transform);
    vec3_t transformedCenter = mulMV3(transform, mesh->center);
    float transformedBoundsRadius = mesh->boundsRadius * object->scale;
    for (int p = 0; p < camera.numPlanes; p++) {
        float distance = distancePlaneV3(camera.planes[p], transformedCenter);
        if (distance < -transformedBoundsRadius) {
            DEBUG_PRINT("DEBUG: Clipped an object fully outside of the clipping volume\n");
            return;
        }
    }

    // If the object is not discarded, transform and project all vertices
    vec3_t *transformed = (vec3_t*) malloc(mesh->numVertices * sizeof(vec3_t));
    vec3_t *camTransformed = (vec3_t*) malloc(mesh->numVertices * sizeof(vec3_t));
    point_t *projected = (point_t*) malloc(mesh->numVertices * sizeof(point_t));
    vec3_t *transformedNormals = (vec3_t*) malloc(mesh->numNormals * sizeof(vec3_t));
    if (projected == NULL || transformed == NULL || camTransformed == NULL || transformedNormals == NULL) {
        fprintf(stderr, "ERROR: Transformed vertices/normals memory couldn't be allocated.\n");
        exit(-1);
    }

    for (int i = 0; i < mesh->numVertices; i++) {
        transformed[i] = mulMV3(object->transform, mesh->vertices[i]);
        camTransformed[i] = mulMV3(camera.transform, transformed[i]);
        projected[i] = projectVertex(camTransformed[i], game->canvas, camera);
    }

    for (int i = 0; i < mesh->numNormals; i++) {
        transformedNormals[i] = mulMV3(object->transform, mesh->normals[i]);
    }

    // Cull, shade and draw each triangle
    for (int i = 0; i < mesh->numTriangles; i++) {
        triangle_t triangle = mesh->triangles[i];

        bool discarded = false;

        // Backface culling
        point_t p0 = projected[triangle.v0];
        point_t p1 = projected[triangle.v1];
        point_t p2 = projected[triangle.v2];
        int area = edgeCross(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
        if (area < 0 && game->backfaceCulling) {
            discarded = true;
        }

        // Fustrum culling
        // TODO: Add config for this
        for (int p = 0; !discarded && p < camera.numPlanes; p++) {
            plane_t plane = camera.planes[p];
            if (distancePlaneV3(plane, camTransformed[triangle.v0]) < 0 &&
                distancePlaneV3(plane, camTransformed[triangle.v1]) < 0 &&
                distancePlaneV3(plane, camTransformed[triangle.v2]) < 0) {
                    DEBUG_PRINT("DEBUG: Clipped triangle fully outside of the camera clipping volume\n");
                    discarded = true;
                    break;
            }

            // TODO: Deal with the case where only one or two vertexes of the triangle
            //       are out of the volume. In this case, we should split the triangle
            //       and create new ones.
        }

        // TODO: Check if this is right, but it fixes overflows when we
        //       interpolate zs and we have z close to 0
        // Don't draw if the triangle has any vertice behind the camera
        if (camTransformed[triangle.v0].z < 0 ||
            camTransformed[triangle.v1].z < 0 ||
            camTransformed[triangle.v2].z < 0) {
            DEBUG_PRINT("DEBUG: Clipped triangle with a vertice behind the camera\n");
            discarded = true;
        }

        if (!discarded) {
            float i0 = 1.0;
            float i1 = 1.0;
            float i2 = 1.0;

            // Shading
            if (game->shaded && game->shadingMode == 0) {
                // Flat shading
                vec3_t v0 = transformed[triangle.v0];
                vec3_t v1 = transformed[triangle.v1];
                vec3_t v2 = transformed[triangle.v2];
                vec3_t normal = triangleNormal(v0, v1, v2);
                float invMag = 1.0f / magnitude(normal);
                vec3_t center = {(v0.x + v1.x + v2.x)/3.0f, (v0.y + v1.y + v2.y)/3.0f, (v0.z + v1.z + v2.z)/3.0f};
                float intensity = shadeVertex(center, normal, invMag, materials[triangle.materialIndex].specularExponent, game);
                i0 = intensity;
                i1 = intensity;
                i2 = intensity;
            } else if (game->shaded && game->shadingMode == 1) {
                // Gouraud shading
                float specularExponent = 900.0f;
                if (mesh->numMaterials != 0) {
                    specularExponent = materials[triangle.materialIndex].specularExponent;
                }

                i0 = shadeVertex(transformed[triangle.v0], transformedNormals[triangle.n0], mesh->invMagnitudeNormals[triangle.n0], specularExponent, game);
                i1 = shadeVertex(transformed[triangle.v1], transformedNormals[triangle.n1], mesh->invMagnitudeNormals[triangle.n1], specularExponent, game);
                i2 = shadeVertex(transformed[triangle.v2], transformedNormals[triangle.n2], mesh->invMagnitudeNormals[triangle.n2], specularExponent, game);
            }

            // Drawing
            if (game->drawWire) {
                drawTriangleWireframe(p0.x, p1.x, p2.x,
                                      p0.y, p1.y, p2.y,
                                      materials[triangle.materialIndex].diffuseColor,
                                      game->canvas);
            }
            
            if (game->drawFilled) {
                vec3_t t0 = {0};
                vec3_t t1 = {0};
                vec3_t t2 = {0};
                int textureWidth = 0;
                int textureHeight = 0;
                uint32_t* texture = NULL;
                if (mesh->numTextureCoords > 0) {
                    t0 = mesh->textureCoords[triangle.t0];
                    t1 = mesh->textureCoords[triangle.t1];
                    t2 = mesh->textureCoords[triangle.t2];
                    textureWidth = materials[triangle.materialIndex].textureWidth;
                    textureHeight = materials[triangle.materialIndex].textureHeight;
                    texture = materials[triangle.materialIndex].texture;
                }
                
                color_t color = COLOR_WHITE;
                float specularExponent = 900.0f;

                if (mesh->numMaterials != 0) {
                    color = materials[triangle.materialIndex].diffuseColor;
                    specularExponent = materials[triangle.materialIndex].specularExponent;
                }
                
                drawTriangleFilled(p0.x, p1.x, p2.x,
                                   p0.y, p1.y, p2.y,
                                   p0.invz, p1.invz, p2.invz,
                                   i0, i1, i2,
                                   transformedNormals[triangle.n0], transformedNormals[triangle.n1], transformedNormals[triangle.n2],
                                   t0, t1, t2,
                                   color, color, color,
                                   specularExponent,
                                   texture, textureWidth, textureHeight,
                                   invCameraTransform,
                                   area,
                                   game);
            }
        }
    }

    free(transformed);
    free(camTransformed);
    free(projected);
    free(transformedNormals);
}

void drawObjects(game_state_t* game) {
    for (int i = 0; i < game->numObjects; i++) {
        drawObject(&game->objects[i], game);
    }
}

// TODO: This is a huge hack, we should have a proper way to solve shading in this case
void drawLights(game_state_t* game) {
    for (int i = 0; i < game->numPointLights; i++) {
        bool shadedBackup = game->shaded;
        game->shaded = false; 
        drawObject(&game->pointLightObjects[i], game);
        game->shaded = shadedBackup;
    }
}

void updateCameraPosition(game_state_t* game) {
    camera_t *camera = &game->camera;
    const uint8_t* keys = game->keys;
    camera_t newCamera;
    vec3_t localTranslation = {0};
    vec3_t newTranslation = camera->translation;
    mat4x4_t newRotation = camera->rotation;
    float elapsedTime = game->elapsedTime / 1000.0f;
    float movementSpeed = camera->movementSpeed * elapsedTime;
    float turningSpeed = camera->turningSpeed * elapsedTime;
    

    if (keys[SDL_SCANCODE_A]) {
        localTranslation.x = -movementSpeed;
    }
    
    if (keys[SDL_SCANCODE_D]) {
        localTranslation.x = movementSpeed;
    }

    if (keys[SDL_SCANCODE_PAGEDOWN]) {
        localTranslation.y = -movementSpeed;
    }
    
    if (keys[SDL_SCANCODE_PAGEUP]) {
        localTranslation.y = movementSpeed;
    }

    if (keys[SDL_SCANCODE_S]) {
        localTranslation.z = -movementSpeed;
    }
    
    if (keys[SDL_SCANCODE_W]) {
        localTranslation.z = movementSpeed;
    }

    if (keys[SDL_SCANCODE_RIGHT]) {
        newRotation = mulMM4(rotationY(-turningSpeed), newRotation);
    }

    if (keys[SDL_SCANCODE_LEFT]) {
        newRotation = mulMM4(rotationY(turningSpeed), newRotation);
    }

    if (keys[SDL_SCANCODE_UP]) {
        newRotation = mulMM4(rotationX(turningSpeed), newRotation);
    }

    if (keys[SDL_SCANCODE_DOWN]) {
        newRotation = mulMM4(rotationX(-turningSpeed), newRotation);
    }

    vec3_t globalTranslation = mulMV3(camera->rotation, localTranslation);
    newTranslation.x += globalTranslation.x;
    newTranslation.y += globalTranslation.y;
    newTranslation.z += globalTranslation.z;

    free(camera->planes);
    *camera  = makeCamera(
        newTranslation,
        newRotation,
        camera->viewportWidth,
        camera->viewportHeight,
        camera->viewportDistance,
        camera->movementSpeed,
        camera->turningSpeed
    );
}

object3D_t rotateObjectY(object3D_t object, float degrees) {
    return makeObject(object.mesh, object.translation, object.scale, mulMM4(rotationY(degrees), object.rotation));
}

void animateObjects(game_state_t* game) {
    float degrees = game->rotationSpeed * (game->elapsedTime / 1000.0f);
    game->objects[0] = rotateObjectY(game->objects[0], fmod(degrees, 360.0f));
}

game_state_t* init() {
    DEBUG_PRINT("INFO: Initializing SDL\n");
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        fprintf(stderr, "ERROR: error initializing SDL: %s\n", SDL_GetError());
        exit(-1);
    }

    DEBUG_PRINT("INFO: Loading meshes and objects\n");
    int numObjects = 1;
    int numMeshes = 4;
    mesh_t* meshes = (mesh_t*) malloc(numMeshes * sizeof(mesh_t));
    object3D_t *objects = (object3D_t*) malloc(numObjects * sizeof(object3D_t));
    if (meshes == NULL || objects == NULL) {
        fprintf(stderr, "ERROR: 3D objects memory couldn't be allocated.\n");
        exit(-1);
    }
        
    meshes[0] = *loadObjFile("assets/light.obj", false);
    meshes[1] = *loadObjFile("assets/snake/snake.obj", true);
    meshes[2] = *loadObjFile("assets/engineer/engineer.obj", false);
    meshes[3] = *loadObjFile("assets/cube.obj", false);

    objects[0] = makeObject(&meshes[2], (vec3_t) {0, 0, 0}, 1.0 , IDENTITY_M4x4);

    DEBUG_PRINT("INFO: Loading lights\n");
    int numAmbientLights = 1;
    int numDirLights = 1;
    int numPointLights = 1;
    ambient_light_t* ambientLights = (ambient_light_t*) malloc(numAmbientLights * sizeof(ambient_light_t));
    dir_light_t* directionalLights = (dir_light_t*) malloc(numDirLights * sizeof(dir_light_t));
    point_light_t* pointLights = (point_light_t*) malloc(numPointLights * sizeof(point_light_t));
    object3D_t *pointLightObjects = (object3D_t*) malloc(numPointLights * sizeof(object3D_t));
    if (pointLights == NULL || directionalLights == NULL || ambientLights == NULL || pointLightObjects == NULL) {
        fprintf(stderr, "ERROR: Lights memory couldn't be allocated.\n");
        exit(-1);
    }

    ambientLights[0] = (ambient_light_t) {0.4};
    directionalLights[0] = (dir_light_t) {0.0, {0.0, -1.0, 1.0}};
    pointLights[0] = (point_light_t) {0.9, {-0.5, 1.5, -2.0}};

    pointLightObjects[0] = makeObject(&meshes[0], pointLights[0].position, 0.05, IDENTITY_M4x4);


    DEBUG_PRINT("INFO: Initializing game state\n");
    uint32_t *frameBuffer = (uint32_t*) malloc(WIDTH * HEIGHT * sizeof(uint32_t));
    float *depthBuffer = (float*) malloc(WIDTH * HEIGHT * sizeof(float));
    game_state_t* game = (game_state_t*) malloc(sizeof(game_state_t));
    if (game == NULL || frameBuffer == NULL || depthBuffer == NULL) {
        fprintf(stderr, "ERROR: Game state memory and buffers couldn't be allocated.\n");
        exit(-1);
    }

    struct canvas_t canvas = {
        .width = WIDTH,
        .height = HEIGHT,
        .hasDepthBuffer = 1,
        .frameBuffer = frameBuffer,
        .depthBuffer = depthBuffer
    };
    
    game->running = true;
    game->elapsedTime = 0;
    game->lastTime = SDL_GetPerformanceCounter();
    game->window = SDL_CreateWindow("Rasterizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
    game->renderer = SDL_CreateRenderer(game->window, -1, 0);
    game->texture = SDL_CreateTexture(game->renderer, SDL_PIXELFORMAT_BGRA32, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
    game->canvas = canvas;
    game->backgroundColor = (color_t) {0, 0, 0};
    game->drawLights    = true;
    game->draw3DObjects =  true;
    game->draw2DObjects =  false;
    game->diffuseLighting = true;
    game->specularLighting = true;
    game->drawWire = false;
    game->drawFilled = true;
    game->shaded = true;
    game->shadingMode = 1;
    game->backfaceCulling = true;
    game->bilinearFiltering = false;
    game->numMeshes = numMeshes;
    game->meshes = meshes;
    game->numObjects = numObjects;
    game->objects = objects;
    
    // Lights
    game->numAmbientLights = numAmbientLights;
    game->ambientLights = ambientLights;
    game->numDirectionalLights = numDirLights;
    game->directionalLights = directionalLights;
    game->numPointLights = numPointLights;
    game->pointLights = pointLights;
    game->pointLightObjects = pointLightObjects;
    
    game->camera = makeCamera(
        (vec3_t) {0, 0, -5},
        rotationY(0.0f),
        VIEWPORT_WIDTH,
        VIEWPORT_HEIGHT,
        VIEWPORT_DISTANCE,
        5.0f,
        90.0f
    );
    game->rotationSpeed = 15.0f;
    game->keys = SDL_GetKeyboardState(NULL);

    #ifdef DEBUGUI
    DEBUG_PRINT("INFO:  Initializing Dear ImGui\n");
    struct nk_context *ctx = nk_sdl_init(game->window, game->renderer);
    struct nk_font_atlas *atlas;
    struct nk_font_config config = nk_font_config(0);
    struct nk_font *font;
    nk_sdl_font_stash_begin(&atlas);
    font = nk_font_atlas_add_default(atlas, 10.0, &config);
    nk_sdl_font_stash_end();
    nk_style_set_font(ctx, &font->handle);
    game->nuklearContext = ctx;
    game->showGUI = true;
    game->toggleGUIKeyPressed = false;
    #endif // DEBUGUI

    return game;
}

void handleEvents(game_state_t* game) {
    DEBUG_PRINT("INFO: Handle events\n");

    #ifdef DEBUGUI
    nk_input_begin(game->nuklearContext);
    #endif // DEBUGUI
    
    while (SDL_PollEvent(&game->event)) {
        #ifdef DEBUGUI
        nk_sdl_handle_event(&game->event);
        #endif // DEBUGUI
        switch (game->event.type) {
        case SDL_QUIT:
            // handling of close button
            DEBUG_PRINT("INFO: Quitting application\n");
            game->running = false;
            break;
        }
    }
}

#ifdef DEBUGUI
void updateDebugUI(game_state_t *game) {
    // Check for toggle key
    if (game->keys[SDL_SCANCODE_SPACE]) {
        if (!game->toggleGUIKeyPressed) {
            game->toggleGUIKeyPressed = true;
            game->showGUI = !game->showGUI; // Toggle state of ImGui display
        }
    } else {
        game->toggleGUIKeyPressed = false;
    }

    if (game->showGUI) {
        DEBUG_PRINT("INFO: Updating GUI\n");
        int row_size = 12;
        struct nk_context *ctx = game->nuklearContext;
        if (nk_begin(ctx, "Settings", nk_rect(0, 0, 220, HEIGHT),
                     NK_WINDOW_BORDER|NK_WINDOW_MOVABLE|NK_WINDOW_CLOSABLE|NK_WINDOW_SCALABLE)) {
            
            nk_layout_row_static(ctx, row_size, 150, 1);
            nk_labelf(ctx, NK_TEXT_LEFT, "FPS: %.2f (%.3lf ms)", floor(1000.0f / game->elapsedTime), game->elapsedTime);

            if (nk_tree_push(ctx, NK_TREE_NODE, "Meshes", NK_MINIMIZED)) {
                for (int i = 0; i < game->numMeshes; i++) {
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_labelf(ctx, NK_TEXT_LEFT, "%s: (vertices %d, triangles: %d, uvs: %d)", game->meshes[i].name, game->meshes[i].numVertices, game->meshes[i].numTriangles, game->meshes[i].numTextureCoords);
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Objects", NK_MAXIMIZED)) {
                for (int i = 0; i < game->numObjects; i++) {
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_labelf(ctx, NK_TEXT_LEFT, "%d: %s at (%.0f, %.0f, %.0f)", i, game->objects[i].mesh->name, game->objects[i].translation.x, game->objects[i].translation.y, game->objects[i].translation.z);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "x", -10.0f, &game->objects[i].translation.x, 10.0f, 0.1f, 0.1f);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "y", -10.0f, &game->objects[i].translation.y, 10.0f, 0.1f, 0.1f);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "z", -10.0f, &game->objects[i].translation.z, 10.0f, 0.1f, 0.1f);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "scale", -0.01f, &game->objects[i].scale, 10.0f, 0.01f, 0.01f);
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Lights", NK_MINIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_checkbox_label(ctx, "Shaded", &game->shaded);
                if (game->shaded) {
                    nk_layout_row_dynamic(ctx, row_size, 3);
                    nk_bool isFlat = game->shadingMode == 0;
                    if (nk_radio_label(ctx, "Flat", &isFlat)) {
                        game->shadingMode = 0;
                    };
                    nk_bool isGouraud = game->shadingMode == 1;
                    if (nk_radio_label(ctx, "Gouraud", &isGouraud)) {
                        game->shadingMode = 1;
                    };

                    nk_bool isPhong = game->shadingMode == 2;
                    if (nk_radio_label(ctx, "Phong", &isPhong)) {
                        game->shadingMode = 2;
                    };

                    nk_layout_row_dynamic(ctx, row_size, 2);
                    nk_checkbox_label(ctx, "Difuse", &game->diffuseLighting);
                    nk_checkbox_label(ctx, "Specular", &game->specularLighting);

                    if (nk_tree_push(ctx, NK_TREE_NODE, "Ambient", NK_MAXIMIZED)) {
                        for (int i = 0; i < game->numAmbientLights; i++) {
                            nk_layout_row_dynamic(ctx, row_size * 4, 1);
                            char label[100];
                            sprintf(label, "Ambient Light %d", i);
                            if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "Intensity", 0.0f, &game->ambientLights[i].intensity, 3.0f, 0.1f, 0.1f);
                                nk_group_end(ctx);
                            }
                        }
                        nk_tree_pop(ctx);
                    }

                    if (nk_tree_push(ctx, NK_TREE_NODE, "Directional", NK_MAXIMIZED)) {
                        for (int i = 0; i < game->numDirectionalLights; i++) {
                            nk_layout_row_dynamic(ctx, row_size * 8, 1);
                            char label[100];
                            sprintf(label, "Directional Light %d", i);
                            if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "Intensity", 0.0f, &game->directionalLights[i].intensity, 3.0f, 0.1f, 0.1f);
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "x", -1.0f, &game->directionalLights[i].direction.x, 1.0f, 0.1f, 0.1f);
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "y", -1.0f, &game->directionalLights[i].direction.y, 1.0f, 0.1f, 0.1f);
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "z", -1.0f, &game->directionalLights[i].direction.z, 1.0f, 0.1f, 0.1f);
                                nk_group_end(ctx);
                            }
                        }
                        nk_tree_pop(ctx);
                    }

                    if (nk_tree_push(ctx, NK_TREE_NODE, "Point", NK_MAXIMIZED)) {
                        for (int i = 0; i < game->numPointLights; i++) {
                            nk_layout_row_dynamic(ctx, row_size * 8, 1);
                            char label[100];
                            sprintf(label, "Point Light %d", i);
                            if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "Intensity", 0.0f, &game->pointLights[i].intensity, 3.0f, 0.1f, 0.1f);
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "x", -10.0f, &game->pointLights[i].position.x, 10.0f, 0.1f, 0.1f);
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "y", -10.0f, &game->pointLights[i].position.y, 10.0f, 0.1f, 0.1f);
                                nk_layout_row_dynamic(ctx, row_size, 1);
                                nk_property_float(ctx, "z", -10.0f, &game->pointLights[i].position.z, 10.0f, 0.1f, 0.1f);
                                nk_group_end(ctx);
                            }
                        }
                        nk_tree_pop(ctx);
                    }
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Scene", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Camera: (%.1f, %.1f, %.1f)", game->camera.translation.x, game->camera.translation.y, game->camera.translation.z);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_label(ctx, "What to draw", NK_TEXT_LEFT);
                nk_layout_row_dynamic(ctx, row_size, 3);
                nk_checkbox_label(ctx, "3D Obj", &game->draw3DObjects);
                nk_checkbox_label(ctx, "2D Obj", &game->draw2DObjects);
                nk_checkbox_label(ctx, "Lights", &game->drawLights);
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Render Options", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 2);
                nk_checkbox_label(ctx, "Wireframe", &game->drawWire);
                nk_checkbox_label(ctx, "Filled", &game->drawFilled);
                nk_layout_row_dynamic(ctx, row_size, 2);
                nk_checkbox_label(ctx, "Backface culling", &game->backfaceCulling);
                nk_checkbox_label(ctx, "Bilinear filtering", &game->bilinearFiltering);
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Background", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_label(ctx, "Color", NK_TEXT_LEFT);
                nk_layout_row_dynamic(ctx, row_size * 10, 1);
                struct nk_colorf nkBackgroundColor = nk_color_cf(nk_rgb(game->backgroundColor.r, game->backgroundColor.g, game->backgroundColor.b));
                nk_color_pick(ctx, &nkBackgroundColor, NK_RGBA);
                game->backgroundColor = (color_t) {
                    (uint8_t) (nkBackgroundColor.r * 255.0f),
                    (uint8_t) (nkBackgroundColor.g * 255.0f),
                    (uint8_t) (nkBackgroundColor.b * 255.0f)
                };

                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Animation", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_property_float(ctx, "Rotation speed", -360.0f, &game->rotationSpeed, 360.0f, 0.1f, 0.1f);
                nk_tree_pop(ctx);
            }

        }
        nk_end(game->nuklearContext);
    }
}
#endif // DEBUGUI

void update(game_state_t* game) {
    DEBUG_PRINT("INFO: Update game state\n");
    updateCameraPosition(game);
    animateObjects(game);
}

void render(point_t p0, point_t p1, point_t p2, game_state_t* game) {
    DEBUG_PRINT("INFO: Rendering scene\n");
    canvas_t canvas = game->canvas;

    // Init depthBuffer
    for (int i = 0; i < canvas.width * canvas.height; i++) {
        canvas.depthBuffer[i] = 0.0;
    }

    DEBUG_PRINT("INFO: Drawing background\n");
    uint32_t backgroundColor = colorToUint32(game->backgroundColor);
    for (int i = 0; i < canvas.width * canvas.height; i++) {
        game->canvas.frameBuffer[i] = backgroundColor;
    }

    // Draw 3D Objects
    if (game->draw3DObjects) {
        DEBUG_PRINT("INFO: Drawing 3D Objects\n");
        drawObjects(game);
    }

    // Draw lights
    if (game->drawLights) {
        DEBUG_PRINT("INFO: Drawing lights\n");
        drawLights(game);
    }
    
    // Draw lines
    if (game->draw2DObjects) {
        DEBUG_PRINT("INFO: Drawing 2D Objects\n");
        if (game->drawWire) {
            DEBUG_PRINT("INFO: Drawing wireframe triangle\n");
            drawTriangleWireframe(p0.x, p1.x, p2.x,
                                  p0.y, p1.y, p2.y,
                                  COLOR_GREEN, game->canvas);
        }
        
        if (game->drawFilled == 1) {
            DEBUG_PRINT("INFO: Drawing triangle\n");
            
            int area = edgeCross(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
            drawTriangleFilled(p0.x, p1.x, p2.x,
                               p0.y, p1.y, p2.y,
                               p0.invz, p1.invz, p2.invz,
                               1.0, 1.0, 1.0,
                               (vec3_t) {0, 0, 0}, (vec3_t) {0, 0, 0}, (vec3_t) {0, 0, 0},
                               (vec3_t) {0, 0, 0}, (vec3_t) {0, 0, 0}, (vec3_t) {0, 0, 0},
                               COLOR_RED, COLOR_GREEN, COLOR_BLUE,
                               0.0,
                               NULL, 0, 0,
                               IDENTITY_M4x4,
                               area,
                               game);
        }
    }
    
    DEBUG_PRINT("INFO: Update backbuffer\n");
    SDL_UpdateTexture(game->texture, NULL, game->canvas.frameBuffer, PITCH);
    SDL_RenderCopy(game->renderer, game->texture, NULL, NULL);
}

#ifdef DEBUGUI
void renderDebugUI(game_state_t* game) {
    if (game->showGUI) {
        DEBUG_PRINT("INFO: Rendering GUI\n");
        nk_sdl_render(NK_ANTI_ALIASING_ON);
    }
}
#endif // DEBUGUI 

void destroy(game_state_t* game) {
    free(game->camera.planes);
    free(game->canvas.frameBuffer);
    free(game->canvas.depthBuffer);

    // Free all triangles, vertices and materials from meshes
    for (int i = 0; i < game->numMeshes; i++) {
        free(game->meshes[i].materials);
        free(game->meshes[i].triangles);
        free(game->meshes[i].vertices);
    } 
    
    free(game->meshes);
    free(game->objects);
    SDL_DestroyTexture(game->texture);
    SDL_DestroyWindow(game->window);

    #ifdef DEBUGUI
    nk_sdl_shutdown();
    #endif // DEBUGUI

    SDL_Quit();
}

int main(int argc, char* argv[])
{
    DEBUG_PRINT("INFO: Initializing game objects\n");
    // TODO: Store 2D objects in game state
    point_t p0 = {541, 199, 1.0f / 0.01f};
    point_t p1 = {613, 279, 1.0f / 0.01f};
    point_t p2 = {453, 399, 1.0f / 0.01f};
    
    game_state_t* game = init();

    while (game->running) {
        handleEvents(game);
        
        #ifdef DEBUGUI
        updateDebugUI(game);
        #endif // DEBUGUI

        update(game);
        render(p0, p1, p2, game);
        
        #ifdef DEBUGUI
        renderDebugUI(game);
        #endif // DEBUGUI

        // Present frame and compute FPS
        DEBUG_PRINT("INFO: Present frame\n");
        SDL_RenderPresent(game->renderer);
        uint64_t currentTime = SDL_GetPerformanceCounter();
        game->elapsedTime = 1000.0 * (currentTime - game->lastTime) / SDL_GetPerformanceFrequency();
        DEBUG_PRINT("INFO: Frame rendered in %.00lf ms (%.0f FPS)\n", game->elapsedTime, floor(1000.0f / game->elapsedTime));
        game->lastTime = currentTime;
    }

    DEBUG_PRINT("INFO: Closing\n");
    destroy(game);
    return 0;
}
