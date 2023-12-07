// #define DEBUG
#define DEBUGUI

// Standard
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

// External
#define SDL_MAIN_HANDLED
#include <SDL.h>
#define CGLTF_IMPLEMENTATION
#include "external/cgltf/cgltf.h"
#define STB_IMAGE_IMPLEMENTATION
#include "external/stb/stb_image.h"
#ifdef DEBUGUI
#include "external/imgui/imgui.h"
#include "external/imgui/imgui_impl_sdl2.h"
#include "external/imgui/imgui_impl_sdlrenderer.h"
#endif // DEBUGUI

// Local
#include "color.h"
#include "vectors.h"
#include "object3D.h"
#include "draw.h"

#ifdef DEBUG
#define DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
#define DEBUG_PRINT(...) do {} while (0)
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ROTATION_SPEED 15.0f // degrees per second


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
    uint32_t*     frameBuffer;
    float*        depthBuffer;
    color_t       backgroundColor;
    bool          drawLights;
    bool          draw3DObjects;
    bool          draw2DObjects;
    bool          drawWire;
    bool          drawFilled;
    bool          diffuseLighting;
    bool          specularLighting;
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

    // GLTF
    cgltf_data* gltfData;
    char*       gltfPath;

    // Camera
    camera_t     camera;
    
    // Animation
    float        rotationSpeed;

    // Input
    const uint8_t* keys;

    // GUI
    #ifdef DEBUGUI
    bool         showGUI;
    bool         toggleGUIKeyPressed;
    #endif // DEBUG
} game_state_t;


// TODO: Code here is a bit repeated between directional and point lights. Maybe refactor?
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

// TODO: Move to draw.cpp
int edgeCross(int ax, int ay, int bx, int by, int px, int py) {
  point_t ab = { bx - ax, by - ay };
  point_t ap = { px - ax, py - ay };
  return ab.x * ap.y - ab.y * ap.x;
}

// TODO: Move to draw.cpp
void drawTriangleWireframe(int x0, int x1, int x2,
                           int y0, int y1, int y2,
                           color_t color, uint32_t* frameBuffer) {
    drawLine(x0, x1, y0, y1, color, frameBuffer);
    drawLine(x1, x2, y1, y2, color, frameBuffer);
    drawLine(x2, x0, y2, y0, color, frameBuffer);
}

// TODO: Remove dependency with game_state_t and move to draw.cpp
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
    c0 = mulScalarColor(i0, c0);
    c1 = mulScalarColor(i1, c1);
    c2 = mulScalarColor(i2, c2);
    
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
                if (invz > game->depthBuffer[y * WIDTH + x]) {
                    uint32_t color = colorToUint32(c0); // Fallback in case of no texture and no shading
                    float light = alpha * i0 + beta * i1 + gamma * i2;

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
                            color_t color_tl = colorFromUint32RGBA(texture[floor_v * textureWidth + floor_u]);
                            color_t color_tr = colorFromUint32RGBA(texture[floor_v * textureWidth + next_u]);
                            color_t color_bl = colorFromUint32RGBA(texture[next_v * textureWidth + floor_u]);
                            color_t color_br = colorFromUint32RGBA(texture[next_v * textureWidth + next_u]);
                            color_t color_b = sumColors(mulScalarColor(1 - frac_u, color_bl), mulScalarColor(frac_u, color_br));
                            color_t color_tp = sumColors(mulScalarColor(1 - frac_u, color_tl), mulScalarColor(frac_u, color_tr));
                            color_typed = sumColors(mulScalarColor(1 - frac_v, color_b), mulScalarColor(frac_v, color_tp));
                        } else {
                            int tex_x = MIN(abs((int)((u_over_z/invz) * textureWidth)), textureWidth - 1);
                            int tex_y = MIN(abs((int)((v_over_z/invz) * textureHeight)), textureHeight - 1);
                            uint32_t color_uint32 = texture[tex_y * textureWidth + tex_x];
                            color_typed = colorFromUint32RGBA(color_uint32);
                        }
                        
                        color_t color_shaded = mulScalarColor(light, color_typed);
                        color = colorToUint32(color_shaded);
                    } else {
                        color_t color_typed = {
                            static_cast<uint8_t>(c0.r * alpha + c1.r * beta + c2.r * gamma),
                            static_cast<uint8_t>(c0.g * alpha + c1.g * beta + c2.g * gamma),
                            static_cast<uint8_t>(c0.b * alpha + c1.b * beta + c2.b * gamma)
                        };
                        color_t color_shaded = mulScalarColor(light, color_typed);
                        color = colorToUint32(color_shaded);
                    }

                    drawPixelDepthBuffer(x, y, invz, color, game->depthBuffer, game->frameBuffer);
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
        projected[i] = projectVertex(camTransformed[i]);
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
            color_t color = COLOR_WHITE;
            float specularExponent = 900.0f;

            if (mesh->numMaterials != 0) {
                color = materials[triangle.materialIndex].diffuseColor;
                specularExponent = materials[triangle.materialIndex].specularExponent;
            }

            // Draw wireframe
            if (game->drawWire) {
                drawTriangleWireframe(p0.x, p1.x, p2.x,
                                      p0.y, p1.y, p2.y,
                                      color,
                                      game->frameBuffer);
            }
            
            if (game->drawFilled) {
                // Gouraud shading
                float i0 = shadeVertex(transformed[triangle.v0], transformedNormals[triangle.n0], mesh->invMagnitudeNormals[triangle.n0], specularExponent, game);
                float i1 = shadeVertex(transformed[triangle.v1], transformedNormals[triangle.n1], mesh->invMagnitudeNormals[triangle.n1], specularExponent, game);
                float i2 = shadeVertex(transformed[triangle.v2], transformedNormals[triangle.n2], mesh->invMagnitudeNormals[triangle.n2], specularExponent, game);

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

uint32_t* loadTextureSTB(char* path, int* width, int* height) {
    int channels;
    uint8_t* texture = stbi_load(path, width, height, &channels, 4);
    if (texture == NULL) {
        fprintf(stderr, "ERROR: Couldn't load texture %s\n", path);
        exit(-1);
    }

    printf("DEBUG: Texture %s loaded with size %dx%d\n", path, *width, *height);
    // Print if the texture is 16 or 8 bit depth
    bool is_16_bit = stbi_is_16_bit(path);
    if (is_16_bit) {
        printf("DEBUG: Texture %s is 16 bit depth\n", path);
    } else {
        printf("DEBUG: Texture %s is 8 bit depth\n", path);
    }
    return (uint32_t*) texture;
}

float wrapOver(float input, float max) {
    float modResult = fmod(input, max);
    if (modResult == 0) {
        // Check if the number of wraps is even or odd
        int numWraps = input / max;
        return (numWraps % 2 == 0) ? 0 : max;
    }
    return modResult;
}

void renderPrimitive(cgltf_accessor* indicesAccessor, cgltf_accessor* positionAccessor, cgltf_accessor* normalAccessor, cgltf_accessor* texcoordAccessor, cgltf_material* material, game_state_t* game) {
    camera_t camera = game->camera;

    cgltf_size indicesCount = indicesAccessor->count;
    printf("DEBUG: Rendering %d triangles\n", (int) indicesCount / 3);
    for (int i=0; i<indicesCount; i+=3) {
        // TODO: remove this
        if (i/3 != 3) {
            // continue;
        }

        printf("DEBUG: Rendering triangle %d\n", i/3);

        uint16_t idx2 = cgltf_accessor_read_index(indicesAccessor, i);
        uint16_t idx1 = cgltf_accessor_read_index(indicesAccessor, i + 1);
        uint16_t idx0 = cgltf_accessor_read_index(indicesAccessor, i + 2);
        printf("DEBUG: Indexes %d, %d, %d\n", idx0, idx1, idx2);

        // Vertices
        vec3_t v0;
        vec3_t v1;
        vec3_t v2;
        bool res1 = cgltf_accessor_read_float(positionAccessor, idx0, &v0.x, 3);
        bool res2 = cgltf_accessor_read_float(positionAccessor, idx1, &v1.x, 3);
        bool res3 = cgltf_accessor_read_float(positionAccessor, idx2, &v2.x, 3);
        if (!res1 || !res2 || !res3) {
            fprintf(stderr, "ERROR: Couldn't read vertex data\n");
            exit(-1);
        }

        v0.z = -v0.z;
        v1.z = -v1.z;
        v2.z = -v2.z;

        printf("DEBUG: Vertices\n");
        printf("DEBUG: v0: %f, %f, %f\n", v0.x, v0.y, v0.z);
        printf("DEBUG: v1: %f, %f, %f\n", v1.x, v1.y, v1.z);
        printf("DEBUG: v2: %f, %f, %f\n", v2.x, v2.y, v2.z);
        
        // Normals
        vec3_t n0 = normalize(crossProduct(sub(v1, v0), sub(v2, v0)));;
        vec3_t n1 = n0;
        vec3_t n2 = n0;
        if (normalAccessor != NULL) {
            // TODO: Check type and component type of normal
            bool res1 = cgltf_accessor_read_float(normalAccessor, idx0, &n0.x, 3);
            bool res2 = cgltf_accessor_read_float(normalAccessor, idx1, &n1.x, 3);
            bool res3 = cgltf_accessor_read_float(normalAccessor, idx2, &n2.x, 3);
            if (!res1 || !res2 || !res3) {
                fprintf(stderr, "ERROR: Couldn't read normal data\n");
                exit(-1);
            }

            n0.z = -n0.z;
            n1.z = -n1.z;
            n2.z = -n2.z;
        }

        // START HERE:
        //    Texture coords are messed up in the sides of the cube
        printf("DEBUG: Getting textures\n");
        vec3_t t0 = {0};
        vec3_t t1 = {0};
        vec3_t t2 = {0};
        if (texcoordAccessor != NULL) {
            bool res1 = cgltf_accessor_read_float(texcoordAccessor, idx0, &t0.x, 2);
            bool res2 = cgltf_accessor_read_float(texcoordAccessor, idx1, &t1.x, 2);
            bool res3 = cgltf_accessor_read_float(texcoordAccessor, idx2, &t2.x, 2);

            if (!res1 || !res2 || !res3) {
                fprintf(stderr, "ERROR: Couldn't read texture data\n");
                exit(-1);
            }

            printf("DEBUG: Textures before transform\n");
            printf("DEBUG: t0: %f, %f\n", t0.x, t0.y);
            printf("DEBUG: t1: %f, %f\n", t1.x, t1.y);
            printf("DEBUG: t2: %f, %f\n", t2.x, t2.y);

            t0.x = wrapOver(t0.x, 1.0f);
            t1.x = wrapOver(t1.x, 1.0f);
            t2.x = wrapOver(t2.x, 1.0f);

            t0.y = wrapOver(t0.y, 1.0f);
            t1.y = wrapOver(t1.y, 1.0f);
            t2.y = wrapOver(t2.y, 1.0f);

            // t0.x = 1.0f - t0.x;
            // t1.x = 1.0f - t1.x;
            // t2.x = 1.0f - t2.x;

            // t0.y = 1.0f - t0.y;
            // t1.y = 1.0f - t1.y;
            // t2.y = 1.0f - t2.y;
        }

        printf("DEBUG: Textures\n");
        printf("DEBUG: t0: %f, %f\n", t0.x, t0.y);
        printf("DEBUG: t1: %f, %f\n", t1.x, t1.y);
        printf("DEBUG: t2: %f, %f\n", t2.x, t2.y);

        // Transform and project vertices
        vec3_t transformedV0 = v0; // mulMV3(object->transform, mesh->vertices[i]);
        vec3_t camTransformedV0 = mulMV3(camera.transform, transformedV0);
        point_t p0 = projectVertex(camTransformedV0);

        vec3_t transformedV1 = v1; // mulMV3(object->transform, mesh->vertices[i]);
        vec3_t camTransformedV1 = mulMV3(camera.transform, transformedV1);
        point_t p1 = projectVertex(camTransformedV1);

        vec3_t transformedV2 = v2; // mulMV3(object->transform, mesh->vertices[i]);
        vec3_t camTransformedV2 = mulMV3(camera.transform, transformedV2);
        point_t p2 = projectVertex(camTransformedV2);

        // Transform normals
        vec3_t transformedN0 = n0; // mulMV3(object->transform, mesh->vertices[i]);
        vec3_t transformedN1 = n1; // mulMV3(object->transform, mesh->vertices[i]);
        vec3_t transformedN2 = n2; // mulMV3(object->transform, mesh->vertices[i]);


        // Backface culling
        // printf("DEBUG: Backface culling\n");
        int area = edgeCross(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
        if (area < 0 && game->backfaceCulling) {
            printf("DEBUG: Clipped a triangle with backface culling\n");
            continue;
        }

        // Fustrum culling
        // TODO: Add config for this
        bool discarded = false;
        for (int p = 0; p < camera.numPlanes; p++) {
            plane_t plane = camera.planes[p];
            if (distancePlaneV3(plane, camTransformedV0) < 0 &&
                distancePlaneV3(plane, camTransformedV1) < 0 &&
                distancePlaneV3(plane, camTransformedV2) < 0) {
                    printf("DEBUG: Clipped triangle fully outside of the camera clipping volume\n");
                    discarded = true;
                    break;
            }

            // TODO: Deal with the case where only one or two vertexes of the triangle
            //       are out of the volume. In this case, we should split the triangle
            //       and create new ones.
        }

        if (discarded) {
            continue;
        }

        // Clip the triangle has any vertice behind the camera
        // printf("DEBUG: Camera culling\n");
        // TODO: Check if this is right, but it fixes overflows when we
        //       interpolate zs and we have z close to 0
        if (camTransformedV0.z < 0 ||
            camTransformedV1.z < 0 ||
            camTransformedV2.z < 0) {
            printf("DEBUG: Clipped triangle with a vertice behind the camera\n");
            continue;
        }

        // Material
        printf("DEBUG: Material\n");
        color_t color = COLOR_WHITE;
        if (material->has_pbr_metallic_roughness && material->pbr_metallic_roughness.base_color_factor != NULL) {
            color = colorFromFloats(material->pbr_metallic_roughness.base_color_factor[0],
                                    material->pbr_metallic_roughness.base_color_factor[1],
                                    material->pbr_metallic_roughness.base_color_factor[2]);
        }

        printf("Load texture\n");
        // Texture
        int textureWidth = 0;
        int textureHeight = 0;
        uint32_t* texture = NULL;
        if (material->has_pbr_metallic_roughness && material->pbr_metallic_roughness.base_color_texture.texture != NULL) {
            cgltf_image* image = material->pbr_metallic_roughness.base_color_texture.texture->image;
            // Load texture from PNG based on image uri using upng
            char *texturePath = (char*) malloc(strlen(game->gltfPath) + strlen(image->uri) + 1);
            if (texturePath == NULL) {
                fprintf(stderr, "ERROR: Texture path memory couldn't be allocated.\n");
                exit(-1);
            }
            strcpy(texturePath, game->gltfPath);
            strcat(texturePath, image->uri);

            printf("DEBUG: Loading texture from %s\n", texturePath);
            texture = loadTextureSTB(texturePath, &textureWidth, &textureHeight);
        }

        float specularExponent = 900.0f;

        // printf("DEBUG: Vertex shading\n");
        float i0 = shadeVertex(transformedV0, transformedN0, 1.0/magnitude(n0), specularExponent, game);
        float i1 = shadeVertex(transformedV1, transformedN1, 1.0/magnitude(n1), specularExponent, game);
        float i2 = shadeVertex(transformedV2, transformedN2, 1.0/magnitude(n2), specularExponent, game);

        drawTriangleWireframe(p0.x, p1.x, p2.x,
                              p0.y, p1.y, p2.y,
                              COLOR_WHITE,
                              game->frameBuffer);

        // Drawing
        // printf("DEBUG: Drawing triangle\n");
        drawTriangleFilled(p0.x, p1.x, p2.x,
                           p0.y, p1.y, p2.y,
                           p0.invz, p1.invz, p2.invz,
                           i0, i1, i2,
                           transformedN0, transformedN1, transformedN2,
                           t0, t1, t2,
                           color, color, color,
                           specularExponent,
                           texture, textureWidth, textureHeight,
                           inverseM4(camera.transform),
                           area,
                           game);
    }
}

void renderMesh(cgltf_mesh* mesh, game_state_t* game) {
    if (mesh == NULL) {
        return;
    }

    printf("DEBUG: Rendering mesh %s\n", mesh->name);
    for (int i = 0; i < mesh->primitives_count; i++) {
        cgltf_primitive* primitive = &mesh->primitives[i];
        printf("DEBUG: Rendering primitive %d with %d attributes\n", i, (int) primitive->attributes_count);
        cgltf_accessor* positionAccessor = NULL;
        cgltf_accessor* normalAccessor = NULL;
        cgltf_accessor* texcoordAccessor = NULL;
        cgltf_accessor* indicesAccessor = primitive->indices;
        cgltf_material* material = primitive->material;
        for (int j = 0; j < primitive->attributes_count; j++) {
            printf("DEBUG: Getting attribute %d accessor\n", j);
            cgltf_attribute* attribute = &primitive->attributes[j];
            if (attribute->type == cgltf_attribute_type_position) {
                positionAccessor = attribute->data;
            } else if (attribute->type == cgltf_attribute_type_normal) {
                normalAccessor = attribute->data;
            } else if (attribute->type == cgltf_attribute_type_texcoord) {
                texcoordAccessor = attribute->data;
            }
        }
        renderPrimitive(indicesAccessor, positionAccessor, normalAccessor, texcoordAccessor, material, game);
    }
}

void renderScene(game_state_t* game) {
    printf("DEBUG: Rendering GLTF scene\n");
    cgltf_scene* scene = game->gltfData->scene;
    // Render each node in the scene
    for (int i = 0; i < scene->nodes_count; i++) {
        cgltf_node* node = scene->nodes[i];
        printf("DEBUG: Rendering node %s\n", node->name);
        renderMesh(node->mesh, game);

        for (int i = 0; i < node->children_count; i++) {
            cgltf_node* child = node->children[i];
            renderMesh(child->mesh, game);
        }
    }
}

// TODO: Hide lights?
void drawLights(game_state_t* game) {
    for (int i = 0; i < game->numPointLights; i++) {
        drawObject(&game->pointLightObjects[i], game);
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

    objects[0] = makeObject(&meshes[2], {0, 0, 0}, 1.0 , IDENTITY_M4x4);

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

    ambientLights[0] = {0.4};
    directionalLights[0] = {0.0, {0.0, -1.0, 1.0}};
    pointLights[0] = {0.9, {-0.5, 1.5, -2.0}};


    pointLightObjects[0] = makeObject(&meshes[0], pointLights[0].position, 0.05, IDENTITY_M4x4);

    DEBUG_PRINT("INFO: Loading GLTF scene\n");
    cgltf_data* gltfData = NULL;
    char *gltfPath = "assets/scenes/boxtextured/";
    char *gltfFilename = "BoxTextured.gltf";
    char *gltfFullPath = (char*) malloc(strlen(gltfPath) + strlen(gltfFilename) + 1);
    if (gltfFullPath == NULL) {
        fprintf(stderr, "ERROR: GLTF scene memory couldn't be allocated.\n");
        exit(-1);
    }
    strcpy(gltfFullPath, gltfPath);
    strcat(gltfFullPath, gltfFilename);

    cgltf_options options = {0};
    cgltf_result result = cgltf_parse_file(&options, gltfFullPath, &gltfData);
    if (result != cgltf_result_success) {
        fprintf(stderr, "ERROR: Couldn't load GLTF scene\n");
        exit(-1);
    }
    result = cgltf_load_buffers(&options, gltfData, gltfPath);
    if (result != cgltf_result_success) {
        fprintf(stderr, "ERROR: Couldn't load GLTF scene buffers\n");
        exit(-1);
    }
    result = cgltf_validate(gltfData);
    if (result != cgltf_result_success) {
        fprintf(stderr, "ERROR: Couldn't validate GLTF scene\n");
        exit(-1);
    }

    // Print some GLTF info
    printf("INFO: GLTF scene loaded\n");
    printf("INFO: GLTF file copywright: %s\n", gltfData->asset.copyright);
    printf("INFO: GLTF file has %d scenes\n", (int) gltfData->scenes_count);
    printf("INFO: GLTF default scene is scene number %d (%s)\n", (int) (gltfData->scene - gltfData->scenes), gltfData->scene->name);
    printf("INFO: GLTF file has %d nodes\n", (int) gltfData->nodes_count);
    for (int i = 0; i < gltfData->nodes_count; i++) {
        printf("INFO: GLTF node %d (%s) has %d children\n", i, gltfData->nodes[i].name, (int) gltfData->nodes[i].children_count);
    }
    printf("INFO: GLTF file has %d meshes\n", (int) gltfData->meshes_count);
    printf("INFO: GLTF file has %d materials\n", (int) gltfData->materials_count);
    printf("INFO: GLTF file has %d cameras\n", (int) gltfData->cameras_count);
    printf("INFO: GLTF file has %d lights\n", (int) gltfData->lights_count);
    printf("INFO: GLTF file has %d images\n", (int) gltfData->images_count);
    printf("INFO: GLTF file has %d textures\n", (int) gltfData->textures_count);
    printf("INFO: GLTF file has %d buffers\n", (int) gltfData->buffers_count);
    printf("INFO: GLTF file has %d samplers\n", (int) gltfData->samplers_count);
    printf("INFO: GLTF file has %d skins\n", (int) gltfData->skins_count);
    printf("INFO: GLTF file has %d animations\n", (int) gltfData->animations_count);
    printf("INFO: GLTF file has %d accessors\n", (int) gltfData->accessors_count);
    printf("INFO: GLTF file has %d buffer views\n", (int) gltfData->buffer_views_count);


    DEBUG_PRINT("INFO: Initializing game state\n");
    uint32_t *frameBuffer = (uint32_t*) malloc(WIDTH * HEIGHT * sizeof(uint32_t));
    float *depthBuffer = (float*) malloc(WIDTH * HEIGHT * sizeof(float));
    game_state_t* game = (game_state_t*) malloc(sizeof(game_state_t));
    if (game == NULL || frameBuffer == NULL || depthBuffer == NULL) {
        fprintf(stderr, "ERROR: Game state memory and buffers couldn't be allocated.\n");
        exit(-1);
    }
    
    game->running = true;
    game->elapsedTime = 0;
    game->lastTime = SDL_GetPerformanceCounter();
    game->window = SDL_CreateWindow("Rasterizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
    game->renderer = SDL_CreateRenderer(game->window, -1, 0);
    game->texture = SDL_CreateTexture(game->renderer, SDL_PIXELFORMAT_BGRA32, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
    game->frameBuffer = frameBuffer;
    game->depthBuffer = depthBuffer;
    game->backgroundColor = {0, 0, 0};
    game->drawLights    = true;
    game->draw3DObjects =  false; // TODO: deal with this
    game->draw2DObjects =  false;
    game->diffuseLighting = true;
    game->specularLighting = true;
    game->drawWire = false;
    game->drawFilled = true;
    game->backfaceCulling = true;
    game->bilinearFiltering = false;
    game->numMeshes = numMeshes;
    game->meshes = meshes;
    game->numObjects = numObjects;
    game->objects = objects;

    // GLTF
    game->gltfData = gltfData;
    game->gltfPath = gltfPath;
    
    // Lights
    game->numAmbientLights = numAmbientLights;
    game->ambientLights = ambientLights;
    game->numDirectionalLights = numDirLights;
    game->directionalLights = directionalLights;
    game->numPointLights = numPointLights;
    game->pointLights = pointLights;
    game->pointLightObjects = pointLightObjects;
    
    game->camera = makeCamera(
        {0, 0, -5},
        rotationY(0.0f),
        VIEWPORT_DISTANCE,
        5.0f,
        90.0f
    );
    game->rotationSpeed = 15.0f;
    game->keys = SDL_GetKeyboardState(NULL);

    #ifdef DEBUGUI
    DEBUG_PRINT("INFO:  Initializing Dear ImGui\n");
    ImGui::CreateContext();
    ImGui_ImplSDL2_InitForSDLRenderer(game->window, game->renderer);
    ImGui_ImplSDLRenderer_Init(game->renderer);
    game->showGUI = true;
    game->toggleGUIKeyPressed = false;
    #endif // DEBUGUI

    return game;
}

void handleEvents(game_state_t* game) {
    DEBUG_PRINT("INFO: Handle events\n");    
    while (SDL_PollEvent(&game->event)) {
        #ifdef DEBUGUI
        ImGui_ImplSDL2_ProcessEvent(&game->event);
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
        ImGui_ImplSDLRenderer_NewFrame();
        ImGui_ImplSDL2_NewFrame(game->window);
        ImGui::NewFrame();

        // Create an ImGui interface
        ImGui::Begin("Settings");

        ImGui::Text("FPS: %.2f (%.3lf ms)", floor(1000.0f / game->elapsedTime), game->elapsedTime);

        if (ImGui::CollapsingHeader("Meshes")) {
            for (int i = 0; i < game->numMeshes; i++) {
                ImGui::Text("%s: (vertices %d, triangles: %d, uvs: %d)", game->meshes[i].name, game->meshes[i].numVertices, game->meshes[i].numTriangles, game->meshes[i].numTextureCoords);
            }
        }

        if (ImGui::CollapsingHeader("Objects")) {
            for (int i = 0; i < game->numObjects; i++) {
                ImGui::Text("%d: %s at (%.0f, %.0f, %.0f)", i, game->objects[i].mesh->name, game->objects[i].translation.x, game->objects[i].translation.y, game->objects[i].translation.z);
                ImGui::SliderFloat("x", &game->objects[i].translation.x, -10.0f, 10.0f);
                ImGui::SliderFloat("y", &game->objects[i].translation.y, -10.0f, 10.0f);
                ImGui::SliderFloat("z", &game->objects[i].translation.z, -10.0f, 10.0f);
                ImGui::SliderFloat("scale", &game->objects[i].scale, -0.01f, 10.0f);
            }
        }

        if (ImGui::CollapsingHeader("Lights")) {
            ImGui::Indent(20.0f);

            ImGui::Checkbox("Difuse", &game->diffuseLighting);
            ImGui::SameLine();
            ImGui::Checkbox("Specular", &game->specularLighting);


            if (ImGui::CollapsingHeader("Ambient")) {
                ImGui::Indent(20.0f);
                for (int i = 0; i < game->numAmbientLights; i++) {
                    ImGui::PushID(i);
                    ImGui::Text("Ambient Light %d", i);
                    ImGui::SliderFloat("Intensity##amb", &game->ambientLights[i].intensity, 0.0f, 3.0f);
                    ImGui::PopID();
                }
                ImGui::Unindent(20.0f);
            }

            if (ImGui::CollapsingHeader("Directional")) {
                ImGui::Indent(20.0f);
                for (int i = 0; i < game->numDirectionalLights; i++) {
                    ImGui::PushID(i);
                    ImGui::Text("Directional Light %d", i);
                    ImGui::SliderFloat("Intensity", &game->directionalLights[i].intensity, 0.0f, 3.0f);
                    ImGui::SliderFloat("x", &game->directionalLights[i].direction.x, -1.0f, 1.0f);
                    ImGui::SliderFloat("y", &game->directionalLights[i].direction.y, -1.0f, 1.0f);
                    ImGui::SliderFloat("z", &game->directionalLights[i].direction.z, -1.0f, 1.0f);
                    ImGui::PopID();
                }
                ImGui::Unindent(20.0f);
            }

            if (ImGui::CollapsingHeader("Point")) {
                ImGui::Indent(20.0f);
                for (int i = 0; i < game->numPointLights; i++) {
                    ImGui::PushID(i);
                    ImGui::Text("Point Light %d", i);
                    ImGui::SliderFloat("Intensity##p", &game->pointLights[i].intensity, 0.0f, 3.0f);
                    ImGui::SliderFloat("x##p", &game->pointLights[i].position.x, -10.0f, 10.0f);
                    ImGui::SliderFloat("y##p", &game->pointLights[i].position.y, -10.0f, 10.0f);
                    ImGui::SliderFloat("z##p", &game->pointLights[i].position.z, -10.0f, 10.0f);
                    ImGui::PopID();
                }
                ImGui::Unindent(20.0f);
            }
            ImGui::Unindent(20.0f);
        }
        
        if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("Camera: (%.1f, %.1f, %.1f)", game->camera.translation.x, game->camera.translation.y, game->camera.translation.z);
            ImGui::Text("What to draw");
            ImGui::Checkbox("3D Obj", &game->draw3DObjects);
            ImGui::SameLine();
            ImGui::Checkbox("2D Obj", &game->draw2DObjects);
            ImGui::SameLine();
            ImGui::Checkbox("Light Bulbs", &game->drawLights);
        }

        if (ImGui::CollapsingHeader("Render Options", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Wireframe", &game->drawWire);
            ImGui::SameLine();
            ImGui::Checkbox("Filled", &game->drawFilled);
            ImGui::Checkbox("Culling", &game->backfaceCulling);
            ImGui::SameLine();
            ImGui::Checkbox("Bilinear Filt.", &game->bilinearFiltering);
        }
        
        
        if (ImGui::CollapsingHeader("Background")) { 
            ImVec4 imBackgroundColor = ImVec4(game->backgroundColor.r/255.0f,
                                                game->backgroundColor.g/255.0f,
                                                game->backgroundColor.b/255.0f,
                                                1.00f);
            ImGui::ColorPicker4(
                "Color",
                (float*)&imBackgroundColor,
                ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_NoSidePreview
            );
            game->backgroundColor = {
                static_cast<uint8_t>(imBackgroundColor.x * 255.0f),
                static_cast<uint8_t>(imBackgroundColor.y * 255.0f),
                static_cast<uint8_t>(imBackgroundColor.z * 255.0f)
            };
        }

        if (ImGui::CollapsingHeader("Animation")) {
            ImGui::SliderFloat("Rotation Speed", &game->rotationSpeed, -360.0f, 360.0f, "%.2f degrees/s");
        }
        
        ImGui::End();
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

    // Init depthBuffer
    memset(game->depthBuffer, 0, WIDTH * HEIGHT * sizeof(float));

    // Background
    DEBUG_PRINT("INFO: Drawing background\n");
    uint32_t backgroundColor = colorToUint32(game->backgroundColor);
    memset(game->frameBuffer, backgroundColor, WIDTH * HEIGHT * sizeof(uint32_t));

    // Draw 3D Objects
    if (game->draw3DObjects) {
        DEBUG_PRINT("INFO: Drawing 3D Objects\n");
        drawObjects(game);
    }

    if (game->gltfData != NULL) {
        DEBUG_PRINT("INFO: Drawing GLTF scene\n");
        renderScene(game);
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
                                  COLOR_GREEN, game->frameBuffer);
        }
        
        if (game->drawFilled == 1) {
            DEBUG_PRINT("INFO: Drawing triangle\n");
            
            int area = edgeCross(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
            drawTriangleFilled(p0.x, p1.x, p2.x,
                               p0.y, p1.y, p2.y,
                               p0.invz, p1.invz, p2.invz,
                               1.0, 1.0, 1.0,
                               {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
                               {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
                               COLOR_RED, COLOR_GREEN, COLOR_BLUE,
                               0.0,
                               NULL, 0, 0,
                               IDENTITY_M4x4,
                               area,
                               game);
        }
    }
    
    DEBUG_PRINT("INFO: Update backbuffer\n");
    SDL_UpdateTexture(game->texture, NULL, (void *) game->frameBuffer, PITCH);
    SDL_RenderCopy(game->renderer, game->texture, NULL, NULL);
}

#ifdef DEBUGUI
void renderDebugUI(game_state_t* game) {
    if (game->showGUI) {
        DEBUG_PRINT("INFO: Rendering GUI\n");
        ImGui::Render(); 
        ImGui_ImplSDLRenderer_RenderDrawData(ImGui::GetDrawData());
    }
}
#endif // DEBUGUI 

void destroy(game_state_t* game) {
    free(game->camera.planes);
    free(game->frameBuffer);
    free(game->depthBuffer);

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
    ImGui_ImplSDLRenderer_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();
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