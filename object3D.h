#ifndef OBJECT_3D_H
#define OBJECT_3D_H

#include "color.h"
#include "vectors.h"

#define M_PI 3.14159265358979323846264338327950288

typedef struct triangle_t {
  int     v0, v1, v2;
  int     t0, t1, t2;
  int     n0, n1, n2;
  int     materialIndex;
} triangle_t;

typedef struct material_t {
    char*   name;
    color_t diffuseColor;
    color_t specularColor;
    float   specularExponent;
    int     textureWidth;
    int     textureHeight;
    uint32_t* texture;
} material_t;

typedef struct mesh_t {
    char*       name;
    int         numVertices;
    vec3_t*     vertices;
    int         numTextureCoords;
    vec3_t*     textureCoords;
    int         numNormals;
    vec3_t*     normals;
    float*      invMagnitudeNormals;
    int         numTriangles;
    triangle_t* triangles;
    int         numMaterials;
    material_t* materials;
    vec3_t      center;
    float       boundsRadius;
} mesh_t;

typedef struct object3D_t {
    mesh_t*  mesh;
    vec3_t   translation;
    float    scale;
    mat4x4_t rotation;
    mat4x4_t transform;
} object3D_t;

typedef struct camera_t {
    vec3_t translation;
    mat4x4_t rotation;
    mat4x4_t transform;
    int      numPlanes;
    plane_t* planes;
    float    viewportDistance;
    float    movementSpeed;
    float    turningSpeed;
} camera_t;

typedef struct ambient_light_t {
    float intensity;
} ambient_light_t;

typedef struct dir_light_t {
    float intensity;
    vec3_t direction;
} dir_light_t;

typedef struct point_light_t {
    float intensity;
    vec3_t position;
} point_light_t;

vec3_t triangleNormal(vec3_t v0, vec3_t v1, vec3_t v2);
vec3_t triangleCenter(vec3_t v0, vec3_t v1, vec3_t v2);
mat4x4_t translationToMatrix(vec3_t vector);
mat4x4_t scaleToMatrix(float scale);
mat4x4_t rotationY(float degrees);
mat4x4_t rotationX(float degrees);
object3D_t makeObject(mesh_t *mesh, vec3_t translation, float scale, mat4x4_t rotation);
camera_t makeCamera(vec3_t translation, mat4x4_t rotation,
                    float viewportDist, float movSpeed, float turnSpeed);
vec3_t meshCenter(vec3_t* vertices, int numVertices);
float meshBoundsRadius(vec3_t* vertices, int numVertices, vec3_t center);
mesh_t* loadObjFile(const char* filename, bool flipTextureVertically);
uint32_t* loadTexture(char* filename, int* textureWidth, int* textureHeight);

#endif