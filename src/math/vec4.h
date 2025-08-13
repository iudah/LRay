#ifndef VEC4_H
#define VEC4_H

#define EPSILON (1e-4f)

typedef struct v4 vec4;
typedef struct m4 mat4;

vec4 *make_vec4(vec4 *v, float[4]);
vec4 *vadd(vec4 *vres, vec4 *va, vec4 *vb);
vec4 *vsub(vec4 *vres, vec4 *va, vec4 *vb);
vec4 *vmul(vec4 *vres, vec4 *va, vec4 *vb);
vec4 *vscale(vec4 *vres, float s, vec4 *vb);
float vdot(float *fres, vec4 *va, vec4 *vb);
vec4 *vnorm(vec4 *res, vec4 *v);
float vmag(float *fres, vec4 *va);

mat4 *make_mat4(mat4 *m, float[16]);
mat4 *mtranslate(mat4 *m, float x, float y, float z);
mat4 *mrotate_x(mat4 *m, float angle_deg);
mat4 *mrotate_y(mat4 *m, float angle_deg);
mat4 *mrotate_z(mat4 *m, float angle_deg);

vec4 *vmdot(vec4 *vres, vec4 *v, mat4 *m);
vec4 *v3mdot(vec4 *vres, vec4 *v, mat4 *m);
vec4 *mv3dot(vec4 *vres, mat4 *m, vec4 *v);
vec4 *mvdot(vec4 *vres, mat4 *m, vec4 *v);
mat4 *mdot(mat4 *mres, mat4 *a, mat4 *b);
mat4 *minv(mat4 *mres, mat4 *M);

#endif