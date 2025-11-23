#include "vec4.h"
#include <stdlib.h>
#include "vec4_shared.h"
#include "quaternion.h"

quaternion *make_quaternion(quaternion *q, float f[4])
{
    float s[] = {f[1], f[2], f[3], f[0]};
    return make_vec4(q, s);
}
quaternion *qadd(quaternion *qres, quaternion *qa, quaternion *qb)
{
    return vadd(qres, qa, qb);
}
quaternion *qsub(vec4 *qres, vec4 *qa, vec4 *qb)
{
    return vsub(qres, qa, qb);
}
quaternion *qmul(vec4 *qres, vec4 *qa, vec4 *qb)
{
    float sa = qa->w;
    float sb = qb->w;
    vec4 va;
    vec4 vb;
    float fint;
    make_vec4(&va, (float[]){qa->x, qa->y, qa->z, 0});
    make_vec4(&vb, (float[]){qb->x, qb->y, qb->z, 0});
    float s = sa * sb - vdot(NULL, &va, &vb);

    vec4 vxv, vs, sv, vtemp, vfinal;
    vcross(&vxv, &va, &vb);
    vscale(&vs, sa, &vb);
    vscale(&sv, sb, &va);
    vadd(&vtemp, &vs, &sv);
    vadd(&vfinal, &vxv, &vtemp);

    return make_quaternion(qres, (float[]){s, vfinal.x, vfinal.y, vfinal.z});
}
quaternion *qscale(vec4 *qres, float s, vec4 *qb)
{
    return vscale(qres, s, qb);
}
