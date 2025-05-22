#include "vec4.h"
#include <math.h>
#include <string.h>
#include <zot.h>

#if defined(__ARM_NEON__) || defined(__AVX__) || defined(__SSE__)
#define use_simd_float32
#else
#define use_float32
#endif

#if defined __arm__ && defined __ARM_FP && !defined __LITTLE_ENDIAN__
#error Sorry but I am trying to finish this project and currently only support little endian arm because that is what I use.
#endif

// Core math functionalities and utilities
#if defined(use_bfp16) || defined(use_float16)
#include <arm_neon.h>

#define SIMD_STRIDE 8
#if defined(use_bfp16)
#define LOAD_SIMD vld1q_bf16
#define DUP_N_SIMD vdupq_n_bf16
#define STORE_SIMD vst1q_bf16

#define SIMD_add vaddq_bf16
#define SIMD_subtract vsubq_bf16
#define SIMD_additive_inverse vnegq_bf16

#elif defined(use_float16)
#define LOAD_SIMD vld1q_f16
#define DUP_N_SIMD vdupq_n_f16
#define STORE_SIMD vst1q_f16
#define SIMD_initial_reciprocal vrecpeq_f16
#define SIMD_correction_factor vrecpsq_f16
#define SIMD_type float16x4_t

#define SIMD_add vaddq_f16
#define SIMD_sum vaddvq_f16
#define SIMD_subtract vsubq_f16
#define SIMD_multiply vmulq_f16
#define SIMD_divide vdivq_f16
#endif

#elif defined(use_simd_float32)
#include <arm_neon.h>
// #include <neon2sse.h>

#define SIMD_STRIDE 4
#define LOAD_SIMD vld1q_f32
#define DUP_N_SIMD vdupq_n_f32
#define STORE_SIMD vst1q_f32

#define SIMD_add vaddq_f32
#define SIMD_subtract vsubq_f32
#define SIMD_additive_inverse vnegq_f32

#else
#define SIMD_type zfl
#define SIMD_STRIDE 1
#endif

#if !(__ARM_FP & 2)
/*No hardware floating point support */

#if defined(use_simd_float32) || defined(use_bfp16)
#define SIMD_type float32x4_t
#define SIMD_initial_reciprocal vrecpeq_f32
#define SIMD_correction_factor vrecpsq_f32
#define SIMD_get_high vget_high_f32
#define SIMD_get_low vget_low_f32
#define SIMD_get_lane vget_lane_f32

#define SIMD_divide vdivq_f32
#define SIMD_multiply_by_scalar vmulq_n_f32
#define SIMD_multiply vmulq_f32

#define SIMD_add_x2 vadd_f32
#define SIMD_padd_x2 vpadd_f32
#define SIMD_sum vaddvq_f32

#define SIMD_min vminq_f32
#define SIMD_min_x2 vmin_f32
#define SIMD_pmin_x2 vpmin_f32
#define SIMD_reduce_min vminvq_f32

#define SIMD_max vmaxq_f32
#define SIMD_max_x2 vmax_f32
#define SIMD_pmax_x2 vpmax_f32
#define SIMD_reduce_max vmaxvq_f32

#define SIMD_float_to_int_with_shift vcvtq_n_s32_f32
#define SIMD_float_to_int vcvtq_s32_f32
#endif

#ifndef use_float32
#define __ai static __inline__ __attribute__((__always_inline__, __nodebug__))

__ai __attribute__((target("neon"))) SIMD_type SIMD_divide(SIMD_type dividend,
                                                           SIMD_type divisor) {
  /*determine an initial estimate of reciprocal of divisor.*/
  auto initial_reciprocal = SIMD_initial_reciprocal(divisor);
  auto correction_factor = SIMD_correction_factor(divisor, initial_reciprocal);
  initial_reciprocal = SIMD_multiply(initial_reciprocal, correction_factor);
  correction_factor = SIMD_correction_factor(divisor, initial_reciprocal);
  initial_reciprocal = SIMD_multiply(initial_reciprocal, correction_factor);

  return SIMD_multiply(dividend, initial_reciprocal);
}
#endif

#if defined(use_float16)
#define float_type zfl
#else
#define float_type float
#endif

#ifndef use_float32

__ai __attribute__((target("neon"))) float_type SIMD_sum(SIMD_type a) {
  auto sum = SIMD_add_x2(SIMD_get_high(a), SIMD_get_high(a));
  sum = SIMD_padd_x2(sum, sum);

  return SIMD_get_lane(sum, 0);
}

__ai __attribute__((target("neon"))) float_type SIMD_reduce_min(SIMD_type a) {
  auto min = SIMD_min_x2(SIMD_get_high(a), SIMD_get_high(a));
  min = SIMD_pmin_x2(min, min);
  return SIMD_get_lane(min, 0);
}

__ai __attribute__((target("neon"))) float_type SIMD_reduce_max(SIMD_type a) {
  auto max = SIMD_max_x2(SIMD_get_high(a), SIMD_get_high(a));
  max = SIMD_pmin_x2(max, max);
  return SIMD_get_lane(max, 0);
}

__ai __attribute__((target("neon"))) SIMD_type SIMD_round(SIMD_type a) {
  // https://stackoverflow.com/a/69770515
  auto a_as_int = SIMD_float_to_int_with_shift(a, 1);
  auto arithmetic_shift_right =
      vsraq_n_s32(a_as_int, a_as_int, 31); // account for negative rounding
  auto floor_round = vrshrq_n_s32(arithmetic_shift_right, 1);
  return vcvtq_f32_s32(floor_round);
}

#define LN2_MUL_INV                                                            \
  1.442695040888963407359924681001892137426645954152985934135449406931109219181185079885526622893506344f
#define LN2                                                                    \
  0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875f

__ai __attribute__((target("neon"))) SIMD_type SIMD_exp(SIMD_type a) {
  auto a_over_ln2 = SIMD_multiply_by_scalar(a, LN2_MUL_INV);

  auto n = SIMD_round(a_over_ln2);

  auto r = SIMD_subtract(a, SIMD_multiply_by_scalar(n, LN2));
  auto r2 = SIMD_multiply(r, r);
  auto r3 = SIMD_multiply(r2, r);
  auto r4 = SIMD_multiply(r3, r);
  auto r5 = SIMD_multiply(r4, r);
  auto r6 = SIMD_multiply(r5, r);

  auto exp_r_poly = SIMD_add(
      DUP_N_SIMD(1.f),
      SIMD_add(
          r,
          SIMD_add(
              SIMD_multiply_by_scalar(r2, 0.5f),
              SIMD_add(
                  SIMD_multiply_by_scalar(
                      r3, 0.16666666666666666666666666666666666666666666667f),
                  SIMD_add(
                      SIMD_multiply_by_scalar(
                          r4, 0.0416666666666666666666666666666666666667f),
                      SIMD_add(
                          SIMD_multiply_by_scalar(
                              r5,
                              0.0083333333333333333333333333333333333333333333333333333333f),
                          SIMD_multiply_by_scalar(
                              r6,
                              0.001388888888888888888888888888888888888888888889)))))));

  // Convert n (float) to integer for bit-level manipulation
  auto n_int = SIMD_float_to_int(n);

  // Compute 2^n by constructing the floating-point exponent:
  // In IEEE 754 single precision, the exponent bias is 127.
  int32x4_t exp_int = vshlq_n_s32(vaddq_s32(n_int, vdupq_n_s32(127)), 23);
  float32x4_t pow2n = vreinterpretq_f32_s32(exp_int);

  // exp(x) ≈ 2^n * poly(r)
  return vmulq_f32(pow2n, exp_r_poly);
}

#endif

#if defined(use_bfp16)
#undef SIMD_divide
#define SIMD_divide vdivq_bf16

#undef SIMD_multiply
#define SIMD_multiply vmulq_bf16

#undef SIMD_type
#define SIMD_type bfloat16x8_t

#undef SIMD_sum
#define SIMD_sum vaddvq_bf16

#undef SIMD_min
#define SIMD_min vminq_bf16

#undef SIMD_max
#define SIMD_max vmaxq_bf16

#undef SIMD_reduce_min
#define SIMD_reduce_min vminvq_bf16

#undef SIMD_reduce_max
#define SIMD_reduce_max vmaxvq_bf16

__ai __attribute__((target("neon"))) SIMD_type simd_operate(
    SIMD_type a, SIMD_type b, float32x4_t(operator)(float32x4_t, float32x4_t)) {
  auto a_low = vcvtq_low_f32_bf16(a);
  auto a_high = vcvtq_high_f32_bf16(a);

  auto b_low = vcvtq_low_f32_bf16(b);
  auto b_high = vcvtq_high_f32_bf16(b);

  auto result_low = operator(a_low, b_low);
  auto result_high = operator(a_high, b_high);

  return vcombine_bf16(vcvt_bf16_f32(result_low), vcvt_bf16_f32(result_high));
}

__ai __attribute__((target("neon"))) SIMD_type SIMD_add(SIMD_type a,
                                                        SIMD_type b) {
  return simd_operate(a, b, vaddq_f32);
}
__ai __attribute__((target("neon"))) SIMD_type SIMD_subtract(SIMD_type a,
                                                             SIMD_type b) {
  return simd_operate(a, b, vsubq_f32);
}
__ai __attribute__((target("neon"))) SIMD_type SIMD_multiply(SIMD_type a,
                                                             SIMD_type b) {
  return simd_operate(a, b, vmulq_f32);
}
__ai __attribute__((target("neon"))) SIMD_type SIMD_divide(SIMD_type a,
                                                           SIMD_type b) {
  return simd_operate(a, b, vdivq_f32);
}
__ai __attribute__((target("neon"))) SIMD_type SIMD_min(SIMD_type a,
                                                        SIMD_type b) {
  return simd_operate(a, b, vminq_f32);
}
__ai __attribute__((target("neon"))) SIMD_type SIMD_max(SIMD_type a,
                                                        SIMD_type b) {
  return simd_operate(a, b, vmaxq_f32);
}

__ai __attribute__((target("neon"))) zfl SIMD_sum(SIMD_type a) {
  auto a_low = vcvtq_low_f32_bf16(a);
  auto a_high = vcvtq_high_f32_bf16(a);

  return (zfl)(vaddvq_f32(a_low) + vaddvq_f32(a_high));
}
__ai __attribute__((target("neon"))) zfl SIMD_reduce_max(SIMD_type a) {
  auto a_low = vcvtq_low_f32_bf16(a);
  auto a_high = vcvtq_high_f32_bf16(a);

  return (zfl)fmaxf(vmaxvq_f32(a_low), vmaxvq_f32(a_high));
}
__ai __attribute__((target("neon"))) zfl SIMD_reduce_min(SIMD_type a) {
  auto a_low = vcvtq_low_f32_bf16(a);
  auto a_high = vcvtq_high_f32_bf16(a);

  return (zfl)fminf(vminvq_f32(a_low), vminvq_f32(a_high));
}

#endif
#endif

#ifdef use_float32
#define SIMD_initial_reciprocal NULL
#define SIMD_correction_factor NULL
#define SIMD_get_high NULL
#define SIMD_get_low NULL
#define SIMD_get_lane NULL

#define SIMD_add NULL
#define SIMD_subtract NULL
#define SIMD_additive_inverse NULL

#define SIMD_divide NULL
#define SIMD_multiply NULL

#define SIMD_add_x2 NULL
#define SIMD_padd_x2 NULL
#define SIMD_sum NULL

#define SIMD_min NULL
#define SIMD_min_x2 NULL
#define SIMD_pmin_x2 NULL
#define SIMD_reduce_min NULL

#define SIMD_max NULL
#define SIMD_max_x2 NULL
#define SIMD_pmax_x2 NULL
#define SIMD_reduce_max NULL

#endif

struct v4 {
  float x, y, z, w;
};

struct m4 {
  vec4 x, y, z, t;
};

vec4 *make_vec4(vec4 *vec, float *n) {
  vec4 *v = vec ? vec : zmalloc(sizeof *v);
  if (n) {
    v->x = n[0];
    v->y = n[1];
    v->z = n[2];
    v->w = n[3];
  } else {
    memset(v, 0, sizeof(vec4));
  }

  return v;
}

vec4 *vadd(vec4 *vres, vec4 *va, vec4 *vb) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  void *v1 = va, *v2 = vb;
  auto a_simd_vector = LOAD_SIMD(((v1)));
  auto b_simd_vector = LOAD_SIMD(((v2)));
  auto simd_result = vaddq_f32(a_simd_vector, b_simd_vector);

  // vec4 *vres = zmalloc(sizeof(vec4));
  STORE_SIMD((void *)vres, simd_result);

  return vres;
}

vec4 *vsub(vec4 *vres, vec4 *va, vec4 *vb) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  void *v1 = va, *v2 = vb;
  auto a_simd_vector = LOAD_SIMD(((v1)));
  auto b_simd_vector = LOAD_SIMD(((v2)));
  auto simd_result = vsubq_f32(a_simd_vector, b_simd_vector);

  // vec4 *vres = zmalloc(sizeof(vec4));
  STORE_SIMD((void *)vres, simd_result);

  return vres;
}

vec4 *vmul(vec4 *vres, vec4 *va, vec4 *vb) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  void *v1 = va, *v2 = vb;
  auto a_simd_vector = LOAD_SIMD(((v1)));
  auto b_simd_vector = LOAD_SIMD(((v2)));
  auto simd_result = vmulq_f32(a_simd_vector, b_simd_vector);

  // vec4 *vres = zmalloc(sizeof(vec4));
  STORE_SIMD((void *)vres, simd_result);

  return vres;
}

vec4 *vscale(vec4 *vres, float s, vec4 *vb) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  void *v2 = vb;
  auto a_simd_vector = DUP_N_SIMD(s);
  auto b_simd_vector = LOAD_SIMD(((v2)));
  auto simd_result = vmulq_f32(a_simd_vector, b_simd_vector);

  // vec4 *vres = zmalloc(sizeof(vec4));
  STORE_SIMD((void *)vres, simd_result);

  return vres;
}

float vdot(float *fres, vec4 *va, vec4 *vb) {
  float f;
  fres = fres ? fres : &f;

  void *v1 = va, *v2 = vb;
  auto a_simd_vector = LOAD_SIMD(((v1)));
  auto b_simd_vector = LOAD_SIMD(((v2)));
  auto prod_result = vmulq_f32(a_simd_vector, b_simd_vector);

  float prodf[4];
  vst1q_f32(prodf, prod_result);

  return *fres = prodf[0] + prodf[1] + prodf[2] + prodf[3];
  //  vaddvq_f32(prod_result);
}

vec4 *vcross(vec4 *vres, vec4 *va, vec4 *vb) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  void *v1 = va, *v2 = vb;
  auto a_simd_vector = LOAD_SIMD(((v1))); // [a0, a1, a2, ?]
  auto b_simd_vector = LOAD_SIMD(((v2))); // [b0, b1, b2, ?]

  // Compute cross product components:
  // result[0] = a[1]*b[2] - a[2]*b[1]
  // result[1] = a[2]*b[0] - a[0]*b[2]
  // result[2] = a[0]*b[1] - a[1]*b[0]
  auto a_lo = vget_low_f32(a_simd_vector);  // [a0, a1]
  auto a_hi = vget_high_f32(a_simd_vector); // [a2, ?]
  auto b_lo = vget_low_f32(b_simd_vector);  // [b0, b1]
  auto b_hi = vget_high_f32(b_simd_vector); // [b2, ?]

  auto cross_xy = vsub_f32(vmul_f32(vrev64_f32(a_hi), b_lo), // [a2*b0, ?*b1]
                           vmul_f32(a_lo, vrev64_f32(b_hi))  // [a0*b2, a1*?]
  );

  auto cross_z = vsub_f32(vmul_f32(a_lo, vrev64_f32(b_lo)), // [a0*b1, a1*b0]
                          vmul_f32(vrev64_f32(a_lo), b_lo)  // [a1*b0, a0*b1]
  );

  // Combine results
  auto cross = vcombine_f32(vrev64_f32(cross_xy), // [result[1], result[0]]
                            cross_z               // [result[2], ?]
  );

  // vec4 *vres = zmalloc(sizeof(vec4));
  STORE_SIMD((void *)vres, cross);
  return vres;
}

float vmag(float *fres, vec4 *va) {
  float f;
  fres = fres ? fres : &f;

  void *v = va;
  auto v_simd_vector = LOAD_SIMD(((v)));
  auto v_sqr = vmulq_f32(v_simd_vector, v_simd_vector);

  float sqr_f[4];
  vst1q_f32(sqr_f, v_sqr);

  return *fres = sqrtf(sqr_f[0] + sqr_f[1] + sqr_f[2] + sqr_f[3]);
}

vec4 *vnorm(vec4 *vres, vec4 *va) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  void *v = va;
  auto v_simd_vector = LOAD_SIMD(((v)));
  float mag;
  auto mag_simd_vector = DUP_N_SIMD(((vmag(&mag, v))));
  auto simd_result = SIMD_divide(v_simd_vector, mag_simd_vector);

  // vec4 *vres = zmalloc(sizeof(vec4));
  STORE_SIMD((void *)vres, simd_result);

  return vres;
}

mat4 *make_mat4(mat4 *mres, float n[16]) {
  mres = mres ? mres : zcalloc(1, sizeof *mres);

  if (n) {
    memcpy(mres, n, 16 * sizeof *n);
  } else {
    mat4 m = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    memcpy(mres, &m, sizeof m);
  }
  return mres;
}

vec4 *vmdot(vec4 *vres, vec4 *v, mat4 *m) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  auto res_x = SIMD_add(
      SIMD_add(SIMD_multiply(DUP_N_SIMD(v->x), LOAD_SIMD((void *)&m->x)),
               SIMD_multiply(DUP_N_SIMD(v->y), LOAD_SIMD((void *)&m->y))),
      SIMD_add(SIMD_multiply(DUP_N_SIMD(v->z), LOAD_SIMD((void *)&m->z)),
               SIMD_multiply(DUP_N_SIMD(v->w), LOAD_SIMD((void *)&m->t))));
  STORE_SIMD(&vres->x, res_x);

  return vres;
}

vec4 *v3mdot(vec4 *vres, vec4 *v, mat4 *m) {
  float w = vres ? vres->w : 0;

  vres = vmdot(vres, v, m);
  vres->w = w;

  return vres;
}

vec4 *mv3dot(vec4 *vres, mat4 *m, vec4 *v) {
  vres = vres ? vres : zcalloc(1, sizeof(struct v4));

  vres->x = vdot(&vres->x, &m->x, v);
  vres->y = vdot(&vres->y, &m->y, v);
  vres->z = vdot(&vres->z, &m->z, v);

  return vres;
}

vec4 *mvdot(vec4 *vres, mat4 *m, vec4 *v) {
  vres = mv3dot(vres, m, v);
  vres->w = vdot(&vres->w, &m->t, v);

  return vres;
}

mat4 *mdot(mat4 *mres, mat4 *a, mat4 *b) {
  mres = mres ? mres : zcalloc(1, sizeof(struct m4));

  auto res_x = SIMD_add(
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->x.x), LOAD_SIMD((void *)&b->x)),
               SIMD_multiply(DUP_N_SIMD(a->x.y), LOAD_SIMD((void *)&b->y))),
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->x.z), LOAD_SIMD((void *)&b->z)),
               SIMD_multiply(DUP_N_SIMD(a->x.w), LOAD_SIMD((void *)&b->t))));
  STORE_SIMD((void *)&mres->x, res_x);

  auto res_y = SIMD_add(
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->y.x), LOAD_SIMD((void *)&b->x)),
               SIMD_multiply(DUP_N_SIMD(a->y.y), LOAD_SIMD((void *)&b->y))),
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->y.z), LOAD_SIMD((void *)&b->z)),
               SIMD_multiply(DUP_N_SIMD(a->y.w), LOAD_SIMD((void *)&b->t))));
  STORE_SIMD((void *)&mres->y, res_y);

  auto res_z = SIMD_add(
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->z.x), LOAD_SIMD((void *)&b->x)),
               SIMD_multiply(DUP_N_SIMD(a->z.y), LOAD_SIMD((void *)&b->y))),
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->z.z), LOAD_SIMD((void *)&b->z)),
               SIMD_multiply(DUP_N_SIMD(a->z.w), LOAD_SIMD((void *)&b->t))));
  STORE_SIMD((void *)&mres->z, res_z);

  auto res_t = SIMD_add(
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->t.x), LOAD_SIMD((void *)&b->x)),
               SIMD_multiply(DUP_N_SIMD(a->t.y), LOAD_SIMD((void *)&b->y))),
      SIMD_add(SIMD_multiply(DUP_N_SIMD(a->t.z), LOAD_SIMD((void *)&b->z)),
               SIMD_multiply(DUP_N_SIMD(a->t.w), LOAD_SIMD((void *)&b->t))));
  STORE_SIMD((void *)&mres->t, res_t);

  return mres;
}

mat4 *minv(mat4 *mres, mat4 *M) {
  mres = mres ? mres : zcalloc(1, sizeof(struct m4));

  float a = M->x.x, b = M->x.y, c = M->x.z, d = M->x.w; // row 1

  float e = M->y.x, f = M->y.y, g = M->y.z, h = M->y.w; // row 2

  float i = M->z.x, j = M->z.y, k = M->z.z, l = M->z.w; // row 3

  float m = M->t.x, n = M->t.y, o = M->t.z, p = M->t.w; // row 4

  float                      //
      kp_lo = k * p - l * o, //
      jp_ln = j * p - l * n, //

      jo_kn = j * o - k * n, //
      io_km = i * o - k * m, //

      ip_lm = i * p - l * m, //
      in_jm = i * n - j * m, //

      cp_do = c * p - d * o, //
      bo_cn = b * o - c * n, //

      bp_dn = b * p - d * n, //
      cl_dk = c * l - d * k, //

      ap_dm = a * p - d * m, //
      ao_cm = a * o - c * m, //

      al_di = a * l - d * i, //
      ak_ci = a * k - c * i, //

      an_bm = a * n - b * m, //
      aj_bi = a * j - b * i, //

      bl_dj = b * l - d * j, bk_cj = b * k - c * j;

  float D = (b * e * (-kp_lo) + c * e * (jp_ln) + d * e * (-jo_kn) +
             f * (a * (kp_lo) + c * (-ip_lm) + d * (io_km)) +
             g * (a * (-jp_ln) + b * (ip_lm) + d * (-in_jm)) +
             h * (a * (jo_kn) + b * (-io_km) + c * (in_jm)));

  if (fabsf(D) < 1e-8) {
    return NULL;
  }

  mat4 inv = {(f * (kp_lo) + g * (-jp_ln) + h * (jo_kn)) / D,
              -(b * (kp_lo) + c * (-jp_ln) + d * (jo_kn)) / D,
              (f * (-cp_do) + g * (bp_dn) + h * (-bo_cn)) / D,
              -(f * (-cl_dk) + g * (bl_dj) + h * (-bk_cj)) / D,

              -(e * (kp_lo) + g * (-ip_lm) + h * (io_km)) / D,
              (a * (kp_lo) + c * (-ip_lm) + d * (io_km)) / D,
              -(-c * e * p + d * e * o + g * (ap_dm) + h * (-ao_cm)) / D,
              (-c * e * l + d * e * k + g * (al_di) + h * (-ak_ci)) / D,

              (e * (jp_ln) + f * (-ip_lm) + h * (in_jm)) / D,
              -(a * (jp_ln) + b * (-ip_lm) + d * (in_jm)) / D,
              (-b * e * p + d * e * n + f * (ap_dm) + h * (-an_bm)) / D,
              -(-b * e * l + d * e * j + f * (al_di) + h * (-aj_bi)) / D,

              -(e * (jo_kn) + f * (-io_km) + g * (in_jm)) / D,
              (a * (jo_kn) + b * (-io_km) + c * (in_jm)) / D,
              -(-b * e * o + c * e * n + f * (ao_cm) + g * (-an_bm)) / D,
              (-b * e * k + c * e * j + f * (ak_ci) + g * (-aj_bi)) / D};

  memcpy(mres, &inv, sizeof(inv));

  return mres;
}

mat4 *mtranslate(mat4 *mres, float x, float y, float z) {
  mres = mres ? mres : zcalloc(1, sizeof(struct m4));

  mres->t.x += x;
  mres->t.y += y;
  mres->t.z += z;

  return mres;
}

mat4 *mrotate_x(mat4 *mres, float angle) {
  mres = mres ? mres : zcalloc(1, sizeof(struct m4));

  float a = angle * M_PI / 180.;
  mat4 rot = {
      {1, 0, 0, 0},              //
      {0, cosf(a), -sinf(a), 0}, //
      {0, sinf(a), cosf(a), 0},  //
      {0, 0, 0, 1},              //
  };
  return mdot(mres, mres, &rot);
}

mat4 *mrotate_y(mat4 *mres, float angle) {
  mres = mres ? mres : zcalloc(1, sizeof(struct m4));

  float a = angle * M_PI / 180.;
  mat4 rot = {
      {cosf(a), 0, sinf(a), 0},  //
      {0, 1, 0, 0},              //
      {-sinf(a), 0, cosf(a), 0}, //
      {0, 0, 0, 1},              //
  };
  return mdot(mres, mres, &rot);
}

mat4 *mrotate_z(mat4 *mres, float angle) {
  mres = mres ? mres : zcalloc(1, sizeof(struct m4));

  float a = angle * M_PI / 180.;
  mat4 rot = {
      {cosf(a), -sinf(a), 0, 0}, //
      {sinf(a), cosf(a), 0, 0},  //
      {0, 0, 1, 0},              //
      {0, 0, 0, 1},              //
  };
  return mdot(mres, mres, &rot);
}
