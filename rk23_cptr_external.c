#include <math.h>

typedef struct {
  int accepted;
  int rejected;
  int nfev;
} rk_stats;

typedef void (*rk23_rhs_fn)(int neqn, double t, const double *y, double *ydot, const void *ctx);

static double weighted_norm(const int n, const double *const y_old, const double *const y_next,
                            const double *const err_vec, const double *const a_tol,
                            const double r_tol) {
  double sum_sq = 0.0;

  for (int i = 0; i < n; ++i) {
    const double scale_i = a_tol[i] + fmax(fabs(y_old[i]), fabs(y_next[i])) * r_tol;
    const double scaled_err = err_vec[i] / scale_i;
    sum_sq += scaled_err * scaled_err;
  }

  return sqrt(sum_sq / (double)n);
}

static void rk23_step(const int n, const rk23_rhs_fn fun, const double t_cur, const double dt,
                      const double *const y_cur, const double *const k1, double *const k2,
                      double *const k3, double *const k4, double *const tmp,
                      const double *const a_tol, const double r_tol, const void *const ctx,
                      double *const y_next, double *const err_val) {
  for (int i = 0; i < n; ++i) {
    tmp[i] = y_cur[i] + dt * 0.5 * k1[i];
  }
  fun(n, t_cur + 0.5 * dt, tmp, k2, ctx);

  for (int i = 0; i < n; ++i) {
    tmp[i] = y_cur[i] + dt * 0.75 * k2[i];
  }
  fun(n, t_cur + 0.75 * dt, tmp, k3, ctx);

  for (int i = 0; i < n; ++i) {
    y_next[i] = y_cur[i] + dt * ((2.0 / 9.0) * k1[i] + (1.0 / 3.0) * k2[i] + (4.0 / 9.0) * k3[i]);
  }
  fun(n, t_cur + dt, y_next, k4, ctx);

  for (int i = 0; i < n; ++i) {
    tmp[i] = dt * ((-5.0 / 72.0) * k1[i] + (1.0 / 12.0) * k2[i] + (1.0 / 9.0) * k3[i] -
                   (1.0 / 8.0) * k4[i]);
  }
  *err_val = weighted_norm(n, y_cur, y_next, tmp, a_tol, r_tol);
}

void rk23_cptr_external(const int *const neqn, const rk23_rhs_fn fun, double *const t, double *const y,
                        const double *const tend, double *const h, const double *const atol,
                        const double *const rtol, double *const work, const void *const ctx,
                        int *const idid, rk_stats *const stats) {
  const int n = *neqn;
  double y_new[n];
  double err;
  double fac;
  double *const k1 = work;
  double *const k2 = work + n;
  double *const k3 = work + 2 * n;
  double *const k4 = work + 3 * n;
  double *const tmp = work + 4 * n;

  *idid = 0;
  stats->accepted = 0;
  stats->rejected = 0;
  stats->nfev = 0;

  fun(n, *t, y, k1, ctx);
  stats->nfev += 1;

  while (*t < *tend) {
    if (*t + *h > *tend) {
      *h = *tend - *t;
    }

    int step_rejected = 0;

    for (;;) {
      if (*h <= fabs(nextafter(*t, INFINITY) - *t)) {
        *idid = -1;
        return;
      }

      rk23_step(n, fun, *t, *h, y, k1, k2, k3, k4, tmp, atol, *rtol, ctx, y_new, &err);
      stats->nfev += 3;
      fac = 0.9 * pow(1.0 / fmax(err, 1.0e-10), 1.0 / 3.0);

      if (err < 1.0) {
        if (step_rejected) {
          fac = fmin(1.0, fac);
        }

        *t += *h;
        for (int i = 0; i < n; ++i) {
          y[i] = y_new[i];
          k1[i] = k4[i];
        }
        *h *= fmax(0.2, fmin(5.0, fac));
        stats->accepted += 1;
        break;
      }

      *h *= fmax(0.2, fmin(5.0, fac));
      step_rejected = 1;
      stats->rejected += 1;
    }
  }
}
