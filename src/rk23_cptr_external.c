#include <stdbool.h>
#include <math.h>

typedef struct {
  int accepted;
  int rejected;
  int nfev;
} rk_stats;

#define CTRL_I 1
#define CTRL_PI 2

typedef void (*rk23_rhs_fn)(int neqn, double t, const double y[],
                            double ydot[], void *ctx);

static double weighted_norm(int n, const double y_old[restrict static n],
                            const double y_next[restrict static n],
                            const double err_vec[restrict static n],
                            const double a_tol[restrict static n], double r_tol) {
  double sum_sq = 0.0;

  for (int i = 0; i < n; ++i) {
    const double scale_i = a_tol[i] + fmax(fabs(y_old[i]), fabs(y_next[i])) * r_tol;
    const double scaled_err = err_vec[i] / scale_i;
    sum_sq += scaled_err * scaled_err;
  }

  return sqrt(sum_sq / (double)n);
}

static void rk23_step(int n, rk23_rhs_fn fun, double t_cur, double dt,
                      const double y_cur[restrict static n],
                      const double k1[restrict static n], double k2[restrict static n],
                      double k3[restrict static n], double k4[restrict static n],
                      double tmp[restrict static n],
                      const double a_tol[restrict static n], double r_tol,
                      void *ctx, double y_next[restrict static n], double *err_val) {
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

void rk23_cptr_external(int neqn, rk23_rhs_fn fun, double *t,
                        double y[restrict static neqn], double tend, double *h,
                        const double atol[restrict static neqn], double rtol,
                        double work[restrict static 5 * neqn], void *ctx, int *idid,
                        rk_stats *stats, int controller_kind) {
  int n = neqn;
  double t_cur = *t;
  double h_cur = *h;
  double y_new[n];
  double err;
  double fac;
  double err_prev = 1.0;
  bool have_prev = false;
  double *restrict k1 = &work[0 * n];
  double *restrict k2 = &work[1 * n];
  double *restrict k3 = &work[2 * n];
  double *restrict k4 = &work[3 * n];
  double *restrict tmp = &work[4 * n];

  *idid = 0;
  *stats = (rk_stats){};

  fun(n, t_cur, y, k1, ctx);
  stats->nfev += 1;

  while (t_cur < tend) {
    if (t_cur + h_cur > tend) {
      h_cur = tend - t_cur;
    }

    bool step_rejected = false;

    do {
      const double ulp = nextafter(t_cur, INFINITY) - t_cur;
      if (fabs(h_cur) <= 10.0 * ulp) {
        *idid = -1;
        goto finish;
      }

      rk23_step(n, fun, t_cur, h_cur, y, k1, k2, k3, k4, tmp, atol, rtol, ctx, y_new, &err);
      stats->nfev += 3;
      const double err_now = fmax(err, 1.0e-10);
      if (controller_kind == CTRL_PI && have_prev) {
        fac = 0.9 * pow(err_now, -0.7 / 3.0) * pow(fmax(err_prev, 1.0e-10), -0.4 / 3.0);
      } else {
        fac = 0.9 * pow(err_now, -1.0 / 3.0);
      }

      if (err < 1.0) {
        if (step_rejected) {
          fac = fmin(1.0, fac);
        }

        t_cur += h_cur;
        for (int i = 0; i < n; ++i) {
          y[i] = y_new[i];
          k1[i] = k4[i];
        }
        h_cur *= fmax(0.2, fmin(5.0, fac));
        if (controller_kind == CTRL_PI) {
          err_prev = err_now;
          have_prev = true;
        }
        stats->accepted += 1;
        break;
      }

      h_cur *= fmax(0.2, fmin(5.0, fac));
      step_rejected = true;
      if (controller_kind == CTRL_PI) {
        have_prev = false;
      }
      stats->rejected += 1;
    } while (true);
  }

finish:
  *t = t_cur;
  *h = h_cur;
}
