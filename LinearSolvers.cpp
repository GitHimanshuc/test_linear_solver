#include "LinearSolvers.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "VectorOps.hpp"

namespace {

// Helper functions for GMRES
void GeneratePlaneRotation(const double dx,
                           const double dy,
                           double &cs_var,
                           double &sn_var) {
  if (dy == 0.0) {
    cs_var = 1.0;
    sn_var = 0.0;
  } else if (std::abs(dy) > std::abs(dx)) {
    const double temp = dx / dy;
    sn_var = 1.0 / std::sqrt(1.0 + temp*temp);
    cs_var = temp * sn_var;
  } else {
    const double temp = dy / dx;
    cs_var = 1.0 / std::sqrt(1.0 + temp*temp);
    sn_var = temp * cs_var;
  }
}

void ApplyPlaneRotation(double &dx,
                        double &dy,
                        const double cs_var,
                        const double sn_var) {
  const double temp = cs_var * dx + sn_var * dy;
  dy = -sn_var * dx + cs_var * dy;
  dx = temp;
}

// Update projects x_vec out of the Krylov subspace so that we can form a
// new one and restart the iterations there.
// This function will fail if we use size_t instead of int
// because of the i-- and j-- in the for loops.
void Update(std::vector<double> &x_vec,
            const int k,
            const std::vector<double> &h,
            const std::vector<double> &s_vec,
            const std::vector<std::vector<double>> &v_vec,
            const int m) {
  std::vector<double> y = s_vec;

  for (int i = k; i >= 0; i--) {
    y[i] /= h[m * i + i];
    for (int j = i - 1; j >= 0; j--) {
      y[j] -= h[m * j + i] * y[i];
    }
  }

  for (int j = 0; j <= k; j++) {
    x_vec = Sum(x_vec, ScalarMultiplication(y[j], v_vec[j]));
  }
}

}  // namespace

std::vector<double> linear_solvers::GMRES(
    const std::vector<double> &rhs_vector,
    const std::vector<double> &initial_guess,
    matrix_wrapper &matrix_action) {
  rel_error_ =
      std::numeric_limits<double>::signaling_NaN();  // Resets the value for
                                                     // each run
  const size_t n = rhs_vector.size();
  const size_t m = restart_algorithm_after_this_iteration_;
  std::vector<double> x(n, 0.0), w(n, 0.0), r(n, 0.0);
  x = initial_guess;

  double normb =
      std::sqrt(NormSquared(matrix_action.inv_preconditioner(rhs_vector)));
  r = matrix_action.inv_preconditioner(
      Difference(rhs_vector, matrix_action.mat_vec_operation(x)));
  double beta = std::sqrt(NormSquared(r));

  if (normb == 0.0) normb = 1;

  if ((rel_error_ = std::sqrt(NormSquared(r)) / normb) <= tol_) {
    info_ = 0;
    total_iterations_ = 1;
    return x;
  }

  // The name of these variables were taken from the references mentioned in the
  // LinearSolvers.hpp file
  std::vector<double> s(m + 1, 0.0), cs(m + 1, 0.0), sn(m + 1, 0.0);
  std::vector<std::vector<double>> v(m + 1, std::vector<double>(n, 0.0));
  std::vector<double> H((m + 1) * m, 0.0);

  for (size_t j = 1; j < max_iter_;) {
    v[0] = ScalarMultiplication(1.0 / beta, r);
    std::fill(s.begin(), s.end(), 0.0);
    s[0] = beta;

    for (size_t i = 0; i < m && j <= max_iter_; i++, j++) {
      w = matrix_action.inv_preconditioner(
          matrix_action.mat_vec_operation(v[i]));
      for (size_t k = 0; k <= i; k++) {
        H[m * k + i] = DotProduct(w, v[k]);
        w = Difference(w, ScalarMultiplication(H[m * k + i], v[k]));
      }
      H[m * (i + 1) + i] = std::sqrt(NormSquared(w));
      v[i + 1] = ScalarMultiplication(1.0 / H[m * (i + 1) + i], w);

      for (size_t k = 0; k < i; k++) {
        ApplyPlaneRotation(H[m * k + i], H[m * (k + 1) + i], cs[k], sn[k]);
      }

      GeneratePlaneRotation(H[m * i + i], H[m * (i + 1) + i], cs[i], sn[i]);
      ApplyPlaneRotation(H[m * i + i], H[m * (i + 1) + i], cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);

      if ((rel_error_ = std::abs(s[i + 1]) / normb) <= tol_) {
        // Note that the rel_error_ that we calculated just now is the error in
        // the Krylov subspace. Mathematically this error should be same as the
        // error outside the Krylov subspace(which is the error we are
        // interested in). But, the process of projecting the vector out of the
        // Krylov subspace is complicated enough that numerical errors build up
        // and the mathematical result about the error being equal no longer
        // holds. This only becomes an issue when we have very low tolerances,
        // in my testing anything above 1E-14 is usually fine. For low
        // tolerances the simplest way to get around this is to compare the
        // error in our solution outside the Krylov subspace and take a few more
        // iterations if needed.

        // Copy the value of x in case the we need another iteration in this
        // Krylov subspace. Update will project x out of the subspace and we do
        // not want to get out until the error outside the subspace is smaller
        // than the tolerance.
        auto temp_x = x;
        Update(temp_x, static_cast<int>(i), H, s, v, static_cast<int>(m));

        // Now calculate the error outside the Krylov subspace
        rel_error_ =
            std::sqrt(NormSquared(matrix_action.inv_preconditioner(Difference(
                rhs_vector, matrix_action.mat_vec_operation(temp_x)))))
            / normb;

        // Terminated only if this error outside the Krylov subspace is smaller
        // than the tolerance
        if (rel_error_ <= tol_) {
          total_iterations_ = j;
          info_ = 0;
          return temp_x;
        }
      }
    }

    Update(x, static_cast<int>(m - 1), H, s, v, static_cast<int>(m));
    r = matrix_action.inv_preconditioner(
        Difference(rhs_vector, matrix_action.mat_vec_operation(x)));
    beta = std::sqrt(NormSquared(r));
    if ((rel_error_ = beta / normb) <= tol_) {
      total_iterations_ = j;
      info_ = 0;
      return x;
    }
  }

  total_iterations_ = max_iter_;
  info_ = 1;
  return x;
}

std::vector<double> linear_solvers::BiCGSTAB(
    const std::vector<double> &rhs_vector,
    const std::vector<double> &initial_guess,
    matrix_wrapper &matrix_action) {
  rel_error_ =
      std::numeric_limits<double>::signaling_NaN();  // Resets the value for
                                                     // each run
  size_t n = rhs_vector.size();

  // These variables names are taken from the reference mentioned in the
  // LinearSolvers.hpp file.
  double rhoj1 = 1.0, rhoj = 1.0, wj1 = 1.0, wj = 1.0, b = 1.0, a = 1.0;
  std::vector<double> rb0(n, 0.0), rj(n, 0.0), rj1(n, 0.0), pj(n, 0.0),
      pj1(n, 0.0), vj(n, 0.0), vj1(n, 0.0), y(n, 0.0), s(n, 0.0), z(n, 0.0),
      t(n, 0.0), xj(n, 0.0), xj1(n, 0.0);

  xj = initial_guess;

  rj = Difference(rhs_vector, matrix_action.mat_vec_operation(xj));
  rb0 = rj;
  for (size_t i = 0; i < max_iter_; i++) {
    rhoj1 = DotProduct(rb0, rj);

    // If this happens then the BiCGSTAB algorithm has failed
    // The real check is suppposed to be against 0, but we are using 1E-75
    // because in practice we do not want it to become so small that we get
    // overflow errors etc.
    if (std::abs(rhoj1) < 1E-75) {
      info_ = 2;
      total_iterations_ = i;
      return xj1;
    }

    b = (rhoj1 / rhoj) * (a / wj);
    pj1 = Sum(
        rj,
        ScalarMultiplication(b, Difference(pj, ScalarMultiplication(wj, vj))));

    y = matrix_action.inv_preconditioner(pj1);
    vj1 = matrix_action.mat_vec_operation(y);

    a = rhoj1 / DotProduct(rb0, vj1);

    s = Difference(rj, ScalarMultiplication(a, vj1));

    z = matrix_action.inv_preconditioner(s);
    t = matrix_action.mat_vec_operation(z);

    wj1 = DotProduct(t, s) / DotProduct(t, t);

    xj1 =
        Sum(xj, Sum(ScalarMultiplication(a, y), ScalarMultiplication(wj1, z)));

    rj = Difference(rhs_vector, matrix_action.mat_vec_operation(xj1));

    rel_error_ = std::sqrt(std::abs(DotProduct(rj, rj))
                           / std::abs(DotProduct(rhs_vector, rhs_vector)));
    if (rel_error_ < tol_) {
      info_ = 0;
      total_iterations_ = i;
      return xj1;
    }

    if ((i + 1) % restart_algorithm_after_this_iteration_ == 0) {
      rhoj1 = 1.0;
      rhoj = 1.0;
      wj1 = 1.0;
      wj = 1.0;
      b = 1.0;
      a = 1.0;
      rj = Difference(rhs_vector, matrix_action.mat_vec_operation(xj));
      rb0 = rj;
      std::fill(pj.begin(), pj.end(), 0.0);
      std::fill(vj.begin(), vj.end(), 0.0);

    } else {
      rj1 = Difference(s, ScalarMultiplication(wj1, t));
      rj = rj1;
      xj = xj1;
      rhoj = rhoj1;
      wj = wj1;
      pj = pj1;
      vj = vj1;
    }
  }

  info_ = 1;
  return xj1;
}
