#pragma once

#include <cstddef>
#include <limits>
#include <vector>
/**

This header file contains the linear_solvers class for the GMRES/BiCGSTAB solver
and a matrix wrapper class required by the linear solvers. To use the solvers
one has to inherit the matrix_wrapper class and define two funtions
mat_vec_operations and inv_preconditioner. The function mat_vec_operation takes
a vector as an input, it should return the action of the matrix on the input as
a vector. This allows us to work with just the action of the matrix in a
completely matrix free manner. In the same way the function inv_preconditioner
should return the action of the inverse of the preconditioner matrix.

//---------EXAMPLE-----------//
Utils/Math/Tests/TestLinearSolvers contains example usage.

//--------ADDITIONAL INFO-----------//
final_error() : final relative error after GMRES/BiCGSTAB ends
number_of_iterations_taken() : number of iterations taken by the GMRES/BiCGSTAB
convergence_info(): is 0 if GMRES/BiCGSTAB converges;
                  : is 1 if GMRES/BiCGSTAB fails to converge in max_iter_
                    iterations,
                  : is 2 if BiCGSTAB fails (This is different from reaching the
                    max iterations refer to
                    https://doi.org/10.1137/1.9781611971538
                    (SIAM template book page 24).


//-----------REFERENCES------------//
Main GMRES algorithm is taken from
https://www.netlib.org/templates/cpp/gmres.h which is based on
https://doi.org/10.1137/1.9781611971538 (SIAM template book page 18).

BiCGSTAB algorithm is taken from https://doi.org/10.1137/0913035.
Both left and right preconditionings are possible, this implementations uses
right preconditioning. Which is equivalent to setting K_{1} equal to identity in
the above reference.

**/

// This class defines the action of the matrix and the preconditioner.
class matrix_wrapper {
 public:
  // mat_vec_operation should return the matrix times the input vector
  virtual std::vector<double> mat_vec_operation(
      const std::vector<double>& input) = 0;
  // inv_preconditioner should return the inverse of the preconditioner times
  // the input vector
  virtual std::vector<double> inv_preconditioner(
      const std::vector<double>& input) = 0;
};

class linear_solvers {
 public:
  linear_solvers(const double tol,
                 const size_t max_iter,
                 const size_t restart_algorithm_after_this_iteration)
      : tol_(tol),
        max_iter_(max_iter),
        restart_algorithm_after_this_iteration_(
            restart_algorithm_after_this_iteration) {}

  int convergence_info() { return info_; }
  size_t number_of_iterations_taken() { return total_iterations_; }
  double final_error() { return rel_error_; }

  // Linear solvers
  std::vector<double> GMRES(const std::vector<double>& rhs_vector,
                            const std::vector<double>& initial_guess,
                            matrix_wrapper& matrix_action);
  std::vector<double> BiCGSTAB(const std::vector<double>& rhs_vector,
                               const std::vector<double>& initial_guess,
                               matrix_wrapper& matrix_action);

 private:
  const double tol_;  // tolerance, relative error will be smaller than this
  const size_t max_iter_;
  const size_t
      restart_algorithm_after_this_iteration_;  // number of iterations after
                                                // which GMRES(m)/BiCGSTAB is
                                                // restarted
  int info_;                                    // Convergence information
  size_t total_iterations_ = 0;  // number of iterations taken to converge
  double rel_error_ =
      std::numeric_limits<double>::signaling_NaN();  // Relative error
};
