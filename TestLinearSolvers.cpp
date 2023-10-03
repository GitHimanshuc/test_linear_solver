#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include "LinearSolvers.hpp"
#include "VectorOps.hpp"

// This class inherits the matrix_wrapper class and defines the
// mat_vec_operation and inv_preconditioner functions.
// mat_vec_operation and inv_preconditioner functions can be anything as long as
// they return a correct sized vector.
class my_matrix_wrapper : public matrix_wrapper {
 public:
  std::vector<double> matrix_;
  const size_t number_of_rows_;

  my_matrix_wrapper(size_t system_size, std::vector<double> matrix)
      : number_of_rows_(system_size) {
    matrix_ = matrix;
  }

  std::vector<double> mat_vec_operation(const std::vector<double> &input) {
    std::vector<double> output(input.size(), 0.0);

    for (size_t i = 0; i < input.size(); i++) {
      for (size_t j = 0; j < input.size(); j++) {
        output[i] += matrix_[number_of_rows_ * i + j] * input[j];
      }
    }

    return output;
  }
  std::vector<double> inv_preconditioner(const std::vector<double> &input) {
    return input;  // Identity preconditioner
  }
};

// Defining a helper function for checking the results
std::vector<double> matrix_vector_product(const std::vector<double> &input,
                                          const std::vector<double> &matrix) {
  std::vector<double> output(input.size(), 0.0);

  for (size_t i = 0; i < input.size(); i++) {
    for (size_t j = 0; j < input.size(); j++) {
      output[i] += matrix[input.size() * i + j] * input[j];
    }
  }
  return output;
}

int main() {
  const size_t system_size = 100;
  const size_t maximum_iterations = system_size;
  // We do not want restarting so we will set the
  // 'restart_algorithm_after_this_iteration' to maximum_iterations
  const size_t restart_algorithm_after_this_iteration = maximum_iterations;
  const double tolerance = 1E-15;

  std::vector<double> rhs_vector(system_size, 0.0);
  std::vector<double> initial_guess(system_size, 0.0);
  std::vector<double> matrix(system_size * system_size, 0.0);

  for (size_t i = 0; i < system_size; i++) {
    rhs_vector[i] = static_cast<double>(i + 1);
    matrix[system_size * i + i] +=
        static_cast<double>(i + 2);  // Creating the matrix
  }

  my_matrix_wrapper matrix_operations = my_matrix_wrapper(system_size, matrix);
  linear_solvers solvers = linear_solvers(
      tolerance, maximum_iterations, restart_algorithm_after_this_iteration);

  // Testing GMRES
  std::vector<double> ans_gmres =
      solvers.GMRES(rhs_vector, initial_guess, matrix_operations);
  if (solvers.convergence_info() == 0) {
    std::cout << "GMRES converged in " << solvers.number_of_iterations_taken()
              << " iterations. The final error was : " << solvers.final_error()
              << std::endl;
  } else {
    std::cout << "GMRES failed to converge in "
              << solvers.number_of_iterations_taken()
              << " iterations. The final error was : " << solvers.final_error()
              << std::endl;
  }

  // Testing BiCGSTAB
  std::vector<double> ans_bicgstab =
      solvers.BiCGSTAB(rhs_vector, initial_guess, matrix_operations);
  if (solvers.convergence_info() == 0) {
    std::cout << "BiCGSTAB converged in "
              << solvers.number_of_iterations_taken()
              << " iterations. The final error was : " << solvers.final_error()
              << std::endl;
  } else {
    std::cout << "BiCGSTAB failed to converge in "
              << solvers.number_of_iterations_taken()
              << " iterations. The final error was : " << solvers.final_error()
              << std::endl;
  }

  // Checking the results
  // Note that both solvers have inbuilt function final_error() that gives the
  // final relative error. Here we are calculating errors without using the
  // inbuilt function for testing.
  std::vector<double> residue_gmres =
      Difference(matrix_vector_product(ans_gmres, matrix), rhs_vector);
  std::vector<double> residue_bicgstab =
      Difference(matrix_vector_product(ans_bicgstab, matrix), rhs_vector);
  double relative_error_gmres =
      std::sqrt(DotProduct(residue_gmres, residue_gmres)
                / DotProduct(rhs_vector, rhs_vector));
  double relative_error_bicgstab =
      std::sqrt(DotProduct(residue_bicgstab, residue_bicgstab)
                / DotProduct(rhs_vector, rhs_vector));

  std::cout <<
          "GMRES solution error: " << relative_error_gmres << "\n";
  std::cout <<
          "BiCGSTAB solution error: " << relative_error_bicgstab << "\n";

  // return UtilsForTesting::NumberOfTestsFailed() > 0;
  return 0;
}
