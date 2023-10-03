/// \file
/// Defines general vector operations

#ifndef VectorOps_hpp
#define VectorOps_hpp


/// Compute the dot (inner) product of two vectors
///
/// This will return the quantity A . B, computed from the input vectors,
/// which can be any type of container implementing the STL interface.
template <class T, class U>
double DotProduct(const T& A, const U& B) {
  double result = 0;

  auto Ait = A.begin();
  for (auto Bit = B.begin();
       Ait != A.end();
       ++Ait, ++Bit) {
    result += (*Ait)*(*Bit);
  }

  return result;
}

/// Compute the squared norm (magnitude) of a vector
///
/// This will return the quantity ||A||^2, computed from the input vector,
/// which can be any type of container implementing the STL interface.
template <class T>
double NormSquared(const T& A)
{
  return DotProduct(A, A);
}

/// Compute the cross product of two vectors
///
/// This will return the quantity A x B, computed from the input vectors,
/// which can be any type of container implementing the STL interface.
/// The resulting vector is returned, in contrast to the version in
/// CrossProduct.hpp, which allows return value optimization, and is more
/// convenient.
// template <class T>
// T CrossProduct(const T& A, const T& B) {
//   T result = A;

//   for (CyclicIterator E(3); E; ++E)
//   {
//     result[E(0)] = A[E(1)] * B[E(2)] - B[E(1)] * A[E(2)];
//   }

//   return result;
// }

/// Compute the difference between two vectors
///
/// This will return the quantity A - B, computed from the input vectors,
/// which can be any type of container implementing the STL interface.
template <class T>
T Difference(const T& A, const T& B) {
  T result = A;

  auto rit = result.begin();
  for (auto Bit = B.begin();
       rit != result.end();
       ++rit, ++Bit)
  {
    *rit -= *Bit;
  }

  return result;
}

/// Compute the sum of two vectors
///
/// This will return the quantity A + B, computed from the input vectors,
/// which can be any type of container implementing the STL interface.
template <class T>
T Sum(const T& A, const T& B) {
  T result = A;

  auto rit = result.begin();
  for (auto Bit = B.begin();
       rit != result.end();
       ++rit, ++Bit)
  {
    *rit += *Bit;
  }

  return result;
}

/// Compute the product of a scalar and a vector
///
/// This will return the quantity a*B, computed from the input scalar(a) and
/// vector(B), B can be any type of container implementing the STL interface.
template <class Scalar, class T>
T ScalarMultiplication(const Scalar& a, const T& B) {
  T result = B;

  for (auto rit = result.begin(); rit != result.end(); ++rit) {
    *rit *= a;
  }

  return result;
}

#endif