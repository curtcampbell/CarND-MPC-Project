#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

typedef Eigen::Matrix<double, 6, 1> TStateVector;
typedef Eigen::Matrix<double, 4, 1> TCoeffVector;

typedef CPPAD_TESTVECTOR(double) TDvector;
typedef CppAD::ipopt::solve_result<TDvector> TIpOptSolution;

using namespace std;


class MPC_Solution{
public:
  typedef CppAD::ipopt::solve_result<TDvector> TIpOptSolution;

  MPC_Solution(TIpOptSolution mpc_solution):
    mpc_solution_(mpc_solution)
  {}

  double GetX(size_t time_step);
  double GetY(size_t time_step);
  double GetSteerAngle(size_t time_step);
  double GetThrottle(size_t time_step);

private:
  TIpOptSolution mpc_solution_;
};


class MPC {
 public:
  MPC(size_t time_steps, double delta_time);

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  MPC_Solution Solve(const TStateVector& state, const TCoeffVector& coeffs);

 private:
  size_t time_steps_;
  double dt_;

};

#endif /* MPC_H */
