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

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lives easier.
struct VariableOffsets {

  size_t x_offset;
  size_t y_offset;
  size_t psi_offset;
  size_t v_offset;
  size_t cte_offset;
  size_t eps_offset;
  size_t delta_offset;
  size_t a_offset;

  VariableOffsets(size_t num_time_steps) {
    x_offset = 0;
    y_offset = x_offset + num_time_steps;
    psi_offset = y_offset + num_time_steps;
    v_offset = psi_offset + num_time_steps;
    cte_offset = v_offset + num_time_steps;
    eps_offset = cte_offset + num_time_steps;
    delta_offset = eps_offset + num_time_steps;
    a_offset = delta_offset + num_time_steps - 1;
  }

};

template<typename TNum>
TNum deriv(const TCoeffVector& coeffs, TNum x);

template<typename TNum>
TNum poly_eval(const TCoeffVector& coeffs, TNum x);


class MPC_Solution{
public:
  typedef CppAD::ipopt::solve_result<TDvector> TIpOptSolution;

  MPC_Solution(TIpOptSolution mpc_solution, const VariableOffsets& offsets):
    mpc_solution_(mpc_solution),
    offsets_(offsets)
  {}

  double GetX(size_t time_step);
  double GetY(size_t time_step);
  double GetSteerAngle(size_t time_step);
  double GetThrottle(size_t time_step);

private:
  TIpOptSolution mpc_solution_;
  const VariableOffsets offsets_;
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
  VariableOffsets offsets_;

};

#endif /* MPC_H */
