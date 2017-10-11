#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;


//Define some constants in our system.

const size_t num_actuators = 2;
const size_t num_constraints = 6;


// NOTE: feel free to play around with this
// or do something completely different
double ref_v = 25;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

//Specialization for double
template<>
double deriv(const TCoeffVector& coeffs, double x) {
  double result = 0.0;
  for (int i = 1; i < coeffs.SizeAtCompileTime; i++) {
    result += coeffs[i]* i * std::pow(x, i-1);
  }
  return result;
}

template<>
CppAD::AD<double>  deriv(const TCoeffVector& coeffs, CppAD::AD<double>  x) {
  CppAD::AD<double>  result = 0.0;
  for (int i = 1; i < coeffs.SizeAtCompileTime; i++) {
    result += coeffs[i]* i * CppAD::pow(x, i-1);
  }
  return result;
}

template<>
double poly_eval(const TCoeffVector& coeffs, double x){
  double result = 0.0;
  for (int i = 0; i < coeffs.SizeAtCompileTime; i++) {
    result += coeffs[i] * std::pow(x, i);
  }
  return result;
}

template<>
CppAD::AD<double> poly_eval(const TCoeffVector& coeffs, CppAD::AD<double> x){
  CppAD::AD<double> result = 0.0;
  for (int i = 0; i < coeffs.SizeAtCompileTime; i++) {
    result += coeffs[i] * CppAD::pow(x, i);
  }
  return result;
}



class FG_eval {
 public:
  // Fitted polynomial coefficients
  TCoeffVector coeffs;
  const VariableOffsets& offsets_;
  const size_t num_steps_;
  const double dt_;

  FG_eval(TCoeffVector coeffs,
          const VariableOffsets& offsets,
          size_t num_steps,
          double dt):
    coeffs(coeffs),
    offsets_(offsets),
    num_steps_(num_steps),
    dt_(dt)
  {
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars){

    // `fg` is a vector of the cost constraints,
    // `vars` is a vector of variable values (state & actuators)

	  using namespace CppAD;

  	// The cost is stored is the first element of `fg`.
    // Any additions to the cost are added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    //Cost with respect to the reference state
    for (size_t t = 0; t < num_steps_; ++t) {
      fg[0] += pow(vars[offsets_.cte_offset + t], 2);
      fg[0] += pow(vars[offsets_.eps_offset + t], 2);
      fg[0] += pow(vars[offsets_.v_offset + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (size_t t = 0; t < num_steps_ - 1; ++t) {
      fg[0] += pow(vars[offsets_.delta_offset + t], 2);
      fg[0] += pow(vars[offsets_.a_offset + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (size_t t = 0; t < num_steps_ - 2; ++t) {
      fg[0] += pow(vars[offsets_.delta_offset + t + 1] - vars[offsets_.delta_offset + t], 2);
      fg[0] += pow(vars[offsets_.a_offset + t + 1] - vars[offsets_.a_offset + t], 2);
    }

    // Model constraints.
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + offsets_.x_offset] = vars[offsets_.x_offset];
    fg[1 + offsets_.y_offset] = vars[offsets_.y_offset];
    fg[1 + offsets_.psi_offset] = vars[offsets_.psi_offset];
    fg[1 + offsets_.v_offset] = vars[offsets_.v_offset];
    fg[1 + offsets_.cte_offset] = vars[offsets_.cte_offset];
    fg[1 + offsets_.eps_offset] = vars[offsets_.eps_offset];

    // The rest of the constraints
    for (size_t t = 1; t < num_steps_; ++t) {

      AD<double> x1 = vars[offsets_.x_offset + t];
      AD<double> y1 = vars[offsets_.y_offset + t];
      AD<double> psi1 = vars[offsets_.psi_offset + t];
      AD<double> v1 = vars[offsets_.v_offset + t];
      AD<double> cte1 = vars[offsets_.cte_offset + t];
      AD<double> epsi1 = vars[offsets_.eps_offset + t];

      // The state at time t.
      AD<double> x0 = vars[offsets_.x_offset + t - 1];
      AD<double> y0 = vars[offsets_.y_offset + t - 1];
      AD<double> psi0 = vars[offsets_.psi_offset + t - 1];
      AD<double> v0 = vars[offsets_.v_offset + t - 1];
      AD<double> cte0 = vars[offsets_.cte_offset + t - 1];
      AD<double> epsi0 = vars[offsets_.eps_offset + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[offsets_.delta_offset + t - 1];
      AD<double> a0 = vars[offsets_.a_offset + t - 1];

      AD<double> f0 = poly_eval(coeffs, x0);
      AD<double> psides0 = CppAD::atan(deriv(coeffs, x0));

      //Setup some constraints
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.
      fg[1 + offsets_.x_offset   + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt_);
      fg[1 + offsets_.y_offset   + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt_);
      fg[1 + offsets_.psi_offset + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt_);
      fg[1 + offsets_.v_offset   + t] = v1 - (v0 + a0 * dt_);
      fg[1 + offsets_.cte_offset + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt_));
      fg[1 + offsets_.eps_offset + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt_);
    }
  }
};


//
// MPC_Solution definitions
//
double MPC_Solution::GetX(size_t time_step) {
  return mpc_solution_.x[offsets_.x_offset + time_step + 1];
}

double MPC_Solution::GetY(size_t time_step) {
  return mpc_solution_.x[offsets_.y_offset + time_step + 1];
}

double MPC_Solution::GetSteerAngle(size_t time_step) {
  return mpc_solution_.x[offsets_.psi_offset + time_step + 1];
}

double MPC_Solution::GetThrottle(size_t time_step) {
  return mpc_solution_.x[offsets_.a_offset+ time_step + 1];
}

//
// MPC class definition implementation.
//
MPC::MPC(size_t time_steps, double delta_time):
    time_steps_(time_steps),
    dt_(delta_time),
    offsets_(time_steps)
{
}

MPC::~MPC() {}

MPC_Solution MPC::Solve(const TStateVector& state, const TCoeffVector& coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];


  //Set sizes for used for computation.
  size_t n_vars = state.SizeAtCompileTime * time_steps_ + num_actuators * (time_steps_ - 1);
  size_t n_constraints = time_steps_ * num_constraints;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[offsets_.x_offset] = x;
  vars[offsets_.y_offset] = y;
  vars[offsets_.psi_offset] = psi;
  vars[offsets_.v_offset] = v;
  vars[offsets_.cte_offset] = cte;
  vars[offsets_.eps_offset] = epsi;


  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < offsets_.delta_offset; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (size_t i = offsets_.delta_offset; i < offsets_.a_offset; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (size_t i = offsets_.a_offset; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[offsets_.x_offset] = x;
  constraints_lowerbound[offsets_.y_offset] = y;
  constraints_lowerbound[offsets_.psi_offset] = psi;
  constraints_lowerbound[offsets_.v_offset] = v;
  constraints_lowerbound[offsets_.cte_offset] = cte;
  constraints_lowerbound[offsets_.eps_offset] = epsi;

  constraints_upperbound[offsets_.x_offset] = x;
  constraints_upperbound[offsets_.y_offset] = y;
  constraints_upperbound[offsets_.psi_offset] = psi;
  constraints_upperbound[offsets_.v_offset] = v;
  constraints_upperbound[offsets_.cte_offset] = cte;
  constraints_upperbound[offsets_.eps_offset] = epsi;


  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, offsets_, time_steps_, dt_);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  MPC_Solution::TIpOptSolution solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  return MPC_Solution(solution, offsets_);
}
