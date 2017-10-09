#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;


//Define some constants in our system.

////number of timesteps
const size_t N = 20;
//Duration of eath timestep
const double dt = 0.05;

const size_t num_actuators = 2;
const size_t num_constraints = 6;


// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_offset = 0;
size_t y_offset = x_offset + N;
size_t psi_offset = y_offset + N;
size_t v_offset = psi_offset + N;
size_t cte_offset = v_offset + N;
size_t eps_offset = cte_offset + N;
size_t delta_offset = eps_offset + N;
size_t a_offset = delta_offset + N - 1;

// NOTE: feel free to play around with this
// or do something completely different
double ref_v = 40;


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


AD<double> poly_eval(const TCoeffVector& coeffs, AD<double> x){
  AD<double> result = 0.0;
  for (int i = 0; i < coeffs.SizeAtCompileTime; i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}


class FG_eval {
 public:
  // Fitted polynomial coefficients
  TCoeffVector coeffs;
  FG_eval(TCoeffVector coeffs):
    coeffs(coeffs)
  {

  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars){

    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

	  using namespace CppAD;

  	// The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    // TODO: Define the cost related the reference state and
    // any anything you think may be beneficial.

    //Cost with respect to the reference state
    for (size_t t = 0; t < N; ++t) {
      fg[0] += pow(vars[cte_offset + t], 2);
      fg[0] += pow(vars[eps_offset + t], 2);
      fg[0] += pow(vars[v_offset + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (size_t t = 0; t < N - 1; ++t) {
      fg[0] += pow(vars[delta_offset + t], 2);
      fg[0] += pow(vars[a_offset + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (size_t t = 0; t < N - 2; ++t) {
      fg[0] += pow(vars[delta_offset + t + 1] - vars[delta_offset + t], 2);
      fg[0] += pow(vars[a_offset + t + 1] - vars[a_offset + t], 2);
    }
    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_offset] = vars[x_offset];
    fg[1 + y_offset] = vars[y_offset];
    fg[1 + psi_offset] = vars[psi_offset];
    fg[1 + v_offset] = vars[v_offset];
    fg[1 + cte_offset] = vars[cte_offset];
    fg[1 + eps_offset] = vars[eps_offset];

    // The rest of the constraints
    for (size_t t = 1; t < N; ++t) {

      AD<double> x1 = vars[x_offset + t];
      AD<double> y1 = vars[y_offset + t];
      AD<double> psi1 = vars[psi_offset + t];
      AD<double> v1 = vars[v_offset + t];
      AD<double> cte1 = vars[cte_offset + t];
      AD<double> epsi1 = vars[eps_offset + t];

      // The state at time t.
      AD<double> x0 = vars[x_offset + t - 1];
      AD<double> y0 = vars[y_offset + t - 1];
      AD<double> psi0 = vars[psi_offset + t - 1];
      AD<double> v0 = vars[v_offset + t - 1];
      AD<double> cte0 = vars[cte_offset + t - 1];
      AD<double> epsi0 = vars[eps_offset + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_offset + t - 1];
      AD<double> a0 = vars[a_offset + t - 1];

      AD<double> f0 = poly_eval(coeffs, x0);
      AD<double> psides0 = CppAD::atan(coeffs[1]);

      //Setup some constraints
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.
      fg[1 + x_offset + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_offset + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_offset + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_offset + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_offset + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + eps_offset + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};


//
// MPC_Solution definitions
//
double MPC_Solution::GetX(size_t time_step) {
  return mpc_solution_.x[x_offset + time_step + 1];
}

double MPC_Solution::GetY(size_t time_step) {
  return mpc_solution_.x[y_offset + time_step + 1];
}

double MPC_Solution::GetSteerAngle(size_t time_step) {
  return mpc_solution_.x[psi_offset + time_step + 1];
}

double MPC_Solution::GetThrottle(size_t time_step) {
  return mpc_solution_.x[a_offset+ time_step + 1];
}

//
// MPC class definition implementation.
//
MPC::MPC(size_t time_steps, double delta_time):
    time_steps_(time_steps),
    dt_(delta_time)
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
  size_t n_vars = state.SizeAtCompileTime * N + num_actuators * (N - 1);
  size_t n_constraints = N * num_constraints;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[x_offset] = x;
  vars[y_offset] = y;
  vars[psi_offset] = psi;
  vars[v_offset] = v;
  vars[cte_offset] = cte;
  vars[eps_offset] = epsi;


  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < delta_offset; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (size_t i = delta_offset; i < a_offset; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (size_t i = a_offset; i < n_vars; i++) {
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

  constraints_lowerbound[x_offset] = x;
  constraints_lowerbound[y_offset] = y;
  constraints_lowerbound[psi_offset] = psi;
  constraints_lowerbound[v_offset] = v;
  constraints_lowerbound[cte_offset] = cte;
  constraints_lowerbound[eps_offset] = epsi;

  constraints_upperbound[x_offset] = x;
  constraints_upperbound[y_offset] = y;
  constraints_upperbound[psi_offset] = psi;
  constraints_upperbound[v_offset] = v;
  constraints_upperbound[cte_offset] = cte;
  constraints_upperbound[eps_offset] = epsi;


  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

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

  return MPC_Solution(solution);
}
