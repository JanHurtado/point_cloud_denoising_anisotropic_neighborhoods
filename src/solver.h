#ifndef SOLVER_H
#define SOLVER_H

#include <ilcplex/cplex.h>
#include "point_cloud.h"

void quadratic_programming_solver(	int num_vars, PCDMatrix & linear_term, PCDMatrix & quadratic_term, PCDMatrix & constraint_coefficients_diag,
									double single_constraint, vector<double> & solution);

#endif // SOLVER_H