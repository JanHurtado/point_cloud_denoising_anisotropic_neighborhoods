#include "solver.h"

void quadratic_programming_solver(int num_vars, PCDMatrix & linear_term, PCDMatrix & quadratic_term, PCDMatrix & constraint_coefficients_diag, double single_constraint, vector<double> & solution)
{
	if (num_vars > 0)
	{
		char     probname[] = "quadratic programming solver";
		int      numcols = num_vars;
		int      numrows = 1;
		int      objsen = 1;
		PCDNumber   *obj = new PCDNumber[numcols];
		PCDNumber   *rhs = new PCDNumber[numrows];
		char     *sense = new char[numrows];
		int      *matbeg = new int[numcols];
		int      *matcnt = new int[numcols];
		int      *matind = new int[numcols];
		PCDNumber   *matval = new PCDNumber[numcols];
		PCDNumber   *lb = new PCDNumber[numcols];
		PCDNumber   *ub = new PCDNumber[numcols];
		int      *qmatbeg = new int[numcols];
		int      *qmatcnt = new int[numcols];
		int      *qmatind = new int[numcols * numcols];
		PCDNumber   *qmatval = new PCDNumber[numcols * numcols];

		int      solstat;
		PCDNumber   objval;
		PCDNumber * x = new PCDNumber[numcols];
		PCDNumber * pi = new PCDNumber[numrows];
		PCDNumber * slack = new PCDNumber[numrows];
		PCDNumber * dj = new PCDNumber[numcols];


		CPXENVptr     env = NULL;
		CPXLPptr      lp = NULL;
		int           status;
		int           cur_numrows, cur_numcols;
		int i, j;

		env = CPXopenCPLEX(&status);

		if (env == NULL) {
			char  errmsg[CPXMESSAGEBUFSIZE];
			fprintf(stderr, "Could not open CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
			goto TERMINATE;
		}

		status = CPXsetintparam(env, CPXPARAM_OptimalityTarget, CPX_OPTIMALITYTARGET_FIRSTORDER);
		if (status) {
			fprintf(stderr,
				"Failure to turn on screen indicator, error %d.\n", status);
			goto TERMINATE;
		}

		rhs[0] = single_constraint;
		sense[0] = 'E';
		int k = 0;
		for (i = 0; i < num_vars; i++)
		{
			lb[i] = 0.0f;
			ub[i] = 1.0f;
			obj[i] = linear_term(i, 0);
			matbeg[i] = i;
			matcnt[i] = 1;
			matind[i] = 0;
			matval[i] = constraint_coefficients_diag(i, i);
			qmatbeg[i] = i * num_vars;
			qmatcnt[i] = num_vars;
			for (j = 0; j < num_vars; j++)
			{
				qmatind[k] = j;
				qmatval[k] = quadratic_term(j, i);
				k++;
			}

		}

		lp = CPXcreateprob(env, &status, probname);

		if (lp == NULL) {
			fprintf(stderr, "Failed to create problem.\n");
			goto TERMINATE;
		}

		status = CPXcopylp(env, lp, numcols, numrows, objsen, obj, rhs,
			sense, matbeg, matcnt, matind, matval,
			lb, ub, NULL);

		if (status) {
			fprintf(stderr, "Failed to copy problem data.\n");
			goto TERMINATE;
		}

		status = CPXcopyquad(env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
		if (status) {
			fprintf(stderr, "Failed to copy quadratic matrix.\n");
			cout << "status " << status << endl;
			goto TERMINATE;
		}

		status = CPXqpopt(env, lp);

		if (status) {
			fprintf(stderr, "Failed to optimize QP.\n");
			goto TERMINATE;
		}

		status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
		if (status) {
			fprintf(stderr, "Failed to obtain solution.\n");
			goto TERMINATE;
		}

		cur_numrows = CPXgetnumrows(env, lp);
		cur_numcols = CPXgetnumcols(env, lp);

		solution.clear();
		solution.resize(cur_numcols);
		for (j = 0; j < cur_numcols; j++) {
			cur_numcols;
			solution[j] = x[j];
		}

	TERMINATE:

		if (lp != NULL) {
			status = CPXfreeprob(env, &lp);
			if (status) {
				fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
			}
		}

		if (env != NULL) {
			status = CPXcloseCPLEX(&env);

			if (status) {
				char  errmsg[CPXMESSAGEBUFSIZE];
				fprintf(stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring(env, status, errmsg);
				fprintf(stderr, "%s", errmsg);
			}
		}

		delete[] obj;
		delete[] rhs;
		delete[] sense;
		delete[] matbeg;
		delete[] matcnt;
		delete[] matval;
		delete[] lb;
		delete[] ub;
		delete[] qmatbeg;
		delete[] qmatcnt;
		delete[] qmatind;
		delete[] qmatval;

		delete[] x;
		delete[] pi;
		delete[] slack;
		delete[] dj;
	}
}
