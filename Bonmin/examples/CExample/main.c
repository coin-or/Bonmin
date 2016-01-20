/* (C) Copyright Shaurya Sharma
   All Rights Reserved.
   This code is published under the Eclipse Public License.

   Authors :
   S. Sharma

   Date :  07.06.2015
*/

#include "MyBonmin.h"

#include <float.h>
#include <stdlib.h>
#include <stdio.h>

enum {NUM_VARS = 4, NUM_CONSTRAINTS = 3};

Index n = NUM_VARS;
Index m = NUM_CONSTRAINTS;
Index nnz_jac_g = 7;
Index nnz_h_lag = 2;
Index index_style = FORTRAN_STYLE;
VariableTypeC variable_types[NUM_VARS] = {BINARY, CONTINUOUS, CONTINUOUS, INTEGER};
LinearityTypeC variable_linearity_types[NUM_VARS] = {LINEAR, NON_LINEAR, NON_LINEAR, LINEAR};
LinearityTypeC constraint_linearity_types[NUM_CONSTRAINTS] = {NON_LINEAR, LINEAR, LINEAR};
Number x_L[NUM_VARS] = {0., 0., 0., 0};
Number x_U[NUM_VARS] = {1., DBL_MAX, DBL_MAX, 5};
Number g_L[NUM_CONSTRAINTS] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
Number g_U[NUM_CONSTRAINTS] = {1./4., 0, 2};
Number starting_point[NUM_VARS] = {0., 0., 0., 0};


int main(int argc, char const *argv[])
{
    BonminProblem bonmin_problem = NULL;
    Number *obj_val = NULL;
    Number *mult_g = (Number*)malloc(m*sizeof(Number));
    Number *mult_x_L = (Number*)malloc(n*sizeof(Number));
    Number *mult_x_U = (Number*)malloc(n*sizeof(Number));
    UserDataPtr user_data = NULL;
    Int result = 0;

    bonmin_problem = CreateBonminProblem(n, x_L, x_U, m, g_L, g_U, nnz_jac_g, nnz_h_lag,
                        index_style, &eval_f, &eval_g, &eval_grad_f,
                        &eval_jac_g, &eval_h, variable_types, variable_linearity_types,
                        constraint_linearity_types, NULL, NULL);
    AddBonminStrOption(bonmin_problem, "bonmin.algorithm", "B-OA");


    result = BonminSolve(bonmin_problem, starting_point, NULL, obj_val,
			 mult_g, mult_x_L, mult_x_U, user_data);

    FreeBonminProblem(bonmin_problem);
    printf("Finished: result = %d\n", result);
    return 0;
}
