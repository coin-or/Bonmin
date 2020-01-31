/* (C) Copyright Shaurya Sharma
   All Rights Reserved.
   This code is published under the Eclipse Public License.

   Authors :
   S. Sharma

   Date :  07.06.2015
*/

#include "BonStdCInterface.h"

Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nnz_jac_g,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nnz_h_lag, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);
