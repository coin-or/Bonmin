/* (C) Copyright Shaurya Sharma
   All Rights Reserved.
   This code is published under the Eclipse Public License.

   Authors :
   S. Sharma

   Date :  07.06.2015
*/

#include "MyBonmin.h"

#include <assert.h>
#include <stdlib.h>

Bool eval_f(Index n, Number* x, Bool new_x,
                            Number* obj_value, UserDataPtr user_data)
{
    assert(n==4);
    *obj_value = - x[0] - x[1] - x[2];
    return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                             Number* grad_f, UserDataPtr user_data)
{
    assert(n==4);
    grad_f[0] = -1.;
    grad_f[1] = -1.;
    grad_f[2] = -1.;
    grad_f[3] = 0.;
    return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
                        Index m, Number* g, UserDataPtr user_data)
{
    assert(n==4);
    assert(m==3);

    g[0] = (x[1] - 1./2.)*(x[1] - 1./2.) + (x[2] - 1./2.)*(x[2] - 1./2.);
    g[1] = x[0] - x[1];
    g[2] = x[0] + x[2] + x[3];

    return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                            Index m, Index nnz_jac_g,
                            Index *iRow, Index *jCol, Number *values,
                            UserDataPtr user_data)
{
    assert(n==4);
    assert(nnz_jac_g == 7);
    if(values == NULL)
    {
        iRow[0] = 2;
        jCol[0] = 1;

        iRow[1] = 3;
        jCol[1] = 1;

        iRow[2] = 1;
        jCol[2] = 2;

        iRow[3] = 2;
        jCol[3] = 2;

        iRow[4] = 1;
        jCol[4] = 3;

        iRow[5] = 3;
        jCol[5] = 3;

        iRow[6] = 3;
        jCol[6] = 4;
    }
    else
    {
        values[0] = 1.;
        values[1] = 1;
        values[2] = 2*x[1] - 1;
        values[3] = -1.;
        values[4] = 2*x[2] - 1;
        values[5] = 1.;
        values[6] = 1.;

    }
    return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
                        Index m, Number *lambda, Bool new_lambda,
                        Index nnz_h_lag, Index *iRow, Index *jCol,
                        Number *values, UserDataPtr user_data)
{
    assert (n==4);
    assert (m==3);
    assert(nnz_h_lag==2);
    if(values==NULL)
    {
        iRow[0] = 2;
        jCol[0] = 2;

        iRow[1] = 3;
        jCol[1] = 3;
    }
    else
    {
        values[0] = 2*lambda[0];
        values[1] = 2*lambda[0];
    }
    return TRUE;
}

