/* (C) Copyright Shaurya Sharma
   All Rights Reserved.
   This code is published under the Eclipse Public License.

   Authors :
   S. Sharma

   Date :  07.06.2015
*/

#ifndef __BONSTDCINTERFACE_H__
#define __BONSTDCINTERFACE_H__

#ifdef __cplusplus
extern "C"
{
#endif

/** Type for all number.
 * 
 * We need to make sure that this is
 * identical with what is defined in Ipopt/src/Common/IpTypes.hpp.
 */
typedef double Number;

/** Type for all incides.
 * 
 * We need to make sure that this is
 * identical with what is defined in Ipopt/src/Common/IpTypes.hpp.
 */
typedef int Index;

/** Type for all integers.
 * 
 * We need to make sure that this is
 * identical with what is defined in Ipopt/src/Common/IpTypes.hpp.
 */
typedef int Int;

/* TODO: avoid duplication of enums */
typedef enum
{
  CONTINUOUS,
  BINARY,
  INTEGER
} VariableTypeC;

typedef enum
{
  LINEAR,
  NON_LINEAR
} LinearityTypeC;

typedef enum
{ C_STYLE = 0,
  FORTRAN_STYLE = 1
} IndexStyleC;

struct BonminBranchingInfo;
struct BonminSosInfo;
struct BonminProblemInfo;
typedef struct BonminProblemInfo* BonminProblem;

typedef int Bool;

#ifndef TRUE
# define TRUE (1)
#endif
#ifndef FALSE
# define FALSE (0)
#endif

/** A pointer for anything that is to be passed between the called
 *  and individual callback function
 */
typedef void * UserDataPtr;

/** Type defining the callback function for evaluating the value of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation.
 */
typedef Bool (*Eval_F_CB)(Index n, Number* x, Bool new_x,
                      Number* obj_value, UserDataPtr user_data);

/** Type defining the callback function for evaluating the gradient of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation.
 */
typedef Bool (*Eval_Grad_F_CB)(Index n, Number* x, Bool new_x,
                           Number* grad_f, UserDataPtr user_data);

/** Type defining the callback function for evaluating the value of
 *  the constraint functions.  Return value should be set to false if
 *  there was a problem doing the evaluation.
 */
typedef Bool (*Eval_G_CB)(Index n, Number* x, Bool new_x,
                      Index m, Number* g, UserDataPtr user_data);

/** Type defining the callback function for evaluating the Jacobian of
 *  the constrant functions.  Return value should be set to false if
 *  there was a problem doing the evaluation.
 */
typedef Bool (*Eval_Jac_G_CB)(Index n, Number *x, Bool new_x,
                          Index m, Index nnz_jac_g,
                          Index *iRow, Index *jCol, Number *values,
                          UserDataPtr user_data);

/** Type defining the callback function for evaluating the Hessian of
 *  the Lagrangian function.  Return value should be set to false if
 *  there was a problem doing the evaluation.
 */
typedef Bool (*Eval_H_CB)(Index n, Number *x, Bool new_x, Number obj_factor,
                      Index m, Number *lambda, Bool new_lambda,
                      Index nnz_h_lag, Index *iRow, Index *jCol,
                      Number *values, UserDataPtr user_data);

/** Type defining the callback function for giving intermediate
 *  execution control to the user.  If set, it is called once per
 *  iteration, providing the user with some information on the state
 *  of the optimization.  This can be used to print some
 *  user-defined output.  It also gives the user a way to terminate
 *  the optimization prematurely.  If this method returns false,
 *  Ipopt will terminate the optimization.
 */
typedef Bool (*Intermediate_CB)(Index alg_mod, /* 0 is regular, 1 is resto */
            Index iter_count, Number obj_value,
            Number inf_pr, Number inf_du,
            Number mu, Number d_norm,
            Number regularization_size,
            Number alpha_du, Number alpha_pr,
            Index ls_trials, UserDataPtr user_data);

BonminProblem CreateBonminProblem(
          Index n
          , Number* x_L
          , Number* x_U
          , Index m
          , Number* g_L
          , Number* g_U
          , Index nnz_jac_g
          , Index nnz_h_lag
          , Index index_style
          , Eval_F_CB eval_f
          , Eval_G_CB eval_g
          , Eval_Grad_F_CB eval_grad_f
          , Eval_Jac_G_CB eval_jac_g
          , Eval_H_CB eval_h
          , VariableTypeC* var_types
          , LinearityTypeC* var_linearity_types
          , LinearityTypeC* constraint_linearity_types
          , struct BonminSosInfo* sos_info
          , struct BonminBranchingInfo* branching_info );


/** Method for freeing a previously created BonminProblem.
 * After freeing an BonminProblem, it cannot be used anymore.
 */
void FreeBonminProblem(BonminProblem bonmin_problem);


/** Function for adding a string option.
 * Returns FALSE the option could not be set (e.g., if keyword is unknown).
 */
Bool AddBonminStrOption(BonminProblem bonmin_problem, char* keyword, char* val);

/** Function for adding a Number option.
 * Returns FALSE the option could not be set (e.g., if keyword is unknown).
 */
Bool AddBonminNumOption(BonminProblem bonmin_problem, char* keyword, Number val);

Bool ReadBonminOptFile( BonminProblem bonmin_problem, char* file_path );

Bool ReadBonminOptString( BonminProblem bonmin_problem, char* option );

/** Function for adding an Int option.
 * Returns FALSE the option could not be set (e.g., if keyword is unknown).
 */
Bool AddBonminIntOption(BonminProblem bonmin_problem, char* keyword, Int val);

/** Function for opening an output file for a given name with given printlevel.
 * Returns false, if there was a problem opening the file.
 */
Bool OpenBonminOutputFile(BonminProblem bonmin_problem, char* file_name,
                     Int print_level);

/** Optional function for setting scaling parameter for the NLP.
 *  This corresponds to the get_scaling_parameters method in TNLP.
 *  If the pointers x_scaling or g_scaling are NULL, then no scaling
 *  for x resp. g is done.
 */
Bool SetBonminProblemScaling(BonminProblem bonmin_problem,
            Number obj_scaling,
            Number* x_scaling,
            Number* g_scaling);

Int BonminSolve(
    BonminProblem bonmin_problem
    , Number* x
    , Number* g
    , Number* obj_val
    , Number* mult_g
    , Number* mult_x_L
    , Number* mult_x_U
    , UserDataPtr user_data );

#ifdef __cplusplus
}
#endif

#endif
