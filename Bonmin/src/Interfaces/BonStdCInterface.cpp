// (C) Copyright Shaurya Sharma
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// S. Sharma
//
// Date :  07.06.2015

#include "BonStdCInterface.h"
#include "BonStdInterfaceTMINLP.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptApplication.hpp"
#include "BonBonminSetup.hpp"
#include "BonCbc.hpp"

using namespace Bonmin;
using namespace Ipopt;


//Note BonminSosInfo and BonminBranchingInfo are copied from BonTMINLP.hpp. This is not ideal
struct BonminSosInfo
{
/** Number of SOS constraints.*/
        int num;

/** Type of sos. At present Only type '1' SOS are supported by Cbc*/
        char * types;

/** priorities of sos constraints.*/
        int * priorities;

/** \name Sparse storage of the elements of the SOS constraints.*/
/** @{ */
/** Total number of non zeroes in SOS constraints.*/
        int numNz;
/** For 0 <= i < nums, start[i] gives the indice of indices and weights arrays at which the description of constraints i begins..*/
        int * starts;
/** indices of elements belonging to the SOS.*/
        int * indices;
/** weights of the elements of the SOS.*/
        double * weights;
};

struct BonminBranchingInfo
{
    /**number of variables*/
    int size;
    /** User set priorities on variables. */
    int * priorities;
    /** User set preferered branching direction. */
    int * branchingDirections;
    /** User set up pseudo costs.*/
    double * upPsCosts;
    /** User set down pseudo costs.*/
    double * downPsCosts;
};

struct BonminProblemInfo
{
    Index n;
    Number* x_L;
    Number* x_U;
    Index m;
    Number* g_L;
    Number* g_U;
    Index nele_jac;
    Index nele_hess;
    Index index_style;
    Eval_F_CB eval_f;
    Eval_G_CB eval_g;
    Eval_Grad_F_CB eval_grad_f;
    Eval_Jac_G_CB eval_jac_g;
    Eval_H_CB eval_h;
    VariableTypeC* var_types;
    LinearityTypeC* var_linearity_types;
    LinearityTypeC* constraint_linearity_types;
    BonminBranchingInfo* branch_info;
    BonminSosInfo* sos_info;
    Intermediate_CB* intermediate_cb;
    BonminSetup bonmin_setup;
    Number obj_scaling;
    Number* x_scaling;
    Number* g_scaling;
};

BonminProblem CreateBonminProblem(
                Index n
                , Number* x_L
                , Number* x_U
                , Index m
                , Number* g_L
                , Number* g_U
                , Index nele_jac
                , Index nele_hess
                , Index index_style
                , Eval_F_CB eval_f
                , Eval_G_CB eval_g
                , Eval_Grad_F_CB eval_grad_f
                , Eval_Jac_G_CB eval_jac_g
                , Eval_H_CB eval_h
                , VariableTypeC* var_types
                , LinearityTypeC* var_linearity_types
                , LinearityTypeC* constraint_linearity_types
                , BonminSosInfo* sos_info
                , BonminBranchingInfo* branch_info )
{
    if ( n<1 || m<0 || !x_L || !x_U || (m>0 && (!g_L || !g_U)) ||
            (m==0 && nele_jac != 0) || (m>0 && nele_jac < 1) || nele_hess < 0 ||
            !eval_f || !eval_grad_f || (m>0 && (!eval_g || !eval_jac_g)) ||
            !var_types || !var_linearity_types || (m>0 && !constraint_linearity_types) )
    {
        return NULL;
    }
    BonminProblem retval = new BonminProblemInfo;

    retval->n = n;
    retval->x_L = new Number[n];

    for (Index i=0; i<n; i++)
    {
        retval->x_L[i] = x_L[i];
    }
    retval->x_U = new Number[n];
    for (Index i=0; i<n; i++)
    {
        retval->x_U[i] = x_U[i];
    }

    retval->m = m;
    if (m>0)
    {
        retval->g_L = new Number[m];
        for (Index i=0; i<m; i++)
        {
            retval->g_L[i] = g_L[i];
        }
        retval->g_U = new Number[m];
        for (Index i=0; i<m; i++)
        {
            retval->g_U[i] = g_U[i];
        }
    }
    else
    {
        retval->g_L = NULL;
        retval->g_U = NULL;
    }

    retval->var_types = new VariableTypeC[n];
    retval->var_linearity_types = new LinearityTypeC[n];
    retval->constraint_linearity_types = new LinearityTypeC[m];
    for (int i = 0; i < n; ++i)
    {
        retval->var_types[i] = var_types[i];
        retval->var_linearity_types[i] = var_linearity_types[i];
    }
    for (int i = 0; i < m; ++i)
    {
        retval->constraint_linearity_types[i] = constraint_linearity_types[i];
    }

    retval->nele_jac = nele_jac;
    retval->nele_hess = nele_hess;
    retval->index_style = index_style;
    retval->eval_f = eval_f;
    retval->eval_g = eval_g;
    retval->eval_grad_f = eval_grad_f;
    retval->eval_jac_g = eval_jac_g;
    retval->eval_h = eval_h;
    retval->intermediate_cb = NULL;

    retval->obj_scaling = 1;
    retval->x_scaling = NULL;
    retval->g_scaling = NULL;

    retval->bonmin_setup.initializeOptionsAndJournalist();

    return retval;
}

void FreeBonminProblem(BonminProblem bonmin_problem)
{
    if (bonmin_problem == NULL)
    {
        return;
    }

    delete [] bonmin_problem->x_L;
    delete [] bonmin_problem->x_U;
    delete [] bonmin_problem->var_types;
    delete [] bonmin_problem->var_linearity_types;
    if (bonmin_problem->m>0)
    {
        delete [] bonmin_problem->constraint_linearity_types;
        delete [] bonmin_problem->g_L;
        delete [] bonmin_problem->g_U;
    }

    if(bonmin_problem->x_scaling != NULL)
    {
        delete [] bonmin_problem->x_scaling;
    }
    if(bonmin_problem->g_scaling != NULL)
    {
        delete [] bonmin_problem->g_scaling;
    }

    delete bonmin_problem;
}

Bool AddBonminStrOption(BonminProblem bonmin_problem, char* keyword, char* val)
{
    std::string tag( keyword );
    std::string value( val );
    return (Bool) bonmin_problem->bonmin_setup.options()->SetStringValue( tag, value );
}

Bool AddBonminNumOption(BonminProblem bonmin_problem, char* keyword, Number val)
{
    std::string tag(keyword);
    Number value=val;
    return (Bool) bonmin_problem->bonmin_setup.options()->SetNumericValue( tag, value );
}

Bool AddBonminIntOption(BonminProblem bonmin_problem, char* keyword, Int val)
{
    std::string tag(keyword);
    Index value=val;
    return (Bool) bonmin_problem->bonmin_setup.options()->SetIntegerValue( tag, value );
}
Bool ReadBonminOptString( BonminProblem bonmin_problem, char* option )
{
    std::string opt_string( option );
    bonmin_problem->bonmin_setup.readOptionsString( opt_string );
    return TRUE;
}

Bool ReadBonminOptFile( BonminProblem bonmin_problem, char* file_path )
{
    std::string file_path_str( file_path );
    if( file_path_str.compare( "" ) == 0 )
    {
        bonmin_problem->bonmin_setup.readOptionsFile();
    }
    else
    {
        bonmin_problem->bonmin_setup.readOptionsFile(file_path_str);
    }
    return TRUE;
}

Bool OpenBonminOutputFile(BonminProblem bonmin_problem, char* file_name, Int print_level)
{
    // TODO
    // std::string name(file_name);
    // EJournalLevel level = EJournalLevel(print_level);
    // return (Bool) bonmin_problem->bonmin_setup.options()->OpenOutputFile(name, level);
    return ( Bool )true;
}

Bool SetBonminProblemScaling(BonminProblem bonmin_problem,
              Number obj_scaling,
              Number* x_scaling,
              Number* g_scaling)
{
    bonmin_problem->obj_scaling = obj_scaling;
    if (x_scaling)
    {
        if (!bonmin_problem->x_scaling)
        {
            bonmin_problem->x_scaling = new Number[bonmin_problem->n];
        }
        for (Index i=0; i<bonmin_problem->n; i++)
        {
            bonmin_problem->x_scaling[i] = x_scaling[i];
        }
    }
    else
    {
        delete [] bonmin_problem->x_scaling;
        bonmin_problem->x_scaling = NULL;
    }
    if (g_scaling)
    {
        if (!bonmin_problem->g_scaling)
        {
            bonmin_problem->g_scaling = new Number[bonmin_problem->m];
        }
        for (Index i=0; i<bonmin_problem->m; i++)
        {
            bonmin_problem->g_scaling[i] = g_scaling[i];
        }
    }
    else
    {
        delete [] bonmin_problem->g_scaling;
        bonmin_problem->g_scaling = NULL;
    }

    return (Bool)true;
}

Int BonminSolve(
  BonminProblem bonmin_problem
, Number* x
, Number* g
, Number* obj_val
, Number* mult_g
, Number* mult_x_L
, Number* mult_x_U
, UserDataPtr user_data)
{
    // TODO return a meaningful status code
    // Initialize and process options
    // ApplicationReturnStatus retval = ipopt_problem->app->Initialize();
    // if (retval!=Solve_Succeeded) {
        // return (::ApplicationReturnStatus) retval;
    // }

    if (!x) {
        return 0;
        // ipopt_problem->app->Jnlst()->Printf(J_ERROR, J_MAIN,
                                            // "Error: Array x with starting point information is NULL.");
        // return (::ApplicationReturnStatus) Invalid_Problem_Definition;
    }

    // Copy the starting point information
    Number* start_x = new Number[bonmin_problem->n];
    for (Index i=0; i<bonmin_problem->n; i++)
    {
        start_x[i] = x[i];
    }
    Number* start_lam = NULL;
    if (mult_g)
    {
        start_lam = new Number[bonmin_problem->m];
        for (Index i=0; i<bonmin_problem->m; i++)
        {
            start_lam[i] = mult_g[i];
        }
    }
    Number* start_z_L = NULL;
    if (mult_x_L)
    {
        start_z_L = new Number[bonmin_problem->n];
        for (Index i=0; i<bonmin_problem->n; i++)
        {
            start_z_L[i] = mult_x_L[i];
        }
    }
    Number* start_z_U = NULL;
    if (mult_x_U)
    {
        start_z_U = new Number[bonmin_problem->n];
        for (Index i=0; i<bonmin_problem->n; i++)
        {
            start_z_U[i] = mult_x_U[i];
        }
    }


    // convert from C type enum to C++ enum
    TMINLP::VariableType* variable_types = new TMINLP::VariableType[bonmin_problem->n];
    TNLP::LinearityType* variable_linearity_types = new TNLP::LinearityType[bonmin_problem->n];
    for (Index i=0; i<bonmin_problem->n; i++)
    {
        variable_types[i] = static_cast<TMINLP::VariableType>( (int)bonmin_problem->var_types[i] ); // yikes
        variable_linearity_types[i] = static_cast<TNLP::LinearityType>( (int)bonmin_problem->var_linearity_types[i] ); // yikes
    }

    TNLP::LinearityType* constraint_linearity_types = NULL;
    if( bonmin_problem->m > 0 )
    {
        constraint_linearity_types = new TNLP::LinearityType[bonmin_problem->m];
        for (Index i=0; i<bonmin_problem->m; i++)
        {
            constraint_linearity_types[i] = static_cast<TNLP::LinearityType>( (int)bonmin_problem->constraint_linearity_types[i] ); // yikes
        }
    }


    SmartPtr<TMINLP> tminlp;
    SmartPtr<StdInterfaceTMINLP> interfaceTMINLP;
    try
    {
        interfaceTMINLP = new StdInterfaceTMINLP(bonmin_problem->n, bonmin_problem->x_L,
                                        bonmin_problem->x_U, bonmin_problem->m,
                                        bonmin_problem->g_L, bonmin_problem->g_U,
                                        bonmin_problem->nele_jac,
                                        bonmin_problem->nele_hess,
                                        bonmin_problem->index_style,
                                        start_x, start_lam, start_z_L, start_z_U,
                                        bonmin_problem->eval_f, bonmin_problem->eval_g,
                                        bonmin_problem->eval_grad_f,
                                        bonmin_problem->eval_jac_g,
                                        bonmin_problem->eval_h,
                                        variable_types,
                                        variable_linearity_types,
                                        constraint_linearity_types,
                                        // bonmin_problem->intermediate_cb,
                                        x, mult_x_L, mult_x_U, g, mult_g,
                                        obj_val, user_data,
                                        bonmin_problem->obj_scaling,
                                        bonmin_problem->x_scaling,
                                        bonmin_problem->g_scaling);

        if( bonmin_problem->sos_info != NULL )
        {
            interfaceTMINLP->getSosInfo()->num                       = bonmin_problem->sos_info->num;
            interfaceTMINLP->getSosInfo()->types                     = bonmin_problem->sos_info->types;
            interfaceTMINLP->getSosInfo()->priorities                = bonmin_problem->sos_info->priorities;
            interfaceTMINLP->getSosInfo()->numNz                     = bonmin_problem->sos_info->numNz;
            interfaceTMINLP->getSosInfo()->starts                    = bonmin_problem->sos_info->starts;
            interfaceTMINLP->getSosInfo()->indices                   = bonmin_problem->sos_info->indices;
            interfaceTMINLP->getSosInfo()->weights                   = bonmin_problem->sos_info->weights;
        }

        if ( bonmin_problem->branch_info != NULL )
        {
            interfaceTMINLP->getBranchingInfo()->size                = bonmin_problem->branch_info->size;
            interfaceTMINLP->getBranchingInfo()->priorities          = bonmin_problem->branch_info->priorities;
            interfaceTMINLP->getBranchingInfo()->branchingDirections = bonmin_problem->branch_info->branchingDirections;
            interfaceTMINLP->getBranchingInfo()->upPsCosts           = bonmin_problem->branch_info->upPsCosts;
            interfaceTMINLP->getBranchingInfo()->downPsCosts         = bonmin_problem->branch_info->downPsCosts;

        }


        tminlp = interfaceTMINLP;

        bonmin_problem->bonmin_setup.initialize( GetRawPtr( tminlp ) );
        Bab branch_and_bound;
        branch_and_bound( bonmin_problem->bonmin_setup );
    }
    catch(TNLPSolver::UnsolvedError *E)
    {
        std::cerr<<"Bonmin has failed to solve a problem"<<std::endl;
    }
    catch(OsiTMINLPInterface::SimpleError &E)
    {
        std::cerr<<E.className()<<"::"<<E.methodName()
        <<std::endl
        <<E.message()<<std::endl;
    }
    catch(CoinError &E)
    {
        std::cerr<<E.className()<<"::"<<E.methodName()
        <<std::endl
        <<E.message()<<std::endl;
    }
    
    return 0;
}
