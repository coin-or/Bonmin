
// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

//Test MySVN
#include "SepaHeuristicInnerApproximation.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"
#include "BonSubMipSolver.hpp"
#include "BonCbcLpStrategy.hpp"

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif

#include "OsiClpSolverInterface.hpp"

#include "OsiAuxInfo.hpp"

#include "CoinTime.hpp"

#include <fstream>

#include <iomanip>

#include "CoinHelperFunctions.hpp"

#define DEBUG_BON_HEURISTIC

namespace Sepa{

HeuristicInnerApproximation::HeuristicInnerApproximation(Bonmin::BonminSetup * setup) :
  CbcHeuristic(), setup_(setup), howOften_(100), mip_(NULL),
        nbAp_(50) {
  Initialize(setup);
}

HeuristicInnerApproximation::HeuristicInnerApproximation(
    const HeuristicInnerApproximation &copy) :
    CbcHeuristic(copy),
        setup_(copy.setup_), 
        howOften_(copy.howOften_), 
        mip_(new Bonmin::SubMipSolver(*copy.mip_)),
        nbAp_(copy.nbAp_) {
}

HeuristicInnerApproximation &
HeuristicInnerApproximation::operator=(const HeuristicInnerApproximation& rhs) {
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    setup_ = rhs.setup_;
    howOften_ = rhs.howOften_;
                nbAp_ = rhs.nbAp_;
    delete mip_;
    if (rhs.mip_)
      mip_ = new Bonmin::SubMipSolver(*rhs.mip_);
  }
  return *this;
}

void HeuristicInnerApproximation::registerOptions(Ipopt::SmartPtr<
    Bonmin::RegisteredOptions> roptions) {
  roptions->SetRegisteringCategory("Initial Approximations descriptions",
      Bonmin::RegisteredOptions::UndocumentedCategory);
  roptions->AddStringOption2("heuristic_inner_approximation",
      "if yes runs the InnerApproximation heuristic", "yes", "no",
      "don't run it", "yes", "runs the heuristic", "");

        roptions->setOptionExtraInfo("heuristic_inner_approximation", 63);
}

void
HeuristicInnerApproximation::Initialize(Bonmin::BonminSetup * b) {

   delete mip_;
   mip_ = new Bonmin::SubMipSolver (*b, "inner_approximation");
   b->options()->GetIntegerValue("number_approximations_initial_outer",
       nbAp_, b->prefix());
}

HeuristicInnerApproximation::~HeuristicInnerApproximation() {
delete mip_;
}

/** Returns a feasible solution to the MINLP 
 * The heuristic constructs a MIP based approximating all univariate functions appearing in nonlinear constraints 
 * The linear approximation is obtained by adding inner chords linking pairs of points until covering the range of each variable **/
int
HeuristicInnerApproximation::solution(double &solutionValue, double *betterSolution)
{
if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
if ((model_->getNodeCount()%howOften_)!=0||model_->getCurrentPassNumber()>1)
return 0;

int returnCode = 0; // 0 means it didn't find a feasible solution

Bonmin::OsiTMINLPInterface * nlp = NULL;
if(setup_->getAlgorithm() == Bonmin::B_BB)
nlp = dynamic_cast<Bonmin::OsiTMINLPInterface *>(model_->solver()->clone());
else
nlp = dynamic_cast<Bonmin::OsiTMINLPInterface *>(setup_->nonlinearSolver()->clone());

Bonmin::TMINLP2TNLP* minlp = nlp->problem();
// set tolerances

//double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

int numberColumns;
int numberRows;
int nnz_jac_g;
int nnz_h_lag;
Ipopt::TNLP::IndexStyleEnum index_style;
minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
    nnz_h_lag, index_style);

//const Bonmin::TMINLP::VariableType* variableType = minlp->var_types();

const double* x_sol = minlp->x_sol();

double* newSolution = new double [numberColumns];
memcpy(newSolution,x_sol,numberColumns*sizeof(double));
double* new_g_sol = new double [numberRows];

bool feasible = true;
// load the problem to OSI
#ifdef DEBUG_BON_HEURISTIC
std::cout << "Loading the problem to OSI\n";
#endif
OsiSolverInterface *si = mip_->solver(); // the MIP solver

bool delete_si = false;
if(si == NULL) {
  si = new OsiClpSolverInterface;
  mip_->setLpSolver(si);
  delete_si = true;
}
CoinMessageHandler * handler = model_->messageHandler()->clone();
si->passInMessageHandler(handler);
si->messageHandler()->setLogLevel(2);
#ifdef DEBUG_BON_HEURISTIC
std::cout << "Loading problem into si\n";
#endif
extractInnerApproximation(*nlp, *si, newSolution, true); // Call the function construncting the inner approximation description 
#ifdef DEBUG_BON_HEURISTIC
std::cout << "problem loaded\n";
std::cout << "**** Running optimization ****\n";
#endif
mip_->optimize(DBL_MAX, 2, 180); // Optimize the MIP
#ifdef DEBUG_BON_HEURISTIC
std::cout << "Optimization finished\n";
#endif
if(mip_->getLastSolution()) { // if the MIP solver returns a feasible solution
  const double* solution = mip_->getLastSolution();
  for (size_t iLCol=0;iLCol<numberColumns;iLCol++) {
    newSolution[iLCol] = solution[iLCol];
  }
}
else
feasible = false;

if(delete_si) {
  delete si;
}
delete handler;

// const double* x_l = minlp->x_l();
// const double* x_u = minlp->x_u();
// const double* g_l = minlp->g_l();
// const double* g_u = minlp->g_u();
// double primalTolerance = 1.0e-6;

#if 0 // Set to 1 if you need to test the feasibility of the returned solution
vector<Ipopt::TNLP::LinearityType>  constTypes(numberRows);
minlp->get_constraints_linearity(numberRows, constTypes());
feasible = true;
for (int iColumn=0;iColumn<numberColumns;iColumn++) {
  double value=newSolution[iColumn];
  if(value - x_l[iColumn] < -1e-8|| value - x_u[iColumn] > 1e-8) {
#ifdef DEBUG_BON_HEURISTIC
    cout<<"It should be infeasible because: "<<endl;
    cout<<"x_l["<<iColumn<<"]= "<<x_l[iColumn]<<" "
    <<"x_sol["<<iColumn<<"]= "<<value<<" "
    <<"x_u["<<iColumn<<"]= "<<x_u[iColumn]<<endl;
#endif
    feasible = false;
    break;
  }
  if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
    if (fabs(floor(value+0.5)-value)>integerTolerance) {
#ifdef DEBUG_BON_HEURISTIC
      cout<<"It should be infeasible because: "<<endl;
      cout<<"x["<<iColumn<<"]= "<<value<<" Should be integer!"<<"\\Integer tolerance = "<<integerTolerance<<"\n";
#endif
      feasible = false;
      break;
    }
  }
}
minlp->eval_g(numberColumns, newSolution, true,
    numberRows, new_g_sol);
for(int iRow=0; iRow<numberRows; iRow++) {
  if(new_g_sol[iRow]<g_l[iRow]-primalTolerance ||
      new_g_sol[iRow]>g_u[iRow]+primalTolerance) {
    if(minlp->optimization_status() != SUCCESS) {
      feasible = false;
      break;
    } else {
#ifdef DEBUG_BON_HEURISTIC
      cout<<"It should be infeasible because: "<<endl;
      cout<<"g_l["<<iRow<<"]= "<<g_l[iRow]<<" "
      <<"g_sol["<<iRow<<"]= "<<new_g_sol[iRow]<<" "
      <<"g_u["<<iRow<<"]= "<<g_u[iRow]<<endl;
      cout<<"primalTolerance= "<<primalTolerance<<endl;
      if(constTypes[iRow] == TNLP::NON_LINEAR)
      cout<<"nonLinear constraint number "<<iRow<<endl;
#endif
      feasible = false;
    }
  }
}
#ifdef DEBUG_BON_HEURISTIC
cout<<"Every thing is feasible"<<endl;
#endif


//#else // fix integer variables and solve the NLP
#if 0
if(feasible ) {

  std::vector<double> memLow(numberColumns);
  std::vector<double> memUpp(numberColumns);
  std::copy(minlp->x_l(), minlp->x_l() + numberColumns, memLow.begin());
  std::copy(minlp->x_u(), minlp->x_u() + numberColumns, memUpp.begin());
  // fix the integer variables and solve the NLP
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
        double value=floor(newSolution[iColumn]+0.5);
        minlp->SetVariableUpperBound(iColumn, value);
        minlp->SetVariableLowerBound(iColumn, value);
  }
  if(feasible) {
    nlp->initialSolve();
    if(minlp->optimization_status() != SUCCESS) {
      feasible = false;
    }
    memcpy(newSolution,minlp->x_sol(),numberColumns*sizeof(double));
  }

 
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
        minlp->SetVariableUpperBound(iColumn, memUpp[iColumn]);
        minlp->SetVariableLowerBound(iColumn, memLow[iColumn]);
    }
  }
}
#endif
#endif

if(feasible) {
  double newSolutionValue;
  minlp->eval_f(numberColumns, newSolution, true, newSolutionValue);
  if(newSolutionValue < solutionValue) {
    memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
    solutionValue = newSolutionValue;
    returnCode = 1;
  }
}

delete [] newSolution;
delete [] new_g_sol;

delete nlp;

#ifdef DEBUG_BON_HEURISTIC
std::cout<<"Inner approximation returnCode = "<<returnCode<<std::endl;
#endif
return returnCode;
}

/** Get an inner-approximation constraint obtained by drawing a chord linking the two given points x and x2. 
 * This only applies to nonlinear constraints featuring univariate functions (f(x) <= y).**/
bool
HeuristicInnerApproximation::getMyInnerApproximation(Bonmin::OsiTMINLPInterface &si, OsiCuts &cs, int ind,
    const double * x, const double * x2) {

  int n, m, nnz_jac_g, nnz_h_lag;
  Ipopt::TNLP::IndexStyleEnum index_style;
  Bonmin::TMINLP2TNLP * problem = si.problem(); 
  problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);


  CoinPackedVector cut;
  double lb = 0;
  double ub = 0;

  double infty = si.getInfinity();

  lb = -infty; // we only compute <= constraints

  double g = 0;
  double g2 = 0;
  double diff = 0;
  double a = 0;
  problem->eval_gi(n, x, 1, ind, g);
  problem->eval_gi(n, x2, 1, ind, g2);
  Bonmin::vector<int> jCol(n);
  int nnz;
  problem->eval_grad_gi(n, x2, 0, ind, nnz, jCol(), NULL);
  Bonmin::vector<double> jValues(nnz);
  problem->eval_grad_gi(n, x2, 0, ind, nnz, NULL, jValues());
  bool add = false;
  //printf("const %i nnz %i\n", ind, nnz);
  for (int i = 0; i < nnz; i++) {
    const int &colIdx = jCol[i];
    if(index_style == Ipopt::TNLP::FORTRAN_STYLE) jCol[i]--;
    diff = x[colIdx] - x2[colIdx];

    if (fabs(diff) >= 1e-8) {
                   a = (g - g2) / diff;
                   cut.insert(colIdx, a);
                   ub = a * x[colIdx] - g;
                   //printf("const %i col %i p[col] %g pp[col] %g g %g g2 %g diff %g\n",ind, colIdx, x[colIdx], x2[colIdx], g, g2, diff);
                   if(a < 1e8)
                     add = true;
    } else {
                  cut.insert(colIdx, jValues[i]);
                  //printf("const %i col %i val %g\n",ind, colIdx, jValues[i]);
    }
  }

  if (add) {

    OsiRowCut newCut;
    newCut.setGloballyValidAsInteger(1);
    newCut.setLb(lb);
    
      //********* Perspective Extension ********//
    int binary_id = 0; // index corresponding to the binary variable activating the corresponding constraint
    const int* ids = problem->get_const_xtra_id(); // vector of indices corresponding to the binary variable activating the corresponding constraint
    // Get the index of the corresponding indicator binary variable
    binary_id = (ids == NULL) ? -1 : ids[ind];
    if(binary_id>0) {// If this hyperplane is a linearization of a disjunctive constraint, we link its righthand side to the corresponding indicator binary variable
        cut.insert(binary_id, -ub); // ∂x ≤ ub => ∂x - ub*z ≤ 0
        newCut.setUb(0);
    }
    else
        newCut.setUb(ub);
    //********* Perspective Extension ********//

    newCut.setRow(cut);
    cs.insert(newCut);
    //newCut.print();
    return true;
  }
  return false;
}

void 
HeuristicInnerApproximation::extractInnerApproximation(Bonmin::OsiTMINLPInterface & nlp, OsiSolverInterface &si,
  const double * x, bool getObj) {
   printf("************  Start extracting inner approx");
   int n;
   int m;
   int nnz_jac_g;
   int nnz_h_lag;
   Ipopt::TNLP::IndexStyleEnum index_style;
   Bonmin::TMINLP2TNLP * problem = nlp.problem(); 
   //Get problem information
   problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
   
   Bonmin::vector<int> jRow(nnz_jac_g);
   Bonmin::vector<int> jCol(nnz_jac_g);
   Bonmin::vector<double> jValues(nnz_jac_g);
   problem->eval_jac_g(n, NULL, 0, m, nnz_jac_g, jRow(), jCol(), NULL);
   if(index_style == Ipopt::TNLP::FORTRAN_STYLE)//put C-style
   {
     for(int i = 0 ; i < nnz_jac_g ; i++){
       jRow[i]--;
       jCol[i]--;
     }
   }
   
   //get Jacobian
   problem->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL,
       jValues());
   
   Bonmin::vector<double> g(m);
   problem->eval_g(n, x, 1, m, g());
   
   Bonmin::vector<int> nonLinear(m);
   //store non linear constraints (which are to be removed from IA)
   int numNonLinear = 0;
   const double * rowLower = nlp.getRowLower();
   const double * rowUpper = nlp.getRowUpper();
   const double * colLower = nlp.getColLower();
   const double * colUpper = nlp.getColUpper();
   assert(m == nlp.getNumRows());
   double infty = si.getInfinity();
   double nlp_infty = nlp.getInfinity();
   Bonmin::vector<Ipopt::TNLP::LinearityType>  constTypes(m);
   Bonmin::vector<Ipopt::TNLP::LinearityType>  varTypes(n);
   problem->get_constraints_linearity(m, constTypes());
   problem->get_variables_linearity(n, varTypes());
   for (int i = 0; i < m; i++) {
     if (constTypes[i] == Ipopt::TNLP::NON_LINEAR) {
       nonLinear[numNonLinear++] = i;
     }
   }
   Bonmin::vector<double> rowLow(m - numNonLinear);
   Bonmin::vector<double> rowUp(m - numNonLinear);
   int ind = 0;
   for (int i = 0; i < m; i++) {
     if (constTypes[i] != Ipopt::TNLP::NON_LINEAR) {
       if (rowLower[i] > -nlp_infty) {
         //   printf("Lower %g ", rowLower[i]);
         rowLow[ind] = (rowLower[i]);
       } else
         rowLow[ind] = -infty;
       if (rowUpper[i] < nlp_infty) {
         //   printf("Upper %g ", rowUpper[i]);
         rowUp[ind] = (rowUpper[i]);
       } else
         rowUp[ind] = infty;
       ind++;
     }
   
   }
   
   CoinPackedMatrix mat(true, jRow(), jCol(), jValues(), nnz_jac_g);
   mat.setDimensions(m, n); // In case matrix was empty, this should be enough
   
   //remove non-linear constraints
   mat.deleteRows(numNonLinear, nonLinear());
   
   int numcols = nlp.getNumCols();
   Bonmin::vector<double> obj(numcols);
   for (int i = 0; i < numcols; i++)
     obj[i] = 0.;
   
   si.loadProblem(mat, nlp.getColLower(), nlp.getColUpper(), 
                  obj(), rowLow(), rowUp());
   const Bonmin::TMINLP::VariableType* variableType = problem->var_types();
   for (int i = 0; i < n; i++) {
     if ((variableType[i] == Bonmin::TMINLP::BINARY) || (variableType[i] == Bonmin::TMINLP::INTEGER))
       si.setInteger(i);
   }
   if (getObj) {
     bool addObjVar = false;
     if (problem->hasLinearObjective()) {
       double zero;
       Bonmin::vector<double> x0(n, 0.);
       problem->eval_f(n, x0(), 1, zero);
       si.setDblParam(OsiObjOffset, -zero);
       //Copy the linear objective and don't create a dummy variable.
       problem->eval_grad_f(n, x, 1, obj());
       si.setObjective(obj());
     } else {
       addObjVar = true;
     }
   
     if (addObjVar) {
       nlp.addObjectiveFunction(si, x);
     }
   }
   
   // Hassan IA initial description
   int InnerDesc = 1;
   if (InnerDesc == 1) {
     OsiCuts cs;
   
     double * p = CoinCopyOfArray(colLower, n);
     double * pp = CoinCopyOfArray(colLower, n);
     double * up = CoinCopyOfArray(colUpper, n);
   
     for (int i = 0; i < n; i++) {
       if (p[i] < -1e3){
          p[i] = pp[i] = -1e3;
       }
       if (up[i] > 1e2){
          up[i] = 1e2;
       }
     } 

     const int& nbAp = nbAp_;
     printf("Generating approximation with %i points.\n", nbAp);
     std::vector<int> nbG(m, 0);// Number of generated points for each nonlinear constraint
   
     std::vector<double> step(n);
     int n_lin = 0;
   
     for (int i = 0; i < n; i++) {
       //if ((variableType[i] == Bonmin::TMINLP::BINARY) || (variableType[i] == Bonmin::TMINLP::INTEGER)) {
       if (varTypes[i] == Ipopt::TNLP::LINEAR) {
         n_lin ++;
         step[i] = 0;
         p[i] = pp[i] = up[i] = 0;
       }
       else {
         step[i] = (up[i] - p[i]) / (nbAp);
       }
     }
     printf("Number of linears %i\n", n_lin);
     Bonmin::vector<double> g_p(m);
     Bonmin::vector<double> g_pp(m);
   
     for (int j = 1; j < nbAp; j++) {
   
       for (int i = 0; i < n; i++) {
         pp[i] += step[i];
       }
   
       problem->eval_g(n, p, 1, m, g_p());
       problem->eval_g(n, pp, 1, m, g_pp());
       double diff = 0;
       for (int i = 0; (i < m ); i++) {
         if (constTypes[i] == Ipopt::TNLP::LINEAR) continue;
         diff = std::abs(g_p[i] - g_pp[i]);
         if (nbG[i] < nbAp - 1) {
           bool status = getMyInnerApproximation(nlp, cs, i, p, pp);// Generate a chord connecting the two points
           if(status == false)
             constTypes[i] = Ipopt::TNLP::LINEAR;
           nbG[i]++;
         }
       }
       std::copy(pp, pp+n, p);
      
     }
   
     for(int i = 0; (i< m); i++) {
         if (constTypes[i] == Ipopt::TNLP::LINEAR) continue;
         getMyInnerApproximation(nlp, cs, i, p, up);// Generate a chord connecting the two points
     }

        delete [] p; 
        delete [] pp;
        delete [] up; 
     si.applyCuts(cs);
   }
   //si.writeLp("IA");
   printf("************  Done extracting inner approx ********");
  }

}

