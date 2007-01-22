// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  07/01/2005

#include "BonAmplInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "BonTMINLP.hpp"
#include "BonAmplTMINLP.hpp"
#include "IpIpoptApplication.hpp"

#include "BonIpoptSolver.hpp"
#include "BonminConfig.h"

#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#endif

#include <cmath>

using namespace Bonmin;
/** Test function for the Osi interface to Ipopt (or any nlp solver). <br>
    If Solver passes all the test then it should have everything needed to be integrated into bonmin. */

void testGetMethods(OsiTMINLPInterface &si)
{
    CoinRelFltEq eq;// to test equality of doubles    
    std::cout<<"Checking get functions"<<std::endl;
      // Problem size
      assert(si.getNumCols()==4);
      assert(si.getNumRows()==3);
      
      //Check bounds on columns
      const double * colLow = si.getColLower();
      assert(eq(colLow[0],0.));
      assert(eq(colLow[1],0.));
      assert(eq(colLow[2],0.));
      assert(eq(colLow[3],0.));
      
      const double * colUp = si.getColUpper();
      assert(colUp[0]>si.getInfinity());
      assert(colUp[1]>si.getInfinity());
      assert(eq(colUp[2],1.));
      assert(eq(colUp[3],5.));      
      //Check bounds on rows
      const double * rowLow = si.getRowLower();
      assert(rowLow[0]<= -si.getInfinity());
      assert(rowLow[1]<= -si.getInfinity());
      assert(rowLow[2]<= -si.getInfinity());
                  
      const double * rowUp = si.getRowUpper();
      assert(eq(rowUp[0], 1./4.));
      assert(eq(rowUp[1], 0.));
      assert(eq(rowUp[2], 2.));

      //check objective sense
      assert(si.getObjSense()==1);
      
      // check variables types
      assert(si.isInteger(0)==0);
      assert(si.isInteger(1)==0);
      assert(si.isInteger(2)==1);
      assert(si.isInteger(3)==1);
      
      assert(si.isContinuous(0)==1);
      assert(si.isContinuous(1)==1);
      assert(si.isContinuous(2)==0);
      assert(si.isContinuous(3)==0);
      
      assert(si.isBinary(0)==0);
      assert(si.isBinary(1)==0);
      assert(si.isBinary(2)==1);
      assert(si.isBinary(3)==0);
      
      assert(si.isIntegerNonBinary(0)==0);
      assert(si.isIntegerNonBinary(1)==0);
      assert(si.isIntegerNonBinary(2)==0);
      assert(si.isIntegerNonBinary(3)==1);
      
      assert(si.isFreeBinary(2)==1);
      si.setColLower(2,1.);
      assert(si.isFreeBinary(2)==0);
      si.setColLower(2,0.);
      
      //assert(si.getInfinity()>1e50);
      std::cout<<"Test passed"<<std::endl;                  
}
void testOptimAndSolutionQuery(OsiTMINLPInterface & si)
{
    CoinRelFltEq eq(1e-07);// to test equality of doubles    
    std::cout<<"Testing optimization methods and solution query"<<std::endl;
    si.initialSolve();
    
    assert(si.isProvenOptimal());
//    assert(si.nCallOptimizeTNLP()==1);
    assert(si.getIterationCount()>0);
    // Optimum of the problem is -( 3/2 + sqrt(5)/2)
    // with x = (1/2 + sqrt(5) y[1]=x and y[2] = 1/2 + sqrt(5)/2
    // (can easily be computed since constraint x-y[1]<=0 imply x = y[1] and the resulting problem has dimension 2
    if(!eq(si.getObjValue(),-( (3./2.) + sqrt(5.)/2.)))
        std::cout<<"Error in objective : "<<fabs(si.getObjValue()+( (3./2.) + sqrt(5.)/2.))<<std::endl;
    
    //Test validity of primal solution
    const double * colsol = si.getColSolution();
    if(!eq(colsol[0],( (1./2.) + 1/sqrt(5.))))
        std::cout<<"Error for y[1]  : "<<fabs(colsol[0]-( (1./2.) + 1/sqrt(5.)))<<std::endl;
    if(!eq(colsol[1],( (1./2.) + 1/(2.*sqrt(5.)))))
        std::cout<<"Error for y[2]  : "<<fabs(colsol[1]-( (1./2.) + 1/(2*sqrt(5.))))<<std::endl;
    if(!eq(colsol[2],( (1./2.) + 1/sqrt(5.))))
        std::cout<<"Error for x  : "<<fabs(colsol[2]-( (1./2.) + 1/sqrt(5.)))<<std::endl;
    //value of z is not tested

    //Test for row activity
    const double * rowAct = si.getRowActivity();
    if(!eq(rowAct[0],1./4.))
        std::cout<<"Error for row activity of c1 : "<<fabs(rowAct[0]-1./4.)<<std::endl;
    if(!eq(rowAct[1],0.))
        std::cout<<"Error for row activity of c2 : "<<fabs(rowAct[1])<<std::endl;
        
     //Check dual values dual for c1 = sqrt(5) c2=1 c3 not tested
     const double * duals = si.getRowPrice();
     if(!eq(duals[0],sqrt(5)))
             std::cout<<"Error dual of c1 : "<<fabs(duals[0]-sqrt(5.))<<std::endl;
     if(!eq(duals[1],1.))
             std::cout<<"Error dual of c2 : "<<fabs(duals[0]-1.)<<std::endl;
             
     std::cout<<"Test passed successfully"<<std::endl;
}

///Test set methods
void testSetMethods(OsiTMINLPInterface &si)
{
    CoinRelFltEq eq(1e-07);// to test equality of doubles    
    si.setColLower(2,1.);
    assert(si.getColLower()[2]==1.);
    si.initialSolve();    
    assert(si.isProvenOptimal());
    assert(eq(si.getColSolution()[2],1));

    CoinWarmStart * ws = si.getWarmStart();
    
    
    si.setColLower(2,0.);
    
    si.setColUpper(2,0.);
    assert(si.getColUpper()[2]==0.);
    si.setWarmStart(ws);

    si.resolve();
    assert(si.isProvenOptimal());
    assert(eq(si.getColSolution()[2],0.));
    
    si.setColUpper(2,1.);
    delete ws;
}

void testOa(Bonmin::AmplInterface &si)
{
        CoinRelFltEq eq(1e-07);// to test equality of doubles    
    OsiClpSolverInterface lp;
    si.extractLinearRelaxation(lp);
    lp.writeMps("toy");
     assert(lp.getNumCols()==5);
      assert(lp.getNumRows()==4);
      //Check bounds on columns
      const double * colLow = lp.getColLower();
      assert(eq(colLow[0],0.));
      assert(eq(colLow[1],0.));
      assert(eq(colLow[2],0.));
      assert(eq(colLow[3],0.));
      
      const double * colUp = lp.getColUpper();
      assert(colUp[0]>=lp.getInfinity());
      assert(colUp[1]>=lp.getInfinity());
      assert(eq(colUp[2],1.));
      assert(eq(colUp[3],5.));      
      //Check bounds on rows
      const double * rowLow = lp.getRowLower();
      std::cout<<rowLow[0]<<"\t"<<lp.getInfinity()<<std::endl;
      assert(rowLow[0]<= -lp.getInfinity());
      assert(rowLow[1]<= -lp.getInfinity());
      assert(rowLow[2]<= -lp.getInfinity());
                  
      const double * rowUp = lp.getRowUpper();
      double sqrt5 = sqrt(5);
      if(!eq(rowUp[0], 1./2. + 3./(2 * sqrt5))){
	double error = fabs(rowUp[0] - 1./2. - 3./(2 * sqrt5));
	std::cout<<"Error in OA for rowUp[0]: "
		 <<error<<std::endl;
      }
      assert(eq(rowUp[1], 0.));
      assert(eq(rowUp[2], 2.));
      assert(eq(rowUp[3], 0.));
      

      //check objective sense
      assert(si.getObjSense()==1);
      
      // check variables types
      assert(si.isInteger(0)==0);
      assert(si.isInteger(1)==0);
      assert(si.isInteger(2)==1);
      assert(si.isInteger(3)==1);
    
       //Now check the full matrix
       const CoinPackedMatrix * mat = lp.getMatrixByCol();
       int  inds[11] = {0, 1, 3, 0, 2, 3, 1, 2, 3, 2, 3};
       double vals[11] = {2. / sqrt(5) , -1., -1., 1./sqrt(5), 1. , -1. , 1. , 1., -1.,1.,-1.};
       assert(mat->getNumElements()==11);
       int k=0;
       for(int i = 0 ; i < si.getNumCols() ; i++)
       {
        for(int j = mat->getVectorStarts()[i] ; j < mat->getVectorStarts()[i] + mat->getVectorLengths()[i] ; j++)
        {
        assert(inds[k]==mat->getIndices()[j]);
        if(!eq(vals[k],mat->getElements()[j])){
	double error = fabs(vals[k] - mat->getElements()[j]);
	std::cout<<"Error in OA for element of constraint matrix "<<k<<": "
		 <<error<<std::endl;
	if(error > 1e-06) throw -1;
      }
          k++;
        }
       }
}

void testFp(Bonmin::AmplInterface &si)
{
        CoinRelFltEq eq(1e-07);// to test equality of doubles
        OsiCuts cuts;
        double x[1] = {0.};
        int ind[1]={1};
        si.getFeasibilityOuterApproximation(1,x,ind,cuts);
        std::cout<<si.getColSolution()[0]<<std::endl;
         std::cout<<si.getColSolution()[1]<<std::endl;
       assert(eq(si.getColSolution()[1],(1./2.)));
}
void interfaceTest(Ipopt::SmartPtr<TNLPSolver> solver)
{
  /**********************************************************************************/
  /*   Test constructors                                                                                                              */
  /**********************************************************************************/
  std::cout<<"Test OsiTMINLPInterface with "
	   <<solver->solverName()<<" solver"<<std::endl;
  // Test usefull constructor
  {
        //read a toy problem and do various tests
//        var x binary;
//        var z integer >= 0 <= 5;
//        var y{1..2} >=0;
//        
//        
//        minimize cost:
//            - x - y[1] - y[2] ;
//            
//        subject to
//            c1: ( y[1] - 1/2 )^2 + (y[2] - 1/2)^2 <= 1/4 ;
//            c2: x - y[1] <= 0 ;
//            c3: x + y[2] + z <= 2;
        
        //Setup Ipopt should be replaced if solver is changed
       const char * args[3] ={"name","mytoy",NULL}; //Ugly, but I don't know how to do differently
       const char ** argv = args;
      AmplInterface si(const_cast<char **&> (argv), solver);
      si.Set_expose_warm_start(1);
    std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"
    <<std::endl<<"Testing usefull constructor"<<std::endl
    <<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      //Start of real tests
      testGetMethods(si);
      testOptimAndSolutionQuery(si);
      testSetMethods(si);
  }
  // Test copy constructor
  {
        //read a toy problem and do various tests
//        var x binary;
//        var z integer >= 0 <= 5;
//        var y{1..2} >=0;
//        
//        
//        minimize cost:
//            - x - y[1] - y[2] ;
//            
//        subject to
//            c1: ( y[1] - 1/2 )^2 + (y[2] - 1/2)^2 <= 1/4 ;
//            c2: x - y[1] <= 0 ;
//            c3: x + y[2] + z <= 2;
        
        //Setup Ipopt should be replaced if solver is changed
      const char * args[3] ={"name","mytoy",NULL}; //Ugly, but I don't know how to do differently
      const char ** argv = args;
      AmplInterface si1(const_cast<char **&> (argv), solver);
      si1.Set_expose_warm_start(1);
      
      OsiTMINLPInterface si(si1);
    std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"
    <<std::endl<<"Testing copy constructor"<<std::endl
    <<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      //Start of real tests
      testGetMethods(si);
      testOptimAndSolutionQuery(si);
      testSetMethods(si);
  }
  
    // Test outer approximation methods
  {
        //Setup Ipopt should be replaced if solver is changed
        const char * args[3] ={"name","mytoy",NULL}; //Ugly, but I don't know how to do differently
        const char ** argv = args;
      Bonmin::AmplInterface si(const_cast<char**&>(argv), solver);
      si.Set_expose_warm_start(1);
      std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"
	       <<std::endl<<"Testing outer approximations related methods"<<std::endl
	       <<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
      testOa(si);
  }
  
  // Test Feasibility Pump methods
//  {
//    //Setup Ipopt should be replaced if solver is changed
//    using namespace Ipopt;
//    SmartPtr<Ipopt::IpoptApplication> app = new Ipopt::IpoptApplication();
//    char * args[3] ={"name","toy3",NULL}; //Ugly, but I don't know how to do differently
//    char ** argv = args;
//    SmartPtr<TMINLP> ampl_tminlp = new AmplTMINLP(ConstPtr(app->Jnlst()),  app->Options(), argv);
//    Bonmin::AmplInterface si(ampl_tminlp);
//    std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"
//	     <<std::endl<<"Testing optimization of some distance over feasible set"<<std::endl
//	     <<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
 //   testFp(si);
//  }
  std::cout<<"All test passed successfully"<<std::endl;
} 

int main()
{
  Ipopt::SmartPtr<IpoptSolver> ipopt_solver = new IpoptSolver;
  interfaceTest(GetRawPtr(ipopt_solver));

#ifdef COIN_HAS_FSQP
  Ipopt::SmartPtr<FilterSolver> filter_solver = new FilterSolver;
  interfaceTest(GetRawPtr(filter_solver));
#endif
  return 0;
}
