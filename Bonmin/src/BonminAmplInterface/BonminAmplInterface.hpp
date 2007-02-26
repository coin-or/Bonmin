#ifndef BonminAmplInterface_H
#define BonminAmplInterface_H
#include "IpoptInterface.hpp"
#include "AmplTMINLP.hpp"

/** Class for providing an Osi interface to Ipopt with an ampl nl file as input. */
class BonminAmplInterface: public IpoptInterface
{
  public:
  /** Default constructor */
  BonminAmplInterface();
  /** Constructor with inputed ampl command line (reads model from nl file)*/ 
  BonminAmplInterface(char **& amplArgs, bool = true);
  /** Copy constructor */
  BonminAmplInterface(const BonminAmplInterface &other);
  /// Clone
  virtual BonminAmplInterface * clone(bool CopyData = true);

  ///Destructor
  virtual ~BonminAmplInterface();

   /**@name Methods to input a problem */
  //@{
  /** Read an ampl . nl file from the given filename */
  virtual void readAmplNlFile(char**& filename,
      std::string* ipopt_file_content =NULL,
      std::string* nl_file_content = NULL,
      bool = true);
  /** write ampl solution file */
  void writeAmplSolFile(std::string message,const double * primalSol = NULL);
  //@}

  /** Fast access to AmplTMINLP */
   const Ipopt::AmplTMINLP * amplModel() const
  {
    return GetRawPtr(amplTminlp_);
  }
 
  protected:
   /** TMINLP problem (the original problem usually an AmplTMINLP).*/
  Ipopt::SmartPtr<Ipopt::AmplTMINLP> amplTminlp_;
};

#endif
