// Copyright (C) 2007, 2008 International Business Machines
// Corporation and others.  All Rights Reserved.
//
// Author:  Andreas Waechter      2007-08-20    IBM
//

#include "BonStrongBranchingSolver.hpp"

namespace Bonmin {

StrongBranchingSolver::StrongBranchingSolver(OsiTMINLPInterface * tminlp_interface)
{
  jnlst_ = tminlp_interface->solver()->journalist();
  DBG_ASSERT(IsValid(jnlst_));
  options_ = tminlp_interface->solver()->options();
  DBG_ASSERT(IsValid(options_));
  reg_options_ = tminlp_interface->solver()->roptions();
  DBG_ASSERT(IsValid(reg_options_));

  options_->GetIntegerValue("bb_log_level", bb_log_level_, tminlp_interface->prefix());
}

StrongBranchingSolver::StrongBranchingSolver(const StrongBranchingSolver & rhs)
{
  jnlst_ = rhs.jnlst_;
  options_ = rhs.options_;
  reg_options_ = rhs.reg_options_;
  bb_log_level_ = rhs.bb_log_level_;
}

StrongBranchingSolver &
StrongBranchingSolver::operator=(const StrongBranchingSolver & rhs)
{
  if (this != &rhs) {
    jnlst_ = rhs.jnlst_;
    options_ = rhs.options_;
    reg_options_ = rhs.reg_options_;
    bb_log_level_ = rhs.bb_log_level_;
  }
  return *this;
}

StrongBranchingSolver::~StrongBranchingSolver ()
{}

}/* Ends Bonmin's namespace.*/
