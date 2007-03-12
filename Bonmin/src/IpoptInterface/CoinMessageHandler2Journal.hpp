//Copyright (C) GAMS Development 2007
// All Rights Reserved.
// This code is published under the Common Public License.
// Authors :
// Stefan Vigerske, Humboldt-University Berlin. 
//
// Date : 04/01/2007 

#ifndef __COINMESSAGEHANDLER2JOURNAL_HPP__
#define __COINMESSAGEHANDLER2JOURNAL_HPP__

//#include "IpoptConfig.h"
#include "IpTypes.hpp"
#include "IpReferenced.hpp"
#include "IpSmartPtr.hpp"
#include "IpJournalist.hpp"
#include "CoinMessageHandler.hpp"

//#ifdef HAVE_CSTDARG
//# include <cstdarg>
//#else
//# ifdef HAVE_STDARG_H
//#  include <stdarg.h>
//# else
//#  error "don't have header file for stdarg"
//# endif
//#endif

#include <string>
//#include <vector>

using namespace Ipopt;

/** An Ipopt Journal that writes to a CoinMessageHandler.
 *  If no CoinMessageHandler is set, it writes to standard output.
 */
class CoinMessageHandler2Journal : public Journal {
public:
  /** Constructor. */
  CoinMessageHandler2Journal(CoinMessageHandler* messagehandler_, const std::string& name, EJournalLevel default_level);

  /** Destructor. */
  virtual ~CoinMessageHandler2Journal();

		 void setMessageHandler(CoinMessageHandler* messagehandler_) {
		 		 messagehandler=messagehandler_;
		 }		 		 		 

protected:
  /** Print to the designated output location */
  virtual void PrintImpl(EJournalCategory category, EJournalLevel level, const char* str);

  /** Printf to the designated output location */
  virtual void PrintfImpl(EJournalCategory category, EJournalLevel level, const char* pformat, va_list ap);

  /** Flush output buffer.*/
  virtual void FlushBufferImpl();

private:
  /** Default Constructor */
  CoinMessageHandler2Journal();

  /** Copy Constructor */
  CoinMessageHandler2Journal(const FileJournal&);

  /** Overloaded Equals Operator */
  void operator=(const CoinMessageHandler2Journal&);

		 /** Message handler to hand on the output to.
		  */
		 CoinMessageHandler* messagehandler;
};

#endif // __COINMESSAGEHANDLER2JOURNAL_HPP__


