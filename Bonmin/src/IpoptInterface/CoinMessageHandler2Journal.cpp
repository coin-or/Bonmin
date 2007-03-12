// Copyright (C) GAMS Development 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Stefan Vigerske, Humboldt-University Berlin. 
//
// Date : 04/01/2007

#include "CoinMessageHandler2Journal.hpp"

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#if defined(_WIN32)
# define vsnprintf _vsnprintf
# define  snprintf  _snprintf
#endif

CoinMessageHandler2Journal::CoinMessageHandler2Journal(CoinMessageHandler* messagehandler_, const std::string& name, EJournalLevel default_level)
: Journal(name, default_level), messagehandler(messagehandler_)
{ }

CoinMessageHandler2Journal::~CoinMessageHandler2Journal() { }

/** Print to the designated output location */
void CoinMessageHandler2Journal::PrintImpl(EJournalCategory category, EJournalLevel level, const char* str) {
  //TODO one should set the detail level of the message according to the journalist print level
 
  if (messagehandler) {
    int length=strlen(str);
    if (!length) return;
      if (length>=1000) {
        // longer then messagehandler buffer size, so we print directly to a file pointer or stdout
        messagehandler->finish();
        if (messagehandler->filePointer())
          fprintf(messagehandler->filePointer(), "%s", str);
        else 
        printf("%s", str);
      }
      else if (str[length-1]=='\n') {
      // if this message ends with a newline, we remove the newline and let the messagehandler finish it
      // the messagehandler will then add the newline again (wonderful complicated, isn't it?)
      const_cast<char*>(str)[length-1]=0;
      *messagehandler << str;
       messagehandler->finish();
       } else if (strlen(messagehandler->messageBuffer())+length<990) {
      // try to cache the message, so we can wait until a newline
      *messagehandler << str;		 		 		 
    } 
    else 
    {
     // message buffer is (almost) full; so we have to print here
     messagehandler->finish();
     *messagehandler << str;
    }
  } 
  else {
     printf("%s", str);
  }
}

/** Printf to the designated output location */
void CoinMessageHandler2Journal::PrintfImpl(EJournalCategory category, EJournalLevel level, const char* pformat, va_list ap) {
  char outBuf[1024];

#ifdef HAVE_VA_COPY
		 va_list apcopy;
		 va_copy(apcopy, ap);
  vsnprintf(outBuf, sizeof(outBuf), pformat, apcopy);
		 va_end(apcopy);
#else
  vsnprintf(outBuf, sizeof(outBuf), pformat, ap);
#endif

		 PrintImpl(category, level, outBuf);
}

/** Flush output buffer.*/
void CoinMessageHandler2Journal::FlushBufferImpl() {
		 if (messagehandler && messagehandler->filePointer())
		 		 fflush(messagehandler->filePointer());		 		 
}

