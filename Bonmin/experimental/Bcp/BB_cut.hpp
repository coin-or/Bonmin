// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#ifndef _BB_CUT_H
#define _BB_CUT_H

#include "BCP_cut.hpp"
#include "BCP_mempool.hpp"
#include "OsiRowCut.hpp"

class BCP_buffer;

/** Simple representation of a cut by storing non zero coefficients only */ 

/****************************************************************************/
class BB_cut : public BCP_cut_algo, public OsiRowCut {

private:

    static BCP_MemPool memPool;

public:

    static inline void * operator new(size_t size) {
	return memPool.alloc(size);
    }

    static inline void operator delete(void *p, size_t size) {
	memPool.free(p, size);
    }

    /// Packing cut to a buffer
    void pack(BCP_buffer& buf) const;

    /**@name Constructors and destructors */
    //@{
    /// Constructor from content of buffer 
    BB_cut(BCP_buffer& buf);

    /// Constructor from an OsiRowCut 
    BB_cut(const OsiRowCut& cut);

    /// Destructor
    ~BB_cut() {}
};

#endif
