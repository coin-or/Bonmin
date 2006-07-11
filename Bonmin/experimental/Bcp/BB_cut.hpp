// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// File from COIN project Examples/Bac/BB_cut.hpp
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

/****************************************************************************/

static inline void
BB_pack_cut(const BCP_cut_algo* cut, BCP_buffer& buf)
{
    const BB_cut* bb_cut = dynamic_cast<const BB_cut*>(cut);
    if (!bb_cut)
	throw BCP_fatal_error("pack_cut_algo() : unknown cut type!\n");
    bb_cut->pack(buf);
}

static inline BCP_cut_algo*
BB_unpack_cut(BCP_buffer& buf)
{
    return new BB_cut(buf);
}

/****************************************************************************/

#endif
