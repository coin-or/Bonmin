// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _BM_VAR_H
#define _BM_VAR_H

#include "BCP_buffer.hpp"
#include "BCP_var.hpp"
#include "BCP_mempool.hpp"

class BCP_buffer;

//#############################################################################
    
class BM_branching_var : public BCP_var_algo {
private:
    static BCP_MemPool memPool;
public:
    static inline void * operator new(size_t size) {
	return memPool.alloc(size);
    }
    static inline void operator delete(void *p, size_t size) {
	memPool.free(p, size);
    }

public:
    int sos_index;
    int split;

public:
    BM_branching_var(int i, int s) : BCP_var_algo(BCP_BinaryVar, 0, 0, 1),
				     sos_index(i), split(s) {}
    BM_branching_var(BCP_buffer& buf)  : BCP_var_algo(BCP_BinaryVar, 0, 0, 1) {
	buf.unpack(sos_index).unpack(split);
    }
    ~BM_branching_var() {}
    void pack(BCP_buffer& buf) const {
	buf.pack(sos_index).pack(split);
    }
};

//#############################################################################
    
static inline void
BM_pack_var(const BCP_var_algo* var, BCP_buffer& buf)
{
    const BM_branching_var* bv = dynamic_cast<const BM_branching_var*>(var);
    if (!bv)
	throw BCP_fatal_error("pack_var_algo() : unknown var type!\n");
    bv->pack(buf);
}

static inline BCP_var_algo*
BM_unpack_var(BCP_buffer& buf)
{
    return new BM_branching_var(buf);
}

//#############################################################################
    
#endif
