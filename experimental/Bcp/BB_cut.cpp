// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "BCP_buffer.hpp"
#include "BB_cut.hpp"

/****************************************************************************/

void
BB_cut::pack(BCP_buffer& buf) const
{
    buf.pack(OsiRowCut::lb())
	.pack(OsiRowCut::ub());
    const CoinPackedVector& v = OsiRowCut::row();
    const int numElem = v.getNumElements();
    buf.pack(v.getIndices(), numElem)
	.pack(v.getElements(), numElem);
}

/****************************************************************************/
BB_cut::BB_cut(BCP_buffer& buf) :
    BCP_cut_algo(-1e40, 1e40), OsiRowCut()
{
    double lb, ub;
    buf.unpack(lb)
	.unpack(ub);
    OsiRowCut::setLb(lb);
    OsiRowCut::setUb(ub);

    int numElem;
    int* indices;
    double* elements;
    buf.unpack(indices, numElem, true)
	.unpack(elements, numElem, true);
    OsiRowCut::setRow(numElem, indices, elements);

    if(numElem > 0) {
	delete[] indices;
	delete[] elements;
    }
}

/****************************************************************************/
BB_cut::BB_cut(const OsiRowCut& cut) :
    BCP_cut_algo(cut.lb(), cut.ub()), OsiRowCut(cut)
{}

/****************************************************************************/
BCP_MemPool BB_cut::memPool(sizeof(BB_cut));
