/****************************************************************
Copyright (C) 1997-2000 Lucent Technologies
Modifications for Coin -  Copyright (C) 2006, International Business Machines Corporation and others.
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

namespace ampl_utils
{
    void
sos_kludge(int nsos, int *sosbeg, double *sosref)
{
  // Adjust sosref if necessary to make monotonic increasing
  int i, j, k;
  double t, t1;
  for(i = j = 0; i++ < nsos; ) {
    k = sosbeg[i];
    t = sosref[j];
    while(++j < k) {
      t1 = sosref[j];
      t += 1e-10;
      if (t1 <= t)
        sosref[j] = t1 = t + 1e-10;
      t = t1;
    }
  }
}
}
