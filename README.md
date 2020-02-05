Bonmin
======

Bonmin (Basic Open-source Nonlinear Mixed INteger programming) is an
open-source code for solving general MINLP (Mixed Integer NonLinear
Programming) problems.
It builds on top of [Cbc](https://github.com/coin-or/Cbc)
and [Ipopt](https://github.com/coin-or/Ipopt).

It is a [COIN-OR](www.coin-or.org) project and licensed under the
EPL (Eclipse Public License). The EPL is a license approved by the OSI (Open Source
Initiative), thus Bonmin is OSI Certified Open Source Software.

This project was initiated in 2004 by IBM and Carnegie Mellon University
as part of a joint effort to study algorithm for MINLP.
You may find additional informations at http://egon.cheme.cmu.edu/ibm/page.htm,
in particular
 - A publicly available library of test instances of Convex MINLPs,
 - Research papers on new algorithmic procedures for MINLP.

Bonmin features several algorithms, including
 - B-BB is a NLP-based branch-and-bound algorithm,
 - B-OA is an outer-approximation decomposition algorithm,
 - B-QG is an implementation of  Quesada and Grossmann's branch-and-cut algorithm,
 - B-Hyb is a hybrid outer-approximation based branch-and-cut algorithm.

Bonmin documentation consists of a users' manual and a reference manual,
which may be found on-line at the project web-site or can be built
from the project source distribution.


BONMIN MAINTAINER : Pierre Bonami

BONMIN WEB-PAGE : https://github.com/coin-or/Bonmin

BONMIN MAILING LIST : http://list.coin-or.org/mailman/listinfo/bonmin

To report bugs, you should create an issue at
https://github.com/coin-or/Bonmin/issues

Bonmin is also available on [NEOS](https://neos-server.org).
