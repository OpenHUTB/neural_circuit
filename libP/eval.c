#include "stdplt.h"

void eval(double x, double y, char *fname, float *pnx, float *pny, char *nfname)
 
/*
 * The point (x,y) in frame 'fname' is evaluated as a position in
 * the frame 'nfname', and returned to the float variables indicated
 * by pnx and pny.
 */
  {
	register FRAME *fp, *nfp;
	register float *ctp;
	float itcat[2][2];
	char *tmpfname;
	float rx, ry;
	float junk;

	tmpfname = fname;
	fp = _find(&tmpfname);
	if (*tmpfname) {
		_err("eval: %s not found\n", _cfname(fname));
	}
	tmpfname = nfname;
	nfp = _find(&tmpfname);
	if (*tmpfname) {
		_err("eval: %s not found\n", _cfname(nfname));
	}
	if (pnx == 0) pnx = &junk;
	if (pny == 0) pny = &junk;
	ctp = (float *)fp->_tcat;
	rx = _T(ctp, x, y) + fp->_xocat - nfp->_xocat;
	ry = _T(ctp, x, y) + fp->_yocat - nfp->_yocat;
	if (_inv2(nfp->_tcat, itcat)) {
		_err("eval: %s has non-invertible transform\n", _cfname(nfname));
	}
	ctp = (float *)itcat;
	*pnx  = *ctp++ * rx;
	*pnx += *ctp++ * ry;
	*pny  = *ctp++ * rx;
	*pny += *ctp++ * ry;
  }
