#include <stdarg.h>
#include "stdplt.h"

#define TWOPI	(2.*PI)
#define CWMIN	16.		/* 64. /* min char width in std pixels */
#define NROT	32768.		/* number of rotation steps */

int _clip();

int textf(const char *fmt, ...)

/* Simulate printf at current pen position. */

  {
	va_list argp;
	register FRAME *dotp;
	register unsigned cw;
	register FILE *rstdplt;
	double crot;
	extern double fmod(double da, double db);

	va_start(argp,fmt);
	dotp = _dotp;
	rstdplt = stdplt;
	if (_clip(&dotp->_xcat,&dotp->_ycat,&dotp->_xcat,&dotp->_ycat, _rwind) < 0)
		return(1);
	if (dotp->_fflags & _DOWIND)
		if (_clip(&dotp->_x, &dotp->_y, &dotp->_x, &dotp->_y, dotp->_wind) < 0)
			return(1);
	if ((dotp->_fflags&_TXTMODE) == 0)
	  {	dotp->_fflags |= _TXTMODE;	/* have to force a move */
		_domove();		/* last action wasn't textf() */
	  }
	if (dotp->_fat != _fat)
	  {	putc(_CMD|'f', rstdplt);
		putwsx(_fat = dotp->_fat, rstdplt);
	  }
	if (dotp->_fflags&_DOTXT)	/* character setting has changed */
	  {	dotp->_fflags &= ~_DOTXT;
		putc(_CMD|'s', rstdplt);
		cw = dotp->_cw/CWMIN + 0.5;
		putc(cw, rstdplt);
		putc(dotp->_font, rstdplt);
		crot = fmod(dotp->_rotcat + dotp->_crot, 2.*PI);
		putwsx((unsigned) ((crot*NROT/TWOPI) + 0.5), rstdplt);
	  }
	putc(_CMD|'t', rstdplt);

	/* _doprnt(fmt, &args, rstdplt);  		/* original */

	/* fprintf(rstdplt,fmt,arg1,arg2,arg3,arg4,arg5,arg6);  /* */

	vfprintf(rstdplt,fmt,argp);  /* */

	putc('\0', rstdplt);
	dotp->_fflags |= _TXTMODE;    /* retain pen posn in successive calls */
	return(0);
  }
