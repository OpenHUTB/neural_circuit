#include "stdplt.h"
#include "string.h"


void _cat(register struct _frame *dotp);
char *gfrname(void);

FRAME *_find(char **afname)
                		/* address of frame name ptr */
/*
 * The character pointer pointed to by afname points to a desired frame name.
 * The character pointer is updated during the search; on return, it is left
 * pointing at the end of the last-matched portion of the name.
 * The return value points to the frame structure corresponding to
 * the matched portion of the name (or "." if no match).
 * The frame returned will be properly concatenated and updated.
 */
  {
	register FRAME *dotp;
	FRAME *dotdotp;
	register char *fname, *mp;

	fname = *afname;		/* desired frame name */
	dotp = (*fname == '/')? (fname++, &_rootx): _dotp;
	while (dotp && *fname)
	  {	if ((mp = _fmatch(".", fname)))
		  {	fname = mp;
			continue;
		  }
		if ((mp = _fmatch("..", fname)))
		  {	fname = mp;
			dotp = dotp->_parent;
			continue;
		  }
		dotdotp = dotp;
		for (dotp = dotp->_child; dotp; dotp = dotp->_sibling)
		  {	if ((mp = _fmatch(dotp->_fname, fname)))
			  {	fname = mp;
				if (dotp->_fflags&_DOCAT)
					_cat(dotp);
				break;
			  }
		  }
	  }
	*afname = fname;
	return(dotp? dotp: dotdotp);
  }

char *gfrname()

{
#define STRBUFSIZ 40

   static char str[STRBUFSIZ],tstr[STRBUFSIZ];
   FRAME *fpnt;

  str[0]=tstr[0]=0;
  for (fpnt=_dotp; fpnt != &_rootx; fpnt=fpnt->_parent) {
     strncpy (str,fpnt->_fname,STRBUFSIZ-1);
     strncat (str,tstr,STRBUFSIZ-1);
     strncpy (tstr,"/",STRBUFSIZ-1);
     strncat (tstr,str,STRBUFSIZ-1);
     strncpy (str,tstr,STRBUFSIZ-1);
  }
return str;
}
