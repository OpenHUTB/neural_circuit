
/* module ncsetvar in program nc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include "nc.h"
#include "ncval.h"
#include "ncio.h"
#include "ncsub.h"
#include "ncsetvar.h"

#define isdigit(c)      ('0' <= (c) && (c) <= '9')

#define MAXVARIABLES 3000
static stype variables[MAXVARIABLES] = {0};
static int nvar = 0;
int help=0;

int varset = 0;
extern stype varval[VVSIZE];

#ifdef __cplusplus
extern "C" {
#endif
	double atof(const char *ptr);
	void exit(int n);
#ifdef __cplusplus
}
#endif

/*----------------------------------------------*/

void cmdlin(int argc, char **argv)

/* Use this procedure to get the commmand line arguments 
   for a program compiled without the nc simulator */

{
   char *cptr;
   char *progname;

 progname = basename(argv[0]);
 do					/* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {
        case 'h': help=1;			/* print help table */
                  break;

        case 'r': argc--; argv++;
		  //drseed = setint(atof(*argv));    /* set random seed */
		  //if (drseed ==0) drseed = 1;
                  break;

        case '-': cptr++;                      /* set variable */
                  varval[varset].name = cptr;  /* remember variable name */
                  argc--; argv++;
                  if (argc) {
		   if (!isdigit(**argv) && ! (**argv == '-' && isdigit(*((*argv)+1)))
		   			&& ! (**argv == '.' && isdigit(*((*argv)+1)))
		   			&& ! (**argv == '-' && *((*argv)+1)=='.' &&
								isdigit(*((*argv)+2)))
			) {
                     varval[varset].cptr = *argv;	/*remember command line str*/
                     varval[varset].type = VSTRING;
                   }
                   else {
                     varval[varset].cptr = *argv;	/*remember command line str*/
                     varval[varset].val  = atof(*argv); /* remember value */
                     varval[varset].type = VNUMBER;
                   }
		  }
                  varset++;
                  if (varset >= VVSIZE) varset = VVSIZE-1;
                  break;

	default:
		ncfprintf (stderr,"%s: unknown switch '%s'\n",progname,*argv);
		exit(1);

        }  /* switch */
      }	  /* if (*cptr=) */
    }  	 /* if (argc) */
  }	/* do */
 while (argc > 0);
}

/*----------------------------------------------*/

void setptr (const char *name, int *iptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VINT;
   variables[nvar].iptr = iptr;
   *iptr = MININT;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

void setptr (const char *name, const char **sptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VSTRING;
   variables[nvar].sptr = (char **)sptr;
   *sptr = NULL;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

void setptr (const char *name, float *fptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VFLOAT;
   variables[nvar].fptr = fptr;
   *fptr = LARGENODE;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

void setptr (const char *name, double *dptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VDOUBLE;
   variables[nvar].dptr = dptr;
   *dptr = LARGENUM;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

/*---------------------------------------------------------*/

void setptrn (const char *name, int *iptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VINT;
   variables[nvar].iptr = iptr;
   // *iptr = MININT;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

void setptrn (const char *name, const char **sptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VSTRING;
   variables[nvar].sptr = (char **)sptr;
   // *sptr = NULL;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

void setptrn (const char *name, float *fptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VFLOAT;
   variables[nvar].fptr = fptr;
   // *fptr = LARGENODE;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

void setptrn (const char *name, double *dptr)
{
   variables[nvar].name = name;
   variables[nvar].type = VDOUBLE;
   variables[nvar].dptr = dptr;
   // *dptr = LARGENUM;
   if (++nvar>=MAXVARIABLES) ncfprintf(stderr,"setptr: too many variables\n");
}

/*---------------------------------------------------------*/

/* Set the pointers to the variables.      */
/* Available types are VINT, VDOUBLE, VFLOAT, VSTRING. */

/* For each new variable that you would like to set from */
/* the command line, add a line in "initptr()" to set a */
/* pointer to it. Then put a call to "setvar()" in your C */
/* code where you would like to have the command line variables set */
/* (by default they are set inside "ncinit()". The pointers must be */
/* set as in this example before "ncinit()" is called. */

/* Example of how to set pointers for command-line variables: */
/* 
void initptr ()
{
        setptr("code",      &codestr);	     // a string variable
        setptr("ntbins",    &ntbins);        // an integer variable
        setptr("file1",     &file1);         // a string variable
        setptr("file2",     &file2);         // a string variable
        setptr("spikthresh",&spikthresh);    // a double variable
        setptr("vstep",     &vstep);	     // a double variable
        setptr("intervals", &intervals);     // a double variable
        setptr("printresp", &printresp);     // an integer variable
        setptr("info",      &info);          // an integer variable
}
*/

/*---------------------------------------------------------*/

double getvarval(const char *name);

void setvariable (stype *var)

/* Procedure to set variables from command line */
/* that have been setup using "setptr()" or "setptrn()" */
/* above. */

{
    int found,i;
    Symbol *s;
    double val;

  for (found=i=0; i<nvar; i++) {
     if (strcmp(variables[i].name,var->name)==0) {
        found = 1;
        break;
     }
  }
    if (found) {
      if (var->type==VSTRING) {				/* command line value is string */
	  switch (variables[i].type) {
	     case VSTRING: *variables[i].sptr = var->cptr; break;

	     case VDOUBLE: 				   /* Here we have variable is number */ 
	     case VFLOAT: 				   /*   but string value on command line */
	     case VINT:  
	          val = getvarval(var->cptr);              /* look in "command line" variable list */
	          if (val==LARGENUM) {
	              s = lookup(var->cptr);               /* look in orig symbol table */
	              if (s) val = s->val;
		      else {
			 ncfprintf(stderr,"# Missing def for var '%.20s' val '%.20s' in command line.\n",
					 		variables[i].name,var->cptr);
			 break;
		      }
		 }
		  switch (variables[i].type) {
	                  case VDOUBLE: *variables[i].dptr = val; break;
	                  case VFLOAT:  *variables[i].fptr = val; break;
	                  case VINT:    *variables[i].iptr = int(val); break;
	         }; break;
	  }
      } else {						/* command line value is number */
	  switch (variables[i].type) {
	     case VSTRING: *variables[i].sptr = var->cptr; break;
	     case VDOUBLE: *variables[i].dptr = var->val; break;
	     case VFLOAT:  *variables[i].fptr = var->val; break;
	     case VINT:    *variables[i].iptr = int(var->val); break;
	  }
      }
      variables[i].vtype = var->type;
      variables[i].cptr = var->cptr;	/* remember orig command-line str */
   }
}

/*---------------------------------------------------------*/

void setexptvar ()

/* procedure to set "expt" variable from command line */

{
    int found,i,j;
    stype *var;
  
  for (j=0; j<varset; j++) {
    var = &varval[j];
    if (strcmp("expt",var->name)==0) {
	setvariable(var);    
	break;
    }
  }
}

/*---------------------------------------------------------*/

double getvarval(const char *name)

/* get value of a variable set on command line */

{
    int found,i;
    double retval;

  for (found=i=0; i<nvar; i++) {
     if (strcmp(variables[i].name,name)==0) {
        found = 1;
        break;
     }
  }
  if (found) {
     if (variables[i].type==VSTRING) retval = LARGENUM;
      else {
       if      (variables[i].type==VDOUBLE) retval = ((double)*variables[i].dptr);
       if      (variables[i].type==VFLOAT) {
 		retval = ((double)*variables[i].fptr);
		if (retval==LARGENODE) retval = LARGENUM;
       }
       else if (variables[i].type==VINT) {
		retval = ((double)*variables[i].iptr);
		if (retval==MININT) retval = LARGENUM;
       }
      }
  }
  else retval = LARGENUM;
  return retval;
}

  
/*---------------------------------------------------------*/

void setvarstr(const char *name)

/* set type of a numeric variable from command line to string */

{
    int found,i;

  for (found=i=0; i<nvar; i++) {
      if (strcmp(variables[i].name,name)==0) {
           found = 1;
           break;
      }
  }
	  /* check to see if command line var should be string */
  if (found) {
     if ((variables[i].vtype == VNUMBER) &&
	  (variables[i].type == VSTRING)) 
            *variables[i].sptr = variables[i].cptr;
  }
}

  
/*---------------------------------------------------------*/

int setvariables()

{
   int i;

   for (i=0; i<varset; i++) {
         setvariable(&varval[i]);
   }
   return i;
}

/*---------------------------------------------------------*/

int setvars(int argc, char **argv) 

/* Set variables in your C program (not for interpreter scripts) */

{

   cmdlin(argc,argv);
   return setvariables();
}

/*---------------------------------------------------------*/

