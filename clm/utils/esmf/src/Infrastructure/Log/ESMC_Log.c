/* $Id$ */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <sys/time.h>

#include <unistd.h>

#include "ESMF_conf.h"

#include "ESMC_Log.h"
#include "ESMC_Error.h"

#ifdef ESMC_HAVE_MPI
#include <mpi.h>
#endif


#ifdef ESMC_HAVE_PTHREADS
#include <pthread.h>
#endif

#ifdef ESMC_HAVE_OMP_THREADS
#include <omp.h> 
#endif


#define MAX_LOGS 10

#include "ESMC_Machine.h"
#include "ESMC_BasicUtil.h"

static char *lognames[] = {
  "INFO",
  "ERROR",
  "TIMING"
};

static char report_matrix[ESMC_LOG_STATES][ESMC_LOG_LEVELS] =
{
  /*            INFO,ERROR,TIMING */
  /* QUIET   */ {  0,    0,     0  },
  /* NORMAL  */ {  1,    1,     1  },
  /* TIMER   */ {  0,    0,     1  },
  /* VERBOSE */ {  1,    1,     1  }
};


/* Static global variables */

static int num_mpi_processes;
static int mpi_rank;

static int last_logstate = -1;

struct timeval start_tv;

/* Non public function prototypes */
int ESMC_LogSetConfigV(ESMC_Log log, ESMC_LogOption option1, va_list args);
int ESMC_LogSetConfigCV(ESMC_Log log, ESMC_LogOption option1, va_list args);

/* These constants are used to find the first string length argument
   in the stdarg list.  They should be smaller than any valid variable
   address and larger than any passed string length.  These constants
   must be used since the length of the format string is not known.  */
#define ESMC_LOG_MIN_ADD 2048
#define ESMC_LOG_MAX_STR 2048

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogGetConvBase"
static char ESMC_LogGetConvBase(char *base_char, char conv_type)
{
  char *syns[] =
  {
    /* First character is the base for all others */
    "diouxX",
    "fFeEgG",
    "s",
    "c",
    NULL
  };
  
  int i;

  for (i = 0; syns[i]; i++)
    {
      if (strchr(syns[i], conv_type))
	{
	  *base_char = syns[i][0];
	  return ESMC_SUCCESS;
	}
    }

  *base_char = '\0';
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogSubChar"
static int ESMC_LogSubChar(char *str, char c)
{
  /* Substitute c for the next two chars, and
     shift string back */
  
  char *p = str;
  
  /* Substitute */
  *p = c;
  /* Advance p to where the \? was, then to next char */
  p += 2;
  for ( ; *p; p++)
    *(p - 1) = *p;
  
  /* Copy zero */
  *(p-1) = *p;

  return ESMC_SUCCESS;
}

#define NUM_INDENT 40
/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogSubCtrls"
static int ESMC_LogSubCtrls(char *str)
{
  /* Substitute \n \t if necessary */
  /* also, substitute '\n' with '\n >' */

  char *p, *pe, *pe_sav;
  typedef enum {gstate = 0, slsh = 1} parse_state;
#ifdef ESMC_SUBSTITUTE_CTRL_CHARS
  parse_state pstate;
  pstate = gstate;
  for (p = str; *p; p++)
    {
      switch(pstate)
	{
	case gstate:
	  if (*p == '\\')
	    pstate = slsh;
	  break;
	case slsh:
	  switch (*p)
	    {
	    case 'n':
	      ESMC_LogSubChar(--p, '\n');
	      break;
	    case 't':
	      ESMC_LogSubChar(--p, '\t');
	      break;

	      /* ...others ... */

	      
	    default:
	      /* False alarm */
	      break;
	    }

	  /* Switch out of state */
	  pstate = gstate;
	  break;
	}
    }
#endif


#ifdef NOTDEF
  pe = strchr(str, '\0');
  for ( p = str; *p; p++)
    {
      if (*p == '\n')
	{
	  for ( pe_sav = (pe + NUM_INDENT); pe > p; pe--)
	    {
	      *(pe + NUM_INDENT) = *pe;
	    }
	  for ( pe = (pe + NUM_INDENT) ; pe > p + 1; pe--)
	    {
	      *pe = ' ';
	    }
	  *pe = '>';
	  pe--;
	  
	  pe = pe_sav;
	}
    }
#endif

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogPrintOut"
static int ESMC_LogPrintOut(int *plen, char *obuf, char len_mod, char conv_type, char *pbuf,
			    void **next_arg, int slen)
{
  char buf[ESMC_LOG_MAX_STR];
  
  if (conv_type == 0)
    {
      sprintf(obuf, pbuf);
      *plen = 0;
      return ESMC_SUCCESS;
    }

  if (slen >=0)
    {
      strncpy(buf, (char*) next_arg, slen);
      buf[slen] = '\0';
      ESMC_LogSubCtrls(buf);
      next_arg = (void*) &buf[0];
      *plen = sprintf(obuf, pbuf, (char*) next_arg);
    }
  else
    {
      /*
      fprintf(stderr, pbuf, *next_arg);
      */
    
    /* Compilers actually parse the format string and
       care about the type sent in to printf */
      switch(conv_type)
	{
	case 'd':
	  if (len_mod == 'l')
	    *plen = sprintf(obuf, pbuf, (long int) *next_arg);
	  else if (len_mod == 'h')
	    *plen = sprintf(obuf, pbuf, *((short int*)next_arg));	    
	  else
	    *plen = sprintf(obuf, pbuf, *((int*)next_arg));
	  break;
	case 'f':
	  if (len_mod == 'l') 
	    *plen = sprintf(obuf, pbuf, *((double*)next_arg));

	  else
	    *plen = sprintf(obuf, pbuf, *((float*)next_arg));	    
	  break;
	}
    }
    return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogOutbuf"
static int ESMC_LogOutbuf(ESMC_Log log, ESMC_LogLevel priority, char *obuf)
{
  char buf[1024];
  char tbuf[1024];
  FILE *fp;
  int node, process, thread;
  char *cptr, *cptr1; 
  struct timeval tv;
  double tstamp;

  ESMC_MachinePInfo(&node, &process, &thread);
  
  if (log->labelio)
    fp = log->log_fp[mpi_rank];
  else
    fp = log->log_fp[0];


  /* This check should be done up front */
  gettimeofday(&tv, NULL);
  tstamp = (tv.tv_sec - start_tv.tv_sec) + 1e-6*(tv.tv_usec - start_tv.tv_usec);

  sprintf(buf, ">> ");
  
  if (log->printtime)
    {
      sprintf(tbuf, "Time=%010.5f,", tstamp);
      strcat(buf, tbuf);
    }
  if (log->printpid)
    {
      sprintf(tbuf, "PID=%d,", process);
      strcat(buf, tbuf);
    }
  if (log->printtid)
    {
      sprintf(tbuf, "TID=%d,", thread);
      strcat(buf, tbuf);
    }
  if (1)
    {
      sprintf(tbuf, "%-6s,", lognames[priority]);
      strcat(buf, tbuf);
    }

  /* Convert the last comma to a colon to end list */
  buf[strlen(buf) - 1] = ':';
  
  if (priority != last_logstate)
    {
      /* Insert an automatic newline and print banner */
      fprintf(fp, "\n");
      last_logstate = priority;
      fprintf(fp, buf);
    }
  
  for (cptr = strchr(obuf, '\n'), cptr1 = obuf; cptr; cptr = strchr(cptr, '\n'))
    {
      *cptr = '\0';
      fprintf(fp, "\n");
      fprintf(fp, buf);
      fprintf(fp, "%s", cptr1);
      cptr++;;
      cptr1 = cptr;
    }
  
  fprintf(fp, "%s", cptr1);

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogPrintC"
int ESMC_LogPrintC(ESMC_Log log, ESMC_LogLevel priority, char *format,
		   va_list args)
{
  char buf[1024];

  ESMC_BasicUtilLockMutex();
  if (!report_matrix[log->logstate][priority])
    {
      /* Do not print message if priority is not right. */
      ESMC_BasicUtilUnlockMutex();
      return ESMC_SUCCESS;
    }
  
  vsprintf(buf, format, args);

  ESMC_LogOutbuf(log, priority, buf);

  ESMC_BasicUtilUnlockMutex();

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogPrintf"
int ESMC_LogPrintf(ESMC_Log log, ESMC_LogLevel priority, char *format,
		      va_list len_args, va_list args)
{
  int i, done;
  unsigned int flen;
  int plen;
  char *p, *p1;
  char base_conv;
  char buf[ESMC_LOG_MAX_STR];
  char pbuf[ESMC_LOG_MAX_STR];
  char obuf[ESMC_LOG_MAX_STR];
  int obuflen;
  
  typedef enum {
    ground = 0,
    percent = 1
  } parse_state;

  parse_state pstate;

  char *lmod = "hl";
  void *next_arg;
  char len_mod, conv_type;
  int slen;

  ESMC_BasicUtilLockMutex();

  if (!report_matrix[log->logstate][priority])
    {
      /* Do not print message if priority is not right. */
      ESMC_BasicUtilUnlockMutex();
      return ESMC_SUCCESS;
    }
  
  /* Initialize output buffer */
  obuf[0] = '\0';
  obuflen = 0;
  
  done = 0;
  i = 0;
  do {
    flen = va_arg(len_args, int);
    /*  printf("arg[%d] = 0x%0x\n", i, flen); */
    
    i++;
    if (flen < ESMC_LOG_MIN_ADD)
      {
	/* printf("Suspect end: flen = %d\n", flen); */
	done = 1;
      }
  } while (!done);

  /* printf("Length of format:%d\n", flen); */

  strncpy(buf, format, flen);

  buf[flen] = '\0';

  ESMC_LogSubCtrls(buf);
  
  /* printf("Buf is:<%s>\n", buf); */

  /* Now on to parsing formats */
  pstate = ground;
  for ( p = buf, p1 = &pbuf[0]; *p; p++, p1++)
  {
    /* Copy char into print buffer */
    *p1 = *p;
    /*
      printf("parsing, pbuf=<%s>, *p=<%c>, *p1=<%c>", pbuf, *p, *p1);
    */
     

    switch (pstate)
      {
      case ground:
	if (*p == '%')
	  {
	    pstate = percent;
	    len_mod = '\0';
	    conv_type = '\0';
	  }
	
	break;

      case percent:
	/* Looking for a lenght, junk, or
	   type */
	if (strchr(lmod, *p))
	  {
	    len_mod = *p;
	  }
	else
	  {
	    ESMC_LogGetConvBase(&base_conv, *p);
	    
            conv_type = base_conv;
	    if (base_conv)
	      {
		*(p1 + 1) = '\0';
		next_arg = va_arg(args, void*);
		if (conv_type == 's')
		  slen = va_arg(len_args, int);
		else
		  slen = -1;
		ESMC_LogPrintOut(&plen, &obuf[obuflen], len_mod, conv_type,
				 pbuf, next_arg, slen);
		obuflen += plen;
		pstate = ground;
		p1 = &pbuf[0];
		
		p1--; /* end of loop will p1++ */
		
	      }
	  }
	
	/* Else we don't care about this char.*/
	
	break;
      }
  }
  
  /* NULL terminate */
  *p1 = *p;
  /* Print end. */
  conv_type = 0;
  ESMC_LogPrintOut(&plen, &obuf[obuflen], 0, conv_type,
		   pbuf, NULL, 0);
  obuflen += plen;

  ESMC_LogOutbuf(log, priority, obuf);

  ESMC_BasicUtilUnlockMutex();

  return ESMC_SUCCESS;
  
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogSetConfig"
int ESMC_LogSetConfig(ESMC_Log log, ESMC_LogOption option1, ...)
{
  va_list args;
  int ret;
  
  va_start(args, option1);
  ret = ESMC_LogSetConfigCV(log, option1, args);
  va_end(args);

  return ret;
    
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogSetConfigCV"
int ESMC_LogSetConfigCV(ESMC_Log log, ESMC_LogOption option1, va_list args)
{
  ESMC_LogOption option;
  ESMC_LogOption ip;
  int ival;

  option = option1;
  while (option != 0)
    {
      ival = va_arg(args, int);
      switch (option)
	{
	case ESMC_LOGCONFIG_PRINTTID:
	  log->printtid = ival;
	  break;
	case ESMC_LOGCONFIG_PRINTPID:
	  log->printpid = ival;
	  break;
	case ESMC_LOGCONFIG_PRINTTIME:
	  log->printtime = ival;
	  break;
	case ESMC_LOGCONFIG_LOGSTATE:
	  ESMC_LogSetState(log, (ESMC_LogState) ival);
	  break;
	default:
	  ESMC_ERRA1(ESMC_ERR_ARG_OUTOFRANGE, 0,
		     "Bad SetConfig option:%d\n", option);	  
	}
      ip = va_arg(args, ESMC_LogOption);
      if (ip)
	option = ip;
      else /* End of list */
	option = 0;
    }
  return ESMC_SUCCESS;
}


/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogSetConfigV"
int ESMC_LogSetConfigV(ESMC_Log log, ESMC_LogOption option1, va_list args)
{
  ESMC_LogOption option;
  ESMC_LogOption *iptr;
  int *ival;

  option = option1;
  while (option != 0)
    {
      ival = va_arg(args, int*);
      switch (option)
	{
	case ESMC_LOGCONFIG_PRINTTID:
	  log->printtid = *ival;
	  break;
	case ESMC_LOGCONFIG_PRINTPID:
	  log->printpid = *ival;
	  break;
	case ESMC_LOGCONFIG_PRINTTIME:
	  log->printtime = *ival;
	  break;
	case ESMC_LOGCONFIG_LOGSTATE:
	  ESMC_LogSetState(log, (ESMC_LogState) *ival);
	  break;
	default:
	  ESMC_ERRA1(ESMC_ERR_ARG_OUTOFRANGE, 0,
		     "Bad SetConfig option:%d\n", option);	  
	}
      iptr = va_arg(args, ESMC_LogOption*);
      if (iptr)
	option = *iptr;
      else /* End of list */
	option = 0;
    }
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogPrint"
int ESMC_LogPrint(ESMC_Log log, ESMC_LogLevel priority, char *format, ...)
{
  va_list args;
  
  va_start(args, format);
  
  ESMC_LogPrintC(log, priority, format, args);
  
  va_end(args);

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogLogNew"
int ESMC_LogNew(ESMC_Log *log, char *logname, ESMC_LogState logstate, int labelio)
{
  int rc;

  *log = (ESMC_Log) malloc(sizeof(ESMC_LogClass));

  rc = ESMC_LogConstruct(*log, logname, logstate, labelio);

  if (rc != ESMC_SUCCESS)
  {
     free(*log);
     *log = 0;
     return rc;
  }

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogConstruct"
int ESMC_LogConstruct(ESMC_Log log, char *logname, ESMC_LogState logstate, int labelio)
{
  char *c;
  int i;
  char fname[1024];

  gettimeofday(&start_tv, NULL);
  
#ifdef ESMC_HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#else
  num_mpi_processes = 1;
#endif

  log->logstate = logstate;
  log->labelio = labelio;

  /* Set defaults for unspecified options */
  log->printtid = 1;
  log->printpid = 1;
  log->printtime = 1;
 
  /* Set up the file descriptor(s) */
  if (labelio)
    {
      log->log_fp = (FILE**) malloc(sizeof(FILE*) * num_mpi_processes);
      log->num_fp = num_mpi_processes;
      for (i = 0; i < num_mpi_processes; i++)
	{
	  sprintf(fname, "%s.%d", logname, i);
	  log->log_fp[i] = fopen(fname, "a");
	  if (!log->log_fp[i])
	    {
	      ESMC_ERRA1(ESMC_ERR_FILE_OPEN, 0, "Could Not open Logfile:%s\n", fname);
	    }
	}
    }
  else
    {
      log->log_fp = (FILE**) malloc(sizeof(FILE*));
      log->num_fp = 1;
      /* If logfile, use that name otherwise just pipe to stdout */
      
      log->log_fp[0] = fopen(logname, "a");
      
      if (!log->log_fp[0])
	{
	  ESMC_ERRA1(ESMC_ERR_FILE_OPEN, 0, "Could Not open Logfile:%s\n", fname);
	}
    }

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogDelete"
int ESMC_LogDelete(ESMC_Log log)
{
  int rc;

  rc = ESMC_LogDestruct(log);

  free(log);

  return rc;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogDestruct"
int ESMC_LogDestruct(ESMC_Log log)     
{
  int i;
  /* Close all logfiles */

  for (i = 0; i < log->num_fp; i++)
    {
      fclose(log->log_fp[i]);
    }

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogSetState"
int ESMC_LogSetState(ESMC_Log log, ESMC_LogState logstate)
{
  log->logstate = logstate;

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_LogFlush"
int ESMC_LogFlush(ESMC_Log log)
{
  int i;

  for (i = 0; i < log->num_fp; i++)
    fflush(log->log_fp[i]);

  return ESMC_SUCCESS;
}

