/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#include <stdarg.h>

#include <string.h>

#include "ESMF_conf.h"

#include "ESMC_Timer.h"
#include "ESMC_Error.h"



static ESMC_Log STDLog = 0;
int ESMC_TimerInitV(char *name, ESMC_TimerOption option1, va_list args);
int ESMC_TimerInitCV(char *name, ESMC_TimerOption option1, va_list args);


#ifdef ESMC_HAVE_MPI
#include <mpi.h>
#endif


/*
** Required OMP calls not available on AIX so use PTHREADS
*/


/* Translate the ESMF defines into Rosinski's defines */

#ifdef ESMC_HAVE_OMP_THREADS
#define THREADED_OMP
#endif

#ifdef ESMC_HAVE_PTHREADS
#define THREADED_PTHREADS
#endif

#ifdef ESMC_HAVE_PCL
#define HAVE_PCL
#endif

/*
** Threading currently doesn't work on SUN so don't enable pthreads or omp
*/

#if ( defined THREADED_OMP )
#include <omp.h>
#elif ( defined THREADED_PTHREADS )
#include <pthread.h>
#endif

#ifdef HAVE_PCL
#include <pcl.h>
#else

/*
** Dummy up pcl stuff if library unavailable
*/   

typedef int PCL_CNT_TYPE;
typedef int PCL_FP_CNT_TYPE;
typedef int PCL_DESCR_TYPE;
#define PCL_MODE_USER       -1
#define PCL_L1DCACHE_MISS   -1
#define PCL_L2CACHE_MISS    -1
#define PCL_CYCLES          -1
#define PCL_ELAPSED_CYCLES  -1
#define PCL_FP_INSTR        -1
#define PCL_LOADSTORE_INSTR -1
#define PCL_INSTR           -1
#define PCL_STALL           -1
#define PCL_COUNTER_MAX      1
#define PCL_SUCCESS          0
#endif

#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif

#define STRMATCH(X,Y) (strcmp((X),(Y)) == 0)

#define AMBIGUOUS -1
#define MAX_THREADS 128

typedef enum {false = 0, true = 1} Boolean;

/*
** User specifiable options.  The values must match their counterparts in header.inc
** Also, we must have pcl_start < all valid pcl values < pcl_end.
** To add a new PCL counter: 
** 1) add the new entry to OptionName below.
** 2) add the appropriate array entry for possible_event[] to t_initialize.c.
** 3) add the appropriate code to the "switch" construct in t_initialize.c
*/

typedef enum {
  usrsys               = 1,
  wall                 = 2,
  pcl_start            = 3,   /* bogus entry delimits start of PCL stuff */
#ifdef HAVE_PCL
  pcl_l1dcache_miss    = 4,
  pcl_l2cache_miss     = 5,
  pcl_cycles           = 6,
  pcl_elapsed_cycles   = 7,
  pcl_fp_instr         = 8,
  pcl_loadstore_instr  = 9,
  pcl_instr            = 10,
  pcl_stall            = 11,
#endif
  pcl_end              = 12   /* bogus entry delimits end of PCL stuff */
} OptionName;

struct Event {
  OptionName name;
  char string[9];
  int index;
};

struct PossibleEvent {
  OptionName name;
  Boolean enabled;
  char string[10];
};

struct node {
  char name[ESMC_PROFILER_MAX_CHARS+1];
  
  int indent_level;        /* indentation level of timer */

  long last_utime;         /* user time from most recent call */
  long last_stime;         /* system time from most recent call */
  long last_wtime_sec;     /* wallclock seconds from most recent call */
  long last_wtime_usec;    /* wallclock microseconds from most recent call */

  long accum_utime;        /* accumulated user time */
  long accum_stime;        /* accumulated system time */
  long accum_wtime_sec;    /* accumulated wallclock seconds */
  long accum_wtime_usec;   /* accumulated wallclock microseconds */

  float max_wtime;         /* maximum wallclock time for each start-stop */
  float min_wtime;         /* minimum wallclock time for each start-stop */

  PCL_CNT_TYPE last_pcl_result[PCL_COUNTER_MAX];
  PCL_CNT_TYPE accum_pcl_result[PCL_COUNTER_MAX];

  Boolean onflg;           /* true => timer is currently on */
  long count;              /* number of calls to t_start for this timer */

  struct node *next;       /* next timer in the linked list */
};

/*
** Globals
*/

extern struct node **timers;
extern struct node **last;
extern struct Options options;
extern long ticks_per_sec;
extern int numthreads;
extern int *max_indent_level;
extern float *overhead;
extern PCL_CNT_TYPE *overhead_pcl;
extern Boolean t_initialized;
extern Boolean wallenabled;
extern Boolean usrsysenabled;
extern struct PossibleEvent possible_event[];

/*
** Needed by PCL library: otherwise unused
*/

extern int counter_list[PCL_COUNTER_MAX];
extern int ncounter;     /* number of counters */
extern int cycles;       /* index of cycle counter */
extern int instr;        /* index of instruction counter */
extern int fp_instr;     /* index of fp instruction counter */
extern int l2cache_miss; /* index of l2cache miss */
extern int jump;
extern PCL_DESCR_TYPE *descr;
extern int nevent;
extern struct Event **event;
extern Boolean pclenabled;
extern Boolean pcl_cyclesenabled;
extern int pcl_cyclesindex;
extern int npossible;


static int max_seen_thread = 0;


struct Stats {		   
  float usr;	   /* user CPU time */
  float sys;	   /* system CPU time */
  float usrsys;	   /* usr + sys */
  float elapse;	   /* elapsed time */
  float max_wtime; /* max elapsed wallclock time per call */
  float min_wtime; /* min elapsed wallclock time per call */

  long count;	   /* number of invocations of this timer */

  PCL_CNT_TYPE pcl_result[PCL_COUNTER_MAX];
};


/* Function prototypes (make static to avoid polluting namespace */

static int t_error (const char *fmt, ...);
static int get_cpustamp (long *usr, long *sys);
static int get_thread_num ();
static int lock_mutex ();
static int unlock_mutex ();
static int t_initialize ();
static char *t_pclstr (int code);
static int t_pr (ESMC_Log log, char *name, int procid);
static void fillstats (struct Stats *stats, struct node *ptr);
static void print_stats_line (ESMC_Log log, struct Stats *stats);
static void print_header (ESMC_Log log, int indent_level);
static int t_reset ();
static int t_setoption (OptionName option, Boolean val);
static int t_stamp (double *wall, double *usr, double *sys);
static int t_start (char *name);
static int t_stop (char *name);

#ifndef HAVE_PCL
static int PCLread (PCL_DESCR_TYPE descr, PCL_CNT_TYPE *i, PCL_CNT_TYPE *j, int k);
#endif


/*
** t_error: error return routine to print a message and return a failure
** value.
**
** Input arguments:
**   fmt: format string
**   variable list of additional arguments for vfprintf
**
** Return value: -1 (failure)
*/

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "t_error"
static int t_error (const char *fmt, ...)
{
  va_list args;
  char buf[1024];

  buf[0] ='\0';
#if ( ! defined DISABLE_TIMERS )
  va_start (args, fmt);

  if (fmt != NULL)
    (void) vsprintf (buf, fmt, args);

  va_end (args);

#endif
  
  ESMC_ERRA(ESMC_ERR_LIB, 0, buf);
}


/*
** get_cpustamp: Invoke the proper system timer and return stats.
**
** Output arguments:
**   usr: user time (usec if USE_GETRUSAGE is defined, ticks otherwise)
**   sys: system time (usec if USE_GETRUSAGE is defined, ticks otherwise)
**
** Return value: 0 (success)
*/

static int get_cpustamp (long *usr, long *sys)
{
  struct tms buf;

  /*
  ** Throw away the wallclock time from times: use gettimeofday instead
  */
  
  (void) times (&buf);
  *usr = buf.tms_utime;
  *sys = buf.tms_stime;

  return 0;
}


/*
** get_thread_num: Obtain logical thread number of calling thread.  If new
** thread, adjust global variables.
*/

static int get_thread_num ()
{
  int mythread = 0 ;  /* return value: default zero for non-threaded case */

  int proc, thread, node;

  ESMC_MachinePInfo(&node, &proc, &thread);

  mythread = thread;

  if (thread > max_seen_thread)
    max_seen_thread = thread;
  
  return mythread;
}

/*
** lock_mutex: lock a mutex for entry into a critical region
*/

/* 
** Array (1 per thread) of linked lists of timers, and last timer in each list
*/

struct node **timers = NULL;
struct node **last = NULL;

long ticks_per_sec;

/*
** Define lock arrays depending upon the type of threading done
*/

float *overhead;                   /* wallclock estimate of timer overhead */
int *max_indent_level;             /* maximum indentation level */
int numthreads            = 1;     /* number of threads.  1 is for no threading */
Boolean t_initialized     = false; /* whether t_initialize has been called */
Boolean wallenabled       = false; /* wallclock timer stats enabled */
Boolean usrsysenabled     = false; /* usr & sys timer stats enabled */
Boolean pclenabled        = false; /* enable PCL library */     
Boolean pcl_cyclesenabled = false; /* enable PCL cycle count */
int pcl_cyclesindex       = -1;    /* index for PCL cycle count */

struct PossibleEvent possible_event[] = {
  {usrsys,               true,  "Usr Sys   "},
  {wall,                 true,  "Wallclock "},
#ifdef HAVE_PCL
  {pcl_start,            false, "          "},  /* bracket PCL entries */
  {pcl_l1dcache_miss,    false, "l1 D miss "},
  {pcl_l2cache_miss,     false, "L2 miss   "},
  {pcl_cycles,           false, "Cycles    "},
  {pcl_elapsed_cycles,   false, "E-Cycles  "},
  {pcl_fp_instr,         false, "FP instr  "},
  {pcl_loadstore_instr,  false, "L/S instr "},
  {pcl_instr,            false, "Instruct  "},
  {pcl_stall,            false, "Stall     "},
  {pcl_end,              false, "          "},  /* bracket PCL entries */
#endif
};

struct Event **event = NULL;
int nevent = 0;
int npossible = sizeof (possible_event) / sizeof (struct PossibleEvent);

/*
** Needed by PCL library: otherwise unused
*/

PCL_DESCR_TYPE *descr;
int counter_list[PCL_COUNTER_MAX];
int ncounter = 0;                  /* number of PCL counters */
PCL_CNT_TYPE *overhead_pcl;        /* overhead counter (cycles) */




/*
** t_initialize (): Initialization routine must be called from single-threaded
**   region before any other timing routines may be called.  The need for this
**   routine could be eliminated if not targetting timing library for threaded
**   capability. 
**
** return value: 0 (success) or -1 (failure)
*/

static int t_initialize ()
{
  int n;             /* index */
  int nbytes;        /* number of bytes for malloc */
  int ret;           /* return code */

/*
** Determine number of ticks per second for conversion use by other t_pr(), t_stamp()
*/

  if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
    return t_error ("t_initialize: token _SC_CLK_TCK is not defined\n");

#if ( ! defined DISABLE_TIMERS )

  if (t_initialized)
    return t_error ("t_initialize has already been called\n");

  /*
  ** Allocate space for global arrays
  */

#if ( defined THREADED_OMP )
  numthreads = omp_get_max_threads();
#elif ( defined THREADED_PTHREADS )
  numthreads = MAX_THREADS;
#endif
  
  nbytes = numthreads * sizeof (struct node *);
  if ((timers = (struct node **) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  if ((last = (struct node **) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  nbytes = numthreads * sizeof (float);
  if ((overhead = (float *) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  nbytes = numthreads * sizeof (PCL_CNT_TYPE);
  if ((overhead_pcl = (PCL_CNT_TYPE *) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  nbytes = numthreads * sizeof (int);
  if ((max_indent_level = (int *) malloc (nbytes)) == 0)
    return t_error ("malloc failure for %d items\n", numthreads);

  /*
  ** Initialize array values
  */

  for (n = 0; n < numthreads; n++) {
    timers[n] = 0;
    last[n] = 0;
    overhead[n] = 0.;
    overhead_pcl[n] = 0;
    max_indent_level[n] = 0;
  }



  if (get_thread_num () > 0) 
    return t_error ("t_initialize: should only be called by master thread\n");

  for (n = 0; n < npossible; n++) {
    if (possible_event[n].enabled) {
      if (possible_event[n].name == usrsys)
	usrsysenabled = true;

      if (possible_event[n].name == wall)
	wallenabled = true;

      if ((event = realloc (event, (nevent+1) * sizeof (struct Event *))) == NULL)
	return t_error ("realloc failure\n");

      if ((event[nevent] = malloc (sizeof (struct Event))) == NULL)
	return t_error ("realloc failure\n");

      event[nevent]->name = possible_event[n].name;
      strcpy (event[nevent]->string, possible_event[n].string);

#ifdef HAVE_PCL

      /*
      ** Set up PCL stuff based on what t_setoption has provided.
      */

      if (event[nevent]->name > pcl_start && event[nevent]->name < pcl_end) {
	pclenabled = true;
	event[nevent]->index = ncounter;

	switch (possible_event[n].name) {

	case pcl_l1dcache_miss:
	  counter_list[ncounter++] = PCL_L1DCACHE_MISS;
	  break;
	  
	case pcl_l2cache_miss: 
	  counter_list[ncounter++] = PCL_L2CACHE_MISS;
	  break;
	  
	case pcl_cycles: 
	  pcl_cyclesindex = ncounter;
	  pcl_cyclesenabled = true;
	  counter_list[ncounter++] = PCL_CYCLES;
	  break;

	case pcl_elapsed_cycles: 
	  counter_list[ncounter++] = PCL_ELAPSED_CYCLES;
	  break;

	case pcl_fp_instr: 
	  counter_list[ncounter++] = PCL_FP_INSTR;
	  break;

	case pcl_loadstore_instr: 
	  counter_list[ncounter++] = PCL_LOADSTORE_INSTR;
	  break;

	case pcl_instr: 
	  counter_list[ncounter++] = PCL_INSTR;
	  break;

	case pcl_stall: 
	  counter_list[ncounter++] = PCL_STALL;
	  break;
	
	default:
	  break;

	}
      }
#endif
      ++nevent;
    }
  }

#ifdef HAVE_PCL

  if (ncounter > 0) {
    int thread;         /* thread number */

    nbytes = numthreads * sizeof (PCL_DESCR_TYPE);
    if ((descr = (PCL_DESCR_TYPE *) malloc (nbytes)) == 0)
      return t_error ("malloc failure: %d items\n", numthreads);

    /*
    ** PCLinit must be called on a per-thread basis.  Therefore must make the call here
    ** rather than in t_initialize.  null timer list flags not initialized.
    ** Also, the critical section is necessary because PCLstart appears not to be
    ** thread-safe.
    */

#pragma omp parallel for
    
    for (thread = 0; thread < numthreads; thread++) {

      unsigned int flags;           /* mode flags needed by PCL */

#pragma omp critical

      {
	if ((ret = PCLinit (&descr[thread])) != PCL_SUCCESS)
	  return t_error ("unable to allocate PCL handle for thread %d. %s\n",
			  thread, t_pclstr (ret));

	/*
	** Always count user mode only
	*/
      
	flags = PCL_MODE_USER;

	if ((ret = PCLquery (descr[thread], counter_list, ncounter, flags)) != PCL_SUCCESS)
	  return t_error ("Bad return from PCLquery thread %d: %s\n", thread, t_pclstr (ret));

	if ((ret = PCLstart (descr[thread], counter_list, ncounter, flags)) != PCL_SUCCESS)
	  return t_error ("PCLstart failed thread=%d: %s\n", thread, t_pclstr (ret));
      }
    }
  }
#endif
  t_initialized = true;
#endif

  return 0;
}


static char *t_pclstr (int code)
{

#if ( defined DISABLE_TIMERS )
  return "";
#endif

#ifdef HAVE_PCL
  switch (code) {

  case PCL_SUCCESS: 
    return "Success";
    
  case PCL_NOT_SUPPORTED:
    return "Event not supported";
    
  case PCL_TOO_MANY_EVENTS:
    return "Too many events";
    
  case PCL_TOO_MANY_NESTINGS:
    return "More nesting levels than allowed";
    
  case PCL_ILL_NESTING:
    return "Bad nesting";
    
  case PCL_ILL_EVENT:
    return "Illegal event identifier";
    
  case PCL_MODE_NOT_SUPPORTED:
    return "Mode not supported";
    
  case PCL_FAILURE:
    return "Failure for unspecified reason";
    
  default:
    return "Unknown error code";
    
  }
#endif
    return "Unknown error code";
}

      
  


/*
** t_pr: print stats for all known timers to a file
**
** Input arguments:
**   procid: Designed for MPI jobs to give a unique output file name, 
**     normally the MPI logical process id.
**
** Return value: 0 (success) or -1 (failure)
*/

static int mhz;            /* processor clock rate (from pcl library) */

static int t_pr (ESMC_Log log, char *name, int procid)
{
  char outfile[11];        /* name of output file: timing.xxx */
			   
  int indent;              /* index over number indentation levels */
  int thread;              /* thread index */
  int ilstart;             /* index over indentation level */
  int n;

  double gttdcost;         /* cost of a single gettimeofday call */
  double deltat;           /* time difference between 2 calls to gettimeofday */
  struct Stats stats;      /* per timer stats */
  struct Stats threadsum;  /* timer sum over threads */
			   
  struct node *ptr, *tptr; /* linked list pointers */
			   
  FILE *fp;                /* output file pointer */

  struct timeval tp1, tp2; /* input to gettimeofday() */
  struct tms buf;          /* input to times() */

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if ( ! t_initialized)
    return t_error ("t_pr: t_initialize has not been called\n");

  /*
  ** Only allow the master thread to print stats
  */

  if (get_thread_num () != 0)
    return 0;

  /*
  ** Estimate wallclock timer overhead: 4 times the cost of a call to gettimeofday
  ** (since each start/stop pair involves 4 such calls).
  */

  gettimeofday (&tp1, NULL);
  gettimeofday (&tp2, NULL);
  gttdcost = 1.e6*(tp2.tv_sec  - tp1.tv_sec) + (tp2.tv_usec - tp1.tv_usec);

  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Wallclock timer cost est.: %8.3f usec per start/stop pair\n", 
	      gttdcost*4.);

  /*
  ** CPU cost estimate: 2 times the cost of a call to times().  Subtract the
  ** cost of a single gettimeofday call to improve the estimate.
  */

  if (usrsysenabled) {

    gettimeofday (&tp1, NULL);
    (void) times (&buf);
    gettimeofday (&tp2, NULL);

    deltat = 1.e6*(tp2.tv_sec  - tp1.tv_sec) + (tp2.tv_usec - tp1.tv_usec);
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "CPU timer cost est.:       %8.3f usec per start/stop pair\n", 
		2.*deltat - gttdcost);
  }
    
  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "CPU accumulation interval is %g seconds\n",
           1./(float) ticks_per_sec);

#ifdef HAVE_PCL
  mhz = PCL_determine_mhz_rate();
  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Clock speed is %d MHz\n", mhz);
#endif

  for (thread = 0; thread < numthreads; thread++) {

    /*
    ** Only print heading for threads that have 1 or more items to report
    */

    if (timers[thread] == NULL) continue;
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "\nStats for thread %d:\n", thread);

    print_header (log, max_indent_level[thread]);

    for (ptr = timers[thread]; ptr != NULL; ptr = ptr->next) {
      if (ptr->onflg) {

	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Timer %s was not off.  No stats will be printed\n",
		 ptr->name);

      } else {

	fillstats (&stats, ptr);

	/*
	** Print stats indented.  If indent_level = AMBIGUOUS (a negative 
	** number) the name will be printed with no indentation.
	*/

	for (indent = 0; indent < ptr->indent_level; indent++)
	  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "  ");
	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%-15s", ptr->name);

	/*
	** If indent_level = AMBIGUOUS (a negative number) we want to loop 
	** from 0
	*/

	ilstart = MAX (0, ptr->indent_level);
	for (indent = ilstart; indent < max_indent_level[thread]; indent++)
	  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "  ");

	print_stats_line (log, &stats);
      }
    }

    if (usrsysenabled || wallenabled)
      ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "\nTIMER OVERHEAD (wallclock seconds) = %12.6f\n", 
	       overhead[thread]);

    if (pcl_cyclesenabled)
      ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "TIMER OVERHEAD (cycles) = %12.6e\n", 
	       (double) overhead_pcl[thread]);
  }

  /*
  ** Print a vertical summary if data exist for more than 1 thread.  The "2"
  ** passed to print_header is so we'll get an extra 4 spaces of indentation
  ** due to the thread number appended to the timer name.
  */

  if (numthreads > 0 && timers[1] != NULL) {
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "\nSame stats sorted by timer with thread number appended:\n");
    print_header (log, 2);

    /*
    ** Print stats for slave threads that match master
    */

    for (ptr = timers[0]; ptr != NULL; ptr = ptr->next) {

      char name[20];
      Boolean found = false;

      /*
      ** Don't bother printing summation stats when only the master thread
      ** invoked the timer
      */

      for (thread = 1; thread < numthreads; thread++)
	for (tptr = timers[thread]; tptr != NULL; tptr = tptr->next) {
	  if (STRMATCH (ptr->name, tptr->name))
	    found = true;
	}
      if ( ! found) continue;

      /*
      ** Initialize stats which sum over threads
      */

      memset (&threadsum, 0, sizeof (threadsum));

      if ( ! ptr->onflg) {
	fillstats (&stats, ptr);
	strcpy (name, ptr->name);
	strcat (name, ".0");
	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%-19s", name);
	print_stats_line (log, &stats);
	threadsum = stats;
      }

      /*
      ** loop over slave threads, printing stats for each and accumulating
      ** sum over threads when the name matches
      */

      for (thread = 1; thread < numthreads; thread++) {
	for (tptr = timers[thread]; tptr != NULL; tptr = tptr->next) {
	  if (STRMATCH (ptr->name, tptr->name)) {
	    if ( ! tptr->onflg) {
	      char num[5];

	      fillstats (&stats, tptr);
	      strcpy (name, tptr->name);
	      sprintf (num, ".%-3d", thread);
	      strcat (name, num);
	      ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%-19s", name);
	      print_stats_line (log, &stats);

	      threadsum.usr      += stats.usr;
	      threadsum.sys      += stats.sys;
	      threadsum.usrsys   += stats.usrsys;
	      threadsum.elapse   += stats.elapse;
	      threadsum.max_wtime = MAX (threadsum.max_wtime, stats.max_wtime);
	      threadsum.min_wtime = MIN (threadsum.min_wtime, stats.min_wtime);
	      threadsum.count    += stats.count;

	      for (n = 0; n < ncounter; n++)
		threadsum.pcl_result[n] += stats.pcl_result[n];
	    }
	    break; /* Go to the next thread */
	  }        /* if (STRMATCH (ptr->name, tptr->name) */
	}          /* loop thru linked list of timers for this thread */
      }            /* loop over slave threads */

      strcpy (name, ptr->name);
      strcat (name, ".sum");
      ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%-19s", name);
      print_stats_line (log, &threadsum);
      ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "\n");

    } /* loop through master timers */

    for (thread = 0; thread < max_seen_thread + 1; thread++) {
      if (usrsysenabled || wallenabled)
	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "OVERHEAD.%-3d (wallclock seconds) = %12.6f\n", 
		 thread, overhead[thread]);

      if (pcl_cyclesenabled)
	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "OVERHEAD.%-3d (cycles) = %12.6e\n", 
		 thread, (double) overhead_pcl[thread]);
    }

  } /* if (numthreads > 0 && timers[1] != NULL */

  return 0;
}

static void fillstats(struct Stats *stats, struct node *ptr)
{
  int n;

  stats->usr       = ptr->accum_utime / (float) ticks_per_sec;
  stats->sys       = ptr->accum_stime / (float) ticks_per_sec;
  stats->usrsys    = stats->usr + stats->sys;
  stats->elapse    = ptr->accum_wtime_sec + 1.e-6 * ptr->accum_wtime_usec;
  stats->max_wtime = ptr->max_wtime;
  stats->min_wtime = ptr->min_wtime;
  stats->count     = ptr->count;

  for (n = 0; n < ncounter; n++)
    stats->pcl_result[n] = ptr->accum_pcl_result[n];
}

static void print_stats_line(ESMC_Log log, struct Stats *stats)
{
  int index;
  int n;
  long long cycles;
  long long instr;
  long long flops;
  long long loadstore;
  long long l2cache;
  long long jump;

  float mflops;
  float ipc;
  float memfp;

  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%9ld ", stats->count);

  if (usrsysenabled)
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%9.3f %9.3f %9.3f ", stats->usr, stats->sys, stats->usrsys);

  if (wallenabled)
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%9.3f %9.3f %9.3f ", stats->elapse, stats->max_wtime, 
	     stats->min_wtime);
  
  for (n = 0; n < nevent; n++) {
    if (event[n]->name > pcl_start && event[n]->name < pcl_end) {
      index = event[n]->index;
      if (stats->pcl_result[index] > 1.e6)
	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%9.3e ", (double) stats->pcl_result[index]);
      else
	ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "%9ld ", (long) stats->pcl_result[index]);
    }
  }

  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "\n");
}

static void print_header (ESMC_Log log, int indent_level)
{
  int i;
  int n;

  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Name           ");
  for (i = 0; i < indent_level; i++)
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "  ");
  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Called    ");

  if (usrsysenabled)
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Usr       Sys       Usr+Sys   ");

  if (wallenabled)
    ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "Wallclock Max       Min       ");

  for (n = 0; n < nevent; n++)
    if (event[n]->name > pcl_start && event[n]->name <= pcl_end)
      ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, event[n]->string);

  ESMC_LogPrint(log, ESMC_LOGLEVEL_TIMER, "\n");
}


/*
** t_reset: reset all known timers to 0
**
** Return value: 0 (success) or -1 (failure)
*/

static int t_reset ()
{
  int n;             /* index over threads */
  struct node *ptr;  /* linked list index */

#if ( ! defined DISABLE_TIMERS )
  if ( ! t_initialized)
    return t_error ("t_reset: t_initialize has not been called\n");

  /*
  ** Only allow the master thread to reset timers
  */

  if (get_thread_num () != 0)
    return 0;

  for (n = 0; n < numthreads; n++) {
    for (ptr = timers[n]; ptr != NULL; ptr = ptr->next) {
      memset (timers[n], 0, sizeof (struct node));
      printf ("Reset accumulators for timer %s to zero\n", ptr->name);
    }
  }
#endif
  return 0;
}



/*
** t_setoption: set option value to true or false.
**
** Input arguments:
**   option: option to be set
**   val:    value to which option should be set
**
** Return value: 0 (success) or -1 (failure)
*/

static int t_setoption (OptionName option, Boolean val)
{
  int n;

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if (t_initialized)
    return (t_error ("t_setoption: Options must be set BEFORE t_initialize\n"));

  for (n = 0; n < npossible; n++) {
    if (possible_event[n].name == option) {
      possible_event[n].enabled = val;

      if (val)
	printf ("t_setoption: option enabled:  %s\n", possible_event[n].string);
      else
	printf ("t_setoption: option disabled: %s\n", possible_event[n].string);

      return 0;
    }
  }

  return (t_error ("t_setoption: Option with enum index %d not available\n",
		     option));
}

/*
** t_stamp: Compute timestamp of usr, sys, and wallclock time (seconds)
**
** Output arguments:
**   wall: wallclock
**   usr:  user time
**   sys:  system time
**
** Return value: 0 (success) or -1 (failure)
*/

static int t_stamp (double *wall, double *usr, double *sys)
{
  struct timeval tp;         /* argument to gettimeofday */
  struct tms buf;            /* argument to times */

  *usr = 0;
  *sys = 0;

  if (times (&buf) == -1)
    return t_error ("t_stamp: times() failed. Timing bogus\n");

  *usr = buf.tms_utime / (double) ticks_per_sec;
  *sys = buf.tms_stime / (double) ticks_per_sec;

  gettimeofday (&tp, NULL);
  *wall = tp.tv_sec + 1.e-6*tp.tv_usec;

  return 0;
}


/*
** t_start.c: start a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

static int t_start (char *name)        /* timer name */
{
  struct timeval tp1, tp2;      /* argument to gettimeofday */
  struct node *ptr;             /* linked list pointer */

  int numchars;                 /* number of characters in timer */
  int mythread;                 /* thread index (of this thread) */
  int indent_level = 0;         /* required indentation level for this timer */
  int ret;                      /* return code */

  PCL_CNT_TYPE i_pcl_result1[PCL_COUNTER_MAX];     /* init. output fm PCLread */
  PCL_CNT_TYPE i_pcl_result2[PCL_COUNTER_MAX];     /* final output fm PCLread */
  PCL_FP_CNT_TYPE fp_pcl_result[PCL_COUNTER_MAX];  /* required by PCLread */

  /*
  ** 1st system timer call is solely for overhead timing
  */

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if (wallenabled)
    gettimeofday (&tp1, NULL);

  if ( ! t_initialized)
    return t_error ("t_start: t_initialize has not been called\n");

  if ((mythread = get_thread_num ()) < 0)
    return t_error ("t_start\n");

  if (ncounter > 0) {
    ret = PCLread (descr[mythread], i_pcl_result1, fp_pcl_result, ncounter);
    if (ret != PCL_SUCCESS)
      return t_error ("t_start: error from PCLread: %s\n", t_pclstr (ret));
  }

  /*
  ** Look for the requested timer in the current list.  For those which don't
  ** match but are currently active, increase the indentation level by 1
  */

  for (ptr = timers[mythread]; ptr != NULL && ! STRMATCH (name, ptr->name); 
       ptr = ptr->next) {

    if (ptr->onflg) 
      indent_level++;
  }

  if (indent_level > max_indent_level[mythread])
    max_indent_level[mythread] = indent_level;
    
  /* 
  ** If a new thing is being timed, add a new link and initialize 
  */

  if (ptr == NULL) {

    if ((ptr = (struct node *) malloc (sizeof (struct node))) == NULL)
      return (t_error ("t_start: malloc failed\n"));

    memset (ptr, 0, sizeof (struct node));
    ptr->indent_level = indent_level;
    ptr->next = NULL;

    if (timers[mythread] == NULL)
      timers[mythread] = ptr;
    else
      last[mythread]->next = ptr;

    last[mythread] = ptr;

    /* 
    ** Truncate input name if longer than ESMC_PROFILER_MAX_CHARS characters 
    */

    numchars = MIN (strlen (name), ESMC_PROFILER_MAX_CHARS);
    strncpy (ptr->name, name, numchars);
    ptr->name[numchars] = '\0';

  } else {

    /*
    ** If computed indentation level is different than before or was
    ** already ambiguous, reset to ambiguous flag value.  This will likely
    ** happen any time the thing being timed is called from more than 1
    ** branch in the call tree.
    */

    if (ptr->indent_level != indent_level) 
      ptr->indent_level = AMBIGUOUS;

    if (ptr->onflg)
      return t_error ("t_start thread %d: timer %s was already on: "
		      "not restarting.\n", mythread, ptr->name);
  }

  ptr->onflg = true;

  if (usrsysenabled)
    if (get_cpustamp (&ptr->last_utime, &ptr->last_stime) < 0)
      return t_error ("t_start: get_cpustamp error");

  /*
  ** The 2nd system timer call is used both for overhead estimation and
  ** the input timer
  */

  if (wallenabled) {

    gettimeofday (&tp2, NULL);
    ptr->last_wtime_sec  = tp2.tv_sec;
    ptr->last_wtime_usec = tp2.tv_usec;
    overhead[mythread] +=       (tp2.tv_sec  - tp1.tv_sec) + 
                          1.e-6*(tp2.tv_usec - tp1.tv_usec);
  }

  if (ncounter > 0) {
    int n;
    int index;

    ret = PCLread (descr[mythread], i_pcl_result2, fp_pcl_result, ncounter); 
    if (ret != PCL_SUCCESS)
      return t_error ("t_start: error from PCLread: %s\n", t_pclstr (ret));

    for (n = 0; n < ncounter; n++) {
      ptr->last_pcl_result[n] = i_pcl_result2[n];
    }

    if (pcl_cyclesenabled) {
      index = pcl_cyclesindex;
      overhead_pcl[mythread] += i_pcl_result2[index] - i_pcl_result1[index];
    }
  } 

  return (0);
}

/*
** This stub should never actually be called
*/

#ifndef HAVE_PCL
static int PCLread (PCL_DESCR_TYPE descr, PCL_CNT_TYPE *i, PCL_CNT_TYPE *j, int k)
{
  return t_error ("PCLread called when library not there\n");
}
#endif


/*
** t_stop: stop a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

static int t_stop (char *name)
{
  long delta_wtime_sec;     /* wallclock change fm t_start() to t_stop() */    
  long delta_wtime_usec;    /* wallclock change fm t_start() to t_stop() */
  float delta_wtime;        /* floating point wallclock change */
  struct timeval tp1, tp2;  /* argument to gettimeofday() */
  struct node *ptr;         /* linked list pointer */

  int mythread;             /* thread number for this process */
  int ret;                  /* return code */

  long usr;
  long sys;

  PCL_CNT_TYPE i_pcl_result1[PCL_COUNTER_MAX];     /* init. output fm PCLread */
  PCL_CNT_TYPE i_pcl_result2[PCL_COUNTER_MAX];     /* final output fm PCLread */
  PCL_FP_CNT_TYPE fp_pcl_result[PCL_COUNTER_MAX];  /* required by PCLread */

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if ( ! t_initialized)
    return t_error ("t_stop: t_initialize has not been called\n");

  /*
  ** The 1st system timer call is used both for overhead estimation and
  ** the input timer
  */

  if (wallenabled)
    gettimeofday (&tp1, NULL);

  if (usrsysenabled && get_cpustamp (&usr, &sys) < 0)
    return t_error (NULL);

  if ((mythread = get_thread_num ()) < 0)
    return t_error ("t_stop\n");

  if (ncounter > 0) {
    ret = PCLread (descr[mythread], i_pcl_result1, fp_pcl_result, ncounter);
    if (ret != PCL_SUCCESS)
      return t_error ("t_stop: error from PCLread: %s\n", t_pclstr (ret));
  }
  
  for (ptr = timers[mythread]; ptr != NULL && ! STRMATCH (name, ptr->name); 
       ptr = ptr->next);

  if (ptr == NULL) 
    return t_error ("t_stop: timer for %s had not been started.\n", name);

  if ( ! ptr->onflg )
    return t_error ("t_stop: timer %s was already off.\n",ptr->name);

  ptr->onflg = false;
  ptr->count++;

  /*
  ** 1st timer stoppage: set max and min to computed values.  Otherwise apply
  ** max or min function
  */

  if (wallenabled) {

    delta_wtime_sec  = tp1.tv_sec  - ptr->last_wtime_sec;
    delta_wtime_usec = tp1.tv_usec - ptr->last_wtime_usec;
    delta_wtime = delta_wtime_sec + 1.e-6*delta_wtime_usec;

    if (ptr->count == 1) {
      ptr->max_wtime = delta_wtime;
      ptr->min_wtime = delta_wtime;
      
    } else {
      
      ptr->max_wtime = MAX (ptr->max_wtime, delta_wtime);
      ptr->min_wtime = MIN (ptr->min_wtime, delta_wtime);
    }

    ptr->accum_wtime_sec  += delta_wtime_sec;
    ptr->accum_wtime_usec += delta_wtime_usec;

    /*
    ** Adjust accumulated wallclock values to guard against overflow in the
    ** microsecond accumulator.
    */

    if (ptr->accum_wtime_usec > 1000000) {
      ptr->accum_wtime_usec -= 1000000;
      ptr->accum_wtime_sec  += 1;
      
    } else if (ptr->accum_wtime_usec < -1000000) {
      
      ptr->accum_wtime_usec += 1000000;
      ptr->accum_wtime_sec  -= 1;
    }

    ptr->last_wtime_sec  = tp1.tv_sec;
    ptr->last_wtime_usec = tp1.tv_usec;

    /*
    ** 2nd system timer call is solely for overhead timing
    */

    gettimeofday (&tp2, NULL);
    overhead[mythread] +=       (tp2.tv_sec  - tp1.tv_sec) + 
                          1.e-6*(tp2.tv_usec - tp1.tv_usec);
  }

  if (usrsysenabled) {
    ptr->accum_utime += usr - ptr->last_utime;
    ptr->accum_stime += sys - ptr->last_stime;
    ptr->last_utime   = usr;
    ptr->last_stime   = sys;
  }

  if (ncounter > 0) {
    int n;
    PCL_CNT_TYPE delta;
    int index;

    for (n = 0; n < ncounter; n++) {
      delta = i_pcl_result1[n] - ptr->last_pcl_result[n];

      /*
      ** Accumulate results only for positive delta
      */

      if (delta < 0) 
	printf ("t_stop: negative delta => probable counter overflow. "
		"Skipping accumulation this round\n"
		"%ld - %ld = %ld\n", (long) i_pcl_result1[n], 
		                     (long) ptr->last_pcl_result[n],
		                     (long) delta);
      else
	ptr->accum_pcl_result[n] += delta;

      ptr->last_pcl_result[n] = i_pcl_result1[n];
    }

    /*
    ** Overhead estimate.  Currently no check for negative delta
    */

    ret = PCLread (descr[mythread], i_pcl_result2, fp_pcl_result, ncounter);
    if (ret != PCL_SUCCESS)
      return t_error ("t_stop: error from PCLread: %s\n", t_pclstr (ret));

    if (pcl_cyclesenabled) {
      index = pcl_cyclesindex;
      overhead_pcl[mythread] += i_pcl_result2[index] - i_pcl_result1[index];
    }
  }

  return 0;
}


/************************* New Interfaces ******************************/
/* Error handling is done by intercepting the t_error function.  This 
   allows for the minimal change in the existing timer library.      */

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerInit"
int ESMC_TimerInit(char *name, ESMC_TimerOption option1, ...)
{
  /* This function uses var args, so it will have to deal with
     the string length as it picks off args. */ 
  va_list args;
  int ret;

  va_start(args, option1);

  ret = ESMC_TimerInitCV(name, option1, args);
  
  va_end(args);

  return ret;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerInitCV"
int ESMC_TimerInitCV(char *name, ESMC_TimerOption option1, va_list args)
{
  int option;
  int iopt;
  static int firstime = 1;
  int ret;

  if (firstime)
    {
      firstime = 0;
    }
  else
    {
      printf("Multiple timers not yet implemented\n");
      return;
    }
  
  /* Repeatedly call t_set_option */
  option = option1;
  while (option != 0)
    {
      switch(option)
	{
	default:
	  iopt = va_arg(args, int);

	  t_setoption((OptionName) option, (Boolean) iopt);
	  break;
	}
      
      /* Get next option */
      iopt = va_arg(args, int);
      if (iopt)
	option = iopt;
      else /* This says that no option can be zero */
	option = 0;
    }

  /* Call the timer libary init */
  t_initialize();

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerTimerInitV"
int ESMC_TimerInitV(char *name, ESMC_TimerOption option1, va_list args)
{
  ESMC_TimerOption option;
  ESMC_TimerOption *iptr;
  static int firstime = 1;
  int ret;

  if (firstime)
    {
      firstime = 0;
    }
  else
    {
      printf("Multiple timers not yet implemented\n");
      return;
    }
  
  /* Repeatedly call t_set_option */
  option = option1;
  while (option != 0)
    {
      switch(option)
	{
	default:
	  iptr = va_arg(args, ESMC_TimerOption*);

	  t_setoption((OptionName) option, (Boolean) *iptr);
	  break;
	}
      
      /* Get next option */
      iptr = va_arg(args, ESMC_TimerOption*);
      if (iptr)
	option = *iptr;
      else /* This says that no option can be zero */
	option = 0;
    }

  /* Call the timer libary init */
  t_initialize();

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerStart"
int ESMC_TimerStart(char *name)
{

  printf("nameBuf:<%s>\n", name);
  t_start(name);
  
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerStop"
int ESMC_TimerStop(char *name)
{
  printf("nameBuf:<%s>\n", name);
  t_stop(name);
  
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerStamp"
int ESMC_TimerStamp(double *wall, double *user, double *sys)
{
  t_stamp(wall, user, sys);
  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerPrint"
int ESMC_TimerPrint(char *name, ESMC_Log log)
{
  int flag;
  int node, process, thread;
  float tfloat;
 
  ESMC_MachinePInfo(&node, &process, &thread);
  
  if (log)
     t_pr(log, name, process);
  else
     t_pr(STDLog, name, process);

  return ESMC_SUCCESS;
}

/*--------------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "ESMC_TimerSetSTDLog"
int ESMC_TimerSetSTDLog(ESMC_Log log)
{
  STDLog = log;

  return ESMC_SUCCESS;
}
