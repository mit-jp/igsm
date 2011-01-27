# $Id$ 
ALL: lib 

CFLAGS    =
FFLAGS    =
SOURCEC	  = 
SOURCEF	  = 
SOURCEH	  = 
OBJSC     = 
OBJSF     = 
TEXFILES  = 
MANSEC	  = 
LIBBASE	  = 
DIRS	  = 
LOCDIR	  = src/Infrastructure/Timer/doc

DVIFILES  = \
		${ESMC_DOCDIR}/Timer_desdoc.dvi \
		${ESMC_DOCDIR}/Timer_reqdoc.dvi \
		${ESMC_DOCDIR}/Timer_refdoc.dvi

PDFFILES  = $(addsuffix .pdf, $(basename $(DVIFILES)))
HTMLFILES = \
		${ESMC_DOCDIR}/Timer_desdoc \
		${ESMC_DOCDIR}/Timer_reqdoc \
		${ESMC_DOCDIR}/Timer_refdoc


include ${ESMF_DIR}/build/${ESMF_ARCH}/base

# Dvi files
${ESMC_DOCDIR}/Timer_desdoc.dvi:
	${DO_LATEX} Timer des
	mv -f Timer_desdoc.dvi ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Timer_reqdoc.dvi:
	${DO_LATEX} Timer req
	mv -f Timer_reqdoc.dvi ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Timer_refdoc.dvi:
	${DO_LATEX} Timer ref
	mv -f Timer_refdoc.dvi ${ESMC_DOCDIR}


# Pdf files
${ESMC_DOCDIR}/Timer_desdoc.pdf: ${ESMC_DOCDIR}/Timer_desdoc.dvi
	dvipdf ${ESMC_DOCDIR}/Timer_desdoc.dvi ${ESMC_DOCDIR}/Timer_desdoc.pdf

${ESMC_DOCDIR}/Timer_reqdoc.pdf: ${ESMC_DOCDIR}/Timer_reqdoc.dvi
	dvipdf ${ESMC_DOCDIR}/Timer_reqdoc.dvi ${ESMC_DOCDIR}/Timer_reqdoc.pdf

${ESMC_DOCDIR}/Timer_refdoc.pdf: ${ESMC_DOCDIR}/Timer_refdoc.dvi
	dvipdf ${ESMC_DOCDIR}/Timer_refdoc.dvi ${ESMC_DOCDIR}/Timer_refdoc.pdf

# Html Files
${ESMC_DOCDIR}/Timer_desdoc: 
	${DOL2H} Timer des
	mv -f Timer_desdoc ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Timer_reqdoc: 
	${DOL2H} Timer req
	mv -f Timer_reqdoc ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Timer_refdoc: 
	${DOL2H} Timer ref
	mv -f Timer_refdoc ${ESMC_DOCDIR}
