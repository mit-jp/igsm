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
LOCDIR	  = src/Infrastructure/Log/doc

DVIFILES  = \
		${ESMC_DOCDIR}/Log_desdoc.dvi \
		${ESMC_DOCDIR}/Log_reqdoc.dvi \
		${ESMC_DOCDIR}/Log_refdoc.dvi

PDFFILES  = $(addsuffix .pdf, $(basename $(DVIFILES)))

HTMLFILES = \
		${ESMC_DOCDIR}/Log_desdoc \
		${ESMC_DOCDIR}/Log_reqdoc \
		${ESMC_DOCDIR}/Log_refdoc

include ${ESMF_DIR}/build/${ESMF_ARCH}/base

# Dvi files
${ESMC_DOCDIR}/Log_desdoc.dvi:
	${DO_LATEX} Log des
	mv -f Log_desdoc.dvi ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Log_reqdoc.dvi:
	${DO_LATEX} Log req
	mv -f Log_reqdoc.dvi ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Log_refdoc.dvi:
	${DO_LATEX} Log ref
	mv -f Log_refdoc.dvi ${ESMC_DOCDIR}


# Pdf files
${ESMC_DOCDIR}/Log_desdoc.pdf: ${ESMC_DOCDIR}/Log_desdoc.dvi
	dvipdf ${ESMC_DOCDIR}/Log_desdoc.dvi ${ESMC_DOCDIR}/Log_desdoc.pdf

${ESMC_DOCDIR}/Log_reqdoc.pdf: ${ESMC_DOCDIR}/Log_reqdoc.dvi
	dvipdf ${ESMC_DOCDIR}/Log_reqdoc.dvi ${ESMC_DOCDIR}/Log_reqdoc.pdf

${ESMC_DOCDIR}/Log_refdoc.pdf: ${ESMC_DOCDIR}/Log_refdoc.dvi
	dvipdf ${ESMC_DOCDIR}/Log_refdoc.dvi ${ESMC_DOCDIR}/Log_refdoc.pdf

# Html Files
${ESMC_DOCDIR}/Log_desdoc: 
	${DOL2H} Log des
	mv -f Log_desdoc ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Log_reqdoc: 
	${DOL2H} Log req
	mv -f Log_reqdoc ${ESMC_DOCDIR}

${ESMC_DOCDIR}/Log_refdoc: 
	${DOL2H} Log ref
	mv -f Log_refdoc ${ESMC_DOCDIR}
