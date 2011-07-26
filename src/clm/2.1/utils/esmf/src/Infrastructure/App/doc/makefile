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
LOCDIR	  = src/Infrastructure/App/doc

DVIFILES  = \
		${ESMC_DOCDIR}/App_desdoc.dvi \
		${ESMC_DOCDIR}/App_refdoc.dvi

#		${ESMC_DOCDIR}/App_reqdoc.dvi 

PDFFILES  = $(addsuffix .pdf, $(basename $(DVIFILES)))

HTMLFILES = \
		${ESMC_DOCDIR}/App_desdoc \
		${ESMC_DOCDIR}/App_refdoc

#		${ESMC_DOCDIR}/App_reqdoc 

include ${ESMF_DIR}/build/${ESMF_ARCH}/base

# Dvi files
${ESMC_DOCDIR}/App_desdoc.dvi:
	${DO_LATEX} App des
	mv -f App_desdoc.dvi ${ESMC_DOCDIR}

${ESMC_DOCDIR}/App_reqdoc.dvi:
	${DO_LATEX} App req
	mv -f App_reqdoc.dvi ${ESMC_DOCDIR}

${ESMC_DOCDIR}/App_refdoc.dvi:
	${DO_LATEX} App ref
	mv -f App_refdoc.dvi ${ESMC_DOCDIR}


# Pdf files
${ESMC_DOCDIR}/App_desdoc.pdf: ${ESMC_DOCDIR}/App_desdoc.dvi
	dvipdf ${ESMC_DOCDIR}/App_desdoc.dvi ${ESMC_DOCDIR}/App_desdoc.pdf

${ESMC_DOCDIR}/App_reqdoc.pdf: ${ESMC_DOCDIR}/App_reqdoc.dvi
	dvipdf ${ESMC_DOCDIR}/App_reqdoc.dvi ${ESMC_DOCDIR}/App_reqdoc.pdf

${ESMC_DOCDIR}/App_refdoc.pdf: ${ESMC_DOCDIR}/App_refdoc.dvi
	dvipdf ${ESMC_DOCDIR}/App_refdoc.dvi ${ESMC_DOCDIR}/App_refdoc.pdf

# Html Files
${ESMC_DOCDIR}/App_desdoc: 
	${DOL2H} App des
	mv -f App_desdoc ${ESMC_DOCDIR}

${ESMC_DOCDIR}/App_reqdoc: 
	${DOL2H} App req
	mv -f App_reqdoc ${ESMC_DOCDIR}

${ESMC_DOCDIR}/App_refdoc: 
	${DOL2H} App ref
	mv -f App_refdoc ${ESMC_DOCDIR}
