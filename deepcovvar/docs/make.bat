@ECHO OFF

REM Makefile for Sphinx documentation
REM

REM You can set these variables from the command line, and also
REM from the environment for the first two.
SET SPHINXOPTS=
SET SPHINXBUILD=sphinx-build
SET SOURCEDIR=source
SET BUILDDIR=build

REM Put it first so that "make" without argument is like "make help".
IF "%1" == "help" (
	%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%
	GOTO :end
)

REM Catch-all target: route all unknown targets to Sphinx using the new
REM "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%
GOTO :end

:end



