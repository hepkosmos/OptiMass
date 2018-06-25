#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.63])
AC_INIT([alm_base], [1.0], [hepkosmos@ibs.re.kr])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config/config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_PREFIX_DEFAULT([${PWD%/*}])

AM_INIT_AUTOMAKE([-Wall foreign])
AC_PROG_CXX
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX(11, noext, mandatory)
AC_PROG_INSTALL
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

AC_CHECK_FUNCS([pow sqrt],,AC_MSG_ERROR([Required math functions not found]))
AC_C_INLINE
AC_TYPE_SIZE_T

# Debug mode
AC_ARG_ENABLE([debug],
[  --enable-debug          turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

# ROOT Installation Check & Get-Substitute Env. Variables
ROOT_PATH(,AC_MSG_NOTICE([ROOT $ROOTVERSION is found.]),AC_MSG_ERROR(ROOT is not found.))

# ROOT MINUIT2 Activation Check
ROOT_FEATURE([minuit2],[AC_MSG_NOTICE([MINUIT2 is enabled.])],[AC_MSG_ERROR(["MINUIT2 in your ROOT is NOT enabled.])])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
