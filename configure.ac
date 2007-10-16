AC_PREREQ(2.50)
AC_INIT([DESCRIPTION])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`${R_HOME}/bin/R CMD config CXX`

AC_PROG_CXX

if test "${GXX}" = yes; then
  gxx_version=`${CXX} -v 2>&1 | grep "^.*g.. version" | \
     sed -e 's/^.*g.. version *//'`
  case ${gxx_version} in
    1.*|2.*|3.*)
      AC_MSG_WARN([Only g++ version 4.0 or greater can be used with MCMCpack.])
      AC_MSG_ERROR([Please use a different compiler.])   
   ;;
  esac
fi

AC_CHECK_HEADERS(ieeefp.h, [MV_HAVE_IEEEFP_H="-DHAVE_IEEEFP_H"], [MV_HAVE_IEEFP_H=""])
AC_CHECK_FUNCS(trunc, [MV_HAVE_TRUNC="-DHAVE_TRUNC"], [MV_HAVE_TRUNC=""])
AC_SUBST(MV_HAVE_IEEEFP_H)
AC_SUBST(MV_HAVE_TRUNC)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
