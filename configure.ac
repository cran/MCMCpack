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
    1.*|2.*)
      AC_MSG_WARN([Only g++ version 3.0 or greater can be used with MCMCpack.])
      AC_MSG_ERROR([Please use a different compiler.])   
   ;;
  esac
fi
