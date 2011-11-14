/@[a-zA-Z_][a-zA-Z_0-9]*@/!b end
s,@SHELL@,|#_!!_#|/bin/sh,g
s,@PATH_SEPARATOR@,|#_!!_#|:,g
s,@PACKAGE_NAME@,|#_!!_#|,g
s,@PACKAGE_TARNAME@,|#_!!_#|,g
s,@PACKAGE_VERSION@,|#_!!_#|,g
s,@PACKAGE_STRING@,|#_!!_#|,g
s,@PACKAGE_BUGREPORT@,|#_!!_#|,g
s,@exec_prefix@,|#_!!_#|${prefix},g
s,@prefix@,|#_!!_#|/usr/local,g
s,@program_transform_name@,|#_!!_#|s\,x\,x\,,g
s,@bindir@,|#_!!_#|${exec_prefix}/bin,g
s,@sbindir@,|#_!!_#|${exec_prefix}/sbin,g
s,@libexecdir@,|#_!!_#|${exec_prefix}/libexec,g
s,@datarootdir@,|#_!!_#|${prefix}/share,g
s,@datadir@,|#_!!_#|${datarootdir},g
s,@sysconfdir@,|#_!!_#|${prefix}/etc,g
s,@sharedstatedir@,|#_!!_#|${prefix}/com,g
s,@localstatedir@,|#_!!_#|${prefix}/var,g
s,@includedir@,|#_!!_#|${prefix}/include,g
s,@oldincludedir@,|#_!!_#|/usr/include,g
s,@docdir@,|#_!!_#|${datarootdir}/doc/${PACKAGE},g
s,@infodir@,|#_!!_#|${datarootdir}/info,g
s,@htmldir@,|#_!!_#|${docdir},g
s,@dvidir@,|#_!!_#|${docdir},g
s,@pdfdir@,|#_!!_#|${docdir},g
s,@psdir@,|#_!!_#|${docdir},g
s,@libdir@,|#_!!_#|${exec_prefix}/lib,g
s,@localedir@,|#_!!_#|${datarootdir}/locale,g
s,@mandir@,|#_!!_#|${datarootdir}/man,g
s,@DEFS@,|#_!!_#|-DPACKAGE_NAME=\\"\\" -DPACKAGE_TARNAME=\\"\\" -DPACKAGE_VERSION=\\"\\" -DPACKAGE_STRING=\\"\\" -DPACKAGE_BUGREPORT=\\"\\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_TRUNC=1,g
s,@ECHO_C@,|#_!!_#|\\c,g
s,@ECHO_N@,|#_!!_#|,g
s,@ECHO_T@,|#_!!_#|,g
s,@LIBS@,|#_!!_#|,g
s,@build_alias@,|#_!!_#|,g
s,@host_alias@,|#_!!_#|,g
s,@target_alias@,|#_!!_#|,g
s,@CXX@,|#_!!_#|g++ -arch x86_64,g
s,@CXXFLAGS@,|#_!!_#|-g -O2,g
s,@LDFLAGS@,|#_!!_#|,g
s,@CPPFLAGS@,|#_!!_#|,g
s,@ac_ct_CXX@,|#_!!_#|,g
s,@EXEEXT@,|#_!!_#|,g
s,@OBJEXT@,|#_!!_#|o,g
s,@CC@,|#_!!_#|gcc,g
s,@CFLAGS@,|#_!!_#|-g -O2,g
s,@ac_ct_CC@,|#_!!_#|gcc,g
s,@CPP@,|#_!!_#|gcc -E,g
s,@GREP@,|#_!!_#|/usr/bin/grep,g
s,@EGREP@,|#_!!_#|/usr/bin/grep -E,g
s,@MV_HAVE_IEEEFP_H@,|#_!!_#|,g
s,@MV_HAVE_TRUNC@,|#_!!_#|-DHAVE_TRUNC,g
s,@LIBOBJS@,|#_!!_#|,g
s,@LTLIBOBJS@,|#_!!_#|,g
:end
s/|#_!!_#|//g
