check_include_files("assert.h" HAVE_ASSERT_H)
check_include_files("float.h" HAVE_FLOAT_H)
check_include_files("dlfcn.h" HAVE_DLFCN_H)
check_include_files("inttypes.h" HAVE_INTTYPES_H)
check_include_files("limits.h" HAVE_LIMITS_H)
check_include_files("math.h" HAVE_MATH_H)
check_include_files("memory.h" HAVE_MEMORY_H)
check_include_files("stddef.h" HAVE_STDDEF_H)
check_include_files("stdint.h" HAVE_STDINT_H)
check_include_files("stdlib.h" HAVE_STDLIB_H)
check_include_files("strings.h" HAVE_STRINGS_H)
check_include_files("string.h" HAVE_STRING_H)
check_include_files("unistd.h" HAVE_UNISTD_H)
check_include_files("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_files("sys/stat.h" HAVE_SYS_STAT_H)
check_include_files("sys/times.h" HAVE_SYS_TIMES_H)
check_include_files("sys/time.h" HAVE_SYS_TIME_H)
check_include_files("sys/types.h" HAVE_SYS_TYPES_H)
check_include_files("sys/wait.h" HAVE_SYS_WAIT_H)
check_function_exists("gethostname" HAVE_GETHOSTNAME)
check_function_exists("getrlimit" HAVE_GETRLIMIT)
check_function_exists("getrusage" HAVE_GETRUSAGE)
check_library_exists("m" "pow" "math.h" HAVE_POW)
check_library_exists("m" "powl" "math.h" HAVE_POWL)
check_library_exists("m" "sqrt" "math.h" HAVE_SQRT)
check_function_exists("strchr" HAVE_STRCHR)
check_function_exists("strstr" HAVE_STRSTR)
check_function_exists("sysconf" HAVE_SYSCONF)
check_type_size("int" SIZEOF_INT)
check_type_size("long" SIZEOF_LONG)
check_type_size("long double" SIZEOF_LONG_DOUBLE)
check_type_size("void*" SIZEOF_VOID_P)
set(CMAKE_EXTRA_INCLUDE_FILES "stddef.h")
check_type_size("ptrdiff_t" PTRDIFF_T)
set(CMAKE_EXTRA_INCLUDE_FILES)
test_big_endian(WORDS_BIGENDIAN)
set(CUDD_VERSION "3.0.0")

configure_file(cudd_config.h.in cudd_config.h @ONLY)

add_library(cudd STATIC
        cplusplus/cuddObj.cc
        cplusplus/cuddObj.hh
        cudd/cuddAddAbs.c
        cudd/cuddAddApply.c
        cudd/cuddAddFind.c
        cudd/cuddAddInv.c
        cudd/cuddAddIte.c
        cudd/cuddAddNeg.c
        cudd/cuddAddWalsh.c
        cudd/cuddAndAbs.c
        cudd/cuddAnneal.c
        cudd/cuddApa.c
        cudd/cuddAPI.c
        cudd/cuddApprox.c
        cudd/cuddBddAbs.c
        cudd/cuddBddCorr.c
        cudd/cuddBddIte.c
        cudd/cuddBridge.c
        cudd/cuddCache.c
        cudd/cuddCheck.c
        cudd/cuddClip.c
        cudd/cuddCof.c
        cudd/cuddCompose.c
        cudd/cuddDecomp.c
        cudd/cuddEssent.c
        cudd/cuddExact.c
        cudd/cuddExport.c
        cudd/cuddGenCof.c
        cudd/cuddGenetic.c
        cudd/cuddGroup.c
        cudd/cuddHarwell.c
        cudd/cuddInit.c
        cudd/cuddInteract.c
        cudd/cuddLCache.c
        cudd/cuddLevelQ.c
        cudd/cuddLinear.c
        cudd/cuddLiteral.c
        cudd/cuddMatMult.c
        cudd/cuddPriority.c
        cudd/cuddRead.c
        cudd/cuddRef.c
        cudd/cuddReorder.c
        cudd/cuddSat.c
        cudd/cuddSign.c
        cudd/cuddSolve.c
        cudd/cuddSplit.c
        cudd/cuddSubsetHB.c
        cudd/cuddSubsetSP.c
        cudd/cuddSymmetry.c
        cudd/cuddTable.c
        cudd/cuddUtil.c
        cudd/cuddWindow.c
        cudd/cuddZddCount.c
        cudd/cuddZddFuncs.c
        cudd/cuddZddGroup.c
        cudd/cuddZddIsop.c
        cudd/cuddZddLin.c
        cudd/cuddZddMisc.c
        cudd/cuddZddPort.c
        cudd/cuddZddReord.c
        cudd/cuddZddSetop.c
        cudd/cuddZddSymm.c
        cudd/cuddZddUtil.c
        dddmp/dddmpBinary.c
        dddmp/dddmpConvert.c
        dddmp/dddmpDbg.c
        dddmp/dddmpLoad.c
        dddmp/dddmpLoadCnf.c
        dddmp/dddmpNodeAdd.c
        dddmp/dddmpNodeBdd.c
        dddmp/dddmpNodeCnf.c
        dddmp/dddmpStoreAdd.c
        dddmp/dddmpStoreBdd.c
        dddmp/dddmpStoreCnf.c
        dddmp/dddmpStoreMisc.c
        dddmp/dddmpUtil.c
        epd/epd.c
        mtr/mtrBasic.c
        mtr/mtrGroup.c
        st/st.c
        util/cpu_stats.c
        util/cpu_time.c
        util/cstringstream.c
        util/datalimit.c
        util/pathsearch.c
        util/pipefork.c
        util/prtime.c
        util/safe_mem.c
        util/strsav.c
        util/texpand.c
        util/ucbqsort.c
        cudd/cudd.h
        cudd/cuddInt.h
        dddmp/dddmp.h
        dddmp/dddmpInt.h
        epd/epd.h
        epd/epdInt.h
        mtr/mtr.h
        mtr/mtrInt.h
        # nanotrav/bnet.h
        # nanotrav/ntr.h
        st/st.h
        util/cstringstream.h
        util/util.h
    )
target_include_directories(cudd PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cudd"
                                       "${CMAKE_CURRENT_SOURCE_DIR}/st"
                                       "${CMAKE_CURRENT_SOURCE_DIR}/mtr"
                                       "${CMAKE_CURRENT_SOURCE_DIR}/epd"
                                       "${CMAKE_CURRENT_SOURCE_DIR}/util"
                                       "${CMAKE_CURRENT_SOURCE_DIR}/cplusplus"
                                       "${CMAKE_CURRENT_BINARY_DIR}")
