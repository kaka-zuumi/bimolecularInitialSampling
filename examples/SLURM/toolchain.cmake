# Compilers
#set(CMAKE_C_COMPILER   cmake)
set(CMAKE_C_COMPILER   cc)
set(CMAKE_CXX_COMPILER c++)
set(MPI_C_COMPILER    /opt/apps/software/mpi/OpenMPI/4.1.4-GCC-11.3.0/bin/mpicc)
set(MPI_CXX_COMPILER   /opt/apps/software/mpi/OpenMPI/4.1.4-GCC-11.3.0/bin/mpiCC)

# Token for private repos
set(CMAIZE_GITHUB_TOKEN xxxxYourGithubTokenHerexxxx)

# Options
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(BUILD_SHARED_LIBS TRUE)
set(BUILD_TESTING TRUE)


set(CMAKE_PREFIX_PATH
    /opt/apps/software/numlib/OpenBLAS/0.3.20-GCC-11.3.0
    /opt/apps/software/mpi/OpenMPI/4.1.4-GCC-11.3.0
    /mnt/lustre/koa/koastore/rsun_group/ruisun/sourcecode/libint/libint-2.6.0
    /mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/nwchemGCCtest1/nwchem-ivybridge-try1/nwchem/bin/LINUX64
)

set(CMAKE_CXX_STANDARD 17)

# BLAS/LAPACK
set(ENABLE_SCALAPACK ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMPI_SKIP_MPICXX")

# Kazuumi addition
set(NWX_MODULE_DIRECTORY /home/kazuumi/rsun_koastore/kazuumiTest1/nwchemexTEST1/nwchemex_gcc_gccNWC_automated_try7/NWChemEx-modules)


