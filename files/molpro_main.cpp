#include "molpro_config.h"
#include "molpro_sha1.h"
#include "molpro.h"
#include "molpro_registry.h"
#include <mpi.h>
#include "ppidd.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <chrono>
#include <ctime>
#ifdef __MINGW32__
#include <winsock2.h>
#endif

#ifdef HAVE_GPERFTOOLS_PROFILER_H
#include <gperftools/profiler.h>
#endif


#define PROGRAM_MOLPRO 0
#define PROGRAM_REGISTRY 1
#define PROGRAM_VERSION 2

extern "C" void declare_version(const char* version);

int main(int argc, char* argv[]) {
 bool debug=false;
 std::string sha1=MOLPRO_SHA1;

 /* With these options the user can modify which implementation of PPIDD
    is used. There is no check as to whether their choice is sensible,
    or even whether it will work at all */
 int ppidd_impl=PPIDD_IMPL_DEFAULT;
 bool enable_gperftools = false;
 for (int i=1; i < argc; ++i) {
  if (strcmp(argv[i],"--ppidd-ga-mpi") == 0) ppidd_impl=PPIDD_IMPL_GA_MPI;
  if (strcmp(argv[i],"--ppidd-mpi2") == 0) ppidd_impl=PPIDD_IMPL_MPI2;
  if (strcmp(argv[i],"--ppidd-no-mpi") == 0) ppidd_impl=PPIDD_IMPL_NO_MPI;
  if (strcmp(argv[i],"--enable-gperftools") == 0) enable_gperftools = true;
 }

 PPIDD_Initialize(&argc,&argv,ppidd_impl);
 int np = PPIDD_Size();
 int me = PPIDD_Rank();
 MPI_Comm mpicomm=MPI_Comm_f2c(PPIDD_Worker_comm());
 if (!getenv("LANG")) putenv(strdup("LANG=C")); /* bug2832 */
 debug = (debug && me == 0) ? true : false;

 declare_version(PACKAGE_VERSION);

 if (debug) std::cout << "argv[0] = " << argv[0] << '\n';

 const char path_sep =
#ifdef _WIN32
                       '\\';
#else
                       '/';
#endif
 std::string libmol="";
 std::string vdate="";
 if (me == 0) { /* maybe only master has accurate argv, broadcast libmol later */

#ifdef __MINGW32__
  libmol.resize(1024);
  GetModuleFileNameA(NULL,&libmol[0],libmol.length());
  libmol.resize(libmol.find('\0'));
#else
  char* real_path=realpath(argv[0],nullptr);
  if (real_path == nullptr) {
   std::cout << "[main] Error: cannot determine libmol\n" ;
   MPI_Abort(MPI_COMM_WORLD,1);
  }
  else {
   libmol=real_path;
   free(real_path);
  }
#endif

#ifdef __cpp_lib_filesystem
  #include <filesystem>
  auto executable = std::filesystem::path{libmol};
  auto last_write_time = std::filesystem::last_write_time(executable);
  auto cftime = std::chrono::system_clock::to_time_t(std::chrono::time_point_cast<std::chrono::system_clock::duration>(
      last_write_time - decltype(last_write_time)::clock::now() + std::chrono::system_clock::now()));
  vdate=std::asctime(std::localtime(&cftime));
#elif __cpp_lib_experimental_filesystem
  #include <experimental/filesystem>
  auto executable = std::experimental/filesystem::path{libmol};
  auto last_write_time = std::experimental/filesystem::last_write_time(executable);
  auto cftime = std::chrono::system_clock::to_time_t(std::chrono::time_point_cast<std::chrono::system_clock::duration>(
      last_write_time - decltype(last_write_time)::clock::now() + std::chrono::system_clock::now()));
  vdate=std::asctime(std::localtime(&cftime));
#endif

  const char* molpro_prefix=getenv("MOLPRO_PREFIX");
  if (molpro_prefix == nullptr) {
   libmol=libmol.substr(0,libmol.find_last_of(path_sep));
   libmol=libmol.substr(0,libmol.find_last_of(path_sep));
  } else libmol=molpro_prefix;
  libmol+=path_sep;
  libmol+="lib";
  libmol+=path_sep;
  if (debug) std::cout << "libmol=" << libmol << '\n';
 }

 int program=PROGRAM_MOLPRO;
 if (me == 0) {
  if (argc > 1 && strcmp(argv[1],"--registry") == 0) program=PROGRAM_REGISTRY;
  for (int i=1; i < argc; ++i) if (strcmp(argv[i],"--version") == 0) program=PROGRAM_VERSION;
 }
 MPI_Bcast(&program,1,MPI_INT,0,mpicomm);

 switch (program) {

  case PROGRAM_REGISTRY:
   if (me == 0) registry_main(argc-1,argv+1,libmol);
   MPI_Barrier(mpicomm);
   break;

  case PROGRAM_VERSION:
   if (me == 0) std::cout << PACKAGE_VERSION << "\n" << sha1 << "\n";
   break;

  case PROGRAM_MOLPRO:
   std::string argstr = "";
   for (int i=1; i < argc; ++i) {
    if (debug) std::cout << "argv[" << i << "] = " << argv[i] << '\n';
    /* we've already parsed the PPIDD related options */
    if (strcmp(argv[i],"--ppidd-ga-mpi") == 0) continue;
    if (strcmp(argv[i],"--ppidd-mpi2") == 0) continue;
    if (strcmp(argv[i],"--ppidd-no-mpi") == 0) continue;
    if (strcmp(argv[i],"--enable-gperftools") == 0) continue;
    if (argstr.length() != 0 ) argstr += " ";
    argstr += argv[i];
   }
   if (debug) std::cout << "argstr  = " << argstr << '\n';

#ifdef HAVE_GPERFTOOLS_PROFILER_H
   if (enable_gperftools) {
    std::string gprof = "molpro.gprof." + std::to_string(me);
    ProfilerStart(gprof.c_str());
   }
#endif

   molpro(libmol.c_str(),libmol.length(),argstr.c_str(),argstr.length(),np,me,ppidd_impl,sha1.c_str(),vdate.c_str());

#ifdef HAVE_GPERFTOOLS_PROFILER_H
   if (enable_gperftools) {
    ProfilerStop();
   }
#endif

 }

 PPIDD_Finalize();

 return 0;
}
