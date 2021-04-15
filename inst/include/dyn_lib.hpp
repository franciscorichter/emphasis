#ifndef EMPAHSIS_DYN_LIB_HPP_INCLUDED
#define EMPAHSIS_DYN_LIB_HPP_INCLUDED

#if defined(_WIN32)
#if !defined(WIN32_LEAN_AND_MEAN)
# define WIN32_LEAN_AND_MEAN
#endif
# include <Windows.h>
#else
# include <dlfcn.h>
#endif
#include <string>
#include "emphasis.hpp"


namespace dll {

#if defined(_WIN32)
  class dynlib
  {
  public:
    dynlib(const std::string& DLL)
    {
      hModule_ = LoadLibrary(DLL.c_str());
      if (NULL == hModule_) {
        throw emphasis::emphasis_error("Unable to load dynamic library");
      }
    }

    ~dynlib()
    {
      if (hModule_) {
        FreeLibrary(hModule_);
      }
    }

    template <typename FPTR>
    FPTR get_address(const char* fname, bool optional)
    {
      FPTR fptr = reinterpret_cast<FPTR>(GetProcAddress(hModule_, fname));
      if (!optional && (nullptr == fptr)) throw emphasis::emphasis_error("Can't load function address");
      return fptr;
    }

  private:
    HMODULE hModule_ = NULL;
  };

#else // _WIN32

  class dynlib
  {
  public:
    dynlib(const std::string& DLL)
    {
      hModule_ = dlopen(DLL.c_str(), RTLD_LAZY);
      if (nullptr == hModule_) {
        throw emphasis::emphasis_error("Unable to load dynamic library");
      }
    }

    ~dynlib()
    {
      if (hModule_) {
        dlclose(hModule_);
      }
    }

    template <typename FPTR>
    FPTR get_address(const char* fname, bool optional)
    {
      FPTR fptr = reinterpret_cast<FPTR>(dlsym(hModule_, fname));
      if (!optional && (nullptr == fptr)) throw emphasis::emphasis_error("Can't load function address");
      return fptr;
    }

  private:
    void* hModule_ = nullptr;
  };
#endif

}

#endif
