/* Thin platform-abstraction shim for dlopen/dlsym/dlclose.
   Called from Fortran via iso_c_binding (bind(c, name="rkb_d*")). */

#ifdef _WIN32
#  include <windows.h>
typedef FARPROC rkb_proc_t;
static void* do_open(const char* path)            { return (void*)LoadLibraryA(path); }
static rkb_proc_t do_sym(void* h, const char* s)  { return GetProcAddress((HMODULE)h, s); }
static void do_close(void* h)                     { FreeLibrary((HMODULE)h); }
#else
#  include <dlfcn.h>
typedef void (*rkb_proc_t)(void);
static void* do_open(const char* path)            { return dlopen(path, RTLD_NOW); }
static rkb_proc_t do_sym(void* h, const char* s)  { return (rkb_proc_t)dlsym(h, s); }
static void do_close(void* h)                     { dlclose(h); }
#endif

/* Open a shared library; returns NULL on failure. */
void* rkb_dlopen(const char* path) { return do_open(path); }

/* Look up a symbol; returns the function pointer as a generic function pointer
   so Fortran can receive it as type(c_funptr). */
rkb_proc_t rkb_dlsym(void* handle, const char* sym) { return do_sym(handle, sym); }

/* Close a previously opened library. */
void rkb_dlclose(void* handle) { do_close(handle); }
