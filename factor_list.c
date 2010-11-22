# 1 "pollard.c"
# 1 "C:\\Users\\tokko\\Documents\\My Dropbox\\tokko\\avalg\\project 2\\Integer-factorization\\Integer-factorization//"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "pollard.c"
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 1 3
# 15 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/_mingw.h" 1 3
# 31 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/_mingw.h" 3
       
# 32 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/_mingw.h" 3
# 16 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 2 3





# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 1 3 4
# 211 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 3 4
typedef unsigned int size_t;
# 323 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 3 4
typedef short unsigned int wchar_t;
# 22 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 2 3
# 71 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
extern int _argc;
extern char** _argv;




extern int* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p___argc(void);
extern char*** __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p___argv(void);
extern wchar_t*** __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p___wargv(void);
# 112 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
   extern __attribute__ ((__dllimport__)) int __mb_cur_max;
# 137 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
 int* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _errno(void);


 int* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __doserrno(void);
# 149 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
  extern char *** __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p__environ(void);
  extern wchar_t *** __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p__wenviron(void);
# 172 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
  extern __attribute__ ((__dllimport__)) int _sys_nerr;
# 196 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
extern __attribute__ ((__dllimport__)) char* _sys_errlist[];
# 209 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
extern unsigned __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) int* __p__osver(void);
extern unsigned __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) int* __p__winver(void);
extern unsigned __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) int* __p__winmajor(void);
extern unsigned __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) int* __p__winminor(void);







extern __attribute__ ((__dllimport__)) unsigned int _osver;
extern __attribute__ ((__dllimport__)) unsigned int _winver;
extern __attribute__ ((__dllimport__)) unsigned int _winmajor;
extern __attribute__ ((__dllimport__)) unsigned int _winminor;
# 260 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
 char** __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p__pgmptr(void);

 wchar_t** __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __p__wpgmptr(void);
# 293 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
extern __attribute__ ((__dllimport__)) int _fmode;
# 303 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
 double __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) atof (const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) atoi (const char*);
 long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) atol (const char*);

 double __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wtof (const wchar_t *);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wtoi (const wchar_t *);
 long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wtol (const wchar_t *);


double __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __strtod (const char*, char**);



static

__inline__ double __attribute__((__cdecl__)) __attribute__ ((__nothrow__))
strtod (const char* __restrict__ __nptr, char** __restrict__ __endptr)
{
  return __strtod(__nptr, __endptr);
}
float __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) strtof (const char * __restrict__, char ** __restrict__);
long double __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) strtold (const char * __restrict__, char ** __restrict__);




 long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) strtol (const char*, char**, int);
 unsigned long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) strtoul (const char*, char**, int);



 long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wcstol (const wchar_t*, wchar_t**, int);
 unsigned long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wcstoul (const wchar_t*, wchar_t**, int);
 double __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wcstod (const wchar_t*, wchar_t**);

float __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wcstof( const wchar_t * __restrict__, wchar_t ** __restrict__);
long double __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wcstold (const wchar_t * __restrict__, wchar_t ** __restrict__);


 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wgetenv(const wchar_t*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wputenv(const wchar_t*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wsearchenv(const wchar_t*, const wchar_t*, wchar_t*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wsystem(const wchar_t*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wmakepath(wchar_t*, const wchar_t*, const wchar_t*, const wchar_t*, const wchar_t*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wsplitpath (const wchar_t*, wchar_t*, wchar_t*, wchar_t*, wchar_t*);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wfullpath (wchar_t*, const wchar_t*, size_t);




 size_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wcstombs (char*, const wchar_t*, size_t);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wctomb (char*, wchar_t);

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) mblen (const char*, size_t);
 size_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) mbstowcs (wchar_t*, const char*, size_t);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) mbtowc (wchar_t*, const char*, size_t);

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) rand (void);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) srand (unsigned int);

 void* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) calloc (size_t, size_t) __attribute__ ((__malloc__));
 void* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) malloc (size_t) __attribute__ ((__malloc__));
 void* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) realloc (void*, size_t);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) free (void*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) abort (void) __attribute__ ((__noreturn__));
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) exit (int) __attribute__ ((__noreturn__));


int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) atexit (void (*)(void));

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) system (const char*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) getenv (const char*);


 void* __attribute__((__cdecl__)) bsearch (const void*, const void*, size_t, size_t,
          int (*)(const void*, const void*));
 void __attribute__((__cdecl__)) qsort(void*, size_t, size_t,
      int (*)(const void*, const void*));

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) abs (int) __attribute__ ((__const__));
 long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) labs (long) __attribute__ ((__const__));
# 393 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
typedef struct { int quot, rem; } div_t;
typedef struct { long quot, rem; } ldiv_t;

 div_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) div (int, int) __attribute__ ((__const__));
 ldiv_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ldiv (long, long) __attribute__ ((__const__));







 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _beep (unsigned int, unsigned int) __attribute__ ((__deprecated__));

 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _seterrormode (int) __attribute__ ((__deprecated__));
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _sleep (unsigned long) __attribute__ ((__deprecated__));

 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _exit (int) __attribute__ ((__noreturn__));



typedef int (* _onexit_t)(void);
_onexit_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _onexit( _onexit_t );

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _putenv (const char*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _searchenv (const char*, const char*, char*);

 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ecvt (double, int, int*, int*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fcvt (double, int, int*, int*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _gcvt (double, int, char*);

 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _makepath (char*, const char*, const char*, const char*, const char*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _splitpath (const char*, char*, char*, char*, char*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fullpath (char*, const char*, size_t);

 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _itoa (int, char*, int);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ltoa (long, char*, int);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ultoa(unsigned long, char*, int);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _itow (int, wchar_t*, int);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ltow (long, wchar_t*, int);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ultow (unsigned long, wchar_t*, int);


 long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _atoi64(const char *);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _i64toa(long long, char *, int);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ui64toa(unsigned long long, char *, int);
 long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wtoi64(const wchar_t *);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _i64tow(long long, wchar_t *, int);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _ui64tow(unsigned long long, wchar_t *, int);

 unsigned int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _rotl(unsigned int, int) __attribute__ ((__const__));
 unsigned int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _rotr(unsigned int, int) __attribute__ ((__const__));
 unsigned long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _lrotl(unsigned long, int) __attribute__ ((__const__));
 unsigned long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _lrotr(unsigned long, int) __attribute__ ((__const__));

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _set_error_mode (int);
# 485 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) putenv (const char*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) searchenv (const char*, const char*, char*);

 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) itoa (int, char*, int);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ltoa (long, char*, int);


 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ecvt (double, int, int*, int*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fcvt (double, int, int*, int*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) gcvt (double, int, char*);
# 505 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdlib.h" 3
void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _Exit(int) __attribute__ ((__noreturn__));

extern inline __attribute__((__gnu_inline__)) void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _Exit(int __status)
 { _exit (__status); }


typedef struct { long long quot, rem; } lldiv_t;

lldiv_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) lldiv (long long, long long) __attribute__ ((__const__));

long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) llabs(long long);

extern inline __attribute__((__gnu_inline__)) long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) llabs(long long _j)
  {return (_j >= 0 ? _j : -_j);}


long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) strtoll (const char* __restrict__, char** __restrict, int);
unsigned long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) strtoull (const char* __restrict__, char** __restrict__, int);


long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) atoll (const char *);


long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wtoll (const wchar_t *);
char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) lltoa (long long, char *, int);
char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ulltoa (unsigned long long , char *, int);
wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) lltow (long long, wchar_t *, int);
wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ulltow (unsigned long long, wchar_t *, int);



extern inline __attribute__((__gnu_inline__)) long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) atoll (const char * _c)
 { return _atoi64 (_c); }
extern inline __attribute__((__gnu_inline__)) char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) lltoa (long long _n, char * _c, int _i)
 { return _i64toa (_n, _c, _i); }
extern inline __attribute__((__gnu_inline__)) char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ulltoa (unsigned long long _n, char * _c, int _i)
 { return _ui64toa (_n, _c, _i); }
extern inline __attribute__((__gnu_inline__)) long long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wtoll (const wchar_t * _w)
  { return _wtoi64 (_w); }
extern inline __attribute__((__gnu_inline__)) wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) lltow (long long _n, wchar_t * _w, int _i)
 { return _i64tow (_n, _w, _i); }
extern inline __attribute__((__gnu_inline__)) wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ulltow (unsigned long long _n, wchar_t * _w, int _i)
 { return _ui64tow (_n, _w, _i); }
# 2 "pollard.c" 2
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 1 3
# 26 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 1 3 4
# 352 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 3 4
typedef short unsigned int wint_t;
# 27 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 2 3

# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stdarg.h" 1 3 4
# 40 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stdarg.h" 3 4
typedef __builtin_va_list __gnuc_va_list;
# 29 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 2 3
# 129 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
typedef struct _iobuf
{
 char* _ptr;
 int _cnt;
 char* _base;
 int _flag;
 int _file;
 int _charbuf;
 int _bufsiz;
 char* _tmpfname;
} FILE;
# 154 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
extern __attribute__ ((__dllimport__)) FILE _iob[];
# 169 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fopen (const char*, const char*);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) freopen (const char*, const char*, FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fflush (FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fclose (FILE*);

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) remove (const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) rename (const char*, const char*);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) tmpfile (void);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) tmpnam (char*);


 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _tempnam (const char*, const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _rmtmp(void);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _unlink (const char*);


 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) tempnam (const char*, const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) rmtmp(void);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) unlink (const char*);



 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) setvbuf (FILE*, char*, int, size_t);

 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) setbuf (FILE*, char*);
# 204 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_fprintf(FILE*, const char*, ...);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_printf(const char*, ...);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_sprintf(char*, const char*, ...);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_snprintf(char*, size_t, const char*, ...);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_vfprintf(FILE*, const char*, __gnuc_va_list);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_vprintf(const char*, __gnuc_va_list);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_vsprintf(char*, const char*, __gnuc_va_list);
extern int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __mingw_vsnprintf(char*, size_t, const char*, __gnuc_va_list);
# 293 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fprintf (FILE*, const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) printf (const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) sprintf (char*, const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vfprintf (FILE*, const char*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vprintf (const char*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vsprintf (char*, const char*, __gnuc_va_list);
# 308 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __msvcrt_fprintf(FILE*, const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __msvcrt_printf(const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __msvcrt_sprintf(char*, const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __msvcrt_vfprintf(FILE*, const char*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __msvcrt_vprintf(const char*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) __msvcrt_vsprintf(char*, const char*, __gnuc_va_list);





 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _snprintf (char*, size_t, const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _vsnprintf (char*, size_t, const char*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _vscprintf (const char*, __gnuc_va_list);
# 331 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) snprintf (char *, size_t, const char *, ...);
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vsnprintf (char *, size_t, const char *, __gnuc_va_list);

int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vscanf (const char * __restrict__, __gnuc_va_list);
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vfscanf (FILE * __restrict__, const char * __restrict__,
       __gnuc_va_list);
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vsscanf (const char * __restrict__,
       const char * __restrict__, __gnuc_va_list);







 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fscanf (FILE*, const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) scanf (const char*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) sscanf (const char*, const char*, ...);




 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgetc (FILE*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgets (char*, int, FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fputc (int, FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fputs (const char*, FILE*);
 char* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) gets (char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) puts (const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ungetc (int, FILE*);







 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _filbuf (FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _flsbuf (int, FILE*);



extern inline __attribute__((__gnu_inline__)) int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) getc (FILE* __F)
{
  return (--__F->_cnt >= 0)
    ? (int) (unsigned char) *__F->_ptr++
    : _filbuf (__F);
}

extern inline __attribute__((__gnu_inline__)) int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) putc (int __c, FILE* __F)
{
  return (--__F->_cnt >= 0)
    ? (int) (unsigned char) (*__F->_ptr++ = (char)__c)
    : _flsbuf (__c, __F);
}

extern inline __attribute__((__gnu_inline__)) int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) getchar (void)
{
  return (--(&_iob[0])->_cnt >= 0)
    ? (int) (unsigned char) *(&_iob[0])->_ptr++
    : _filbuf ((&_iob[0]));
}

extern inline __attribute__((__gnu_inline__)) int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) putchar(int __c)
{
  return (--(&_iob[1])->_cnt >= 0)
    ? (int) (unsigned char) (*(&_iob[1])->_ptr++ = (char)__c)
    : _flsbuf (__c, (&_iob[1]));}
# 412 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 size_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fread (void*, size_t, size_t, FILE*);
 size_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fwrite (const void*, size_t, size_t, FILE*);





 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fseek (FILE*, long, int);
 long __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ftell (FILE*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) rewind (FILE*);
# 455 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
typedef long long fpos_t;




 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgetpos (FILE*, fpos_t*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fsetpos (FILE*, const fpos_t*);





 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) feof (FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ferror (FILE*);
# 480 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) clearerr (FILE*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) perror (const char*);






 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _popen (const char*, const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _pclose (FILE*);


 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) popen (const char*, const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) pclose (FILE*);





 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _flushall (void);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fgetchar (void);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fputchar (int);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fdopen (int, const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fileno (FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fcloseall (void);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fsopen (const char*, const char*, int);

 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _getmaxstdio (void);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _setmaxstdio (int);
# 522 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgetchar (void);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fputchar (int);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fdopen (int, const char*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fileno (FILE*);
# 534 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/sys/types.h" 1 3
# 21 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/sys/types.h" 3
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 1 3 4
# 149 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 3 4
typedef int ptrdiff_t;
# 22 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/sys/types.h" 2 3





typedef long __time32_t;




typedef long long __time64_t;
# 45 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/sys/types.h" 3
typedef __time32_t time_t;






typedef long _off_t;


typedef _off_t off_t;







typedef unsigned int _dev_t;





typedef _dev_t dev_t;






typedef short _ino_t;


typedef _ino_t ino_t;






typedef int _pid_t;


typedef _pid_t pid_t;






typedef unsigned short _mode_t;


typedef _mode_t mode_t;






typedef int _sigset_t;


typedef _sigset_t sigset_t;





typedef long _ssize_t;


typedef _ssize_t ssize_t;





typedef long long fpos64_t;




typedef long long off64_t;



typedef unsigned int useconds_t;
# 535 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 2 3
extern inline __attribute__((__gnu_inline__)) FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fopen64 (const char* filename, const char* mode)
{
  return fopen (filename, mode);
}

int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fseeko64 (FILE*, off64_t, int);






extern inline __attribute__((__gnu_inline__)) off64_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ftello64 (FILE * stream)
{
  fpos_t pos;
  if (fgetpos(stream, &pos))
    return -1LL;
  else
   return ((off64_t) pos);
}
# 563 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fwprintf (FILE*, const wchar_t*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wprintf (const wchar_t*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _snwprintf (wchar_t*, size_t, const wchar_t*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vfwprintf (FILE*, const wchar_t*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vwprintf (const wchar_t*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _vsnwprintf (wchar_t*, size_t, const wchar_t*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _vscwprintf (const wchar_t*, __gnuc_va_list);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fwscanf (FILE*, const wchar_t*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wscanf (const wchar_t*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) swscanf (const wchar_t*, const wchar_t*, ...);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgetwc (FILE*);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fputwc (wchar_t, FILE*);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) ungetwc (wchar_t, FILE*);



 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) swprintf (wchar_t*, const wchar_t*, ...);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vswprintf (wchar_t*, const wchar_t*, __gnuc_va_list);



 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgetws (wchar_t*, int, FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fputws (const wchar_t*, FILE*);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) getwc (FILE*);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) getwchar (void);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _getws (wchar_t*);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) putwc (wint_t, FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _putws (const wchar_t*);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) putwchar (wint_t);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wfdopen(int, const wchar_t *);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wfopen (const wchar_t*, const wchar_t*);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wfreopen (const wchar_t*, const wchar_t*, FILE*);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wfsopen (const wchar_t*, const wchar_t*, int);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wtmpnam (wchar_t*);
 wchar_t* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wtempnam (const wchar_t*, const wchar_t*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wrename (const wchar_t*, const wchar_t*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wremove (const wchar_t*);
 void __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wperror (const wchar_t*);
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _wpopen (const wchar_t*, const wchar_t*);



int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) snwprintf (wchar_t* s, size_t n, const wchar_t* format, ...);
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vsnwprintf (wchar_t* s, size_t n, const wchar_t* format, __gnuc_va_list arg);

extern inline __attribute__((__gnu_inline__)) int __attribute__((__cdecl__)) __attribute__ ((__nothrow__))
vsnwprintf (wchar_t* s, size_t n, const wchar_t* format, __gnuc_va_list arg)
  { return _vsnwprintf ( s, n, format, arg);}

int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vwscanf (const wchar_t * __restrict__, __gnuc_va_list);
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vfwscanf (FILE * __restrict__,
         const wchar_t * __restrict__, __gnuc_va_list);
int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) vswscanf (const wchar_t * __restrict__,
         const wchar_t * __restrict__, __gnuc_va_list);
# 625 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/stdio.h" 3
 FILE* __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) wpopen (const wchar_t*, const wchar_t*);






 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fgetwchar (void);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _fputwchar (wint_t);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _getw (FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) _putw (int, FILE*);


 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fgetwchar (void);
 wint_t __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) fputwchar (wint_t);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) getw (FILE*);
 int __attribute__((__cdecl__)) __attribute__ ((__nothrow__)) putw (int, FILE*);
# 3 "pollard.c" 2
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 1 3
# 52 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
# 1 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/include/stddef.h" 1 3 4
# 53 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 2 3
# 193 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
typedef unsigned long int mp_limb_t;
typedef long int mp_limb_signed_t;


typedef unsigned long int mp_bitcnt_t;




typedef struct
{
  int _mp_alloc;

  int _mp_size;


  mp_limb_t *_mp_d;
} __mpz_struct;




typedef __mpz_struct MP_INT;
typedef __mpz_struct mpz_t[1];

typedef mp_limb_t * mp_ptr;
typedef const mp_limb_t * mp_srcptr;







typedef long int mp_size_t;
typedef long int mp_exp_t;


typedef struct
{
  __mpz_struct _mp_num;
  __mpz_struct _mp_den;
} __mpq_struct;

typedef __mpq_struct MP_RAT;
typedef __mpq_struct mpq_t[1];

typedef struct
{
  int _mp_prec;



  int _mp_size;


  mp_exp_t _mp_exp;
  mp_limb_t *_mp_d;
} __mpf_struct;


typedef __mpf_struct mpf_t[1];


typedef enum
{
  GMP_RAND_ALG_DEFAULT = 0,
  GMP_RAND_ALG_LC = GMP_RAND_ALG_DEFAULT
} gmp_randalg_t;


typedef struct
{
  mpz_t _mp_seed;
  gmp_randalg_t _mp_alg;
  union {
    void *_mp_lc;
  } _mp_algdata;
} __gmp_randstate_struct;
typedef __gmp_randstate_struct gmp_randstate_t[1];



typedef const __mpz_struct *mpz_srcptr;
typedef __mpz_struct *mpz_ptr;
typedef const __mpf_struct *mpf_srcptr;
typedef __mpf_struct *mpf_ptr;
typedef const __mpq_struct *mpq_srcptr;
typedef __mpq_struct *mpq_ptr;
# 542 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
 void __gmp_set_memory_functions (void *(*) (size_t), void *(*) (void *, size_t, size_t), void (*) (void *, size_t))

                                                   ;


 void __gmp_get_memory_functions (void *(**) (size_t), void *(**) (void *, size_t, size_t), void (**) (void *, size_t))

                                                                                ;


 extern const int __gmp_bits_per_limb;


 extern int __gmp_errno;


 extern const char * const __gmp_version;






 void __gmp_randinit (gmp_randstate_t, gmp_randalg_t, ...);


 void __gmp_randinit_default (gmp_randstate_t);


 void __gmp_randinit_lc_2exp (gmp_randstate_t, mpz_srcptr, unsigned long int, mp_bitcnt_t)

                          ;


 int __gmp_randinit_lc_2exp_size (gmp_randstate_t, mp_bitcnt_t);


 void __gmp_randinit_mt (gmp_randstate_t);


 void __gmp_randinit_set (gmp_randstate_t, const __gmp_randstate_struct *);


 void __gmp_randseed (gmp_randstate_t, mpz_srcptr);


 void __gmp_randseed_ui (gmp_randstate_t, unsigned long int);


 void __gmp_randclear (gmp_randstate_t);


 unsigned long __gmp_urandomb_ui (gmp_randstate_t, unsigned long);


 unsigned long __gmp_urandomm_ui (gmp_randstate_t, unsigned long);





 int __gmp_asprintf (char **, const char *, ...);



 int __gmp_fprintf (FILE *, const char *, ...);
# 621 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
 int __gmp_printf (const char *, ...);


 int __gmp_snprintf (char *, size_t, const char *, ...);


 int __gmp_sprintf (char *, const char *, ...);
# 659 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
 int __gmp_fscanf (FILE *, const char *, ...);



 int __gmp_scanf (const char *, ...);


 int __gmp_sscanf (const char *, const char *, ...);
# 688 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
 void *__gmpz_realloc (mpz_ptr, mp_size_t);



 void __gmpz_abs (mpz_ptr, mpz_srcptr);



 void __gmpz_add (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_add_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_addmul (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_addmul_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_and (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_array_init (mpz_ptr, mp_size_t, mp_size_t);


 void __gmpz_bin_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_bin_uiui (mpz_ptr, unsigned long int, unsigned long int);


 void __gmpz_cdiv_q (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_cdiv_q_2exp (mpz_ptr, mpz_srcptr, unsigned long);


 unsigned long int __gmpz_cdiv_q_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_cdiv_qr (mpz_ptr, mpz_ptr, mpz_srcptr, mpz_srcptr);


 unsigned long int __gmpz_cdiv_qr_ui (mpz_ptr, mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_cdiv_r (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_cdiv_r_2exp (mpz_ptr, mpz_srcptr, mp_bitcnt_t);


 unsigned long int __gmpz_cdiv_r_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 unsigned long int __gmpz_cdiv_ui (mpz_srcptr, unsigned long int) __attribute__ ((__pure__));


 void __gmpz_clear (mpz_ptr);


 void __gmpz_clears (mpz_ptr, ...);


 void __gmpz_clrbit (mpz_ptr, mp_bitcnt_t);


 int __gmpz_cmp (mpz_srcptr, mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_cmp_d (mpz_srcptr, double) __attribute__ ((__pure__));


 int __gmpz_cmp_si (mpz_srcptr, signed long int) __attribute__ ((__pure__));


 int __gmpz_cmp_ui (mpz_srcptr, unsigned long int) __attribute__ ((__pure__));


 int __gmpz_cmpabs (mpz_srcptr, mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_cmpabs_d (mpz_srcptr, double) __attribute__ ((__pure__));


 int __gmpz_cmpabs_ui (mpz_srcptr, unsigned long int) __attribute__ ((__pure__));


 void __gmpz_com (mpz_ptr, mpz_srcptr);


 void __gmpz_combit (mpz_ptr, mp_bitcnt_t);


 int __gmpz_congruent_p (mpz_srcptr, mpz_srcptr, mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_congruent_2exp_p (mpz_srcptr, mpz_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 int __gmpz_congruent_ui_p (mpz_srcptr, unsigned long, unsigned long) __attribute__ ((__pure__));


 void __gmpz_divexact (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_divexact_ui (mpz_ptr, mpz_srcptr, unsigned long);


 int __gmpz_divisible_p (mpz_srcptr, mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_divisible_ui_p (mpz_srcptr, unsigned long) __attribute__ ((__pure__));


 int __gmpz_divisible_2exp_p (mpz_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 void __gmpz_dump (mpz_srcptr);


 void *__gmpz_export (void *, size_t *, int, size_t, int, size_t, mpz_srcptr);


 void __gmpz_fac_ui (mpz_ptr, unsigned long int);


 void __gmpz_fdiv_q (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_fdiv_q_2exp (mpz_ptr, mpz_srcptr, mp_bitcnt_t);


 unsigned long int __gmpz_fdiv_q_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_fdiv_qr (mpz_ptr, mpz_ptr, mpz_srcptr, mpz_srcptr);


 unsigned long int __gmpz_fdiv_qr_ui (mpz_ptr, mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_fdiv_r (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_fdiv_r_2exp (mpz_ptr, mpz_srcptr, mp_bitcnt_t);


 unsigned long int __gmpz_fdiv_r_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 unsigned long int __gmpz_fdiv_ui (mpz_srcptr, unsigned long int) __attribute__ ((__pure__));


 void __gmpz_fib_ui (mpz_ptr, unsigned long int);


 void __gmpz_fib2_ui (mpz_ptr, mpz_ptr, unsigned long int);


 int __gmpz_fits_sint_p (mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_fits_slong_p (mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_fits_sshort_p (mpz_srcptr) __attribute__ ((__pure__));



 int __gmpz_fits_uint_p (mpz_srcptr) __attribute__ ((__pure__));




 int __gmpz_fits_ulong_p (mpz_srcptr) __attribute__ ((__pure__));




 int __gmpz_fits_ushort_p (mpz_srcptr) __attribute__ ((__pure__));



 void __gmpz_gcd (mpz_ptr, mpz_srcptr, mpz_srcptr);


 unsigned long int __gmpz_gcd_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_gcdext (mpz_ptr, mpz_ptr, mpz_ptr, mpz_srcptr, mpz_srcptr);


 double __gmpz_get_d (mpz_srcptr) __attribute__ ((__pure__));


 double __gmpz_get_d_2exp (signed long int *, mpz_srcptr);


 long int __gmpz_get_si (mpz_srcptr) __attribute__ ((__pure__));


 char *__gmpz_get_str (char *, int, mpz_srcptr);



 unsigned long int __gmpz_get_ui (mpz_srcptr) __attribute__ ((__pure__));




 mp_limb_t __gmpz_getlimbn (mpz_srcptr, mp_size_t) __attribute__ ((__pure__));



 mp_bitcnt_t __gmpz_hamdist (mpz_srcptr, mpz_srcptr) __attribute__ ((__pure__));


 void __gmpz_import (mpz_ptr, size_t, int, size_t, int, size_t, const void *);


 void __gmpz_init (mpz_ptr);


 void __gmpz_init2 (mpz_ptr, mp_bitcnt_t);


 void __gmpz_inits (mpz_ptr, ...);


 void __gmpz_init_set (mpz_ptr, mpz_srcptr);


 void __gmpz_init_set_d (mpz_ptr, double);


 void __gmpz_init_set_si (mpz_ptr, signed long int);


 int __gmpz_init_set_str (mpz_ptr, const char *, int);


 void __gmpz_init_set_ui (mpz_ptr, unsigned long int);



 size_t __gmpz_inp_raw (mpz_ptr, FILE *);




 size_t __gmpz_inp_str (mpz_ptr, FILE *, int);



 int __gmpz_invert (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_ior (mpz_ptr, mpz_srcptr, mpz_srcptr);


 int __gmpz_jacobi (mpz_srcptr, mpz_srcptr) __attribute__ ((__pure__));




 int __gmpz_kronecker_si (mpz_srcptr, long) __attribute__ ((__pure__));


 int __gmpz_kronecker_ui (mpz_srcptr, unsigned long) __attribute__ ((__pure__));


 int __gmpz_si_kronecker (long, mpz_srcptr) __attribute__ ((__pure__));


 int __gmpz_ui_kronecker (unsigned long, mpz_srcptr) __attribute__ ((__pure__));


 void __gmpz_lcm (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_lcm_ui (mpz_ptr, mpz_srcptr, unsigned long);




 void __gmpz_lucnum_ui (mpz_ptr, unsigned long int);


 void __gmpz_lucnum2_ui (mpz_ptr, mpz_ptr, unsigned long int);


 int __gmpz_millerrabin (mpz_srcptr, int) __attribute__ ((__pure__));


 void __gmpz_mod (mpz_ptr, mpz_srcptr, mpz_srcptr);




 void __gmpz_mul (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_mul_2exp (mpz_ptr, mpz_srcptr, mp_bitcnt_t);


 void __gmpz_mul_si (mpz_ptr, mpz_srcptr, long int);


 void __gmpz_mul_ui (mpz_ptr, mpz_srcptr, unsigned long int);



 void __gmpz_neg (mpz_ptr, mpz_srcptr);



 void __gmpz_nextprime (mpz_ptr, mpz_srcptr);



 size_t __gmpz_out_raw (FILE *, mpz_srcptr);




 size_t __gmpz_out_str (FILE *, int, mpz_srcptr);



 int __gmpz_perfect_power_p (mpz_srcptr) __attribute__ ((__pure__));



 int __gmpz_perfect_square_p (mpz_srcptr) __attribute__ ((__pure__));




 mp_bitcnt_t __gmpz_popcount (mpz_srcptr) __attribute__ ((__pure__));



 void __gmpz_pow_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_powm (mpz_ptr, mpz_srcptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_powm_sec (mpz_ptr, mpz_srcptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_powm_ui (mpz_ptr, mpz_srcptr, unsigned long int, mpz_srcptr);


 int __gmpz_probab_prime_p (mpz_srcptr, int) __attribute__ ((__pure__));


 void __gmpz_random (mpz_ptr, mp_size_t);


 void __gmpz_random2 (mpz_ptr, mp_size_t);


 void __gmpz_realloc2 (mpz_ptr, mp_bitcnt_t);


 unsigned long int __gmpz_remove (mpz_ptr, mpz_srcptr, mpz_srcptr);


 int __gmpz_root (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_rootrem (mpz_ptr,mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_rrandomb (mpz_ptr, gmp_randstate_t, mp_bitcnt_t);


 mp_bitcnt_t __gmpz_scan0 (mpz_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 mp_bitcnt_t __gmpz_scan1 (mpz_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 void __gmpz_set (mpz_ptr, mpz_srcptr);


 void __gmpz_set_d (mpz_ptr, double);


 void __gmpz_set_f (mpz_ptr, mpf_srcptr);



 void __gmpz_set_q (mpz_ptr, mpq_srcptr);



 void __gmpz_set_si (mpz_ptr, signed long int);


 int __gmpz_set_str (mpz_ptr, const char *, int);


 void __gmpz_set_ui (mpz_ptr, unsigned long int);


 void __gmpz_setbit (mpz_ptr, mp_bitcnt_t);



 size_t __gmpz_size (mpz_srcptr) __attribute__ ((__pure__));



 size_t __gmpz_sizeinbase (mpz_srcptr, int) __attribute__ ((__pure__));


 void __gmpz_sqrt (mpz_ptr, mpz_srcptr);


 void __gmpz_sqrtrem (mpz_ptr, mpz_ptr, mpz_srcptr);


 void __gmpz_sub (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_sub_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_ui_sub (mpz_ptr, unsigned long int, mpz_srcptr);


 void __gmpz_submul (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_submul_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_swap (mpz_ptr, mpz_ptr) ;


 unsigned long int __gmpz_tdiv_ui (mpz_srcptr, unsigned long int) __attribute__ ((__pure__));


 void __gmpz_tdiv_q (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_tdiv_q_2exp (mpz_ptr, mpz_srcptr, mp_bitcnt_t);


 unsigned long int __gmpz_tdiv_q_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_tdiv_qr (mpz_ptr, mpz_ptr, mpz_srcptr, mpz_srcptr);


 unsigned long int __gmpz_tdiv_qr_ui (mpz_ptr, mpz_ptr, mpz_srcptr, unsigned long int);


 void __gmpz_tdiv_r (mpz_ptr, mpz_srcptr, mpz_srcptr);


 void __gmpz_tdiv_r_2exp (mpz_ptr, mpz_srcptr, mp_bitcnt_t);


 unsigned long int __gmpz_tdiv_r_ui (mpz_ptr, mpz_srcptr, unsigned long int);


 int __gmpz_tstbit (mpz_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 void __gmpz_ui_pow_ui (mpz_ptr, unsigned long int, unsigned long int);


 void __gmpz_urandomb (mpz_ptr, gmp_randstate_t, mp_bitcnt_t);


 void __gmpz_urandomm (mpz_ptr, gmp_randstate_t, mpz_srcptr);



 void __gmpz_xor (mpz_ptr, mpz_srcptr, mpz_srcptr);






 void __gmpq_abs (mpq_ptr, mpq_srcptr);



 void __gmpq_add (mpq_ptr, mpq_srcptr, mpq_srcptr);


 void __gmpq_canonicalize (mpq_ptr);


 void __gmpq_clear (mpq_ptr);


 void __gmpq_clears (mpq_ptr, ...);


 int __gmpq_cmp (mpq_srcptr, mpq_srcptr) __attribute__ ((__pure__));


 int __gmpq_cmp_si (mpq_srcptr, long, unsigned long) __attribute__ ((__pure__));


 int __gmpq_cmp_ui (mpq_srcptr, unsigned long int, unsigned long int) __attribute__ ((__pure__));


 void __gmpq_div (mpq_ptr, mpq_srcptr, mpq_srcptr);


 void __gmpq_div_2exp (mpq_ptr, mpq_srcptr, mp_bitcnt_t);


 int __gmpq_equal (mpq_srcptr, mpq_srcptr) __attribute__ ((__pure__));


 void __gmpq_get_num (mpz_ptr, mpq_srcptr);


 void __gmpq_get_den (mpz_ptr, mpq_srcptr);


 double __gmpq_get_d (mpq_srcptr) __attribute__ ((__pure__));


 char *__gmpq_get_str (char *, int, mpq_srcptr);


 void __gmpq_init (mpq_ptr);


 void __gmpq_inits (mpq_ptr, ...);



 size_t __gmpq_inp_str (mpq_ptr, FILE *, int);



 void __gmpq_inv (mpq_ptr, mpq_srcptr);


 void __gmpq_mul (mpq_ptr, mpq_srcptr, mpq_srcptr);


 void __gmpq_mul_2exp (mpq_ptr, mpq_srcptr, mp_bitcnt_t);



 void __gmpq_neg (mpq_ptr, mpq_srcptr);




 size_t __gmpq_out_str (FILE *, int, mpq_srcptr);



 void __gmpq_set (mpq_ptr, mpq_srcptr);


 void __gmpq_set_d (mpq_ptr, double);


 void __gmpq_set_den (mpq_ptr, mpz_srcptr);


 void __gmpq_set_f (mpq_ptr, mpf_srcptr);


 void __gmpq_set_num (mpq_ptr, mpz_srcptr);


 void __gmpq_set_si (mpq_ptr, signed long int, unsigned long int);


 int __gmpq_set_str (mpq_ptr, const char *, int);


 void __gmpq_set_ui (mpq_ptr, unsigned long int, unsigned long int);


 void __gmpq_set_z (mpq_ptr, mpz_srcptr);


 void __gmpq_sub (mpq_ptr, mpq_srcptr, mpq_srcptr);


 void __gmpq_swap (mpq_ptr, mpq_ptr) ;





 void __gmpf_abs (mpf_ptr, mpf_srcptr);


 void __gmpf_add (mpf_ptr, mpf_srcptr, mpf_srcptr);


 void __gmpf_add_ui (mpf_ptr, mpf_srcptr, unsigned long int);

 void __gmpf_ceil (mpf_ptr, mpf_srcptr);


 void __gmpf_clear (mpf_ptr);


 void __gmpf_clears (mpf_ptr, ...);


 int __gmpf_cmp (mpf_srcptr, mpf_srcptr) __attribute__ ((__pure__));


 int __gmpf_cmp_d (mpf_srcptr, double) __attribute__ ((__pure__));


 int __gmpf_cmp_si (mpf_srcptr, signed long int) __attribute__ ((__pure__));


 int __gmpf_cmp_ui (mpf_srcptr, unsigned long int) __attribute__ ((__pure__));


 void __gmpf_div (mpf_ptr, mpf_srcptr, mpf_srcptr);


 void __gmpf_div_2exp (mpf_ptr, mpf_srcptr, mp_bitcnt_t);


 void __gmpf_div_ui (mpf_ptr, mpf_srcptr, unsigned long int);


 void __gmpf_dump (mpf_srcptr);


 int __gmpf_eq (mpf_srcptr, mpf_srcptr, unsigned long int) __attribute__ ((__pure__));


 int __gmpf_fits_sint_p (mpf_srcptr) __attribute__ ((__pure__));


 int __gmpf_fits_slong_p (mpf_srcptr) __attribute__ ((__pure__));


 int __gmpf_fits_sshort_p (mpf_srcptr) __attribute__ ((__pure__));


 int __gmpf_fits_uint_p (mpf_srcptr) __attribute__ ((__pure__));


 int __gmpf_fits_ulong_p (mpf_srcptr) __attribute__ ((__pure__));


 int __gmpf_fits_ushort_p (mpf_srcptr) __attribute__ ((__pure__));


 void __gmpf_floor (mpf_ptr, mpf_srcptr);


 double __gmpf_get_d (mpf_srcptr) __attribute__ ((__pure__));


 double __gmpf_get_d_2exp (signed long int *, mpf_srcptr);


 mp_bitcnt_t __gmpf_get_default_prec (void) __attribute__ ((__pure__));


 mp_bitcnt_t __gmpf_get_prec (mpf_srcptr) __attribute__ ((__pure__));


 long __gmpf_get_si (mpf_srcptr) __attribute__ ((__pure__));


 char *__gmpf_get_str (char *, mp_exp_t *, int, size_t, mpf_srcptr);


 unsigned long __gmpf_get_ui (mpf_srcptr) __attribute__ ((__pure__));


 void __gmpf_init (mpf_ptr);


 void __gmpf_init2 (mpf_ptr, mp_bitcnt_t);


 void __gmpf_inits (mpf_ptr, ...);


 void __gmpf_init_set (mpf_ptr, mpf_srcptr);


 void __gmpf_init_set_d (mpf_ptr, double);


 void __gmpf_init_set_si (mpf_ptr, signed long int);


 int __gmpf_init_set_str (mpf_ptr, const char *, int);


 void __gmpf_init_set_ui (mpf_ptr, unsigned long int);



 size_t __gmpf_inp_str (mpf_ptr, FILE *, int);



 int __gmpf_integer_p (mpf_srcptr) __attribute__ ((__pure__));


 void __gmpf_mul (mpf_ptr, mpf_srcptr, mpf_srcptr);


 void __gmpf_mul_2exp (mpf_ptr, mpf_srcptr, mp_bitcnt_t);


 void __gmpf_mul_ui (mpf_ptr, mpf_srcptr, unsigned long int);


 void __gmpf_neg (mpf_ptr, mpf_srcptr);



 size_t __gmpf_out_str (FILE *, int, size_t, mpf_srcptr);



 void __gmpf_pow_ui (mpf_ptr, mpf_srcptr, unsigned long int);


 void __gmpf_random2 (mpf_ptr, mp_size_t, mp_exp_t);


 void __gmpf_reldiff (mpf_ptr, mpf_srcptr, mpf_srcptr);


 void __gmpf_set (mpf_ptr, mpf_srcptr);


 void __gmpf_set_d (mpf_ptr, double);


 void __gmpf_set_default_prec (mp_bitcnt_t) ;


 void __gmpf_set_prec (mpf_ptr, mp_bitcnt_t);


 void __gmpf_set_prec_raw (mpf_ptr, mp_bitcnt_t) ;


 void __gmpf_set_q (mpf_ptr, mpq_srcptr);


 void __gmpf_set_si (mpf_ptr, signed long int);


 int __gmpf_set_str (mpf_ptr, const char *, int);


 void __gmpf_set_ui (mpf_ptr, unsigned long int);


 void __gmpf_set_z (mpf_ptr, mpz_srcptr);


 size_t __gmpf_size (mpf_srcptr) __attribute__ ((__pure__));


 void __gmpf_sqrt (mpf_ptr, mpf_srcptr);


 void __gmpf_sqrt_ui (mpf_ptr, unsigned long int);


 void __gmpf_sub (mpf_ptr, mpf_srcptr, mpf_srcptr);


 void __gmpf_sub_ui (mpf_ptr, mpf_srcptr, unsigned long int);


 void __gmpf_swap (mpf_ptr, mpf_ptr) ;


 void __gmpf_trunc (mpf_ptr, mpf_srcptr);


 void __gmpf_ui_div (mpf_ptr, unsigned long int, mpf_srcptr);


 void __gmpf_ui_sub (mpf_ptr, unsigned long int, mpf_srcptr);


 void __gmpf_urandomb (mpf_t, gmp_randstate_t, mp_bitcnt_t);
# 1501 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
 mp_limb_t __gmpn_add (mp_ptr, mp_srcptr, mp_size_t, mp_srcptr,mp_size_t);




 mp_limb_t __gmpn_add_1 (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t) ;



 mp_limb_t __gmpn_add_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);


 mp_limb_t __gmpn_addmul_1 (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);



 int __gmpn_cmp (mp_srcptr, mp_srcptr, mp_size_t) __attribute__ ((__pure__));






 mp_limb_t __gmpn_divexact_by3c (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);





 mp_limb_t __gmpn_divrem (mp_ptr, mp_size_t, mp_ptr, mp_size_t, mp_srcptr, mp_size_t);


 mp_limb_t __gmpn_divrem_1 (mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_limb_t);


 mp_limb_t __gmpn_divrem_2 (mp_ptr, mp_size_t, mp_ptr, mp_size_t, mp_srcptr);


 mp_size_t __gmpn_gcd (mp_ptr, mp_ptr, mp_size_t, mp_ptr, mp_size_t);


 mp_limb_t __gmpn_gcd_1 (mp_srcptr, mp_size_t, mp_limb_t) __attribute__ ((__pure__));


 mp_limb_t __gmpn_gcdext_1 (mp_limb_signed_t *, mp_limb_signed_t *, mp_limb_t, mp_limb_t);


 mp_size_t __gmpn_gcdext (mp_ptr, mp_ptr, mp_size_t *, mp_ptr, mp_size_t, mp_ptr, mp_size_t);


 size_t __gmpn_get_str (unsigned char *, int, mp_ptr, mp_size_t);


 mp_bitcnt_t __gmpn_hamdist (mp_srcptr, mp_srcptr, mp_size_t) __attribute__ ((__pure__));


 mp_limb_t __gmpn_lshift (mp_ptr, mp_srcptr, mp_size_t, unsigned int);


 mp_limb_t __gmpn_mod_1 (mp_srcptr, mp_size_t, mp_limb_t) __attribute__ ((__pure__));


 mp_limb_t __gmpn_mul (mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);


 mp_limb_t __gmpn_mul_1 (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);


 void __gmpn_mul_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);


 void __gmpn_sqr (mp_ptr, mp_srcptr, mp_size_t);



 mp_limb_t __gmpn_neg (mp_ptr, mp_srcptr, mp_size_t);




 void __gmpn_com (mp_ptr, mp_srcptr, mp_size_t);



 int __gmpn_perfect_square_p (mp_srcptr, mp_size_t) __attribute__ ((__pure__));


 int __gmpn_perfect_power_p (mp_srcptr, mp_size_t) __attribute__ ((__pure__));


 mp_bitcnt_t __gmpn_popcount (mp_srcptr, mp_size_t) __attribute__ ((__pure__));


 mp_size_t __gmpn_pow_1 (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t, mp_ptr);



 mp_limb_t __gmpn_preinv_mod_1 (mp_srcptr, mp_size_t, mp_limb_t, mp_limb_t) __attribute__ ((__pure__));


 void __gmpn_random (mp_ptr, mp_size_t);


 void __gmpn_random2 (mp_ptr, mp_size_t);


 mp_limb_t __gmpn_rshift (mp_ptr, mp_srcptr, mp_size_t, unsigned int);


 mp_bitcnt_t __gmpn_scan0 (mp_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 mp_bitcnt_t __gmpn_scan1 (mp_srcptr, mp_bitcnt_t) __attribute__ ((__pure__));


 mp_size_t __gmpn_set_str (mp_ptr, const unsigned char *, size_t, int);


 mp_size_t __gmpn_sqrtrem (mp_ptr, mp_ptr, mp_srcptr, mp_size_t);



 mp_limb_t __gmpn_sub (mp_ptr, mp_srcptr, mp_size_t, mp_srcptr,mp_size_t);




 mp_limb_t __gmpn_sub_1 (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t) ;



 mp_limb_t __gmpn_sub_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);


 mp_limb_t __gmpn_submul_1 (mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);


 void __gmpn_tdiv_qr (mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);


 void __gmpn_and_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_andn_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_nand_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_ior_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_iorn_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_nior_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_xor_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

 void __gmpn_xnor_n (mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);


 void __gmpn_copyi (mp_ptr, mp_srcptr, mp_size_t);

 void __gmpn_copyd (mp_ptr, mp_srcptr, mp_size_t);

 void __gmpn_zero (mp_ptr, mp_size_t);
# 1681 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
extern __inline__ __attribute__ ((__gnu_inline__)) void
__gmpz_abs (mpz_ptr __gmp_w, mpz_srcptr __gmp_u)
{
  if (__gmp_w != __gmp_u)
    __gmpz_set (__gmp_w, __gmp_u);
  __gmp_w->_mp_size = ((__gmp_w->_mp_size) >= 0 ? (__gmp_w->_mp_size) : -(__gmp_w->_mp_size));
}
# 1705 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
extern __inline__ __attribute__ ((__gnu_inline__))

int
__gmpz_fits_uint_p (mpz_srcptr __gmp_z)
{
  mp_size_t __gmp_n = __gmp_z->_mp_size; mp_ptr __gmp_p = __gmp_z->_mp_d; return (__gmp_n == 0 || (__gmp_n == 1 && __gmp_p[0] <= (~ (unsigned) 0)));;
}




extern __inline__ __attribute__ ((__gnu_inline__))

int
__gmpz_fits_ulong_p (mpz_srcptr __gmp_z)
{
  mp_size_t __gmp_n = __gmp_z->_mp_size; mp_ptr __gmp_p = __gmp_z->_mp_d; return (__gmp_n == 0 || (__gmp_n == 1 && __gmp_p[0] <= (~ (unsigned long) 0)));;
}




extern __inline__ __attribute__ ((__gnu_inline__))

int
__gmpz_fits_ushort_p (mpz_srcptr __gmp_z)
{
  mp_size_t __gmp_n = __gmp_z->_mp_size; mp_ptr __gmp_p = __gmp_z->_mp_d; return (__gmp_n == 0 || (__gmp_n == 1 && __gmp_p[0] <= ((unsigned short) ~0)));;
}




extern __inline__ __attribute__ ((__gnu_inline__))

unsigned long
__gmpz_get_ui (mpz_srcptr __gmp_z)
{
  mp_ptr __gmp_p = __gmp_z->_mp_d;
  mp_size_t __gmp_n = __gmp_z->_mp_size;
  mp_limb_t __gmp_l = __gmp_p[0];






  return (__gmp_n != 0 ? __gmp_l : 0);
# 1761 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
}




extern __inline__ __attribute__ ((__gnu_inline__))

mp_limb_t
__gmpz_getlimbn (mpz_srcptr __gmp_z, mp_size_t __gmp_n)
{
  mp_limb_t __gmp_result = 0;
  if (__builtin_expect ((__gmp_n >= 0 && __gmp_n < ((__gmp_z->_mp_size) >= 0 ? (__gmp_z->_mp_size) : -(__gmp_z->_mp_size))) != 0, 1))
    __gmp_result = __gmp_z->_mp_d[__gmp_n];
  return __gmp_result;
}



extern __inline__ __attribute__ ((__gnu_inline__)) void
__gmpz_neg (mpz_ptr __gmp_w, mpz_srcptr __gmp_u)
{
  if (__gmp_w != __gmp_u)
    __gmpz_set (__gmp_w, __gmp_u);
  __gmp_w->_mp_size = - __gmp_w->_mp_size;
}




extern __inline__ __attribute__ ((__gnu_inline__))

int
__gmpz_perfect_square_p (mpz_srcptr __gmp_a)
{
  mp_size_t __gmp_asize;
  int __gmp_result;

  __gmp_asize = __gmp_a->_mp_size;
  __gmp_result = (__gmp_asize >= 0);
  if (__builtin_expect ((__gmp_asize > 0) != 0, 1))
    __gmp_result = __gmpn_perfect_square_p (__gmp_a->_mp_d, __gmp_asize);
  return __gmp_result;
}




extern __inline__ __attribute__ ((__gnu_inline__))

mp_bitcnt_t
__gmpz_popcount (mpz_srcptr __gmp_u)
{
  mp_size_t __gmp_usize;
  mp_bitcnt_t __gmp_result;

  __gmp_usize = __gmp_u->_mp_size;
  __gmp_result = (__gmp_usize < 0 ? (~ (unsigned long) 0) : 0);
  if (__builtin_expect ((__gmp_usize > 0) != 0, 1))
    __gmp_result = __gmpn_popcount (__gmp_u->_mp_d, __gmp_usize);
  return __gmp_result;
}




extern __inline__ __attribute__ ((__gnu_inline__))

void
__gmpz_set_q (mpz_ptr __gmp_w, mpq_srcptr __gmp_u)
{
  __gmpz_tdiv_q (__gmp_w, (&((__gmp_u)->_mp_num)), (&((__gmp_u)->_mp_den)));
}




extern __inline__ __attribute__ ((__gnu_inline__))

size_t
__gmpz_size (mpz_srcptr __gmp_z)
{
  return ((__gmp_z->_mp_size) >= 0 ? (__gmp_z->_mp_size) : -(__gmp_z->_mp_size));
}






extern __inline__ __attribute__ ((__gnu_inline__)) void
__gmpq_abs (mpq_ptr __gmp_w, mpq_srcptr __gmp_u)
{
  if (__gmp_w != __gmp_u)
    __gmpq_set (__gmp_w, __gmp_u);
  __gmp_w->_mp_num._mp_size = ((__gmp_w->_mp_num._mp_size) >= 0 ? (__gmp_w->_mp_num._mp_size) : -(__gmp_w->_mp_num._mp_size));
}



extern __inline__ __attribute__ ((__gnu_inline__)) void
__gmpq_neg (mpq_ptr __gmp_w, mpq_srcptr __gmp_u)
{
  if (__gmp_w != __gmp_u)
    __gmpq_set (__gmp_w, __gmp_u);
  __gmp_w->_mp_num._mp_size = - __gmp_w->_mp_num._mp_size;
}
# 2103 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
extern __inline__ __attribute__ ((__gnu_inline__))

mp_limb_t
__gmpn_add (mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize)
{
  mp_limb_t __gmp_c;
  do { mp_size_t __gmp_i; mp_limb_t __gmp_x; __gmp_i = (__gmp_ysize); if (__gmp_i != 0) { if (__gmpn_add_n (__gmp_wp, __gmp_xp, __gmp_yp, __gmp_i)) { do { if (__gmp_i >= (__gmp_xsize)) { (__gmp_c) = 1; goto __gmp_done; } __gmp_x = (__gmp_xp)[__gmp_i]; } while ((((__gmp_wp)[__gmp_i++] = (__gmp_x + 1) & ((~ ((mp_limb_t) (0))) >> 0)) == 0)); } } if ((__gmp_wp) != (__gmp_xp)) do { mp_size_t __gmp_j; ; for (__gmp_j = (__gmp_i); __gmp_j < (__gmp_xsize); __gmp_j++) (__gmp_wp)[__gmp_j] = (__gmp_xp)[__gmp_j]; } while (0); (__gmp_c) = 0; __gmp_done: ; } while (0);
  return __gmp_c;
}




extern __inline__ __attribute__ ((__gnu_inline__))

mp_limb_t
__gmpn_add_1 (mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n)
{
  mp_limb_t __gmp_c;
  do { mp_size_t __gmp_i; mp_limb_t __gmp_x, __gmp_r; __gmp_x = (__gmp_src)[0]; __gmp_r = __gmp_x + (__gmp_n); (__gmp_dst)[0] = __gmp_r; if (((__gmp_r) < ((__gmp_n)))) { (__gmp_c) = 1; for (__gmp_i = 1; __gmp_i < (__gmp_size);) { __gmp_x = (__gmp_src)[__gmp_i]; __gmp_r = __gmp_x + 1; (__gmp_dst)[__gmp_i] = __gmp_r; ++__gmp_i; if (!((__gmp_r) < (1))) { if ((__gmp_src) != (__gmp_dst)) do { mp_size_t __gmp_j; ; for (__gmp_j = (__gmp_i); __gmp_j < (__gmp_size); __gmp_j++) (__gmp_dst)[__gmp_j] = (__gmp_src)[__gmp_j]; } while (0); (__gmp_c) = 0; break; } } } else { if ((__gmp_src) != (__gmp_dst)) do { mp_size_t __gmp_j; ; for (__gmp_j = (1); __gmp_j < (__gmp_size); __gmp_j++) (__gmp_dst)[__gmp_j] = (__gmp_src)[__gmp_j]; } while (0); (__gmp_c) = 0; } } while (0);
  return __gmp_c;
}




extern __inline__ __attribute__ ((__gnu_inline__))

int
__gmpn_cmp (mp_srcptr __gmp_xp, mp_srcptr __gmp_yp, mp_size_t __gmp_size)
{
  int __gmp_result;
  do { mp_size_t __gmp_i; mp_limb_t __gmp_x, __gmp_y; (__gmp_result) = 0; __gmp_i = (__gmp_size); while (--__gmp_i >= 0) { __gmp_x = (__gmp_xp)[__gmp_i]; __gmp_y = (__gmp_yp)[__gmp_i]; if (__gmp_x != __gmp_y) { (__gmp_result) = (__gmp_x > __gmp_y ? 1 : -1); break; } } } while (0);
  return __gmp_result;
}




extern __inline__ __attribute__ ((__gnu_inline__))

mp_limb_t
__gmpn_sub (mp_ptr __gmp_wp, mp_srcptr __gmp_xp, mp_size_t __gmp_xsize, mp_srcptr __gmp_yp, mp_size_t __gmp_ysize)
{
  mp_limb_t __gmp_c;
  do { mp_size_t __gmp_i; mp_limb_t __gmp_x; __gmp_i = (__gmp_ysize); if (__gmp_i != 0) { if (__gmpn_sub_n (__gmp_wp, __gmp_xp, __gmp_yp, __gmp_i)) { do { if (__gmp_i >= (__gmp_xsize)) { (__gmp_c) = 1; goto __gmp_done; } __gmp_x = (__gmp_xp)[__gmp_i]; } while ((((__gmp_wp)[__gmp_i++] = (__gmp_x - 1) & ((~ ((mp_limb_t) (0))) >> 0)), __gmp_x == 0)); } } if ((__gmp_wp) != (__gmp_xp)) do { mp_size_t __gmp_j; ; for (__gmp_j = (__gmp_i); __gmp_j < (__gmp_xsize); __gmp_j++) (__gmp_wp)[__gmp_j] = (__gmp_xp)[__gmp_j]; } while (0); (__gmp_c) = 0; __gmp_done: ; } while (0);
  return __gmp_c;
}




extern __inline__ __attribute__ ((__gnu_inline__))

mp_limb_t
__gmpn_sub_1 (mp_ptr __gmp_dst, mp_srcptr __gmp_src, mp_size_t __gmp_size, mp_limb_t __gmp_n)
{
  mp_limb_t __gmp_c;
  do { mp_size_t __gmp_i; mp_limb_t __gmp_x, __gmp_r; __gmp_x = (__gmp_src)[0]; __gmp_r = __gmp_x - (__gmp_n); (__gmp_dst)[0] = __gmp_r; if (((__gmp_x) < ((__gmp_n)))) { (__gmp_c) = 1; for (__gmp_i = 1; __gmp_i < (__gmp_size);) { __gmp_x = (__gmp_src)[__gmp_i]; __gmp_r = __gmp_x - 1; (__gmp_dst)[__gmp_i] = __gmp_r; ++__gmp_i; if (!((__gmp_x) < (1))) { if ((__gmp_src) != (__gmp_dst)) do { mp_size_t __gmp_j; ; for (__gmp_j = (__gmp_i); __gmp_j < (__gmp_size); __gmp_j++) (__gmp_dst)[__gmp_j] = (__gmp_src)[__gmp_j]; } while (0); (__gmp_c) = 0; break; } } } else { if ((__gmp_src) != (__gmp_dst)) do { mp_size_t __gmp_j; ; for (__gmp_j = (1); __gmp_j < (__gmp_size); __gmp_j++) (__gmp_dst)[__gmp_j] = (__gmp_src)[__gmp_j]; } while (0); (__gmp_c) = 0; } } while (0);
  return __gmp_c;
}




extern __inline__ __attribute__ ((__gnu_inline__))

mp_limb_t
__gmpn_neg (mp_ptr __gmp_rp, mp_srcptr __gmp_up, mp_size_t __gmp_n)
{
  mp_limb_t __gmp_ul, __gmp_cy;
  __gmp_cy = 0;
  do {
      __gmp_ul = *__gmp_up++;
      *__gmp_rp++ = -__gmp_ul - __gmp_cy;
      __gmp_cy |= __gmp_ul != 0;
  } while (--__gmp_n != 0);
  return __gmp_cy;
}
# 2260 "c:\\mingw\\bin\\../lib/gcc/mingw32/4.5.0/../../../../include/gmp.h" 3
enum
{
  GMP_ERROR_NONE = 0,
  GMP_ERROR_UNSUPPORTED_ARGUMENT = 1,
  GMP_ERROR_DIVISION_BY_ZERO = 2,
  GMP_ERROR_SQRT_OF_NEGATIVE = 4,
  GMP_ERROR_INVALID_ARGUMENT = 8
};
# 4 "pollard.c" 2

# 1 "settings.h" 1
# 6 "pollard.c" 2
# 1 "factor_list.h" 1





typedef struct factor_list
{
 mpz_t * value;
 struct factor_list * next;
} factor_list;

factor_list * factor_list_add(factor_list ** f, mpz_t * v);
# 7 "pollard.c" 2
# 1 "pollard.h" 1






int pollard(factor_list ** f, const mpz_t n);
int rho(mpz_t res, const mpz_t N);
void f(mpz_t result, mpz_t x, const mpz_t N);
int floyd(const mpz_t N, mpz_t divisor);
int brent(const mpz_t N, mpz_t divisor);
# 8 "pollard.c" 2






int pollard(factor_list ** f, const mpz_t n)
{

 if ((__builtin_constant_p (1) && (1) == 0 ? ((n)->_mp_size < 0 ? -1 : (n)->_mp_size > 0) : __gmpz_cmp_ui (n,1)) <= 0)
 {
  return 1;
 }

 else if (__gmpz_probab_prime_p(n, 10))
 {
  mpz_t * v = malloc(sizeof(mpz_t));
  __gmpz_init_set(*v, n);
  factor_list_add(f, v);
  return 1;
 }






 mpz_t divisor; __gmpz_init(divisor);
 mpz_t divend; __gmpz_init(divend);

 if (__gmpz_perfect_square_p(n))
 {
  __gmpz_sqrt(divisor, n);
  __gmpz_set(divend, divisor);
 }
 else
 {
  if (! rho(divisor, n))
   return 0;

  __gmpz_divexact(divend, n, divisor);
 }





 int r = pollard(f, divisor);
 __gmpz_clear(divisor);

 if (! r)
  return 0;

 r = pollard(f, divend);
 __gmpz_clear(divend);

 if (! r)
  return 0;

 return 1;
}

int rho(mpz_t result, const mpz_t N)
{

 if ((! (((N)->_mp_size != 0) & ((int) ((N)->_mp_d[0])))))
 {
  __gmpz_set_ui(result, 2);
  return 1;
 }

 mpz_t divisor; __gmpz_init_set_ui(divisor, 1);


 if(! brent(N, divisor))



 {
  __gmpz_clear(divisor);
  return 0;
 }


 __gmpz_set(result, divisor);
 __gmpz_clear(divisor);
 return 1;
}

int brent(const mpz_t N, mpz_t divisor)
{
 unsigned int iterations = 0;
 unsigned int x0 = 2;
 mpz_t power; __gmpz_init_set_ui(power, 1);
 mpz_t lambda; __gmpz_init_set_ui(lambda, 1);
 mpz_t tortoise; __gmpz_init_set_ui(tortoise, x0);
 mpz_t hare; __gmpz_init(hare); f(hare, tortoise, N);
 mpz_t diff; __gmpz_init(diff);
 int ret = 1;

 while((__builtin_constant_p (1) && (1) == 0 ? ((divisor)->_mp_size < 0 ? -1 : (divisor)->_mp_size > 0) : __gmpz_cmp_ui (divisor,1))==0)
 {
  if(iterations++>200000)
  {



   ret = 0;
   break;
  }
  if(__gmpz_cmp(power, lambda)==0)
  {
   __gmpz_set(tortoise, hare);
   __gmpz_mul_ui(power, power, 2);
   __gmpz_set_ui(lambda, 0);
  }
  f(hare, hare, N);
  __gmpz_add_ui(lambda, lambda, 1);

  __gmpz_sub(diff, tortoise, hare);
  __gmpz_abs(diff, diff);






  __gmpz_gcd(divisor,diff,N);

 }
 __gmpz_clear(diff);
 __gmpz_clear(power);
 __gmpz_clear(lambda);
 __gmpz_clear(tortoise);
 __gmpz_clear(hare);
 return __gmpz_cmp(divisor, N)==0?0:ret;
}

int floyd(const mpz_t N, mpz_t divisor)
{
 mpz_t x; __gmpz_init_set_ui(x, 2);
 mpz_t y; __gmpz_init_set_ui(y, 2);

 unsigned long iterations = 0;
 while((__builtin_constant_p (1) && (1) == 0 ? ((divisor)->_mp_size < 0 ? -1 : (divisor)->_mp_size > 0) : __gmpz_cmp_ui (divisor,1)) == 0)
 {
  if (++iterations == 200000)
  {



   __gmpz_clear(x);
   __gmpz_clear(y);
   return 0;
  }





  f(x, x, N);






  f(y, y, N);
  f(y, y, N);




  mpz_t diff; __gmpz_init(diff);
  __gmpz_sub(diff, x, y);
  __gmpz_abs(diff, diff);




  if ((__builtin_constant_p (0) && (0) == 0 ? ((diff)->_mp_size < 0 ? -1 : (diff)->_mp_size > 0) : __gmpz_cmp_ui (diff,0)) == 0)
  {



   __gmpz_clear(x);
   __gmpz_clear(y);
   return 0;
  }

  __gmpz_gcd(divisor, diff, N);
  __gmpz_clear(diff);



 }

 __gmpz_clear(x);
 __gmpz_clear(y);
 return 1;
}

void f(mpz_t result, mpz_t x, const mpz_t N)
{
 __gmpz_set(result, x);
 __gmpz_mul(result, result, result);
 __gmpz_add_ui(result, result, 1);
 __gmpz_mod(result, result, N);
}
