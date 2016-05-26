// Functions in this file is copy from raster_tuning project.

#define _CRT_SECURE_NO_WARNINGS 1  // Suppress MSVC warnings about strcat

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>

#if defined(_WIN32) || defined(__WIN32__)
#  define _WINDOWS_USE_
#  include <windows.h>
#else  // Linux, no Mac OS support
#  undef _WINDOWS_USE_
#  include <sys/types.h>
#  include <sys/stat.h>    // for stat()
#  include <unistd.h>
#endif

#ifndef _WINDOWS_USE_
#include <sys/file.h>   // for checking file lock state
#include <errno.h>
#include <err.h>
#endif

// Convert string to array, won't change the value of the one not given
// c_str format: value1[@position1] ...   (position start from 1)
// Return the number of numbers in c_str, or the negitave position of error
// Only have simple error checking, use with careful
int Str2Arr(const char *c_str, double *a, int size)
{
  std::istringstream sin(c_str);
  std::string sst;
  int cnt=0, j=0;
  while (sin>>sst) {
    std::string::size_type idx = sst.find('@');
    if (idx == std::string::npos) {              // not found '@'
      std::istringstream numin(sst);
      if (++j>size) return -cnt-1;               // out of range
      numin>>a[j-1];
    } else {
      sst.at(idx) = ' ';                         // apart the two numbers
      std::istringstream numin(sst);
      double r=0;
      numin>>r>>j;
      if (j>size || j<1) return -cnt-1;          // out of range
      a[j-1] = r;
    }
    cnt++;
  }
  return cnt;
}

// Convert the positioned number array string to double array
// Actual works are done in Str2Arr(). This function is slow and ugly,
//  consider rewrite it using C++ string
// pp  : start position in argv[]
// vec : restore the numbers in argv to vec[]
// size: size of vec, usually the number of neurons
int ReadOneLongCmdPara(int argc, char *argv[], int pp, double *vec, int size)
{
  // calculate the length of argv[]
  int argvs_sz = 10;  // 10 is just for safe
  for (int k=pp; k<argc; k++) {
    argvs_sz += (int)strlen(argv[k]);
  }
  char *tmp_str = (char*)calloc(argvs_sz, sizeof(char));
  int q = pp;
  while (++pp < argc) {
    if (argv[pp][0]=='-' &&
        !(('0'<=argv[pp][1] && argv[pp][1]<='9') || argv[pp][1]=='.')) {
      pp--;
      break;
    }
    if (argv[pp][0]=='@' || tmp_str[strlen(tmp_str)-1] == '@') {
      strcat(tmp_str, argv[pp]);
    } else {
      strcat(tmp_str, " ");
      strcat(tmp_str, argv[pp]);
    }
  }
  int rt = Str2Arr(tmp_str, vec, size);
  free(tmp_str);
  if (rt==0) {
    printf("Warning: nothing after \"%s\" (expect numbers)\n", argv[q]);
  } else if (rt<0) {
    printf("Error: Index out of range! See \"%s\"\n", argv[q]);
    return -pp;
  }
  return pp;
}

int CheckDirAndCreate(const std::string &filepath)
{
  std::string path = filepath;
  int l = (int)path.size();
  if (l>=2 && path[0]=='.' && path[1]=='/') {         // extract the path name
    l -= 2;
    for (int i=0; i<l; i++) { path[i] = path[i+2]; }
  }
  while (--l>=0) {
    if (path[l] != '\\' && path[l] != '/')
      path.pop_back();
    else
      break;
  }
  if (l < 0) return 0;
  if (l == 0 && path[0] == '.') return 0;
  if (l>1) path[l] = '\0';                // so that we can save as /foo.txt
  struct stat sb;
  if (stat(path.c_str(), &sb) == -1) {
    std::string cl;
    printf("Seems no the \"%s\" directory. Trying to create one for you...", path.c_str());
#ifdef _WINDOWS_USE_
    cl = "md " + path;
#else
    cl = "mkdir -p " + path;
#endif
    int rt = system(cl.c_str());
    if (rt) {
      printf("\n  (Failed: Return value of mkdir: %d)\n", rt);
    } else {
      printf("done\n");
    }
    return rt;
  } else {
    if ((sb.st_mode & S_IFMT) != S_IFDIR) {
      printf("Error: Seems \"%s\" is not a directory!\n", path.c_str());
      return 1;
    }
  }
  return 0;
}

