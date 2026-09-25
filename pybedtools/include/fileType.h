/*****************************************************************************
  fileType.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the MIT license (as of Jan 2022)
******************************************************************************/
#ifndef FILETYPE_H
#define FILETYPE_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>

#if !defined(S_ISREG) && defined(_S_IFMT) && defined(_S_IFREG)
#define S_ISREG(m) (((m) & _S_IFMT) == _S_IFREG)
#endif
#if !defined(S_ISDIR) && defined(_S_IFMT) && defined(_S_IFDIR)
#define S_ISDIR(m) (((m) & _S_IFMT) == _S_IFDIR)
#endif

using namespace std;

/*****************************************************************************
  Convenience functions to detect whether a given file is
  "regular" and/or "gzipped".

  Kindly contributed by Assaf Gordon.
******************************************************************************/
string string_error(int errnum);
bool isRegularFile(const string& filename);
bool isGzipFile(const string& filename);

#endif /* FILETYPE_H */
