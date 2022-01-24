/*****************************************************************************
  fileType.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the MIT license (as of Jan 2022)
******************************************************************************/

#include "fileType.h"


/*
   returns TRUE if the file is a regular file:
     not a pipe/device.

   This implies that the file can be opened/closed/seek'd multiple times without losing information
 */
bool isRegularFile(const string& filename) {
       struct stat buf ;
       int i;

       i = stat(filename.c_str(), &buf);
       if (i!=0) {
               cerr << "BEDTools Error: can't determine file type of '" << filename << "': " << strerror(errno) << endl;
               return false;
       }
       if (S_ISREG(buf.st_mode))
               return true;

       return false;
}


/*
   returns TRUE if the file has a GZIP header.
   Should only be run on regular files.
 */
bool isGzipFile(const string& filename) {
       //see http://www.gzip.org/zlib/rfc-gzip.html#file-format
       struct  {
               unsigned char id1;
               unsigned char id2;
               unsigned char cm;
       } gzip_header;
       ifstream f(filename.c_str(), ios::in|ios::binary);
       if (!f)
               return false;

       if (!f.read((char*)&gzip_header, sizeof(gzip_header)))
               return false;

       if ( gzip_header.id1 == 0x1f
                       &&
                       gzip_header.id2 == 0x8b
                       &&
                       gzip_header.cm == 8 )
               return true;

       return false;
}
