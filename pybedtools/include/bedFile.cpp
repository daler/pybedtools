/*****************************************************************************
bedFile.cpp

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licensed under the MIT license (as of Jan 2022)
******************************************************************************/
#include "bedFile.h"

/*******************************************
Class methods
*******************************************/

// Constructor
BedFile::BedFile(string bedFile)
: bedFile(bedFile),
_typeIsKnown(false),
_lineNum(0)
{}

// Destructor
BedFile::~BedFile(void) {
}


int BedFile::Open(void) {
    if (bedFile == "stdin") {
        _bedStream = &cin;
    }
    // New method thanks to Assaf Gordon
    else if ((isGzipFile(bedFile) == false) && (isRegularFile(bedFile) == true)) {
       // open an ifstream
        ifstream beds(bedFile.c_str(), ios::in);

        // can we open the file?
        if ( !beds ) {
            cerr << "BEDTools Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            return -1;
        }
        else {
            // if so, close it (this was just a test)
            beds.close();       
            // now set a pointer to the stream so that we
            _bedStream = new ifstream(bedFile.c_str(), ios::in);
        }
    } 
    else if ((isGzipFile(bedFile) == true) && (isRegularFile(bedFile) == true)) {        
        igzstream beds(bedFile.c_str(), ios::in);
        if ( !beds ) {
            cerr << "BEDTools Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            return -1;
        }
        else {
            // if so, close it (this was just a test)
            beds.close();       
            // now set a pointer to the stream so that we
            _bedStream = new igzstream(bedFile.c_str(), ios::in);
        }
    }
    else {
        cerr << "BEDTools Error: Unexpected file type (" << bedFile << "). Exiting!" << endl;
        return -1;
    }
    return 1;
}

// Rewind the pointer back to the beginning of the file
void BedFile::Rewind(void) {
    _bedStream->seekg(0, ios::beg);
}

// Jump to a specific byte in the file
void BedFile::Seek(unsigned long offset) {
    _bedStream->seekg(offset);
}

// Close the BED file
void BedFile::Close(void) {
    if (bedFile != "stdin") delete _bedStream;
}


BED BedFile::GetNextBed() {

    BED bed;
    
    // make sure there are still lines to process.
    // if so, tokenize, validate and return the BED entry.
    if (_bedStream->good()) {
        string bedLine;
        vector<string> bedFields;
        bedFields.reserve(12);

        // parse the bedStream pointer
        getline(*_bedStream, bedLine);
        _lineNum++;

        // split into a string vector.
        Tokenize(bedLine,bedFields);

        // load the BED struct as long as it's a valid BED entry.
        bed.status = parseLine(bed, bedFields);
        bed.fields = bedFields;
        return bed;
    }
    else {
        // default if file is closed or EOF
        bed.status = BED_INVALID;
        return bed;
    }
}

vector<BED> BedFile::FindOverlapsPerBin(const BED &bed, float overlapFraction) {
    vector<BED> hits;

    BIN startBin, endBin;
    startBin = (bed.start   >> _binFirstShift);
    endBin   = ((bed.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::iterator bedItr = bedMap[bed.chrom][j].begin();
            vector<BED>::iterator bedEnd = bedMap[bed.chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size   = (float) bed.end-bed.start;
                int maxStart = max(bed.start, bedItr->start);
                int minEnd   = min(bed.end, bedItr->end);
                int overlap  = minEnd - maxStart;
                float ofrac = (overlap/size);
                // Note: zero-length features with no overlap (e.g., chr1:1-1
                // to chr1:5-500) are ofrac = (-4/0) which in C++ is -inf.
                //
                // Zero-length features with exactly zero overlap (0/0;
                // chr1:1-1 to chr1:1-1) in C++ is -nan.
                //
                // Note that in cbedtools, default overlapFraction is 0 and
                // currently only positive values are supported.
                //
                // Also note that a zero-length feature that overlaps something
                // else will *always* be considered a hit, regardless of
                // overlapFraction.
                if ((ofrac >= overlapFraction) || ( (size == 0) && (overlap == 0)))
                {
                    bedItr->o_start = maxStart;
                    bedItr->o_end   = minEnd;
                    hits.push_back(*bedItr);
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return hits;
}

vector<BED> BedFile::FindOverlapsPerBin(const BED &bed, bool forceStrand, float overlapFraction) {
    vector<BED> hits;

    BIN startBin, endBin;
    startBin = (bed.start >> _binFirstShift);
    endBin = ((bed.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::iterator bedItr = bedMap[bed.chrom][j].begin();
            vector<BED>::iterator bedEnd = bedMap[bed.chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size   = (float) bed.end-bed.start;
                int maxStart = max(bed.start, bedItr->start);
                int minEnd   = min(bed.end, bedItr->end);
                int overlap  = minEnd - maxStart;
                float ofrac = (overlap/size);
                if (
                        (bed.strand == bedItr->strand)
                        && (
                            (ofrac >= overlapFraction)
                            || ((size == 0) && (overlap == 0))
                            )
                   )
                {
                    bedItr->o_start = maxStart;
                    bedItr->o_end   = minEnd;
                    hits.push_back(*bedItr);
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return hits;
}


int BedFile::FindAnyOverlapsPerBin(const BED &bed, float overlapFraction) {
    BIN startBin, endBin;
    startBin = (bed.start   >> _binFirstShift);
    endBin   = ((bed.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[bed.chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[bed.chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size   = (float) bed.end-bed.start;
                int maxStart = max(bed.start, bedItr->start);
                int minEnd   = min(bed.end, bedItr->end);
                int overlap  = minEnd - maxStart;
                float ofrac = (overlap/size);
                if ((ofrac >= overlapFraction) || ( (size == 0) && (overlap == 0)))
                {
                    return 1;
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return 0;
}


int BedFile::FindAnyOverlapsPerBin(const BED &bed, bool forceStrand, float overlapFraction) {
    BIN startBin, endBin;
    startBin = (bed.start >> _binFirstShift);
    endBin = ((bed.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[bed.chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[bed.chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size   = (float) bed.end-bed.start;
                int maxStart = max(bed.start, bedItr->start);
                int minEnd   = min(bed.end, bedItr->end);
                int overlap  = minEnd - maxStart;
                float ofrac = (overlap/size);
                if (
                        (bed.strand == bedItr->strand)
                        && (
                            (ofrac >= overlapFraction)
                            || ((size == 0) && (overlap == 0))
                            )
                   )
                {
                    return 1;
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return 0;
}


int BedFile::CountOverlapsPerBin(const BED &bed, float overlapFraction) {
    BIN startBin, endBin;
    startBin = (bed.start   >> _binFirstShift);
    endBin   = ((bed.end-1) >> _binFirstShift);
    int count = 0;
    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[bed.chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[bed.chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size   = (float) bed.end-bed.start;
                int maxStart = max(bed.start, bedItr->start);
                int minEnd   = min(bed.end, bedItr->end);
                int overlap  = minEnd - maxStart;
                float ofrac = (overlap/size);
                if ((ofrac >= overlapFraction) || ( (size == 0) && (overlap == 0)))
                {
                    count++;
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return count;
}


int BedFile::CountOverlapsPerBin(const BED &bed, bool forceStrand, float overlapFraction) {
    BIN startBin, endBin;
    startBin = (bed.start >> _binFirstShift);
    endBin = ((bed.end-1) >> _binFirstShift);
    int count = 0;
    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[bed.chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[bed.chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size   = (float) bed.end-bed.start;
                int maxStart = max(bed.start, bedItr->start);
                int minEnd   = min(bed.end, bedItr->end);
                int overlap  = minEnd - maxStart;
                float ofrac = (overlap/size);
                if (
                        (bed.strand == bedItr->strand)
                        && (
                            (ofrac >= overlapFraction)
                            || ((size == 0) && (overlap == 0))
                            )
                   )
                {
                    count++;
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return count;
}


void BedFile::setFileType (FileType type) {
    _fileType    = type;
    _typeIsKnown = true;
}


void BedFile::setBedType (int colNums) {
    bedType = colNums;
}

void BedFile::loadBedFileIntoMap() {

    BED bed, nullBed;
    //BedLineStatus bedStatus;

    Open();
    bed = GetNextBed();
    while ( bed.status != BED_INVALID) {
        if (bed.status == BED_VALID) {
            BIN bin = getBin(bed.start, bed.end);
            bedMap[bed.chrom][bin].push_back(bed);
            bed = nullBed;
        }
        bed = GetNextBed();
    }
    Close();
}
