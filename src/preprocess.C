#define PTHREAD
#define MAIN
#define HHBLITS

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <sstream>
#include <vector>
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror(), strerror(errno)
#include <cassert>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::vector;
using std::pair;

extern "C" {
#include <ffindex.h>     // fast index-based database reading
}

#include "cs.h"          // context-specific pseudocounts
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "abstract_state_matrix.h"

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hash.C"          // context-specific pseudocounts
#include "hhdecl.C"      // Constants, global variables, struct Parameters
#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM

#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfullalignment.h" // class FullAlignment
#include "hhhitlist.h"   // class Hit

#include "hhhmm.C"	 // class HMM
#include "hhalignment.C" // class Alignment

#include "preprocess.h"

typedef struct {
    char *name;
    char *subseq;
    int pos;
} Loc;

vector<Loc*> locs[H];
int locs_size = 0;

static int flush() {
    printf("FLUSHING\n");
    char path_buf[16];
    for (int i = 0; i < (1 << (2 * (K - T))); i++) {
        // Flush to filesystem
        sprintf(path_buf, "store/%03x", i);
        FILE *fout = fopen(path_buf, "a");
        if (!fout)
            return 1;
        for (int j = 0; j < (1 << (2 * T)); j++) {
            int index = (i << (2 * T)) + j;
            for (int k = 0; (size_t) k < locs[index].size(); k++) {
                fprintf(fout, "%x %s %s %d\n",
                        index, locs[index][k]->name,
                        locs[index][k]->subseq, locs[index][k]->pos);
            }
        }
        fclose(fout);
    }
    for (int i = 0; i < H; i++) {
        for (int j = 0; (size_t) j < locs[i].size(); j++) {
            delete[] locs[i][j]->name;
            delete[] locs[i][j]->subseq;
            free(locs[i][j]);
        }
        locs[i].clear();
    }
    locs_size = 0;
    return 0;
}

int read_cs4(FILE *fin) {
    static int count = 0;
    if (++count % 1000 == 0)
        printf("read %d cs4s\n", count);

    char name_buf[1 << 16];
    char buf[1 << 16];

    // read name
    if (fgetline(name_buf, sizeof(name_buf), fin) == NULL)
        return 1;
    char *name_ptr = name_buf;
    while (*name_ptr != '|')
        name_ptr++;
    strcpy(name_ptr + 10, ".hhm");

    // read seq
    fgetline(buf, sizeof(buf), fin);
    int seq_size = strlen(buf);
    for (int i = 0; i + K < seq_size; i += JUMP) {
        int h = hash(buf, i);
        if (h == -1)
            continue;
        if (locs[h].size() < (size_t) CHAIN_LIM) {
            // Push data into locs table
            char *name = new char[32];
            strcpy(name, name_ptr + 1);
            char *subseq = new char[RR + K + RR + 2];
            char *ptr = subseq;
            for (int j = (i < RR ? 0 : i - RR);
                    j < seq_size && j < i + K + RR; j++) {
                *ptr++ = buf[j];
            }
            *ptr = '\0';
            Loc *loc = (Loc*) malloc(sizeof(Loc));
            loc->name = name;
            loc->subseq = subseq;
            loc->pos = i;
            locs[h].push_back(loc);
            locs_size++;
            if (locs_size > LOCS_LIM) {
                if (flush())
                    return 1;  // error
            }
        }
    }

    return 0;
}

/*
int read_HMM(FILE *fin) {
    HMM *hmm = new HMM();
    hmm->Read(fin, NULL);

    if (hmm->n_seqs == 0) {
        delete hmm;
        return 1;
    }
    while (fgetc(fin) != '\n');  // read 'Done!' line

    char buf[1 << 16];
    int seq_size = to_cs4(hmm, buf);

    for (int i = 0; i + K < seq_size; i += JUMP) {
        int h = hash(buf, i);
        if (h == -1)
            continue;
        if (locs[h].size() < (size_t) CHAIN_LIM) {
            // Copy correct part of name
            char *ptr = hmm->name;
            while (*ptr != '|')
                ptr++;
            char *name = new char[32];
            strcpy(name, ptr + 1);
            ptr = name;
            while (*ptr != '|')
                ptr++;
            strcpy(ptr, ".hhm");

            // Push data into locs table
            Loc *loc = malloc(sizeof(Loc));
            loc->name = name;
            loc->subseq = subseq;
            loc->pos = i;
            locs[h].push_back(loc);
            locs_size++;
            if (locs_size > LOCS_LIM) {
                if (flush())
                    return 1;  // error
            }
        }
    }

    delete hmm;

    return 0;
}
*/

void parseHMM(FILE *fin) {
    HMM *hmm = new HMM();
    hmm->Read(fin, NULL);
    char buf[1 << 16];
    to_cs4(hmm, buf);
    printf("code: %s\n", buf);
}

int main(int argc, char **argv)
{
    FILE *fin = fopen(argv[1], "r");

    while (read_cs4(fin) == 0);
    if (flush())
        printf("error with flushing\n");

    return 0;
}
