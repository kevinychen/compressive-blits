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
using std::make_pair;

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

const int K = 10;  // k-mer length for hashing
const int H = (1 << (2 * K));  // number of distinct hashes
const int T = 4;  // remainder of k-mer (# bits per file)
const int CHAIN_LIM = 64;  // max locs for a given hash
const int JUMP = 4;  // overlap between different k-mers
const int LOCS_LIM = 1 << 24;  // max size of locs table

vector< pair<int, int> > locs[H];
int line = 0;
int locs_size = 0;

static int ord(char c) {
    return c - 'P';
}

static int hash(char *buf, int i) {
    int h = 0;
    for (int j = i; j < i + K; j++) {
        if (buf[j] == '-')
            return -1;
        h = h * 4 + ord(buf[j]);
    }
    return h;
}

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
                fprintf(fout, "%x %d %d\n", index,
                        locs[index][k].first, locs[index][k].second);
            }
        }
        fclose(fout);
    }
    for (int i = 0; i < H; i++)
        locs[i].clear();
    locs_size = 0;
    return 0;
}

int read_HMM(FILE *fin) {
    HMM *hmm = new HMM();
    hmm->Read(fin, NULL);

    if (hmm->n_seqs == 0) {
        delete hmm;
        return 1;
    }
    while (fgetc(fin) != '\n');  // read 'Done!' line
    line++;

    int seq_size = -1;
    for (int i = 0; i < hmm->n_seqs; i++)
        if (i == hmm->ncons) {
            for (int k = 0; hmm->seq[i][k]; k++)
                if (hmm->seq[i][k] >= 'A' && hmm->seq[i][k] <= 'Z')
                    hmm->seq[i][k] += ('a' - 'A');
        } else {
            int j = 0;
            for (int k = 0; hmm->seq[i][k]; k++)
                if (hmm->seq[i][k] == '-' ||
                        (hmm->seq[i][k] >= 'A' && hmm->seq[i][k] <= 'Z')) {
                    hmm->seq[i][j++] = hmm->seq[i][k];
                }
            hmm->seq[i][j] = '\0';
            if (seq_size == -1)
                seq_size = j;
            else if (seq_size != j)
                throw Exception("Mismatched sizes: %d != %d\n", seq_size, j);
        }

    char buf[1 << 16];
    for (int i = 0; i < seq_size; i++) {
        int freqs[128];
        for (int j = 'A'; j <= 'Z'; j++)
            freqs[j] = 0;
        for (int k = 0; k < hmm->n_seqs; k++)
            freqs[hmm->seq[k][i]]++;
        int max_j = 'A';
        for (int j = 'A'; j <= 'Z'; j++)
            if (freqs[j] > freqs[max_j])
                max_j = j;

        switch (max_j) {
            case 'A':
            case 'I':
            case 'L':
            case 'M':
            case 'V':
                buf[i] = 'P';
                break;
            case 'C':
            case 'F':
            case 'W':
            case 'Y':
                buf[i] = 'Q';
                break;
            case 'D':
            case 'E':
            case 'H':
            case 'K':
            case 'R':
                buf[i] = 'R';
                break;
            case 'G':
            case 'N':
            case 'P':
            case 'Q':
            case 'S':
            case 'T':
                buf[i] = 'S';
                break;
            default:
                printf("char: %c\n", max_j);
                buf[i] = '-';
        }
    }

    for (int i = 0; i + K < seq_size; i += JUMP) {
        int h = hash(buf, i);
        if (h == -1)
            continue;
        if (locs[h].size() < (size_t) CHAIN_LIM) {
            locs[h].push_back( make_pair(line, i) );
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

int main(int argc, char **argv)
{
    FILE *fin = fopen(argv[1], "r");

    while (read_HMM(fin) == 0);
    if (flush())
        printf("error with flushing\n");

    return 0;
}
