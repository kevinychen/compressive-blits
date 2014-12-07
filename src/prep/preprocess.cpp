
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "align.h"

using namespace std;

const size_t K = 8;
const size_t R = 24;

vector<string*> seqs;
vector< pair<size_t, size_t> > locs[1 << (2 * K)];  // size 4^K

static int ord(char c) {
    return c - 'P';
}

static int hash(string s, int i) {
    int h = 0;
    for (int j = i; j < i + K; j++)
        h = h * 4 + ord(s[j]);
    return h;
}

int main(int argc, char *argv[]) {
    ifstream in(argv[1]);  // cs4s file
    string s;
    int line = 0;
    int histo[200]; for (int i = 0; i < 200; i++) histo[i] = 0;
    int num_chars = 0;
    int total_saved = 0;
    while (getline(in, s)) {
        if (s[0] == '#')  // header name
            continue;
        num_chars += s.size();
        seqs.push_back(new string(s));

        for (int i = 0; i + K < s.size(); i++) {
            int h = hash(s, i);
            int best_sw = 0;
            for (int j = 0; j < locs[h].size(); j++) {
                pair<size_t, size_t> loc = locs[h][j];

                const char *c_str1 = seqs[loc.first]->c_str();
                const char *c_str2 = s.c_str();

                size_t start1 = loc.second < R ? 0 : loc.second - R;
                size_t start2 = i < R ? 0 : i - R;
                size_t end1 = min(seqs[loc.first]->size(), loc.second + K + R);
                size_t end2 = min(s.size(), i + K + R);

                int sw = smith_waterman(
                        c_str1 + start1,
                        c_str2 + start2,
                        end1 - start1,
                        end2 - start2);
                if (sw > best_sw)
                    best_sw = sw;
            }
            if (locs[h].size() < 10)
                locs[h].push_back( make_pair(line, i) );
            i += R;
            if (best_sw < 200) histo[best_sw]++;
            if (best_sw > 3 * K)
                total_saved += best_sw / 3 - K;
        }

        line++;
        if (line % 1000 == 0)
            printf("line: %d\n", line);
        if (line == 10000)
            break;
    }
    for (int i = 0; i < 200; i++)
        printf("histo %d %d\n", i, histo[i]);
    printf("num chars: %d\n", num_chars);
    printf("total saved: %d\n", total_saved);
    return 0;
}
