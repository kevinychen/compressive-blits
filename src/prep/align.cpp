#include "align.h"

#include <stdio.h>
#include <cstdlib>

const int GAP_PEN = 4;

const int TB_START = 0;
const int TB_DIAG = 1;
const int TB_DOWN = 2;
const int TB_RIGHT = 3;

static int score(char c1, char c2) {
    return c1 == c2 ? 3 : -2;
}

int smith_waterman(const char *seq1, const char *seq2, int len1, int len2) {
    int F[len1 + 1][len2 + 1];
    int TB[len1 + 1][len2 + 1];

    for (int i = 0; i <= len1; i++) {
        F[i][0] = 0;
        TB[i][0] = TB_START;
    }
    for (int j = 0; j <= len2; j++) {
        F[0][j] = 0;
        TB[0][j] = TB_START;
    }

    int best = 0, best_i = -1, best_j = -1;
    for (int i = 1; i <= len1; i++)
        for (int j = 1; j <= len2; j++) {
            int val;

            F[i][j] = 0;
            TB[i][j] = TB_START;

            val = F[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]);
            if (val > F[i][j]) {
                F[i][j] = val;
                TB[i][j] = TB_DIAG;
            }

            val = F[i - 1][j] - GAP_PEN;
            if (val > F[i][j]) {
                F[i][j] = val;
                TB[i][j] = TB_DOWN;
            }

            val = F[i][j - 1] - GAP_PEN;
            if (val > F[i][j]) {
                F[i][j] = val;
                TB[i][j] = TB_RIGHT;
            }

            if (F[i][j] > best) {
                best = F[i][j];
                best_i = i;
                best_j = j;
            }
        }

    char buf1[64], buf2[64];
    int index1 = sizeof(buf1), index2 = sizeof(buf2);
    buf1[--index1] = buf2[--index2] = '\0';
    while (true) {
        switch (TB[best_i][best_j]) {
            case TB_START:
//                printf("%s\n", buf1 + index1);
//                printf("%s\n", buf2 + index2);
//                printf("score: %d\n", best);
                return best;
            case TB_DIAG:
                buf1[--index1] = seq1[--best_i];
                buf2[--index2] = seq2[--best_j];
                break;
            case TB_DOWN:
                buf1[--index1] = seq1[--best_i];
                buf2[--index2] = '-';
                break;
            case TB_RIGHT:
                buf1[--index1] = '-';
                buf2[--index2] = seq2[--best_j];
                break;
            default:
                printf("Error\n");
                exit(1);
        }
    }
}

/*
int main(int argc, char *argv[]) {
    const char *seq1 = "ORANGETHINGS";
    const char *seq2 = "WEGETSWINGS";
    int sw = smith_waterman(seq1, seq2, 12, 11);
    printf("%d\n", sw);
    return 0;
}
*/
