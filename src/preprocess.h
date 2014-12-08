const int K = 10;  // k-mer length for hashing
const int H = (1 << (2 * K));  // number of distinct hashes
const int T = 4;  // remainder of k-mer (# bits per file)
const int RR = 7;  // window size around target k-mer
const int LL = 32;  // threshhold for window alignment
const int CHAIN_LIM = 64;  // max locs for a given hash
const int JUMP = 2;  // overlap between different k-mers
const int LOCS_LIM = 1 << 24;  // max size of locs table

int ord(char c) {
    return c - 'P';
}

int hash(char *buf, int i) {
    int h = 0;
    for (int j = i; j < i + K; j++) {
        if (buf[j] == '-')
            return -1;
        h = h * 4 + ord(buf[j]);
    }
    return h;
}

int to_cs4(HMM *hmm, char *buf) {
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
    return seq_size;
}


