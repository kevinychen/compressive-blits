// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "alignment-inl.h"
#include "application.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"
#include "library_pseudocounts-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "pssm.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

struct ClusterAppOptions {
    static const int kAssignMatchColsByQuery = -1;

    ClusterAppOptions() { Init(); }
    virtual ~ClusterAppOptions() {}

    // Set csbuild default parameters
    void Init() {
        informat         = "auto";
        outformat        = "seq";
        pc_admix         = 0.30;
        pc_ali           = 4.0;
        pc_engine        = "auto";
        match_assign     = kAssignMatchColsByQuery;
        weight_center    = 1.6;
        weight_decay     = 0.85;
        weight_as        = 1000.0;
        verbose          = true;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (infile.empty()) throw Exception("No input file provided!");
        if (alphabetfile.empty()) throw Exception("No abstract states provided!");
    }

    // The input alignment file with training data.
    string infile;
    // The output file.
    string outfile;
    // The file to which the output should be appended.
    string appendfile;
    // Input file with context profile library or HMM for generating pseudocounts
    string modelfile;
    // Input file with profile library to be used as abstract state alphabet
    string alphabetfile;
    // Input file format
    string informat;
    // Output file format (abstract state sequence or abstract state profile)
    string outformat;
    // Overall pseudocount admixture
    double pc_admix;
    // Constant in pseudocount calculation for alignments
    double pc_ali;
    // Pseudocount engine
    string pc_engine;
    // Match column assignment for FASTA alignments
    int match_assign;
    // Weight of central column in multinomial emission
    double weight_center;
    // Exponential decay of window weights
    double weight_decay;
    // Weight in emission calculation of abstract states
    double weight_as;
    // verbose output
    bool verbose;
};  // ClusterAppOptions


template<class Abc>
class ClusterApp : public Application {
  private:
    // Runs the csbuild application.
    virtual int Run();
    // Parses command line options.
    virtual void ParseOptions(GetOpt_pp& ops);
    // Prints options summary to stream.
    virtual void PrintOptions() const;
    // Prints short application description.
    virtual void PrintBanner() const;
    // Prints usage banner to stream.
    virtual void PrintUsage() const;
    // Profile output helper methods
    void OutputProfile1(const CountProfile<Abc>& prof, FILE *fout) const;
    void OutputProfile2(const CountProfile<Abc>& prof, FILE *fout) const;
    void OutputProfile3(const CountProfile<Abc>& prof, FILE *fout) const;
    // Writes profile to outfile
    void WriteProfile(const CountProfile<Abc>& prof) const;

    // Parameter wrapper
    ClusterAppOptions opts_;
    // Profile library with abstract states
    scoped_ptr< ContextLibrary<Abc> > as_lib_;
    // Profile library for context pseudocounts
    scoped_ptr< ContextLibrary<Abc> > pc_lib_;
    // Pseudocount engine
    scoped_ptr< Pseudocounts<Abc> > pc_;
};  // class ClusterApp



template<class Abc>
void ClusterApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('i', "infile", opts_.infile, opts_.infile);
    ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
    ops >> Option('a', "appendfile", opts_.appendfile, opts_.appendfile);
    ops >> Option('I', "informat", opts_.informat, opts_.informat);
    ops >> Option('O', "outformat", opts_.outformat, opts_.outformat);
    ops >> Option('M', "match-assign", opts_.match_assign, opts_.match_assign);
    ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
    ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
    ops >> Option('A', "alphabet", opts_.alphabetfile, opts_.alphabetfile);
    ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
    ops >> Option('p', "pc-engine", opts_.pc_engine, opts_.pc_engine);
    ops >> Option('w', "weight", opts_.weight_as, opts_.weight_as);
    ops >> Option('v', "verbose", opts_.verbose, opts_.verbose);

    opts_.Validate();

    if (strcmp(opts_.outfile.c_str(), "stdout") == 0)
        opts_.verbose = false;

    if (opts_.outfile.empty() && opts_.appendfile.empty())
        opts_.outfile = GetBasename(opts_.infile, false) + ".as";
    if (opts_.informat == "auto")
        opts_.informat = GetFileExt(opts_.infile);
    if (opts_.pc_engine == "auto" && !opts_.modelfile.empty())
        opts_.pc_engine = GetFileExt(opts_.modelfile);
}

template<class Abc>
void ClusterApp<Abc>::PrintBanner() const {
    fputs("Cluster profile columns", out_);
}

template<class Abc>
void ClusterApp<Abc>::PrintUsage() const {
    fputs("Usage: cluster -i <infile> -A <alphabetlib> [options]\n", out_);
}

template<class Abc>
void ClusterApp<Abc>::PrintOptions() const {
    fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
            "Input file with alignment or sequence");
    fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "Output file for generated abstract state sequence (def: <infile>.as)");
    fprintf(out_, "  %-30s %s\n", "-a, --append <file>", "Append generated abstract state sequence to this file");
    fprintf(out_, "  %-30s %s (def=%s)\n", "-I, --informat prf|seq|fas|...", "Input format: prf, seq, fas, a2m, or a3m", opts_.informat.c_str());
    fprintf(out_, "  %-30s %s (def=%s)\n", "-O, --outformat seq|prf", "Outformat: abstract state sequence or profile", opts_.outformat.c_str());
    fprintf(out_, "  %-30s %s\n", "-M, --match-assign [0:100]", "Make all FASTA columns with less than X% gaps match columns");
    fprintf(out_, "  %-30s %s\n", "", "(def: make columns with residue in first sequence match columns)");
    fprintf(out_, "  %-30s %s (def=off)\n", "-A, --alphabet <file>", "Abstract state alphabet consisting of exactly 219 states");
    fprintf(out_, "  %-30s %s (def=off)\n", "-D, --context-data <file>", "Add context-specific pseudocounts using given context-data");
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]", "Pseudocount admix for context-specific pseudocounts", opts_.pc_admix);
    fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[", "Constant in pseudocount calculation for alignments", opts_.pc_ali);
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-w, --weight [0,inf[", "Weight of abstract state column in emission calculation", opts_.weight_as);
}

template<class Abc>
void ClusterApp<Abc>::OutputProfile1(const CountProfile<Abc> &profile, FILE *fout) const {
    profile.Write(fout);
}

template<class Abc>
void ClusterApp<Abc>::OutputProfile2(const CountProfile<Abc> &profile, FILE *fout) const {
    // Print counts matrix and neff vector as negative logs scaled by 'kScale'
    const Profile<Abc> counts(profile.counts);
    const Vector<double> neff(profile.neff);
    for (size_t i = 0; i < counts.length(); ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a) {
            fprintf(fout, "%d\t", iround((double) counts[i][a] / neff[i] * kScale));
        }
        fprintf(fout, "\n");
    }
}

template<class Abc>
void ClusterApp<Abc>::OutputProfile3(const CountProfile<Abc> &profile, FILE *fout) const {
    // Convert profile to sequence of characters from a K-letter alphabet
    const Profile<Abc> counts(profile.counts);
    const Vector<double> neff(profile.neff);
    char buf[4096];
    for (size_t i = 0; i < counts.length(); ++i) {
        size_t max_a = 0, max_val = 0, val;
        for (size_t a = 0; a < Abc::kSize; ++a) {
            val = iround((double) counts[i][a] / neff[i] * kScale);
            if (val > max_val) {
                max_a = a;
                max_val = val;
            }
        }
        switch (Abc::kIntToChar[max_a]) {
            case 'A':
            case 'I':
            case 'L':
            case 'M':
            case 'V':
//                fputs("A", fout);
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
                printf("char: %c\n", Abc::kIntToChar[max_a]);
                buf[i] = '-';
        }
    }
    buf[counts.length()] = '\0';
    fputs(buf, fout);
}

template<class Abc>
void ClusterApp<Abc>::WriteProfile(const CountProfile<Abc> &profile) const {
    if (!opts_.outfile.empty()) {
        FILE* fout;
        if (strcmp(opts_.outfile.c_str(), "stdout") == 0)
            fout = stdout;
        else
            fout = fopen(opts_.outfile.c_str(), "w");
        if (!fout) throw Exception("Can't write to output file '%s'!", opts_.outfile.c_str());

        OutputProfile3(profile, fout);

        if (opts_.verbose)
            fprintf(out_, "Wrote abstract state sequence to %s\n", opts_.outfile.c_str());
        fclose(fout);
    }
    if (!opts_.appendfile.empty()) {
        throw Exception("not implemented");
    }
}

template<class Abc>
int ClusterApp<Abc>::Run() {
    // Setup pseudocount engine
    if (!opts_.modelfile.empty() && opts_.pc_engine == "lib") {
        if (opts_.verbose)
            fprintf(out_, "Reading context library for pseudocounts from %s ...\n",
                GetBasename(opts_.modelfile).c_str());
        FILE* fin = fopen(opts_.modelfile.c_str(), "r");
        if (!fin)
            throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
        pc_lib_.reset(new ContextLibrary<Abc>(fin));
        TransformToLog(*pc_lib_);
        fclose(fin);
        pc_.reset(new LibraryPseudocounts<Abc>(*pc_lib_, opts_.weight_center,
                                               opts_.weight_decay));

    }

    // Setup abstract state engine
    if (opts_.verbose)
        fprintf(out_, "Reading abstract state alphabet from %s ...\n",
            GetBasename(opts_.alphabetfile).c_str());
    FILE* fin = fopen(opts_.alphabetfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!",
                              opts_.alphabetfile.c_str());
    as_lib_.reset(new ContextLibrary<Abc>(fin));
    if (as_lib_->size() != AS219::kSize)
        throw Exception("Abstract state alphabet should have %zu states but actually "
                        "has %zu states!", AS219::kSize, as_lib_->size());
    if (static_cast<int>(as_lib_->wlen()) != 1)
        throw Exception("Abstract state alphabet should have a window length of %zu "
                        "but actually has %zu!", 1, as_lib_->wlen());
    TransformToLog(*as_lib_);
    fclose(fin);

    string header;
    CountProfile<Abc> profile;  // input profile we want to translate

    if (strcmp(opts_.infile.c_str(), "stdin") == 0)
        fin = stdin;
    else
        fin = fopen(opts_.infile.c_str(), "r");
    if (!fin)
        throw Exception("Unable to read input file '%s'!", opts_.infile.c_str());

    if (opts_.informat == "prf") {  // read count profile from infile
        profile = CountProfile<Abc>(fin);;
        if (profile.name.empty()) header = GetBasename(opts_.infile, false);
        else header = profile.name;

        if (pc_) {
            if (opts_.verbose)
                fprintf(out_, "Adding cs-pseudocounts (admix=%.2f) ...\n", opts_.pc_admix);
            CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
            profile.counts = pc_->AddTo(profile, admix);
            Normalize(profile.counts, profile.neff);
        }

    } else if (opts_.informat == "seq") {  // build profile from sequence
        Sequence<Abc> seq(fin);
        header = seq.header();
        profile = CountProfile<Abc>(seq);

        if (pc_) {
            if (opts_.verbose)
                fprintf(out_, "Adding cs-pseudocounts (admix=%.2f) ...\n", opts_.pc_admix);
            profile.counts = pc_->AddTo(seq, ConstantAdmix(opts_.pc_admix));
        }

    } else {  // build profile from alignment
        AlignmentFormat f = AlignmentFormatFromString(opts_.informat);
        Alignment<Abc> ali(fin, f);
        header = ali.name();

        if (f == FASTA_ALIGNMENT) {
            if (opts_.match_assign == ClusterAppOptions::kAssignMatchColsByQuery)
                ali.AssignMatchColumnsBySequence(0);
            else
                ali.AssignMatchColumnsByGapRule(opts_.match_assign);
        }
        profile = CountProfile<Abc>(ali);

        if (pc_) {
            if (opts_.verbose)
                fprintf(out_, "Adding cs-pseudocounts (admix=%.2f) ...\n", opts_.pc_admix);
            CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
            profile.counts = pc_->AddTo(profile, admix);
            Normalize(profile.counts, profile.neff);
        }
    }
    fclose(fin);  // close input file

    /*
    for (int col = 0; col < 5; col++) {
        printf("Column %d\n", col);
        for (int i = 0; i < Abc::kAny; i++)
            printf("%d: %.4f\n", i, profile.counts[col][i]);
        printf("\n");
    }
    */

    WriteProfile(profile);

    return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
    return cs::ClusterApp<cs::AA>().main(argc, argv, stdout, "cluster");
}
