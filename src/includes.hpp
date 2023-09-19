#include <string>
#include <vector>

// ORF class
class OpenReadingFrame{
    std::string rna_sequence;
    std::string parent_sequence;
    std::string protein_sequence;
    int start;
    int stop;
    int strand;
    int frame;
public:
    OpenReadingFrame();
    void set_props(std::string rna_sequence_, std::string parent_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_);
    void set_start(int start_);
    int get_start();
    void set_stop(int stop_);
    int get_stop();
    void set_strand(int strand_);
    int get_strand();
    void set_frame(int frame_);
    int get_frame();
    std::string get_rna_sequence();
    void set_rna_sequence(std::string rs);
    std::string get_parent_sequence();
    void set_parent_sequence(std::string ps);
    std::string get_protein_sequence();
    void set_protein_sequence(std::string ps);
};

// Sequence class
class Sequence {
    std::string name;
    std::string seq_str;
    std::string type;
public:
    Sequence(std::string &seq_str_);
    void set_name(std::string &name_);
    std::string get_name();
    void set_seq(std::string &seq_str_);
    std::string get_seq();
    std::string get_type();

    // return length
    int length();

    // checking equality
    bool const operator==(const Sequence& other) const;

    // checking <
    bool const operator<(const Sequence& other) const;

    // function to set type
    void set_type(std::string seq_type);

    // find codon
    int find_codon(std::string codon);

    // reverse
    Sequence reverse();

    // complement
    Sequence complement();

    // reverse complement
    Sequence reverse_complement();

    // frame shift
    Sequence frame_shift(int frame);

    // translate
    Sequence translate();

    // write fasta
    void write_fasta(std::string file_name);

    // get orfs
    std::vector<OpenReadingFrame> get_open_reading_frames(unsigned int min_len);
};

class SeqSet {
    std::string name;
    // pointers to sequences!
    std::vector<Sequence> records;
    std::string type;
public:
    SeqSet() {};
    void set_name(std::string &name_) {name = name_;}
    std::string get_name() {return name;}
    int size() { return records.size(); }
    void set_type(std::string seq_type=""){
        if ( records.size() == 0){
            std::cout << "WARNING: there are no sequences in the set to determine type." << std::endl;
            return;
        }
        std::string prev;
        for (auto& t : records){
            t.set_type(seq_type);
            if(prev.length() > 0 && t.get_type().compare(prev) != 0){
                std::cout << "WARNING: set contains sequences of different types! The assigned type is probably wrong." << std::endl;
                prev = t.get_type();
            }
        }
        // set the type based on the first sequence
        // (Assumes all sequences have the same type -- TODO: add check.)
        type = records[0].get_type();
    }
    std::string get_type(){ return type; };
    void set_records(std::vector<Sequence> &records_) {records = records_;
    this->set_type();}
    std::vector<Sequence> get_records() { return records;}

    // function to be able to add sequence to the set
    void add_sequence(Sequence seq){
        // TODO: add check for type
        records.push_back(seq);
        return;
    }

    // dealign
    void dealign(){
        std::string seq;
        for(auto& t : records){
            seq = t.get_seq();
            seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
            t.set_seq(seq);
        }
        return;
    }

    // alphabetize
    void alphabetize(){
        std::sort(records.begin(), records.end());
    }

    // remove duplicates
    void remove_duplicates(){
        if(records.size() < 2){
            return;
        }

        // sort
        this->alphabetize();

        // flag
        std::string flag = "**";

        unsigned int idx=0;
        std::string name = records[0].get_name();

        // label duplicates, append their names to the first occurrence
        for (unsigned int i=1; i < records.size(); i++){
            if (records[i] == records[i-1]){
                name += '_' + records[i].get_name();
                records[i].set_name(flag);
            }
            else{
                if(i - idx > 1){
                    records[idx].set_name(name);
                }
                name = records[i].get_name();
                idx=i;
            }
            std::cout << i << " " << name << std::endl;
        }

        // remove the labeled ones
        records.erase(std::remove_if(
            records.begin(), records.end(), 
            [](auto t){return t.get_name().compare("**")==0;}),
            records.end());
    }

    // add elements from another SeqSet instance
    void add_set(SeqSet* sset){
        for (auto t : sset->get_records()){
            records.push_back(t);
        }
    }

    // write fasta
    void write_fasta(std::string file_name){
        std::ofstream out_file;
        out_file.open(file_name);
        for (auto t : records){
            out_file << ">" << name << std::endl;
            for (int i = 0; i < t.length(); i++){
                if (i != 0 && i % 79 == 0){
                    out_file << std::endl;
                }
                out_file << t.get_seq()[i];
            }
            out_file << std::endl;
        }
        out_file.close();
    }

    // read fasta
    void read_fasta(std::string file_name){
        std::ifstream in_file (file_name);
        std::string line;
        std::string name;
        std::string seq;

        // make sure file is opened
        if (!in_file.is_open()){
            std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
            return;
        }

        int seq_ct = 0;
        while(std::getline(in_file, line)){
            if (line[0] == '>')
            {
                if (seq_ct > 0){
                    Sequence new_seq(seq);
                    new_seq.set_name(name);
                    records.push_back(new_seq);
                    seq = "";
                }
                // get new name
                name = line.substr(1, line.length() - 1);
                seq_ct++;
            }
            else{
                // add new line to seq
                seq += line;
            }
        }

        // add the final one
        Sequence new_seq(seq);
        new_seq.set_name(name);
        records.push_back(new_seq);
    }

    // write phylip
    void write_phylip(std::string file_name, std::string mode){
        // in phylip format, all sequences must be of the same length
        int n_seqs = records.size();
        std::unordered_set<int> n_chars;
        unsigned int max_name_chars = 0;

        for (auto t : records){
            n_chars.insert(t.length());
            // get the longest name for later use
            if(t.get_name().length() > max_name_chars) max_name_chars = t.get_name().length();
        }
        if (n_chars.size() != 1){
            std::cout << "ERROR: sequences must be of the same length." << std::endl;
            return;
        }
        std::ofstream out_file;
        out_file.open(file_name);
        auto num_chars = n_chars.begin();

        // first lines differ depending on whether we want the sequential or interleaved versions
        if (mode.compare("sequential") == 0){
            out_file << n_seqs << "\t" << *num_chars << std::endl;
        }
        else if (mode.compare("interleaved") == 0){
            out_file << n_seqs << "\t" << *num_chars << "\tI" <<std::endl; 
        }
        else{
            std::cout << "ERROR: mode not understood." << std::endl;
            return;
        }

        // sequential is pretty straightforward
        if (mode.compare("sequential") == 0){
            for (auto t : records){
                out_file << t.get_name() << "    " << t.get_seq() << std::endl;
            }
        }
        // interleaved -> write in 6 blocks of 10
        else if (mode.compare("interleaved") == 0){
            int total_lines = *num_chars / 60 + 1;
            for (int i=0; i < total_lines; i++){
                for (auto t : records){
                    if (i == 0){
                        out_file << t.get_name();
                        // proper number of spaces
                        for (unsigned int j = 0; j < (max_name_chars - t.get_name().length() + 2); j++) out_file << " ";
                    }
                    else
                    {
                        // proper number of spaces
                        for (unsigned int j = 0; j < (max_name_chars + 2); j++) out_file << " ";
                    }
                    // write seq_str in chunks
                    int start = i * 60;
                    int len = fmin(60, t.length() - start);
                    std::string relevant = t.get_seq().substr(start, len);
                    int count = 0;
                    for(auto c : relevant){
                        if (count != 0 && count % 10 == 0){
                            out_file << " ";
                        }
                        out_file << c;
                        count++;
                    }
                    out_file << std::endl;
                }
                out_file << std::endl;
            }
        }

        out_file.close();
    }

    // read only names and lengths
    std::unordered_map<std::string, unsigned int> get_names_and_lengths_from_fasta(std::string file_name){
        std::unordered_map<std::string, unsigned int> lengths;

        std::ifstream in_file(file_name);
        std::string line;
        std::string name;
        std::string seq;

        // TODO: refactor -- this is nearly a duplicate of read_fasta

        // make sure file is open
        if(!in_file.is_open()){
            std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
            return lengths;
        }

        int seq_ct = 0;
        while(std::getline(in_file, line)){
            if (line[0] == '>')
            {
                if (seq_ct > 0){
                    lengths[name] = seq.length();
                    seq = "";
                }
                name = line.substr(1, line.length() - 1);
                seq_ct++;
            }
            else{
                seq += line;
            }
        }
        return lengths;
    }

    /*
    TODO:
    - filter support for read_fasta
    - __getitem__
    - __add__
    - get frequencies
    - read fastq
    - read phylip?
    - write fastq?
    - dealign (remove all gaps)
    - filter by frequency ?
    - find kmers
    */
};