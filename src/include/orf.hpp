#include <string>

// ORF class
class OpenReadingFrame{
    std::string rna_sequence;
    std::string protein_sequence;
    int start;
    int stop;
    int strand;
    int frame;
public:
    OpenReadingFrame();
    void set_props(std::string rna_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_);
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
    std::string get_protein_sequence();
    void set_protein_sequence(std::string ps);
};