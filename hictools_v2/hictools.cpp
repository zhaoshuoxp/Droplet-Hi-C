//
// hictools.cpp
// Optimized, Robust, and Auto-detect Version
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "cxstring.hpp" 

using namespace std;

class combine_hic{
private:
    static string clean_readname(const string& raw_name) {
        size_t space_pos = raw_name.find_first_of(" \t");
        string base = (space_pos != string::npos) ? raw_name.substr(0, space_pos) : raw_name;
        
        size_t colon_r_pos = base.rfind(":r");
        if (colon_r_pos != string::npos && colon_r_pos > 0) {
            bool is_barcode = true;
            for (size_t i = colon_r_pos + 2; i < base.length(); ++i) {
                char c = base[i];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                    is_barcode = false; 
                    break;
                }
            }
            if (is_barcode) {
                base = base.substr(0, colon_r_pos);
            }
        }
        return base;
    }

    static bool file_exists(const string& name) {
        if (FILE *file = fopen(name.c_str(), "r")) {
            fclose(file);
            return true;
        }
        return false;
    }

public:
    static void run(string mode, string r2);
};

class convert_hic2{
private:

public:
    static void run(string prefix);
};

void combine_hic::run(string mode, string r2){
    long long total = 0;
    long long pass = 0;
    string cmd;
    string file_r1, file_r2, file_r3;

    if (file_exists(r2 + "_R1.fq.gz")) file_r1 = r2 + "_R1.fq.gz";
    else if (file_exists(r2 + "_R1.fastq.gz")) file_r1 = r2 + "_R1.fastq.gz";
    else { cerr << "Error: Cannot find R1 file (.fq.gz or .fastq.gz) for " << r2 << endl; exit(1); }

    if (file_exists(r2 + "_R2.fq.gz")) file_r2 = r2 + "_R2.fq.gz";
    else if (file_exists(r2 + "_R2.fastq.gz")) file_r2 = r2 + "_R2.fastq.gz";
    else { cerr << "Error: Cannot find R2 file (.fq.gz or .fastq.gz) for " << r2 << endl; exit(1); }

    if (file_exists(r2 + "_R3.fq.gz")) file_r3 = r2 + "_R3.fq.gz";
    else if (file_exists(r2 + "_R3.fastq.gz")) file_r3 = r2 + "_R3.fastq.gz";
    else { cerr << "Error: Cannot find R3 file (.fq.gz or .fastq.gz) for " << r2 << endl; exit(1); }

    cmd = "zcat " + file_r1;
    FILE * red1 = popen(cmd.c_str(), "r");
    if(!red1) { cerr << "Error opening R1: " << cmd << endl; exit(1); }

    cmd = "gzip - > " + r2 + "_R1_combined.fq.gz";
    FILE * outfile1 = popen(cmd.c_str(), "w");

    cmd = "zcat " + file_r2;
    FILE * red2 = popen(cmd.c_str(), "r");
    if(!red2) { cerr << "Error opening R2: " << cmd << endl; exit(1); }

    cmd = "zcat " + file_r3;
    FILE * red3 = popen(cmd.c_str(), "r");
    if(!red3) { cerr << "Error opening R3: " << cmd << endl; exit(1); }

    cmd = "gzip - > " + r2 + "_R3_combined.fq.gz";
    FILE * outfile2 = popen(cmd.c_str(), "w");

    char buffer[4096];
    fqline in_line1;
    fqline in_line2;
    fqline in_line3;
    string line1_raw;
    string readname_base;

    while(fgets(buffer, sizeof(buffer), red1)){
        ++total;
        
        line1_raw = buffer; 
        if(!line1_raw.empty() && line1_raw.back() == '\n') line1_raw.pop_back();
        if(!line1_raw.empty() && line1_raw.back() == '\r') line1_raw.pop_back();
        
        in_line1.read_part_record(red1, line1_raw);
        in_line2.read_full_record(red2);
        in_line3.read_full_record(red3);

        if(mode == "atac"){
            if(in_line2.seq.length() > 16) in_line2.seq.resize(16);
            if(in_line2.qual.length() > 16) in_line2.qual.resize(16);
            
            readname_base = clean_readname(in_line2.readname);

            in_line2.readname = readname_base + ":" + in_line1.seq + ":" + in_line1.qual;
            in_line2.write_record(outfile1); 
            
            in_line2.readname = readname_base + ":" + in_line3.seq + ":" + in_line3.qual;
            in_line2.write_record(outfile2);
            ++pass;

        } else if(mode == "arc"){
            if(in_line2.seq.length() < 24) continue;

            in_line2.seq = in_line2.seq.substr(8, 16);
            in_line2.qual = in_line2.qual.substr(8, 16);

            readname_base = clean_readname(in_line2.readname);

            in_line2.readname = readname_base + ":" + in_line1.seq + ":" + in_line1.qual;
            in_line2.write_record(outfile1);
            
            in_line2.readname = readname_base + ":" + in_line3.seq + ":" + in_line3.qual;
            in_line2.write_record(outfile2);
            ++pass;

        } else if(mode == "rna"){
            if(in_line1.seq.length() < 28) continue;

            string barcode = in_line1.seq.substr(0, 16);
            string umi = in_line1.seq.substr(16, 12);
            
            in_line1.seq = barcode;
            in_line1.qual = in_line1.qual.substr(0, 16);

            readname_base = clean_readname(in_line1.readname);

            in_line1.readname = readname_base + ":" + umi + ":" + in_line3.seq + ":" + in_line3.qual;
            in_line1.write_record(outfile2);
            ++pass;
        }
    }

    pclose(red1);
    pclose(red2);
    pclose(red3);
    pclose(outfile1);
    pclose(outfile2);

    double ratio = 0.0;
    if (total > 0) {
        ratio = ((double)pass / total) * 100.0;
    }

    cout << "==================================================\n(10X) Barcode Locator Report: " << r2 << endl;
    cout << "barcodes modes:\t\t" << mode << endl;
    cout << "# total raw reads:\t\t" << total << endl << "# of full barcoded reads:\t" << pass << endl;
    cout << "% of full length barcode reads:\t" << ratio << "%\n==================================================" << endl << endl;
    return;
}

void convert_hic2::run(string prefix){
    long long total = 0;
    long long pass = 0;
    string cmd = "cat " + prefix; 
    FILE * inbam = popen(cmd.c_str(), "r");
    if(!inbam){ cerr << "Error opening SAM: " << cmd << endl; exit(1); }

    string out_name = prefix.substr(0, prefix.length()-4) + "_cov.fq.gz";
    cmd = "gzip - > " + out_name;
    FILE * fout = popen(cmd.c_str(), "w");

    samline align_line;
    fqline fastq_line;
    char buffer[4096];
    
    while(fgets(buffer, sizeof(buffer), inbam)){
        if(buffer[0] == '@') continue;
        
        string line_str(buffer); 
        ++total;
        align_line.init(line_str);
        
        if(align_line.chr == "*") continue;

        vector<string> tmp = cxstring::split(align_line.readname, ":");
        if (tmp.size() < 3) continue; 

        int seq_idx = -1;
        for (size_t i = 1; i < tmp.size(); ++i) {
            if (tmp[i].length() > 20) {
                bool is_dna = true;
                for (char c : tmp[i]) {
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                        is_dna = false; break;
                    }
                }
                if (is_dna) {
                    seq_idx = i;
                    break;
                }
            }
        }

        if (seq_idx == -1) {
            if (tmp.size() > 7) seq_idx = 7;
            else continue; 
        }

        fastq_line.seq = tmp[seq_idx];
        size_t seq_len = fastq_line.seq.length();
        
        size_t expected_min_len = 2 * seq_len + 2; 
        if (align_line.readname.length() >= expected_min_len) {
            fastq_line.qual = align_line.readname.substr(align_line.readname.length() - seq_len);
            string base_name = align_line.readname.substr(0, align_line.readname.length() - expected_min_len);
            fastq_line.readname = "@" + align_line.chr + ":" + base_name;
        } else {
            fastq_line.qual = fastq_line.seq; 
            fastq_line.readname = "@" + align_line.chr + ":" + tmp[0];
        }
        
        fastq_line.mark = "+";
        fastq_line.write_record(fout);
        ++pass;
    }
    
    pclose(inbam);
    pclose(fout);
    
    string base_name = prefix.substr(0, prefix.length()-4);
    cout << "==================================================\n" << total << " reads processed in " << base_name << endl;
    cout << pass << " mapped reads in " << base_name << "\n==================================================" << endl;
    return;
}

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        cerr << "Usage: hictools [mode] ..." << endl;
        return 1;
    }

    string mod(argv[1]);
    
    if(mod == "combine_hic"){
        if(argc < 4){
            cerr<<"Usage: hictools combine_hic [atac/arc/rna] [prefix]"<<endl;
            return 1;
        }
        string mode = argv[2];
        if((mode != "atac") && (mode != "arc") && (mode != "rna")){
            cerr << "Invalid mode: " << mode << endl;
            return 1;
        }
        combine_hic::run(mode, argv[3]);
        return 0;
    }

    if(mod == "convert_hic2"){
        if(argc < 3){
            cerr << "hictools convert_hic2 [sample_BC.sam]." << endl;
            return 1;
        }
        convert_hic2::run(argv[2]);
        return 0;
    } 
    
    cerr << "Unknown module: " << mod << endl;
    return 1;
}