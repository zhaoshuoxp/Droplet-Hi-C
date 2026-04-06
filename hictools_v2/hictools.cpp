//
// hictools.cpp
// Ultimate Performance Version (kseq.h + string_view + OpenMP + pigz)
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string_view>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <omp.h>


#include "kseq.h"
KSEQ_INIT(int, read)

using namespace std;

bool file_exists(const string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}


FILE* get_reader(const string& path, const string& threads) {
    string cmd = "pigz -p " + threads + " -dc " + path + " 2>/dev/null || zcat " + path;
    FILE* fp = popen(cmd.c_str(), "r");
    if (!fp) { cerr << "Error opening reader for " << path << endl; exit(1); }
    return fp;
}

FILE* get_writer(const string& path, const string& threads) {
    string cmd = "pigz -p " + threads + " -c > " + path + " 2>/dev/null || gzip -c > " + path;
    FILE* fp = popen(cmd.c_str(), "w");
    if (!fp) { cerr << "Error opening writer for " << path << endl; exit(1); }
    return fp;
}

string_view clean_readname(string_view raw_name) {
    auto space_pos = raw_name.find_first_of(" \t");
    string_view base = (space_pos != string_view::npos) ? raw_name.substr(0, space_pos) : raw_name;
    
    auto colon_r_pos = base.rfind(":r");
    if (colon_r_pos != string_view::npos && colon_r_pos > 0) {
        string_view suffix = base.substr(colon_r_pos + 2);
        bool is_barcode = true;
        for (char c : suffix) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                is_barcode = false; 
                break;
            }
        }
        if (is_barcode) {
            return base.substr(0, colon_r_pos);
        }
    }
    return base;
}

void run_combine_hic(const string& mode, const string& r2_prefix, const string& threads) {
    string file_r1, file_r2, file_r3;
    const char* exts[] = {".fq.gz", ".fastq.gz"};
    
    for (const char* ext : exts) {
        if (file_exists(r2_prefix + "_R1" + ext)) file_r1 = r2_prefix + "_R1" + ext;
        if (file_exists(r2_prefix + "_R2" + ext)) file_r2 = r2_prefix + "_R2" + ext;
        if (file_exists(r2_prefix + "_R3" + ext)) file_r3 = r2_prefix + "_R3" + ext;
    }

    if (file_r1.empty() || file_r2.empty() || file_r3.empty()) {
        cerr << "Error: Cannot find matching _R1, _R2, or _R3 (.fq.gz/.fastq.gz) for prefix " << r2_prefix << endl;
        exit(1);
    }

    FILE* red1 = get_reader(file_r1, threads);
    FILE* red2 = get_reader(file_r2, threads);
    FILE* red3 = get_reader(file_r3, threads);

    FILE* out1 = get_writer(r2_prefix + "_R1_combined.fq.gz", threads);
    FILE* out2 = get_writer(r2_prefix + "_R3_combined.fq.gz", threads);

    kseq_t *seq1 = kseq_init(fileno(red1));
    kseq_t *seq2 = kseq_init(fileno(red2));
    kseq_t *seq3 = kseq_init(fileno(red3));

    unsigned long long total = 0, pass = 0;
    
    string out1_buf, out2_buf;
    out1_buf.reserve(4 * 1024 * 1024); // 4MB buffer
    out2_buf.reserve(4 * 1024 * 1024);

    while (kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0 && kseq_read(seq3) >= 0) {
        total++;
        
        string_view n1(seq1->name.s, seq1->name.l);
        string_view n2(seq2->name.s, seq2->name.l);
        
        string_view s1(seq1->seq.s, seq1->seq.l);
        string_view q1(seq1->qual.s, seq1->qual.l);
        string_view s2(seq2->seq.s, seq2->seq.l);
        string_view q2(seq2->qual.s, seq2->qual.l);
        string_view s3(seq3->seq.s, seq3->seq.l);
        string_view q3(seq3->qual.s, seq3->qual.l);

        if (mode == "atac") {
            string_view seq = s2.length() > 16 ? s2.substr(0, 16) : s2;
            string_view qual = q2.length() > 16 ? q2.substr(0, 16) : q2;
            string_view base_name = clean_readname(n2);

            out1_buf.append("@").append(base_name).append(":").append(s1).append(":").append(q1)
                    .append("\n").append(seq).append("\n+\n").append(qual).append("\n");
            
            out2_buf.append("@").append(base_name).append(":").append(s3).append(":").append(q3)
                    .append("\n").append(seq).append("\n+\n").append(qual).append("\n");
            pass++;
        } 
        else if (mode == "arc") {
            if (s2.length() < 24) continue;
            string_view seq = s2.substr(8, 16);
            string_view qual = q2.substr(8, 16);
            string_view base_name = clean_readname(n2);

            out1_buf.append("@").append(base_name).append(":").append(s1).append(":").append(q1)
                    .append("\n").append(seq).append("\n+\n").append(qual).append("\n");
            
            out2_buf.append("@").append(base_name).append(":").append(s3).append(":").append(q3)
                    .append("\n").append(seq).append("\n+\n").append(qual).append("\n");
            pass++;
        } 
        else if (mode == "rna") {
            if (s1.length() < 28) continue;
            string_view barcode = s1.substr(0, 16);
            string_view umi = s1.substr(16, 12);
            string_view qual = q1.substr(0, 16);
            string_view base_name = clean_readname(n1);

            out2_buf.append("@").append(base_name).append(":").append(umi).append(":").append(s3).append(":").append(q3)
                    .append("\n").append(barcode).append("\n+\n").append(qual).append("\n");
            pass++;
        }

        if (out1_buf.size() >= 2 * 1024 * 1024) {
            fwrite(out1_buf.data(), 1, out1_buf.size(), out1);
            out1_buf.clear();
        }
        if (out2_buf.size() >= 2 * 1024 * 1024) {
            fwrite(out2_buf.data(), 1, out2_buf.size(), out2);
            out2_buf.clear();
        }
    }

    if (!out1_buf.empty()) fwrite(out1_buf.data(), 1, out1_buf.size(), out1);
    if (!out2_buf.empty()) fwrite(out2_buf.data(), 1, out2_buf.size(), out2);

    kseq_destroy(seq1); pclose(red1);
    kseq_destroy(seq2); pclose(red2);
    kseq_destroy(seq3); pclose(red3);
    pclose(out1); pclose(out2);

    double ratio = (total > 0) ? ((double)pass / total) * 100.0 : 0.0;
    cout << "==================================================\n(10X) Barcode Locator Report: " << r2_prefix << endl;
    cout << "barcodes modes:\t\t" << mode << endl;
    cout << "# total raw reads:\t\t" << total << "\n# of full barcoded reads:\t" << pass << endl;
    cout << "% of full length barcode reads:\t" << ratio << "%\n==================================================\n\n";
}


void run_convert_hic2(const string& prefix, const string& threads) {

    int num_threads = std::stoi(threads);
    omp_set_num_threads(num_threads);

    FILE* inbam = fopen(prefix.c_str(), "r");
    if (!inbam) { cerr << "Error opening SAM: " << prefix << endl; exit(1); }

    string out_name = (prefix.length() > 4 && prefix.substr(prefix.length() - 4) == ".sam") 
                      ? prefix.substr(0, prefix.length() - 4) + "_cov.fq.gz" 
                      : prefix + "_cov.fq.gz";
                      
    FILE* fout = get_writer(out_name, threads);

    const int BATCH_SIZE = 100000; 
    vector<string> lines(BATCH_SIZE);
    vector<string> out_lines(BATCH_SIZE);
    
    char buffer[4096];
    unsigned long long total = 0, pass = 0;
    int current_batch_size = 0;

    auto process_batch = [&]() {
        #pragma omp parallel for
        for (int i = 0; i < current_batch_size; ++i) {
            out_lines[i].clear();
            string_view line = lines[i];
            
            if (line.empty() || line[0] == '@') continue;
            
            size_t t1 = line.find('\t'); if (t1 == string_view::npos) continue;
            size_t t2 = line.find('\t', t1 + 1); if (t2 == string_view::npos) continue;
            size_t t3 = line.find('\t', t2 + 1); if (t3 == string_view::npos) continue;

            string_view qname = line.substr(0, t1);
            string_view chr = line.substr(t2 + 1, t3 - t2 - 1);
            if (chr == "*") continue;

            vector<string_view> tmp;
            size_t start = 0, end;
            while ((end = qname.find(':', start)) != string_view::npos) {
                tmp.push_back(qname.substr(start, end - start));
                start = end + 1;
            }
            tmp.push_back(qname.substr(start));

            if (tmp.size() < 3) continue;

            int seq_idx = -1;
            for (size_t j = 1; j < tmp.size(); ++j) {
                if (tmp[j].length() > 20) {
                    bool is_dna = true;
                    for (char c : tmp[j]) {
                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') { is_dna = false; break; }
                    }
                    if (is_dna) { seq_idx = j; break; }
                }
            }

            if (seq_idx == -1) {
                if (tmp.size() > 7) seq_idx = 7;
                else continue;
            }

            string_view seq = tmp[seq_idx];
            size_t seq_len = seq.length();
            size_t expected_min_len = 2 * seq_len + 2;

            string final_name, final_qual;
            if (qname.length() >= expected_min_len) {
                final_qual = string(qname.substr(qname.length() - seq_len));
                string_view base_name = qname.substr(0, qname.length() - expected_min_len);
                final_name = "@" + string(chr) + ":" + string(base_name);
            } else {
                final_qual = string(seq);
                final_name = "@" + string(chr) + ":" + string(tmp[0]);
            }

            out_lines[i] = final_name + "\n" + string(seq) + "\n+\n" + final_qual + "\n";
        }

        string batch_out;
        batch_out.reserve(BATCH_SIZE * 200);
        for (int i = 0; i < current_batch_size; ++i) {
            if (!out_lines[i].empty()) {
                batch_out += out_lines[i];
                pass++;
            }
        }
        if (!batch_out.empty()) fwrite(batch_out.data(), 1, batch_out.size(), fout);
    };

    while (fgets(buffer, sizeof(buffer), inbam)) {
        size_t len = strlen(buffer);
        if (len > 0 && buffer[len-1] == '\n') buffer[len-1] = '\0';
        if (len > 1 && buffer[len-2] == '\r') buffer[len-2] = '\0';
        
        lines[current_batch_size++] = buffer;
        if (buffer[0] != '@') total++;
        
        if (current_batch_size == BATCH_SIZE) {
            process_batch();
            current_batch_size = 0;
        }
    }
    if (current_batch_size > 0) process_batch(); 

    fclose(inbam); pclose(fout);

    string base_name = (prefix.length() > 4 && prefix.substr(prefix.length() - 4) == ".sam") ? prefix.substr(0, prefix.length() - 4) : prefix;
    cout << "==================================================\n" << total << " reads processed in " << base_name << endl;
    cout << pass << " mapped reads in " << base_name << "\n==================================================\n";
}

int main(int argc, char *argv[]) {
    if (argc < 2) { cerr << "Usage: hictools [mode] ..." << endl; return 1; }

    string mod(argv[1]);
    
    if (mod == "combine_hic") {
        if (argc < 4) { cerr << "Usage: hictools combine_hic [atac/arc/rna] [prefix] [threads]" << endl; return 1; }
        string mode = argv[2];
        if (mode != "atac" && mode != "arc" && mode != "rna") { cerr << "Invalid mode: " << mode << endl; return 1; }
        string threads = (argc > 4) ? argv[4] : "8";
        run_combine_hic(mode, argv[3], threads);
    } 
    else if (mod == "convert_hic2") {
        if (argc < 3) { cerr << "Usage: hictools convert_hic2 [sample_BC.sam] [threads]" << endl; return 1; }
        string threads = (argc > 3) ? argv[3] : "8";
        run_convert_hic2(argv[2], threads);
    } 
    else { cerr << "Unknown module: " << mod << endl; return 1; }
    
    return 0;
}