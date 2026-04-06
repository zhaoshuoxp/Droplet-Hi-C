//
// hictools.cpp
// All-in-One Ultimate Version
// Includes: combine_hic, convert_hic2, fast_map, end_to_end
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
#include <unordered_map>

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
                is_barcode = false; break;
            }
        }
        if (is_barcode) return base.substr(0, colon_r_pos);
    }
    return base;
}

inline uint64_t encode_seq(string_view seq) {
    uint64_t val = 0;
    for (char c : seq) {
        uint64_t b = 4;
        if (c == 'A' || c == 'a') b = 0;
        else if (c == 'C' || c == 'c') b = 1;
        else if (c == 'G' || c == 'g') b = 2;
        else if (c == 'T' || c == 't') b = 3;
        val = (val << 3) | b;
    }
    return val;
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
        cerr << "Error: Cannot find raw files for " << r2_prefix << endl; exit(1);
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
    out1_buf.reserve(4 * 1024 * 1024);
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

        if (mode == "atac" || mode == "arc") {
            string_view seq, qual;
            if (mode == "atac") {
                seq = s2.length() > 16 ? s2.substr(0, 16) : s2;
                qual = q2.length() > 16 ? q2.substr(0, 16) : q2;
            } else {
                if (s2.length() < 24) continue;
                seq = s2.substr(8, 16);
                qual = q2.substr(8, 16);
            }
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
            fwrite(out1_buf.data(), 1, out1_buf.size(), out1); out1_buf.clear();
        }
        if (out2_buf.size() >= 2 * 1024 * 1024) {
            fwrite(out2_buf.data(), 1, out2_buf.size(), out2); out2_buf.clear();
        }
    }

    if (!out1_buf.empty()) fwrite(out1_buf.data(), 1, out1_buf.size(), out1);
    if (!out2_buf.empty()) fwrite(out2_buf.data(), 1, out2_buf.size(), out2);

    kseq_destroy(seq1); pclose(red1); kseq_destroy(seq2); pclose(red2); kseq_destroy(seq3); pclose(red3);
    pclose(out1); pclose(out2);

    cout << "==================================================\n(10X) Combine Report: " << r2_prefix << "\n# reads: " << total << " | pass: " << pass << "\n==================================================\n";
}

void run_convert_hic2(const string& prefix, const string& threads) {
    omp_set_num_threads(std::stoi(threads));

    FILE* inbam = fopen(prefix.c_str(), "r");
    if (!inbam) { cerr << "Error opening SAM: " << prefix << endl; exit(1); }

    string out_name = (prefix.length() > 4 && prefix.substr(prefix.length() - 4) == ".sam") 
                      ? prefix.substr(0, prefix.length() - 4) + "_cov.fq.gz" : prefix + "_cov.fq.gz";
    FILE* fout = get_writer(out_name, threads);

    const int BATCH_SIZE = 100000; 
    vector<string> lines(BATCH_SIZE);
    vector<string> out_lines(BATCH_SIZE);
    char buffer[4096];
    unsigned long long total = 0, pass = 0;
    int current_batch = 0;

    auto process_batch = [&]() {
        #pragma omp parallel for
        for (int i = 0; i < current_batch; ++i) {
            out_lines[i].clear();
            string_view line = lines[i];
            if (line.empty() || line[0] == '@') continue;
            
            size_t t1 = line.find('\t'), t2 = line.find('\t', t1 + 1), t3 = line.find('\t', t2 + 1);
            if (t1 == string_view::npos || t2 == string_view::npos || t3 == string_view::npos) continue;

            string_view qname = line.substr(0, t1);
            string_view chr = line.substr(t2 + 1, t3 - t2 - 1);
            if (chr == "*") continue;

            vector<string_view> tmp; size_t start = 0, end;
            while ((end = qname.find(':', start)) != string_view::npos) {
                tmp.push_back(qname.substr(start, end - start)); start = end + 1;
            }
            tmp.push_back(qname.substr(start));

            if (tmp.size() < 3) continue;
            int seq_idx = -1;
            for (size_t j = 1; j < tmp.size(); ++j) {
                if (tmp[j].length() > 20 && tmp[j].find_first_not_of("ACGTNacgtn") == string_view::npos) { seq_idx = j; break; }
            }
            if (seq_idx == -1) seq_idx = (tmp.size() > 7) ? 7 : -1;
            if (seq_idx == -1) continue;

            string_view seq = tmp[seq_idx];
            size_t min_len = 2 * seq.length() + 2;

            string fname, fqual;
            if (qname.length() >= min_len) {
                fqual = string(qname.substr(qname.length() - seq.length()));
                fname = "@" + string(chr) + ":" + string(qname.substr(0, qname.length() - min_len));
            } else {
                fqual = string(seq); fname = "@" + string(chr) + ":" + string(tmp[0]);
            }
            out_lines[i] = fname + "\n" + string(seq) + "\n+\n" + fqual + "\n";
        }

        string batch_out; batch_out.reserve(BATCH_SIZE * 200);
        for (int i = 0; i < current_batch; ++i) {
            if (!out_lines[i].empty()) { batch_out += out_lines[i]; pass++; }
        }
        if (!batch_out.empty()) fwrite(batch_out.data(), 1, batch_out.size(), fout);
    };

    while (fgets(buffer, sizeof(buffer), inbam)) {
        size_t len = strlen(buffer);
        if (len > 0 && buffer[len-1] == '\n') buffer[len-1] = '\0';
        if (len > 1 && buffer[len-2] == '\r') buffer[len-2] = '\0';
        lines[current_batch++] = buffer;
        if (buffer[0] != '@') total++;
        if (current_batch == BATCH_SIZE) { process_batch(); current_batch = 0; }
    }
    if (current_batch > 0) process_batch(); 

    fclose(inbam); pclose(fout);
    cout << "==================================================\nConvert Report\nReads: " << total << " | Mapped: " << pass << "\n==================================================\n";
}

void run_fast_map(const string& whitelist_path, const string& input_fq, const string& output_fq, const string& threads) {
    omp_set_num_threads(std::stoi(threads));
    cout << "Building 1-mismatch dictionary from: " << whitelist_path << " ..." << endl;
    
    std::unordered_map<uint64_t, uint32_t> dict; dict.reserve(50000000);
    vector<string> wl_strs; ifstream wl_file(whitelist_path);
    if (!wl_file) { cerr << "Cannot open whitelist!" << endl; exit(1); }
    
    string line; uint32_t idx = 0;
    while (getline(wl_file, line)) {
        if (line.empty() || line[0] == '>') continue; 
        if (line.back() == '\r') line.pop_back();     
        if (line.length() != 16) continue;
        
        wl_strs.push_back(line);
        auto insert_dict = [&](uint64_t key, uint32_t val) {
            auto it = dict.find(key);
            if (it == dict.end()) dict[key] = val; else if (it->second != val) it->second = 0xFFFFFFFF; 
        };
        insert_dict(encode_seq(line), idx);
        for (int i = 0; i < 16; ++i) {
            char orig = line[i]; const char bases[] = {'A', 'C', 'G', 'T', 'N'};
            for (char b : bases) {
                if (b == orig) continue; string mut = line; mut[i] = b; insert_dict(encode_seq(mut), idx);
            }
        }
        idx++;
    }

    FILE* in_fp = get_reader(input_fq, threads); FILE* out_fp = get_writer(output_fq, threads);
    kseq_t* seq = kseq_init(fileno(in_fp));
    struct FastqItem { string name, seq, qual; };
    const int BATCH_SIZE = 100000;
    vector<FastqItem> batch(BATCH_SIZE); vector<string> out_lines(BATCH_SIZE);
    unsigned long long total = 0, pass = 0; int current_batch = 0;

    auto process_batch = [&]() {
        #pragma omp parallel for
        for (int i = 0; i < current_batch; ++i) {
            out_lines[i].clear();
            if (batch[i].seq.length() != 16) continue;
            auto it = dict.find(encode_seq(batch[i].seq));
            if (it != dict.end() && it->second != 0xFFFFFFFF) {
                const string& correct_bc = wl_strs[it->second];
                string_view qname = batch[i].name;
                vector<string_view> tmp; size_t start = 0, end;
                while ((end = qname.find(':', start)) != string_view::npos) {
                    tmp.push_back(qname.substr(start, end - start)); start = end + 1;
                }
                tmp.push_back(qname.substr(start));

                if (tmp.size() < 3) continue;
                int seq_idx = -1;
                for (size_t j = 1; j < tmp.size(); ++j) {
                    if (tmp[j].length() > 20 && tmp[j].find_first_not_of("ACGTNacgtn") == string_view::npos) { seq_idx = j; break; }
                }
                if (seq_idx == -1) seq_idx = (tmp.size() > 7) ? 7 : -1;
                if (seq_idx == -1) continue;

                string_view bio_seq = tmp[seq_idx]; size_t min_len = 2 * bio_seq.length() + 2;
                string fname, fqual;
                if (qname.length() >= min_len) {
                    fqual = string(qname.substr(qname.length() - bio_seq.length()));
                    fname = "@" + correct_bc + ":" + string(qname.substr(0, qname.length() - min_len));
                } else {
                    fqual = string(bio_seq); fname = "@" + correct_bc + ":" + string(tmp[0]);
                }
                out_lines[i] = fname + "\n" + string(bio_seq) + "\n+\n" + fqual + "\n";
            }
        }
        string batch_out; batch_out.reserve(BATCH_SIZE * 200);
        for (int i = 0; i < current_batch; ++i) {
            if (!out_lines[i].empty()) { batch_out += out_lines[i]; pass++; }
        }
        if (!batch_out.empty()) fwrite(batch_out.data(), 1, batch_out.size(), out_fp);
    };

    while (kseq_read(seq) >= 0) {
        batch[current_batch].name = string(seq->name.s, seq->name.l);
        batch[current_batch].seq = string(seq->seq.s, seq->seq.l);
        batch[current_batch].qual = string(seq->qual.s, seq->qual.l);
        current_batch++; total++;
        if (current_batch == BATCH_SIZE) { process_batch(); current_batch = 0; }
    }
    if (current_batch > 0) process_batch();
    kseq_destroy(seq); pclose(in_fp); pclose(out_fp);
    cout << "==================================================\nFast Map Report\nReads: " << total << " | Pass: " << pass << "\n==================================================\n";
}


inline uint64_t encode_rc_seq(string_view seq) {
    uint64_t val = 0;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        uint64_t b = 4; // N
        if (*it == 'A' || *it == 'a') b = 3;      // T
        else if (*it == 'C' || *it == 'c') b = 2; // G
        else if (*it == 'G' || *it == 'g') b = 1; // C
        else if (*it == 'T' || *it == 't') b = 0; // A
        val = (val << 3) | b;
    }
    return val;
}

void run_end_to_end(const string& mode, const string& r2_prefix, const string& whitelist_path, const string& threads) {
    int num_threads = std::stoi(threads);
    omp_set_num_threads(num_threads);

    cout << "Mode: " << mode << " | Loading Whitelist: " << whitelist_path << "..." << endl;
    std::unordered_map<uint64_t, uint32_t> dict;
    dict.reserve(50000000); 
    vector<string> wl_strs;
    ifstream wl_file(whitelist_path);
    if (!wl_file) { cerr << "Cannot open whitelist!" << endl; exit(1); }
    
    string line;
    uint32_t idx = 0;
    while (getline(wl_file, line)) {
        if (line.empty() || line[0] == '>') continue; 
        if (!line.empty() && line.back() == '\r') line.pop_back();     
        if (line.length() != 16) continue;
        
        wl_strs.push_back(line);
        uint64_t exact_key = encode_seq(line);
        auto insert_dict = [&](uint64_t key, uint32_t val) {
            auto it = dict.find(key);
            if (it == dict.end()) dict[key] = val;
            else if (it->second != val) dict[key] = 0xFFFFFFFF; 
        };
        insert_dict(exact_key, idx);
        for (int i = 0; i < 16; ++i) {
            char orig = line[i];
            const char bases[] = {'A', 'C', 'G', 'T', 'N'};
            for (char b : bases) {
                if (b == orig) continue;
                string mut = line; mut[i] = b;
                insert_dict(encode_seq(mut), idx);
            }
        }
        idx++;
    }
    cout << "Dictionary built with " << wl_strs.size() << " barcodes (1-mismatch enabled)." << endl;

    string f1, f2, f3;
    const char* exts[] = {".fq.gz", ".fastq.gz"};
    for (const char* e : exts) {
        if (file_exists(r2_prefix + "_R1" + e)) f1 = r2_prefix + "_R1" + e;
        if (file_exists(r2_prefix + "_R2" + e)) f2 = r2_prefix + "_R2" + e;
        if (file_exists(r2_prefix + "_R3" + e)) f3 = r2_prefix + "_R3" + e;
    }
    if (f1.empty() || f2.empty() || f3.empty()) { cerr << "Input files not found for prefix " << r2_prefix << endl; exit(1); }

    FILE *red1 = get_reader(f1, threads), *red2 = get_reader(f2, threads), *red3 = get_reader(f3, threads);
    FILE *out1 = get_writer(r2_prefix + "_R1_BC_cov.fq.gz", threads);
    FILE *out2 = get_writer(r2_prefix + "_R3_BC_cov.fq.gz", threads);

    kseq_t *s1 = kseq_init(fileno(red1)), *s2 = kseq_init(fileno(red2)), *s3 = kseq_init(fileno(red3));

    const int BATCH_SIZE = 100000;
    struct ReadData { 
        string name, seq1, qual1, seq2, seq3, qual3;
        string out1, out2;
        bool pass = false;
    };
    vector<ReadData> batch(BATCH_SIZE);
    unsigned long long total = 0, total_pass = 0;
    bool done = false;

    cout << "Processing reads with " << num_threads << " threads..." << endl;

    while (!done) {
        int count = 0;
        while (count < BATCH_SIZE && kseq_read(s1) >= 0 && kseq_read(s2) >= 0 && kseq_read(s3) >= 0) {
            ReadData &r = batch[count];
            r.name = s2->name.s;
            r.seq1 = s1->seq.s; r.qual1 = s1->qual.s;
            r.seq2 = s2->seq.s;
            r.seq3 = s3->seq.s; r.qual3 = s3->qual.s;
            r.pass = false;
            r.out1.clear(); r.out2.clear();
            count++;
        }
        if (count == 0) break;
        if (count < BATCH_SIZE) done = true;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < count; ++i) {
            ReadData &r = batch[i];
            string_view raw_bc;
            if (mode == "atac") raw_bc = r.seq2.length() >= 16 ? string_view(r.seq2.data(), 16) : "";
            else if (mode == "arc") raw_bc = r.seq2.length() >= 24 ? string_view(r.seq2.data() + 8, 16) : "";
            else if (mode == "rna") raw_bc = r.seq1.length() >= 16 ? string_view(r.seq1.data(), 16) : "";

            if (!raw_bc.empty()) {
                // 💡 双向嗅探匹配策略
                // 由于原 Bowtie 脚本使用 --nofw，数据主要是反向互补的，所以优先比对 RC
                uint64_t key_rc = encode_rc_seq(raw_bc);
                auto it = dict.find(key_rc);
                
                // 如果 RC 没匹配上（或发生错配碰撞），则降级尝试正向比对
                if (it == dict.end() || it->second == 0xFFFFFFFF) {
                    uint64_t key_fw = encode_seq(raw_bc);
                    it = dict.find(key_fw);
                }

                if (it != dict.end() && it->second != 0xFFFFFFFF) {
                    const string &cbc = wl_strs[it->second];
                    string_view bname = clean_readname(r.name);
                    
                    if (mode == "rna") {
                        string_view umi(r.seq1.data() + 16, 12);
                        r.out2 = "@" + cbc + ":" + string(bname) + ":" + string(umi) + "\n" + r.seq3 + "\n+\n" + r.qual3 + "\n";
                    } else {
                        r.out1 = "@" + cbc + ":" + string(bname) + "\n" + r.seq1 + "\n+\n" + r.qual1 + "\n";
                        r.out2 = "@" + cbc + ":" + string(bname) + "\n" + r.seq3 + "\n+\n" + r.qual3 + "\n";
                    }
                    r.pass = true;
                }
            }
        }

        for (int i = 0; i < count; ++i) {
            total++;
            if (batch[i].pass) {
                if (!batch[i].out1.empty()) fwrite(batch[i].out1.data(), 1, batch[i].out1.size(), out1);
                if (!batch[i].out2.empty()) fwrite(batch[i].out2.data(), 1, batch[i].out2.size(), out2);
                total_pass++;
            }
        }
        
        if (total % 10000000 == 0) {
            cout << "Processed " << total / 1000000 << "M reads. Valid matches so far: " << total_pass << endl;
        }
    }

    kseq_destroy(s1); kseq_destroy(s2); kseq_destroy(s3);
    pclose(red1); pclose(red2); pclose(red3); pclose(out1); pclose(out2);
    
    double ratio = (total > 0) ? ((double)total_pass / total) * 100.0 : 0.0;
    cout << "==================================================\n(10X) End-to-End Processing Report\n";
    cout << "Mode:\t\t\t" << mode << "\n";
    cout << "Total raw reads:\t" << total << "\nValid barcoded reads:\t" << total_pass << "\n";
    cout << "Mapping Rate:\t\t" << ratio << "%\n==================================================\n";
}

int main(int argc, char *argv[]) {
    if (argc < 2) { 
        cerr << "Usage: hictools [module] [options]...\n"
             << "Modules:\n"
             << "  combine_hic   [mode] [prefix] [threads]\n"
             << "  convert_hic2  [sam_file] [threads]\n"
             << "  fast_map      [whitelist.fa] [input_fq] [output_fq] [threads]\n"
             << "  end_to_end    [mode] [prefix] [whitelist.fa] [threads]\n";
        return 1; 
    }

    string mod(argv[1]);
    
    if (mod == "combine_hic") {
        if (argc < 4) { cerr << "Usage: hictools combine_hic [atac/arc/rna] [prefix] [threads]" << endl; return 1; }
        run_combine_hic(argv[2], argv[3], (argc > 4) ? argv[4] : "8");
    } 
    else if (mod == "convert_hic2") {
        if (argc < 3) { cerr << "Usage: hictools convert_hic2 [sam_file] [threads]" << endl; return 1; }
        run_convert_hic2(argv[2], (argc > 3) ? argv[3] : "8");
    }
    else if (mod == "fast_map") {
        if (argc < 5) { cerr << "Usage: hictools fast_map [whitelist.fa] [input.fq.gz] [output.fq.gz] [threads]" << endl; return 1; }
        run_fast_map(argv[2], argv[3], argv[4], (argc > 5) ? argv[5] : "8");
    }
    else if (mod == "end_to_end") {
        if (argc < 5) { cerr << "Usage: hictools end_to_end [atac/arc/rna] [prefix] [whitelist.fa] [threads]" << endl; return 1; }
        run_end_to_end(argv[2], argv[3], argv[4], (argc > 5) ? argv[5] : "8");
    } 
    else { 
        cerr << "Unknown module: " << mod << endl; 
        return 1; 
    }
    
    return 0;
}