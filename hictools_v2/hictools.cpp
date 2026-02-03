//
// hictools.cpp
// Optimized Version
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "cxstring.hpp" 
// 注意：不要 include .cpp 文件，这是不规范的做法，编译时链接即可

using namespace std;

class combine_hic{
private:

public:
    static void run(string mode, string r2, string ty);
};

class convert_hic2{
private:

public:
    static void run(string prefix);
};

void combine_hic::run(string mode, string r2, string ty){
  
    // processing scHiC with 10X barcodes
    int total = 0;
    int pass = 0;
    string s1, s2, s3;
    string cmd;

    // 构造命令字符串
    if(ty == "gz"){
        s1 = "zcat ";
        s2 = "_R1.fq.gz";
    } else if(ty == "bz2"){
        s1 = "bzcat ";
        s2 = "_R1.fastq.bz2";
    }
    cmd = s1 + r2 + s2;
    FILE * red1 = popen(cmd.c_str(), "r");
    if(!red1) { cerr << "Error opening R1: " << cmd << endl; exit(1); }

    // R1 Output
    s1 = "gzip - > ";
    s2 = r2 + "_R1_combined.fq.gz";
    cmd = s1 + s2;
    FILE * outfile1 = popen(cmd.c_str(), "w");

    // Open R2
    if(ty == "gz"){
        s1 = "zcat ";
        s2 = "_R2.fq.gz";
    } else if(ty == "bz2"){
        s1 = "bzcat ";
        s2 = "_R2.fastq.bz2";
    }
    cmd = s1 + r2 + s2;
    FILE * red2 = popen(cmd.c_str(), "r");
    if(!red2) { cerr << "Error opening R2: " << cmd << endl; exit(1); }

    // Open R3
    if(ty == "gz"){
        s1 = "zcat ";
        s2 = "_R3.fq.gz";
    } else if(ty == "bz2"){
        s1 = "bzcat ";
        s2 = "_R3.fastq.bz2";
    }
    cmd = s1 + r2 + s2;
    FILE * red3 = popen(cmd.c_str(), "r");
    if(!red3) { cerr << "Error opening R3: " << cmd << endl; exit(1); }

    // R3 Output
    s1 = "gzip - > ";
    s2 = r2 + "_R3_combined.fq.gz";
    cmd = s1 + s2;
    FILE * outfile2 = popen(cmd.c_str(), "w");

    // 变量预声明，减少循环内内存分配
    char buffer[4096]; // 增大 buffer
    fqline in_line1;
    fqline in_line2;
    fqline in_line3;
    string line1_raw;
    string readname_base;
    size_t space_pos;

    // 主循环
    // 逻辑：手动读取 R1 的第一行(readname)，然后用 fqline 读取剩余部分
    while(fgets(buffer, sizeof(buffer), red1)){
        ++total;
        
        // 处理 R1 Readname
        line1_raw = buffer; 
        if(!line1_raw.empty() && line1_raw.back() == '\n') line1_raw.pop_back();
        if(!line1_raw.empty() && line1_raw.back() == '\r') line1_raw.pop_back();
        
        in_line1.read_part_record(red1, line1_raw);
        in_line2.read_full_record(red2);
        in_line3.read_full_record(red3);

        if(mode == "atac"){
            // ATAC: R2 前 16bp 是 Barcode
            // 优化：直接检查长度
            // if(in_line2.seq.length() != 16) continue; // 原代码逻辑似乎不严格跳过，只截取

            // 优化：使用 resize 避免 substr 产生新对象
            if(in_line2.seq.length() > 16) in_line2.seq.resize(16);
            if(in_line2.qual.length() > 16) in_line2.qual.resize(16);
            
            // 优化：移除 stringstream，使用 find 截取 readname
            space_pos = in_line2.readname.find_first_of(" \t");
            if(space_pos != string::npos) {
                readname_base = in_line2.readname.substr(0, space_pos);
            } else {
                readname_base = in_line2.readname;
            }

            // 构造新的 readname: R2ReadName:R1Seq:R1Qual
            in_line2.readname = readname_base + ":" + in_line1.seq + ":" + in_line1.qual;
            
            in_line2.write_record(outfile1); // R1 combined (实际上内容是 R2 的 Barcode + R1 info)
            
            // Output 2: Readname 相同，Seq/Qual 来自 R3
            in_line2.readname = readname_base + ":" + in_line3.seq + ":" + in_line3.qual;
            // 注意：这里原代码写的是 in_line2.write_record(outfile2)，但此时 in_line2.seq 还是 Barcode
            // 原逻辑是把 Barcode 当作 Seq 写入 R1_combined 和 R3_combined，
            // 并在 Readname 中携带原始信息。这逻辑有点奇怪但保持原样。
            in_line2.write_record(outfile2);
            
            ++pass;

        } else if(mode == "arc"){
            // ARC: 8-24bp 是 Barcode (共16bp)
            // 检查长度是否足够
            if(in_line2.seq.length() < 24) continue;

            // 优化：截取中间部分，这必须创建新字符串或移动
            in_line2.seq = in_line2.seq.substr(8, 16); // 取 8 开始的 16bp
            in_line2.qual = in_line2.qual.substr(8, 16);

            space_pos = in_line2.readname.find_first_of(" \t");
            if(space_pos != string::npos) {
                readname_base = in_line2.readname.substr(0, space_pos);
            } else {
                readname_base = in_line2.readname;
            }

            in_line2.readname = readname_base + ":" + in_line1.seq + ":" + in_line1.qual;
            in_line2.write_record(outfile1);
            
            in_line2.readname = readname_base + ":" + in_line3.seq + ":" + in_line3.qual;
            in_line2.write_record(outfile2);
            
            ++pass;

        } else if(mode == "rna"){
            // RNA: R1 前 16bp Barcode, 16-28bp UMI (12bp)
            if(in_line1.seq.length() < 28) continue; // 安全检查

            string barcode = in_line1.seq.substr(0, 16);
            string umi = in_line1.seq.substr(16, 12);
            
            // 原逻辑修改了 in_line1 的 seq/qual 为 Barcode
            in_line1.seq = barcode;
            in_line1.qual = in_line1.qual.substr(0, 16);

            space_pos = in_line1.readname.find_first_of(" \t");
            if(space_pos != string::npos) {
                readname_base = in_line1.readname.substr(0, space_pos);
            } else {
                readname_base = in_line1.readname;
            }

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

    // 修复 Bug：整数溢出修复，并增加被除数为0的检查
    double ratio = 0.0;
    if (total > 0) {
        ratio = ((double)pass / total) * 100.0;
    }

    cout << "==================================================\n(10X) Barcode Locator Report: " << r2 << endl;
    cout << "barcodes modes:\t\t" << mode << endl;
    cout << "# total raw reads:\t\t" << total << endl << "# of full barcoded reads:\t" << pass << endl;
    // 使用 fixed 和 setprecision 可以让输出更好看，这里保持原样但数值正确
    cout << "% of full length barcode reads:\t" << ratio << "%\n==================================================" << endl << endl;
    return;
}

void convert_hic2::run(string prefix){
    // processing scHiC with 10X barcodes: prefix_R1.combined.fq.gz
    // prefix_BC.sam is the mapped barcodes.
    
    int total = 0;
    int pass = 0;
    string cmd = "cat " + prefix; // prefix 应该是 sam 文件名
    FILE * inbam = popen(cmd.c_str(), "r");
    if(!inbam){ cerr << "Error opening SAM: " << cmd << endl; exit(1); }

    string out_name = prefix.substr(0, prefix.length()-4) + "_cov.fq.gz";
    cmd = "gzip - > " + out_name;
    FILE * fout = popen(cmd.c_str(), "w");

    samline align_line;
    fqline fastq_line;
    char buffer[4096];
    
    while(fgets(buffer, sizeof(buffer), inbam)){
        // 快速跳过 Header
        if(buffer[0] == '@') continue;
        
        // 优化：使用 samline::init 解析
        // 这里需要先把 buffer 转 string 或者修改 samline::init 接受 char*
        // 为了兼容 cxstring.hpp，我们构造 string，但得益于之前的优化，开销已降低
        string line_str(buffer); 
        
        ++total;
        align_line.init(line_str);
        
        if(align_line.chr == "*") continue;

        // 原逻辑解析 readname: Name:Info:Info...
        // 优化：cxstring::split 现在效率很高，可以使用
        vector<string> tmp = cxstring::split(align_line.readname, ":");
        
        // 确保切分后的部分足够多，防止越界崩溃
        if (tmp.size() < 8) {
            // cerr << "Warning: Malformed readname at line " << total << endl;
            continue; 
        }

        // 构造 FASTQ
        // tmp[0] 是原始 Readname (可能含部分)，tmp 后面是追加的信息
        // 原代码逻辑比较硬编码，假设格式固定。
        // tmp[7] 是 seq
        
        // 优化字符串拼接
        fastq_line.readname = "@" + align_line.chr + ":" + tmp[0] + ":" + tmp[1] + ":" + tmp[2] + ":" + tmp[3] + ":" + tmp[4] + ":" + tmp[5] + ":" + tmp[6];
        fastq_line.seq = tmp[7];
        
        // fix bug for NovaSeq (Qual string length match)
        // 原逻辑：从 align_line.readname 的末尾截取和 seq 长度一样的 qual
        if (align_line.readname.length() >= tmp[7].length()) {
             fastq_line.qual = align_line.readname.substr(align_line.readname.length() - tmp[7].length());
        } else {
             fastq_line.qual = tmp[7]; // Fallback, 虽然这意味数据有问题
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
        if(argc < 3){
            cerr<<"hictools combine_hic [atac/arc/rna] [prefix] [gz/bz2]"<<endl;
            return 1;
        }
        string mode = argv[2];
        if((mode != "atac") && (mode != "arc") && (mode != "rna")){
            cerr << "Invalid mode: " << mode << endl;
            cerr << "Need to specify mode (arc or atac or rna)!"<<endl;
            return 1;
        }
        // Default to gz if not specified
        string type = (argc >= 5) ? argv[4] : "gz";
        combine_hic::run(mode, argv[3], type);
        return 0;
    }

    if(mod == "convert_hic2"){
        if(argc < 3){
            cerr << "hictools convert_hic2 [sample_BC.sam]. This will add barcode to the beginning of read name than the end." << endl;
            return 1;
        }
        convert_hic2::run(argv[2]);
        return 0;
    } 
    
    cerr << "Unknown module: " << mod << endl;
    return 1;
}