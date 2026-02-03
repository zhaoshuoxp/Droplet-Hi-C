//
//  cxstring.hpp
//  upstoools
//
//  Optimized Version
//

#ifndef cxstring_hpp
#define cxstring_hpp

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>

using namespace std;

class cxstring {
public:
    // 优化：传参改为 const reference
    static vector<string> split(const string &s, const string &seperator);
    static string chomp(string str); // 传值后内部修改返回，或者改为传引用修改均可，此处保持接口
    static bool is_bam_header(const std::string& buffer);
    static int str2int(const string &str);
    static string int2str(int i);
    static string uc(string input);
    static string reverse_complmentary(string input);
    static int hamming_dist(const string &str1, const string &str2);
};

class samline {
public:
    string readname, chr, cigar, chrnext, seq, qual, other;
    int flag, pos, mapq, posnext, plen;
        
    // 优化：不再依赖 split 产生临时 vector，直接解析
    void init(const string &str);
    
    void empty() {
        readname.clear();
        flag = 0;
        chr.clear();
        pos = 0;
        mapq = 0;
        cigar.clear();
        chrnext.clear();
        posnext = 0;
        plen = 0;
        seq.clear();
        qual.clear();
        other.clear();
    }
    
    void write(FILE * outfile);
};

class fqline {
public:
    string readname, seq, mark, qual;

    // 辅助函数：快速清理读取的一行（去换行符）
    void _assign_chomp(string &dest, const char* buffer) {
        dest = buffer;
        if (!dest.empty() && dest.back() == '\n') dest.pop_back();
        if (!dest.empty() && dest.back() == '\r') dest.pop_back();
    }

    void read_part_record(FILE * infile, const string &rn) {
        readname = rn;
        char buffer[4096]; // 稍微加大 buffer 比较安全
        
        if(fgets(buffer, sizeof(buffer), infile)) _assign_chomp(seq, buffer);
        if(fgets(buffer, sizeof(buffer), infile)) _assign_chomp(mark, buffer);
        if(fgets(buffer, sizeof(buffer), infile)) _assign_chomp(qual, buffer);
    }
    
    void read_full_record(FILE * infile) {
        char buffer[4096];
        
        if(fgets(buffer, sizeof(buffer), infile)) {
             // 模拟原逻辑：stringstream >> readname 会读取到空格为止
             // 这里手动处理：读取整行后截取到第一个空格
             string tmp = buffer;
             size_t space_pos = tmp.find_first_of(" \t\n\r");
             if(space_pos != string::npos) readname = tmp.substr(0, space_pos);
             else readname = tmp;
             // 如果原逻辑不需要严格模拟 >> 的截断行为，直接用 _assign_chomp 会更快
        }
        
        if(fgets(buffer, sizeof(buffer), infile)) _assign_chomp(seq, buffer);
        if(fgets(buffer, sizeof(buffer), infile)) _assign_chomp(mark, buffer);
        if(fgets(buffer, sizeof(buffer), infile)) _assign_chomp(qual, buffer);
    }
    
    void write_record(FILE * outfile) {
        // 使用 fputc/fputs 比构建大字符串稍微快一点点，也更省内存
        fputs(readname.c_str(), outfile); fputc('\n', outfile);
        fputs(seq.c_str(), outfile);      fputc('\n', outfile);
        fputs(mark.c_str(), outfile);     fputc('\n', outfile);
        fputs(qual.c_str(), outfile);     fputc('\n', outfile);
    }
};

#endif /* cxstring_hpp */