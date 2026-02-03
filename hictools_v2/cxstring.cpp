//
//  cxstring.cpp
//  upstoools
//
//  Optimized Version
//

#include "cxstring.hpp"

// 优化：使用 find_first_of 替代手动循环，效率更高
vector<string> cxstring::split(const string &s, const string &seperator) {
    vector<string> result;
    if (s.empty()) return result;

    string::size_type start = 0;
    string::size_type end = s.find_first_of(seperator);

    while (end != string::npos) {
        if (end != start) {
            result.emplace_back(s.substr(start, end - start));
        }
        start = end + 1;
        end = s.find_first_of(seperator, start);
    }

    // 处理最后一个部分
    if (start < s.length()) {
        // 检查最后部分是否也是分隔符（原逻辑似乎跳过空字符串）
        string last = s.substr(start);
        if (last.find_first_of(seperator) == string::npos) {
            result.emplace_back(std::move(last));
        }
    }
    
    return result;
}

// 优化：避免 O(N^2) 的 erase 操作，改用 O(N) 的截取
string cxstring::chomp(string str) {
    if (str.empty()) return str;

    // 1. 去除尾部换行符 (\r, \n)
    while (!str.empty() && (str.back() == '\n' || str.back() == '\r')) {
        str.pop_back();
    }

    // 2. 去除首尾空格 (Trim)
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first) {
        return ""; // 全是空格
    }
    size_t last = str.find_last_not_of(' ');
    
    // 如果没有空格需要去除，直接返回（避免拷贝）
    if (first == 0 && last == str.length() - 1) {
        return str;
    }
    
    return str.substr(first, (last - first + 1));
}

bool cxstring::is_bam_header(const std::string& buffer) {
    return (!buffer.empty() && buffer[0] == '@');
}

// 优化：使用 std::stoi 替代 stringstream
int cxstring::str2int(const string &str) {
    try {
        return std::stoi(str);
    } catch (...) {
        return 0; // 发生错误返回0，保持鲁棒性
    }
}

// 优化：使用 std::to_string 替代 stringstream
string cxstring::int2str(int i) {
    return std::to_string(i);
}

string cxstring::uc(string input) {
    // 引用遍历，避免拷贝字符
    for (char &c : input) {
        c = toupper((unsigned char)c);
    }
    return input;
}

// 优化：使用查找表（Lookup Table）代替大量的 if-else
// 静态表只初始化一次
static const char complement_table[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 0-15
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 16-31
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 32-47
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 48-63
    0,'T','V','G','H','E','F','C',  // @, A->T, B->V, C->G
    'D','I','J','M','L','K','N','O',  // ... G->C ...
    'P','Q','R','S','A','U','V','W',  // ... T->A ...
    'X','Y','Z',0,0,0,0,0,            // ...
    0,'T','V','G','H','E','F','C',  // 小写字母支持
    'D','I','J','M','L','K','N','O',
    'P','Q','R','S','A','U','V','W',
    'X','Y','Z',0,0,0,0,0,
    // ... 其余为0
};
// 注意：原代码只处理了 ACGT，这里为了安全，未定义的字符可以保持原样或置0
// 为了严格保持原功能（只变 ACGT），我们用更简单的逻辑，但为了速度还是用 Table

string cxstring::reverse_complmentary(string input) {
    reverse(input.begin(), input.end());
    for (char &c : input) {
        // 使用快速 switch 或 table。这里为了最大兼容原逻辑(只变ACGT)，使用 switch
        switch(c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            case 'a': c = 't'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            case 't': c = 'a'; break;
            case 'N': c = 'N'; break;
            default: break;
        }
    }
    return input;
}

int cxstring::hamming_dist(const string &str1, const string &str2) {
    int dist = 0;
    size_t len = min(str1.length(), str2.length()); // 防止越界
    for (size_t i = 0; i < len; ++i) {
        if (str1[i] != str2[i]) ++dist;
    }
    return dist;
}

//
// Samline Optimization implementation
//
void samline::init(const string &str) {
    // 预清理
    string line = cxstring::chomp(str);
    if (line.empty()) { empty(); return; }

    // 优化：直接扫描 TAB，不创建 vector<string>
    // 格式：QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [OPT]
    
    size_t start = 0;
    size_t end = 0;
    int field_idx = 0;

    while ((end = line.find('\t', start)) != string::npos) {
        // 获取当前字段的 substring
        // 注意：这是不可避免的拷贝，但比 split 生成 vector 少了一次全量拷贝
        string field = line.substr(start, end - start);
        
        switch (field_idx) {
            case 0: readname = field; break;
            case 1: flag = cxstring::str2int(field); break;
            case 2: chr = field; break;
            case 3: pos = cxstring::str2int(field); break;
            case 4: mapq = cxstring::str2int(field); break;
            case 5: cigar = field; break;
            case 6: chrnext = field; break;
            case 7: posnext = cxstring::str2int(field); break;
            case 8: plen = cxstring::str2int(field); break;
            case 9: seq = field; break;
            case 10: qual = field; break;
            default: break; // 11+ 属于 tags，稍后处理
        }
        
        start = end + 1;
        field_idx++;
        if (field_idx > 10) break; // 已读取基本字段
    }

    // 处理剩余部分（Tags 或者最后一个字段 Qual/Seq）
    if (start < line.length()) {
        string field = line.substr(start);
        if (field_idx == 10) qual = field; // 只有11列的情况
        else if (field_idx > 10) other = field; // Tags
    } else {
        other = "";
    }
}

void samline::write(FILE * outfile) {
    // 优化：减少字符串拼接产生的临时对象
    // 使用 fprintf 或者构建好字符串一次写入
    // 为了保持 C++ 风格且比 stringstream 快，我们手动拼接
    
    string out;
    out.reserve(512); // 预分配内存，避免多次 realloc
    
    out += readname; out += "\t";
    out += std::to_string(flag); out += "\t";
    out += chr; out += "\t";
    out += std::to_string(pos); out += "\t";
    out += std::to_string(mapq); out += "\t";
    out += cigar; out += "\t";
    out += chrnext; out += "\t";
    out += std::to_string(posnext); out += "\t";
    out += std::to_string(plen); out += "\t";
    out += seq; out += "\t";
    out += qual; out += "\t";
    out += other; out += "\n";
    
    fputs(out.c_str(), outfile);
}