//
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "cxstring.hpp"
#include "cxstring.cpp"

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
  
  // processing scHiC with 10X barcodes: R1, R3 are paired-end reads, R2 is 10X index reads, with fixed length 16bp
  // mode: atac (use 10XATAC barcodes, 16bp) or arc (use 10Xmultiome barcodes, 8+16bp)
  int total = 0;
  int pass = 0;
  string s1;
  string s2;
  string s3;
  if(ty == "gz"){
  s1 = "zcat ";
  s2 = r2 + "_R1.fq.gz";
  }
  else if(ty == "bz2"){
    s1 = "bzcat ";
    s2 = r2 + "_R1.fastq.bz2";
  }
  s3 = s1 + s2;
  FILE * red1;
  red1 = popen(s3.c_str(), "r");

  // read R1 for DNA
  //if(mode != "rna"){
  s1 = "gzip - > ";
  s2 = r2 + "_R1_combined.fq.gz";
  s3 = s1 + s2;
  FILE * outfile1;
  outfile1 = popen(s3.c_str(), "w");
  //}
  if(ty == "gz"){
    s1 = "zcat ";
    s2 = r2 + "_R2.fq.gz";
  }
  else if(ty == "bz2"){
    s1 = "bzcat ";
    s2 = r2 + "_R2.fastq.bz2";
  }
  s3 = s1 + s2;
  FILE * red2;
  red2 = popen(s3.c_str(), "r");

  // read R3 for DNA & RNA
  if(ty == "gz"){
  s1 = "zcat ";
  s2 = r2 + "_R3.fq.gz";
  }
  else if(ty == "bz2"){
    s1 = "bzcat ";
    s2 = r2 + "_R3.fastq.bz2";
  }
  s3 = s1 + s2;
  FILE * red3;
  red3 = popen(s3.c_str(), "r");
  s1 = "gzip - > ";
  s2 = r2 + "_R3_combined.fq.gz";
  s3 = s1 + s2;
  FILE * outfile2;
  outfile2 = popen(s3.c_str(), "w");

  char buffer[2000];
  fqline in_line1;
  fqline in_line2;
  fqline in_line3;
  while(fgets(buffer, sizeof(buffer), red1)){
    ++total;
    string line1(buffer);
    line1 = cxstring::chomp(line1);
    in_line1.read_part_record(red1, line1);
    in_line2.read_full_record(red2);
      in_line3.read_full_record(red3);
    if(mode == "atac"){
      string new_seq = in_line2.seq.substr(0,16); 
      // if(in_line2.seq.length()!=16)continue;
      in_line2.seq = new_seq;
      in_line2.qual = in_line2.qual.substr(0,16);
      string a;
      stringstream as;
      as << in_line2.readname;
      as >> a;
      in_line2.readname = a +  ":" + in_line1.seq + ":" + in_line1.qual;
      //in_line.readname = in_line.readname + ":" + umi;
      in_line2.write_record(outfile1);
      in_line2.readname = a +  ":" + in_line3.seq + ":" + in_line3.qual;
      in_line2.write_record(outfile2);
      ++pass;
    } else if(mode == "arc"){
      string new_seq = in_line2.seq.substr(8,24);
      if(in_line2.seq.length()!=24)continue;
      in_line2.seq = new_seq;
      in_line2.qual = in_line2.qual.substr(8,24);
      string a;
      stringstream as;
      as << in_line2.readname;
      as >> a;
      in_line2.readname = a +  ":" + in_line1.seq + ":" + in_line1.qual;
      in_line2.write_record(outfile1);
      in_line2.readname = a +  ":" + in_line3.seq + ":" + in_line3.qual;
      in_line2.write_record(outfile2);
      ++pass;
    } else if(mode == "rna"){
      string new_seq = in_line1.seq.substr(0,16);
      string umi = in_line1.seq.substr(16,28);
      //if(in_line1.seq.length()!=28)continue;
      in_line1.seq = new_seq;
      in_line1.qual = in_line1.qual.substr(0,16);
      string a;
      stringstream as;
      as << in_line1.readname;
      as >> a;
      in_line1.readname = a +  ":" + umi + ":" + in_line3.seq + ":" + in_line3.qual;
      //in_line.readname = in_line.readname + ":" + umi;
      in_line1.write_record(outfile2);
      ++pass;
    }
  }
  pclose(red1);
  pclose(red2);
  pclose(red3);
  pclose(outfile1);
  pclose(outfile2);

  int aaa = pass*10000/total;
  float ratio = (float)aaa / 100;
  cout << "==================================================\n(10X) Barcode Locator Report: " << r2 << endl;
  cout << "barcodes modes:\t\t" << mode << endl;
  cout << "# total raw reads:\t\t" << total << endl << "# of full barcoded reads:\t" << pass << endl;
  cout << "% of full length barcode reads:\t" << ratio << "%\n==================================================" << endl << endl;
  return;
}

void convert_hic2::run(string prefix){
  // processing scHiC with 10X barcodes: prefix_R1.combined.fq.gz, prefix_R3.combined.fq.gz are combined reads. 
  // prefix_BC.sam is the mapped barcodes. Add barcode to HEAD of the readname
  int total = 0;
  int pass = 0;
  string s1 = "cat ";
  string s2 = prefix;
  string s3 = s1 + s2;
  FILE * inbam;
  inbam = popen(s3.c_str(), "r");
  s1 = "gzip - > ";
  s2 = prefix.substr(0, prefix.length()-4) + "_cov.fq.gz";
  s3 = s1 + s2;
  FILE * fout;
  fout = popen(s3.c_str(), "w");
  samline align_line;
  fqline fastq_line;
  char buffer[2000];
  while(fgets(buffer, sizeof(buffer), inbam)){
    string line(buffer);
    line=cxstring::chomp(line);
    if(line.substr(0, 1) == "@"){
      continue;
    }
    ++total;
    align_line.init(line);
    if(align_line.chr == "*")continue;
    vector<string> tmp = cxstring::split(align_line.readname, ":");
    fastq_line.readname = "@" + align_line.chr + ":" + tmp[0] + ":" + tmp[1] + ":" + tmp[2] + ":" + tmp[3] + ":" + tmp[4] + ":" + tmp[5] + ":" + tmp[6];
    fastq_line.seq = tmp[7];
    //fastq_line.qual = tmp[8];
    fastq_line.qual = align_line.readname.substr(align_line.readname.length()-tmp[7].length(), tmp[7].length()); // fix bug for NovaSeq
    fastq_line.mark = "+";
    fastq_line.write_record(fout);
    ++pass;
  }
  pclose(inbam);
  pclose(fout);
  cout << "==================================================\n" << total << " reads processed in " << prefix.substr(0, prefix.length()-4) << endl;
  cout << pass << " mapped reads in " << prefix.substr(0, prefix.length()-4) << "\n==================================================" << endl;
  return;
}


int main(int argc, const char * argv[]) {
  string mod(argv[1]);
  if(mod == "combine_hic"){
    if(argc < 3){
        cerr<<"hictools combine_hic [atac/arc/rna] [prefix] [gz/bz2]"<<endl;
        return 1;
      }
      if((argv[2] != string("atac")) && (argv[2] != string("arc")) && (argv[2] != string("rna"))){
        cerr<<argv[2];
        cerr<<"Need to specify mode (arc or atac or rna)!"<<endl;
        return 1;
      }
      if(argc < 5){
        combine_hic::run(argv[2], argv[3], "gz");
        return 0;
      }
      combine_hic::run(argv[2], argv[3], argv[4]);
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
}
