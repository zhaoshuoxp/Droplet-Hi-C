use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio, exit};

// ==========================================
// 工具函数：检查文件是否存在
// ==========================================
fn file_exists(path: &str) -> bool {
    Path::new(path).exists()
}

// ==========================================
// 工具函数：创建智能解压和压缩的管道
// 自动优先使用多线程的 pigz，失败则降级为 gzip/zcat
// ==========================================
fn get_reader(path: &str, threads: &str) -> BufReader<std::process::ChildStdout> {
    let child = Command::new("pigz")
        .args(["-p", threads, "-dc"]) // 动态传入解压线程
        .arg(path)
        .stdout(Stdio::piped())
        .spawn()
        .unwrap_or_else(|_| {
            Command::new("zcat")
                .arg(path)
                .stdout(Stdio::piped())
                .spawn()
                .unwrap_or_else(|_| panic!("Error: Neither pigz nor zcat is available to read {}", path))
        });
    BufReader::with_capacity(256 * 1024, child.stdout.unwrap())
}

fn get_writer(path: &str, threads: &str) -> BufWriter<std::process::ChildStdin> {
    let out_file = File::create(path).unwrap_or_else(|_| panic!("Cannot create output file: {}", path));
    let child = Command::new("pigz")
        .args(["-p", threads, "-c"]) // 动态传入压缩线程
        .stdin(Stdio::piped())
        .stdout(Stdio::from(out_file.try_clone().unwrap()))
        .spawn()
        .unwrap_or_else(|_| {
            Command::new("gzip")
                .arg("-c")
                .stdin(Stdio::piped())
                .stdout(Stdio::from(out_file))
                .spawn()
                .expect("Error: Neither pigz nor gzip is available")
        });
    BufWriter::with_capacity(256 * 1024, child.stdin.unwrap())
}

// ==========================================
// 核心逻辑 1：清洗 Readname (零内存分配切片)
// ==========================================
fn clean_readname(raw_name: &str) -> &str {
    let trimmed = raw_name.trim_end();
    // 截取到第一个空格
    let base = match trimmed.find(|c| c == ' ' || c == '\t') {
        Some(pos) => &trimmed[..pos],
        None => trimmed,
    };

    // 探测并剥离异常的 10X 附加条码，如 ":rACGT..."
    if let Some(colon_r_pos) = base.rfind(":r") {
        if colon_r_pos > 0 {
            let suffix = &base[colon_r_pos + 2..];
            let is_barcode = suffix.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'N'));
            if is_barcode {
                return &base[..colon_r_pos];
            }
        }
    }
    base
}

// 高效复用内存的 Fastq 记录载体
struct FastqRecord {
    name: String,
    seq: String,
    qual: String,
}

impl FastqRecord {
    fn new() -> Self {
        FastqRecord {
            name: String::with_capacity(256),
            seq: String::with_capacity(256),
            qual: String::with_capacity(256),
        }
    }

    fn read_from<R: BufRead>(&mut self, reader: &mut R) -> std::io::Result<bool> {
        self.name.clear();
        if reader.read_line(&mut self.name)? == 0 {
            return Ok(false);
        }
        self.seq.clear();
        reader.read_line(&mut self.seq)?;
        let mut plus = String::new(); // + 行直接抛弃，节约内存
        reader.read_line(&mut plus)?;
        self.qual.clear();
        reader.read_line(&mut self.qual)?;

        // 去掉末尾的换行符
        if self.name.ends_with('\n') { self.name.pop(); }
        if self.name.ends_with('\r') { self.name.pop(); }
        if self.seq.ends_with('\n') { self.seq.pop(); }
        if self.seq.ends_with('\r') { self.seq.pop(); }
        if self.qual.ends_with('\n') { self.qual.pop(); }
        if self.qual.ends_with('\r') { self.qual.pop(); }

        Ok(true)
    }
}

// ==========================================
// 模块：combine_hic
// ==========================================
fn run_combine_hic(mode: &str, r2_prefix: &str, threads: &str) {
    let mut file_r1 = String::new();
    let mut file_r2 = String::new();
    let mut file_r3 = String::new();

    for ext in &[".fq.gz", ".fastq.gz"] {
        if file_exists(&format!("{}_R1{}", r2_prefix, ext)) { file_r1 = format!("{}_R1{}", r2_prefix, ext); }
        if file_exists(&format!("{}_R2{}", r2_prefix, ext)) { file_r2 = format!("{}_R2{}", r2_prefix, ext); }
        if file_exists(&format!("{}_R3{}", r2_prefix, ext)) { file_r3 = format!("{}_R3{}", r2_prefix, ext); }
    }

    if file_r1.is_empty() || file_r2.is_empty() || file_r3.is_empty() {
        eprintln!("Error: Cannot find matching _R1, _R2, or _R3 (.fq.gz/.fastq.gz) for prefix {}", r2_prefix);
        exit(1);
    }

    // 传入指定的线程数
    let mut red1 = get_reader(&file_r1, threads);
    let mut red2 = get_reader(&file_r2, threads);
    let mut red3 = get_reader(&file_r3, threads);

    let mut outfile1 = get_writer(&format!("{}_R1_combined.fq.gz", r2_prefix), threads);
    let mut outfile2 = get_writer(&format!("{}_R3_combined.fq.gz", r2_prefix), threads);

    let mut total: u64 = 0;
    let mut pass: u64 = 0;

    let mut r1 = FastqRecord::new();
    let mut r2 = FastqRecord::new();
    let mut r3 = FastqRecord::new();

    while r1.read_from(&mut red1).unwrap() {
        r2.read_from(&mut red2).unwrap();
        r3.read_from(&mut red3).unwrap();
        total += 1;

        if mode == "atac" {
            let seq = if r2.seq.len() > 16 { &r2.seq[..16] } else { &r2.seq };
            let qual = if r2.qual.len() > 16 { &r2.qual[..16] } else { &r2.qual };
            let base_name = clean_readname(&r2.name);

            let new_name = format!("{}:{}:{}", base_name, r1.seq, r1.qual);
            write!(outfile1, "{}\n{}\n+\n{}\n", new_name, seq, qual).unwrap();

            let new_name3 = format!("{}:{}:{}", base_name, r3.seq, r3.qual);
            write!(outfile2, "{}\n{}\n+\n{}\n", new_name3, seq, qual).unwrap();
            pass += 1;

        } else if mode == "arc" {
            if r2.seq.len() < 24 { continue; }
            let seq = &r2.seq[8..24];
            let qual = &r2.qual[8..24];
            let base_name = clean_readname(&r2.name);

            let new_name = format!("{}:{}:{}", base_name, r1.seq, r1.qual);
            write!(outfile1, "{}\n{}\n+\n{}\n", new_name, seq, qual).unwrap();

            let new_name3 = format!("{}:{}:{}", base_name, r3.seq, r3.qual);
            write!(outfile2, "{}\n{}\n+\n{}\n", new_name3, seq, qual).unwrap();
            pass += 1;

        } else if mode == "rna" {
            if r1.seq.len() < 28 { continue; }
            let barcode = &r1.seq[0..16];
            let umi = &r1.seq[16..28];
            let qual = &r1.qual[0..16];
            let base_name = clean_readname(&r1.name);

            let new_name = format!("{}:{}:{}:{}", base_name, umi, r3.seq, r3.qual);
            write!(outfile2, "{}\n{}\n+\n{}\n", new_name, barcode, qual).unwrap();
            pass += 1;
        }
    }

    let ratio = if total > 0 { (pass as f64 / total as f64) * 100.0 } else { 0.0 };
    println!("==================================================");
    println!("(10X) Barcode Locator Report: {}", r2_prefix);
    println!("barcodes modes:\t\t{}", mode);
    println!("# total raw reads:\t\t{}", total);
    println!("# of full barcoded reads:\t{}", pass);
    println!("% of full length barcode reads:\t{:.2}%", ratio);
    println!("==================================================\n");
}

// ==========================================
// 模块：convert_hic2
// ==========================================
fn run_convert_hic2(prefix: &str, threads: &str) {
    let mut total: u64 = 0;
    let mut pass: u64 = 0;

    // 直接读取 SAM 文件 (不用 cat 建立子进程，效率更高)
    let file = File::open(prefix).unwrap_or_else(|_| panic!("Error opening SAM: {}", prefix));
    let reader = BufReader::with_capacity(256 * 1024, file);

    let out_name = if prefix.ends_with(".sam") {
        format!("{}_cov.fq.gz", &prefix[..prefix.len() - 4])
    } else {
        format!("{}_cov.fq.gz", prefix)
    };
    
    // 传入指定的线程数用于压缩
    let mut fout = get_writer(&out_name, threads);

    for line_result in reader.lines() {
        let line = line_result.unwrap();
        if line.starts_with('@') || line.is_empty() { continue; }
        total += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }
        
        let qname = fields[0];
        let chr = fields[2];
        if chr == "*" { continue; }

        let tmp: Vec<&str> = qname.split(':').collect();
        if tmp.len() < 3 { continue; }

        // 动态定位纯 DNA 序列位置
        let mut seq_idx: i32 = -1;
        for (i, part) in tmp.iter().enumerate().skip(1) {
            if part.len() > 20 && part.chars().all(|c| matches!(c, 'A'|'C'|'G'|'T'|'N')) {
                seq_idx = i as i32;
                break;
            }
        }

        if seq_idx == -1 {
            if tmp.len() > 7 { seq_idx = 7; } else { continue; }
        }

        let seq = tmp[seq_idx as usize];
        let seq_len = seq.len();
        let expected_min_len = 2 * seq_len + 2;

        let final_name;
        let final_qual;

        if qname.len() >= expected_min_len {
            final_qual = &qname[qname.len() - seq_len..];
            let base_name = &qname[..qname.len() - expected_min_len];
            final_name = format!("@{}:{}", chr, base_name);
        } else {
            final_qual = seq; // Fallback
            final_name = format!("@{}:{}", chr, tmp[0]);
        }

        write!(fout, "{}\n{}\n+\n{}\n", final_name, seq, final_qual).unwrap();
        pass += 1;
    }

    let base_name = if prefix.ends_with(".sam") { &prefix[..prefix.len() - 4] } else { prefix };
    println!("==================================================");
    println!("{} reads processed in {}", total, base_name);
    println!("{} mapped reads in {}", pass, base_name);
    println!("==================================================");
}

// ==========================================
// 主入口
// ==========================================
fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: hictools [mode] ...");
        exit(1);
    }

    let mod_name = &args[1];

    if mod_name == "combine_hic" {
        if args.len() < 4 {
            eprintln!("Usage: hictools combine_hic [atac/arc/rna] [prefix] [threads]");
            exit(1);
        }
        let mode = &args[2];
        if mode != "atac" && mode != "arc" && mode != "rna" {
            eprintln!("Invalid mode: {}", mode);
            exit(1);
        }
        let prefix = &args[3];
        // 如果提供了第4个参数则作为线程数，否则默认 8
        let threads = if args.len() > 4 { &args[4] } else { "8" };
        
        run_combine_hic(mode, prefix, threads);

    } else if mod_name == "convert_hic2" {
        if args.len() < 3 {
            eprintln!("Usage: hictools convert_hic2 [sample_BC.sam] [threads]");
            exit(1);
        }
        // 如果提供了第3个参数则作为线程数，否则默认 8
        let threads = if args.len() > 3 { &args[3] } else { "8" };
        
        run_convert_hic2(&args[2], threads);

    } else {
        eprintln!("Unknown module: {}", mod_name);
        exit(1);
    }
}
