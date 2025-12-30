options(shiny.maxRequestSize = 2 * 1024^3)

source("blockzoom.f.R", local = TRUE)
source("locuszoom+LD.R", local = TRUE)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  shiny, data.table, readr, DT, ggplot2, CMplot, rhdf5,
  ggrepel, scales, ggrastr, rtracklayer, GenomicRanges, GenomeInfoDb,cowplot,jsonlite
)

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
try(rhdf5::h5closeAll(), silent = TRUE)

Sys.setenv(OMP_NUM_THREADS = "4", 
           MKL_NUM_THREADS = "4", 
           OPENBLAS_NUM_THREADS = "4",
           RCPP_PARALLEL_NUM_THREADS = "4")

# HDL
.pick <- function(nms, alts){ z <- intersect(alts, nms); if (length(z)) z[1] else NA_character_ }
.clean_snp <- function(x){
  x <- toupper(trimws(as.character(x)))
  # rsid 原样返回（只去空白）
  is_rsid <- grepl("^RS\\d+$", x)
  y <- x
  # chr:pos(:alleles) → 只保留前两段 "CHR:POS"
  y[!is_rsid] <- sub("^CHR", "", y[!is_rsid])
  y[!is_rsid] <- sub("^([0-9]{1,2}|X|Y|MT|M)[:|_/ ]+([0-9]+).*$", "\\1:\\2", y[!is_rsid], perl=TRUE)
  y
}
.ensure_Z <- function(dt){
  nm0 <- names(dt)
  nm  <- toupper(trimws(sub("^\ufeff","", nm0)))
  data.table::setnames(dt, nm0, nm, skip_absent = TRUE)
  
  if ("Z" %in% names(dt)) return(suppressWarnings(as.numeric(dt$Z)))
  
  b <- intersect(c("BETA","B"), names(dt))[1]
  s <- intersect(c("SE","STDERR","SE_BETA","STDERR_BETA"), names(dt))[1]
  if (!is.na(b) && !is.na(s))
    return(suppressWarnings(as.numeric(dt[[b]])/as.numeric(dt[[s]])))
  
  pcol <- intersect(c("P","PVAL","P_VALUE","PVALUE","NEG_LOG_PVALUE"), names(dt))[1]
  if (!is.na(pcol)) {
    p <- if (pcol == "NEG_LOG_PVALUE")
      10^(-suppressWarnings(as.numeric(dt[[pcol]])))
    else
      suppressWarnings(as.numeric(dt[[pcol]]))
    p[!is.finite(p) | p<=0] <- NA
    pminpos <- suppressWarnings(min(p[p>0], na.rm=TRUE)); if (!is.finite(pminpos)) pminpos <- 1e-12
    p[is.na(p)] <- pminpos/10
    p <- pmax(pmin(p, 1-1e-16), 1e-300)
    return(qnorm(p/2, lower.tail=FALSE))
  }
  stop("sumstats 需要 Z 或 (BETA,SE) 或 P/NEG_LOG_PVALUE")
}
.norm_chr_pos <- function(dt){
  if ("CHR" %in% names(dt)) {
    dt[, CHR := gsub("^chr", "", as.character(CHR), ignore.case = TRUE)]
    dt[tolower(CHR)=="x", CHR := "23"]
    dt[tolower(CHR)=="y", CHR := "24"]
    dt[tolower(CHR) %in% c("mt","m"), CHR := "25"]
    suppressWarnings(dt[, CHR := as.integer(CHR)])
  }
  if ("POS" %in% names(dt)) {
    suppressWarnings(dt[, POS := as.integer(POS)])
  }
  dt[]
}
.h5_ls_safe <- function(file, ...) {
  on.exit(try(rhdf5::h5closeAll(), silent = TRUE), add = TRUE)
  rhdf5::h5ls(file, all = TRUE, ...)
}
.h5_read_safe <- function(file, name, ...) {
  on.exit(try(rhdf5::h5closeAll(), silent = TRUE), add = TRUE)
  rhdf5::h5read(file, name, ...)
}
.H5_GRPS_CACHE <- new.env(parent = emptyenv())
.get_h5_grps <- function(h5file){
  key <- normalizePath(h5file, winslash = "/", mustWork = FALSE)
  if (exists(key, envir = .H5_GRPS_CACHE)) return(get(key, envir = .H5_GRPS_CACHE))
  info   <- data.table::as.data.table(.h5_ls_safe(h5file))
  has_ds <- info[otype == "H5I_DATASET", .(group = as.character(group), name)]
  grps   <- intersect(has_ds[name == "snplist", unique(group)],
                      has_ds[name == "ldblk",   unique(group)])
  assign(key, grps, envir = .H5_GRPS_CACHE)
  grps
}
rank_blocks_topN_by_bed_simple <- function(gwas_dt, bed_file){
  bed <- .read_ld_blocks_bed(bed_file)  # CHR, START, END, BLK
  dt  <- data.table::as.data.table(gwas_dt)
  
  # 只需要 SNP/CHR/POS/P
  dt <- dt[, .(SNP, CHR = as.integer(CHR), POS = as.integer(POS), P = as.numeric(P))]
  dt <- dt[is.finite(CHR) & is.finite(POS) & is.finite(P)]
  if (!nrow(dt)) return(data.table::data.table())
  
  out_list <- vector("list", nrow(bed))
  for (i in seq_len(nrow(bed))) {
    b <- bed[i]
    x <- dt[CHR == b$CHR & POS >= b$START & POS <= b$END]
    if (!nrow(x)) next
    j <- which.min(x$P)
    out_list[[i]] <- data.table::data.table(
      blk = as.character(b$BLK),
      `lead snp` = as.character(x$SNP[j]),
      chr = as.integer(x$CHR[j]),
      pos = as.integer(x$POS[j]),
      p   = as.numeric(x$P[j])
    )
  }
  data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

# tiny helpers
scale_p <- function(p, min_p = 1e-324, max_log10 = 50) {
  p <- ifelse(p < min_p, min_p, p)
  log10p <- -log10(p)
  log10p.max <- max(log10p, na.rm=TRUE)
  log10p <- ifelse(
    is.na(p) | log10p.max < max_log10, 
    log10p, 
    ifelse(log10p <= 10, log10p, 10 + (log10p - 10) * (max_log10 - 10) / (log10p.max - 10))
  )
  10^(-log10p)
}
read_gwas <- function(path, build_disp = c("37","38")) {
  # 0) 规范 build 标记（默认用第一个，即 "37"）
  build_disp <- as.character(build_disp)[1L]
  
  # 1) 读取
  dt <- data.table::fread(
    file = path,
    header = TRUE,
    na.strings = c("NA","NaN","","Inf","-Inf","."),
    data.table = TRUE,
    fill = TRUE
  )
  stopifnot(is.data.frame(dt), nrow(dt) > 0)
  
  # 2) 列名标准化（转大写、去BOM）
  nm0 <- names(dt)
  nmU <- toupper(trimws(sub("^\ufeff","", nm0)))
  data.table::setnames(dt, nm0, nmU, skip_absent = TRUE)
  
  # 通用列名匹配器：忽略大小写和 # 前缀
  nm_nohash <- sub("^#", "", names(dt))
  names_map  <- stats::setNames(names(dt), nm_nohash)
  names(names_map) <- toupper(names(names_map))
  pick1 <- function(alts){
    altsU <- toupper(sub("^#","", alts))
    hit <- names_map[match(altsU, names(names_map))]
    hit <- hit[!is.na(hit)]
    if (length(hit)) hit[1] else NA_character_
  }
  
  # 3) 基础列
  chr_col <- pick1(c("CHR","CHROM","#CHROM","CHROMOSOME"))
  
  # 根据 build_disp 选择位置列候选
  if (build_disp == "37") {
    pos_alts <- c(
      "POS.37", "POS37", "BP_HG19",     # 优先 37 专用列
      "POS","BP","POSITION","BASE_PAIR_LOCATION"  # 通用列作为后备
    )
  } else if (build_disp == "38") {
    pos_alts <- c(
      "POS.38", "POS38", "BP_HG38",     # 优先 38 专用列
      "POS","BP","POSITION","BASE_PAIR_LOCATION","BP_HG19"  # 没有 38 列时退回通用
    )
  } else {
    # 兜底：什么都试一遍
    pos_alts <- c(
      "POS","BP","POSITION","BASE_PAIR_LOCATION",
      "BP_HG19","POS37","POS.37","POS38","POS.38","BP_HG38"
    )
  }
  pos_col <- pick1(pos_alts)
  
  snp_col <- pick1(c("SNP","MARKERNAME","ID","RSID","VARIANT"))
  
  # 4) 缺CHR/POS则尝试从SNP解析（形如 chr1:123456 或 1:123456）
  if (is.na(chr_col) || is.na(pos_col)) {
    if (!is.na(snp_col)) {
      s <- as.character(dt[[snp_col]])
      if (any(grepl(":", s))) {
        sp <- data.table::tstrsplit(gsub("(?i)^CHR","", s), "[:_]", keep = 1:2)
        dt[, CHR := sp[[1]]]
        dt[, POS := suppressWarnings(as.numeric(sp[[2]]))]
        chr_col <- "CHR"; pos_col <- "POS"
      }
    }
  }
  if (is.na(chr_col) || is.na(pos_col)) {
    stop("GWAS 里需要提供染色体与位置列（支持 CHR/CHROM/#CHROM 与 POS/BP）。")
  }
  
  # 5) 统一列名
  data.table::setnames(dt, chr_col, "CHR", skip_absent = TRUE)
  data.table::setnames(dt, pos_col, "POS", skip_absent = TRUE)
  if (!is.na(snp_col)) data.table::setnames(dt, snp_col, "SNP", skip_absent = TRUE)
  
  # 5.5) 规范 CHR/POS
  if ("CHR" %in% names(dt)) {
    dt[, CHR := gsub("(?i)^chr", "", as.character(CHR), perl = TRUE)]
    dt[tolower(CHR) == "x",  CHR := "23"]
    dt[tolower(CHR) == "y",  CHR := "24"]
    dt[tolower(CHR) %in% c("mt","m"), CHR := "25"]
    suppressWarnings(dt[, CHR := as.integer(CHR)])
  }
  if ("POS" %in% names(dt)) {
    suppressWarnings(dt[, POS := as.integer(POS)])
  }
  
  # 6) P 值
  p_col <- pick1(c("P","PVAL","P_VALUE","PVALUE","NEG_LOG_PVALUE"))
  if (is.na(p_col)) {
    stop("缺少 P 值列（P/PVAL/P_VALUE/PVALUE/NEG_LOG_PVALUE 均未找到）")
  } else if (toupper(p_col) == "NEG_LOG_PVALUE") {
    dv <- suppressWarnings(as.numeric(dt[[p_col]]))
    P  <- 10^(-dv)
  } else {
    if (p_col != "P") data.table::setnames(dt, p_col, "P", skip_absent = TRUE)
    P <- suppressWarnings(as.numeric(dt[["P"]]))
  }
  # 清洗 + 下限裁切
  P[!is.finite(P) | P <= 0] <- NA_real_
  P[P >= 1] <- 1 - 1e-16
  P <- pmax(P, 1e-50)
  dt[, P := P]
  
  # 7) 返回（仅保留有效 CHR/POS）
  dt <- dt[is.finite(CHR) & is.finite(POS)]
  dt
}
.infer_neff_from_dt <- function(dt, fallback = 1e5, neff_min = 5e3, neff_max = 2e6){
  if (!is.data.frame(dt) || !nrow(dt)) return(fallback)
  nm <- toupper(names(dt))
  get_num <- function(x) suppressWarnings(as.numeric(x))
  pick1 <- function(keys){
    hit <- match(keys, nm); hit <- hit[!is.na(hit)]
    if (!length(hit)) return(NA_real_)
    v <- get_num(dt[[ hit[1] ]]); v <- v[is.finite(v) & v>0]
    if (!length(v)) return(NA_real_) else median(v)
  }
  # 1) 直接 NEFF
  neff <- pick1(c("NEFF","N_EFF","N_EFFECTIVE","NEFFECTIVE"))
  # 2) 没有就用 N
  if (!is.finite(neff)) neff <- pick1(c("N","N_TOTAL","SAMPLESIZE","SAMPLE_SIZE"))
  # 3) 病例/对照
  if (!is.finite(neff)){
    n1 <- pick1(c("NCASE","N_CASE","CASES"))
    n2 <- pick1(c("NCTRL","N_CONTROL","CONTROLS"))
    if (is.finite(n1) && is.finite(n2)) neff <- 4/(1/n1 + 1/n2)
  }
  if (!is.finite(neff)) neff <- fallback
  # 合理范围裁剪
  neff <- min(max(neff, neff_min), neff_max)
  neff
}
.to_num <- function(x, default){
  if (is.null(x) || length(x) == 0) return(default)
  v <- suppressWarnings(as.numeric(gsub("\\s+", "", as.character(x))))
  if (length(v) == 0 || !is.finite(v)) default else v
}
.parse_chrpos_from_str <- function(x){
  s <- as.character(x)
  s <- gsub(",", "", s)
  # 抓 "chr1:123456", "1_123456", "1 123456", "1:123456_A_G" 等
  m <- regexec("(?i)^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[[:punct:]_ ]+([0-9]+)", s, perl = TRUE)
  mt <- regmatches(s, m)
  chr <- vapply(mt, function(v) if (length(v) >= 3) v[2] else NA_character_, "")
  pos <- vapply(mt, function(v) if (length(v) >= 3) v[3] else NA_character_, "")
  chr[toupper(chr) == "X"]  <- "23"
  chr[toupper(chr) == "Y"]  <- "24"
  chr[toupper(chr) %in% c("MT","M")] <- "25"
  list(CHR = suppressWarnings(as.integer(chr)),
       POS = suppressWarnings(as.integer(pos)))
}
.is_nonempty_tsv <- function(p){
  if (!file.exists(p)) return(FALSE)
  info <- file.info(p)
  !is.na(info$size) && info$size > 0
}

.read_ld_blocks_bed <- function(path_bed){
  stopifnot(file.exists(path_bed))
  bed <- data.table::fread(path_bed, header = FALSE, data.table = TRUE)
  
  # 如果第一行像表头（首列含非数字或含 "chr"），自动当作表头再读一次
  if (nrow(bed) > 0) {
    first1 <- as.character(bed$V1[1])
    looks_header <- grepl("(?i)chr|start|end|block|name", first1) || !suppressWarnings(is.finite(as.numeric(first1)))
    if (looks_header) {
      bed <- data.table::fread(path_bed, header = TRUE, data.table = TRUE)
      # 统一成 V1..Vn 方便处理
      old <- names(bed)
      data.table::setnames(bed, old, paste0("V", seq_along(old)))
    }
  }
  
  if (ncol(bed) < 3) stop("ld block.bed 至少需要 3 列（chr,start,end）")
  
  # 仅重命名前 3 列，避免出现 NA new-name
  data.table::setnames(
    bed,
    old = paste0("V", 1:3),
    new = c("CHR","START0","END0"),
    skip_absent = TRUE
  )
  
  # 第 4 列如果存在则作为块名，否则生成 1..n
  if ("V4" %in% names(bed)) {
    bed[, BLK := as.character(V4)]
  } else if ("name" %in% names(bed)) {           # 兼容少数带 name 列的情况
    bed[, BLK := as.character(name)]
  } else {
    bed[, BLK := as.character(seq_len(.N))]
  }
  
  # 规范 chr：去掉"chr"，X/Y/MT→23/24/25
  bed[, CHR := gsub("(?i)^chr", "", as.character(CHR))]
  bed[tolower(CHR)=="x",  CHR:="23"]
  bed[tolower(CHR)=="y",  CHR:="24"]
  bed[tolower(CHR) %in% c("mt","m"), CHR:="25"]
  suppressWarnings(bed[, CHR := as.integer(CHR)])
  
  # 0-based half-open -> 1-based closed
  bed[, START := as.integer(START0) + 1L]
  bed[, END   := as.integer(END0)]
  bed <- bed[is.finite(CHR) & is.finite(START) & is.finite(END) & END >= START]
  
  bed[, .(CHR, START, END, BLK)]
}
parse_pos_to_bp <- function(x){
  if (is.null(x) || length(x) == 0) return(NA_integer_)
  s <- as.character(x)[1]
  # "12.34 Mb" -> bp
  if (grepl("Mb", s, ignore.case = TRUE)) {
    v <- suppressWarnings(as.numeric(sub(" .*", "", s)))
    return(as.integer(v * 1e6))
  }
  # 纯数字 -> bp
  suppressWarnings(as.integer(s))
}
collapse_blocks_window <- function(
    dt,
    pthr        = 5e-8,
    window_kb   = 500,
    min_h2      = -Inf,
    top_per_chr = Inf,
    top_total   = Inf
){
  x <- data.table::as.data.table(dt)
  x <- x[is.finite(p) & p <= pthr & is.finite(chr) & nzchar(pos)]
  if (!nrow(x)) return(x[0])
  
  # 自动识别 SNP 列
  snp_col <- if ("lead snp" %in% names(x)) "lead snp" else if ("snp" %in% names(x)) "snp" else NA_character_
  if (is.na(snp_col)) stop("collapse_blocks_window: 未找到 SNP 列（支持 'lead snp' 或 'snp'）")
  
  # 解析区间中心
  posbp <- .parse_blk_pos_bp(x$pos)
  x[, `:=`(start_bp = posbp$start_bp, end_bp = posbp$end_bp, center_bp = posbp$center_bp)]
  x <- x[is.finite(center_bp)]
  data.table::setorder(x, p)
  
  keep <- rep(FALSE, nrow(x))
  
  for (i in seq_len(nrow(x))) {
    if (isTRUE(keep[i])) next
    if (is.na(keep[i]))  next
    
    keep[i] <- TRUE
    idx_same <- which(x$chr == x$chr[i] & isFALSE(keep) & !is.na(keep))
    if (length(idx_same)) {
      dist <- abs(x$center_bp[idx_same] - x$center_bp[i])
      keep[idx_same[dist < window_kb*1000]] <- NA
    }
  }
  
  # —— 需要的列：优先 blk_h5，其次 blk —— 
  blk_col <- if ("blk_h5" %in% names(x)) "blk_h5" else if ("blk" %in% names(x)) "blk" else NA_character_
  if (is.na(blk_col)) stop("collapse_blocks_window: 未找到 block 列（blk_h5 或 blk）")
  
  base_cols <- c(blk_col, snp_col, "chr", "pos", "p", "h2", "se")
  base_cols <- intersect(base_cols, names(x))
  
  out <- x[keep %in% TRUE, ..base_cols]
  data.table::setnames(out, snp_col, "lead snp", skip_absent = TRUE)
  
  # ① h2 阈
  if (is.finite(min_h2)) out <- out[is.finite(h2) & h2 >= min_h2]
  # ② 每条染色体上限
  if (is.finite(top_per_chr)) { out[, ord := seq_len(.N), by = chr]; out <- out[ord <= top_per_chr][, ord := NULL] }
  # ③ 全局上限
  if (is.finite(top_total) && nrow(out) > top_total) out <- out[order(p, -h2)][seq_len(top_total)]
  
  # —— 最终排序：chr + blk_col（数字化）——
  suppressWarnings(out[, (blk_col) := as.integer(get(blk_col))])
  if (all(c("chr", blk_col) %in% names(out))) {
    data.table::setorderv(out, c("chr", blk_col))
  } else {
    data.table::setorderv(out, blk_col)
  }
  
  # —— 严格锁定最终 7 列（blk_h5 优先；没有就用 blk）——
  if (blk_col != "blk_h5" && "blk_h5" %in% names(out)) blk_col <- "blk_h5"
  want <- c(blk_col, "lead snp", "chr", "pos", "p", "h2", "se")
  out <- out[, want, with = FALSE]
  
  out[]
  
}

.chr_cache_paths_abs <- function(gwas_name, build_disp, anc, chr_i){
  cache_dir <- file.path(CACHE_ROOT, gwas_name)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- sprintf("%s_%s_chr%02d", build_disp, anc, as.integer(chr_i))
  list(tsv = file.path(cache_dir, paste0(stem, ".tsv")))
}
.all_blocks_raw_path <- function(gwas_name, build_disp, anc){
  cache_dir <- file.path(CACHE_ROOT, gwas_name)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(cache_dir, sprintf("all_blocks_raw_%s_%s.tsv", build_disp, anc))
}
.merge_chr_blocks_raw <- function(gwas_name, build_disp, anc){
  files <- vapply(1:22, function(i){
    .chr_cache_paths_abs(gwas_name, build_disp, anc, i)$tsv
  }, character(1))
  
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (!length(files)) return(invisible(NULL))
  
  dt_all <- data.table::rbindlist(
    lapply(files, function(f) data.table::fread(f, sep = "\t", header = TRUE)),
    fill = TRUE
  )
  # 保险排序（你每个 chr 已经 setorder 过，但合并后再排一次更干净）
  if (all(c("chr","blk_h5") %in% names(dt_all))) data.table::setorder(dt_all, chr, blk_h5)
  
  out_path <- .all_blocks_raw_path(gwas_name, build_disp, anc)
  .atomic_write_tsv(dt_all, out_path)
  invisible(out_path)
}
.MEMO <- new.env(parent = emptyenv())
.read_cached_tsv <- function(path, mapping = NULL, ...) {
  if (!file.exists(path)) return(NULL)
  
  # 如果没有传 mapping（兼容旧调用），则直接读全表
  if (is.null(mapping)) {
    return(data.table::fread(path, header = TRUE, data.table = TRUE))
  }
  
  # 获取用户选定的列
  raw_cols <- as.character(unlist(mapping))
  
  dt <- tryCatch(
    data.table::fread(
      file = path,
      header = TRUE,
      data.table = TRUE,
      showProgress = FALSE,
      select = raw_cols
    ),
    error = function(e) NULL
  )
  
  if (is.null(dt)) return(NULL)
  
  # 统一列名
  data.table::setnames(dt, old = raw_cols, new = names(mapping))
  
  # 清洗
  dt <- dt[!is.na(p) & p > 0]
  if ("chr" %in% names(dt)) dt[, chr := as.character(chr)]
  
  return(dt)
}
.atomic_write_tsv <- function(dt, path){
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- sprintf("%s.tmp_%s_%d", path, Sys.getpid(), sample.int(.Machine$integer.max, 1))
  readr::write_tsv(dt, tmp)
  stopifnot(file.exists(tmp), file.size(tmp) > 0)
  file.rename(tmp, path)
  invisible(TRUE)
}
.file_nonempty <- function(p){
  is.character(p) && length(p)==1 && nzchar(p) &&
    file.exists(p) && isTRUE(try(file.info(p)$size > 0, silent = TRUE))
}

# ===== 绝对路径根 =====
Sys.setenv(SHINY_APP_ROOT = "D:/hello")
APP_ROOT <- tryCatch({
  r <- Sys.getenv("SHINY_APP_ROOT", unset = "")
  if (!nzchar(r)) r <- getwd()
  normalizePath(r, winslash = "/", mustWork = TRUE)
}, error = function(e) normalizePath(getwd(), winslash = "/", mustWork = FALSE))
DATA_GWAS_DIR <- file.path(APP_ROOT, "data", "gwas")
FILES_DIR     <- file.path(APP_ROOT, "files")
LD_BASE_DIR   <- file.path(APP_ROOT, "data", "ldblk")
CACHE_ROOT    <- file.path(APP_ROOT, "cache")
LD_BLOCK_BED <- file.path(FILES_DIR, "ld_block.bed")

IMG_CACHE_ROOT <- CACHE_ROOT
dir.create(CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)

.save_png_atomic <- function(path, width, height, res = 144, draw) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- sprintf("%s.tmp_%s_%d", path, Sys.getpid(), sample.int(.Machine$integer.max, 1))
  ok <- FALSE
  png(tmp, width = width, height = height, res = res)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    if (ok && file.exists(tmp) && file.size(tmp) > 0) {
      file.rename(tmp, path)
    } else {
      unlink(tmp, force = TRUE)
    }
  }, add = TRUE)
  draw()
  ok <- TRUE
  invisible(TRUE)
}

.SNPLIST_CACHE <- new.env(parent = emptyenv())
.get_snplist_from_h5file <- function(h5file){
  key <- normalizePath(h5file, winslash = "/", mustWork = FALSE)
  if (exists(key, envir = .SNPLIST_CACHE, inherits = FALSE)) {
    return(get(key, envir = .SNPLIST_CACHE, inherits = FALSE))
  }
  grps <- .get_h5_grps(h5file)  # 只返回同时具备 snplist/ldblk 的 group :contentReference[oaicite:2]{index=2}
  if (!length(grps)) {
    assign(key, character(0), envir = .SNPLIST_CACHE)
    return(character(0))
  }
  lst <- lapply(grps, function(g) {
    nm <- paste0(g, "/snplist")
    x  <- tryCatch(.h5_read_safe(h5file, nm), error = function(e) NULL)
    as.character(x)
  })
  snps <- unique(unlist(lst, use.names = FALSE))
  snps <- snps[!is.na(snps) & nzchar(snps)]
  assign(key, snps, envir = .SNPLIST_CACHE)
  snps
}
.get_union_snplist_from_h5map <- function(h5map){
  # h5map 可能是 list / named vector；这里尽量兼容
  files <- unique(as.character(unlist(h5map, use.names = FALSE)))
  files <- files[file.exists(files)]
  if (!length(files)) return(character(0))
  unique(unlist(lapply(files, .get_snplist_from_h5file), use.names = FALSE))
}

collapse_to_loci <- function(dt, p_thr = 5e-8, window_kb = 250) {
  # ① 列名统一（兼容小写）
  if ("chr" %in% names(dt)) data.table::setnames(dt, "chr", "CHR", skip_absent = TRUE)
  if ("pos" %in% names(dt)) data.table::setnames(dt, "pos", "POS", skip_absent = TRUE)
  if ("p"   %in% names(dt)) data.table::setnames(dt, "p",   "P",   skip_absent = TRUE)
  
  # ② 兜底：阈值与窗口
  p_thr     <- suppressWarnings(as.numeric(p_thr))
  if (!is.finite(p_thr) || p_thr <= 0) p_thr <- 5e-8
  window_kb <- suppressWarnings(as.integer(window_kb))
  if (!is.finite(window_kb) || window_kb <= 0) window_kb <- 250L
  kbp <- as.integer(window_kb) * 1000L
  
  # ③ 强制数值化（避免字符/空串造成 NA）
  DT <- data.table::as.data.table(dt)
  if ("CHR" %in% names(DT)) DT[, CHR := suppressWarnings(as.integer(CHR))]
  if ("POS" %in% names(DT)) DT[, POS := suppressWarnings(as.integer(POS))]
  if ("P"   %in% names(DT)) {
    # 科学计数/逗号 → 数值；并做你原来的上下限裁剪
    Pnum <- suppressWarnings(as.numeric(gsub(",", "", as.character(DT[["P"]]))))
    Pnum[!is.finite(Pnum) | Pnum <= 0] <- NA_real_
    Pnum[Pnum >= 1] <- 1 - 1e-16
    Pnum <- pmax(Pnum, 1e-50)
    DT[, P := Pnum]
  }
  
  # ④ 过滤：必须是有限 CHR/POS/P，且 P<=阈值
  DT <- DT[is.finite(CHR) & is.finite(POS) & is.finite(P) & P <= p_thr,
           .(SNP, CHR = as.integer(CHR), POS = as.integer(POS),
             P = as.numeric(P),
             H2 = if ("h2" %in% names(DT)) as.numeric(h2) else NA_real_)]
  
  if (!nrow(DT)) return(DT[0])
  
  data.table::setorder(DT, P, CHR, POS)
  
  out <- data.table::data.table(SNP=character(), CHR=integer(),
                                POS=integer(), P=double(), H2=double())
  
  while (nrow(DT)) {
    lead <- DT[1]
    out  <- rbind(out, lead)
    
    left  <- lead$POS - kbp
    right <- lead$POS + kbp
    
    DT <- DT[!(CHR == lead$CHR & POS >= left & POS <= right)]
  }
  
  unique(out, by = c("CHR","POS","SNP"))[]
}
collapse_to_independent_loci <- function(dt, ldb = NULL, pthr = 5e-8, window_kb = 500) {
  collapse_to_loci(dt, p_thr = pthr, window_kb = window_kb)
}

# ----------------------------【Server 主体】------------------------------------
server <- function(input, output, session){
  
  gwas_data_rv <- reactiveVal(NULL)
  genes_gr_rv <- reactiveVal(NULL)
  
  # 修改后的基因数据加载函数：使用 .gtf.gz 文件
  .load_gene_gr <- function(f) {
    if (grepl("\\.gz$", f, ignore.case = TRUE)) {
      # 直接读取 .gtf.gz 文件
      gr <- tryCatch({
        rtracklayer::import(f)
      }, error = function(e) {
        # 如果加载失败，输出错误信息并返回空
        message("Gene annotation loading failed: ", e$message)
        return(NULL)
      })
      
      if (is.null(gr)) {
        message("Gene annotation not ready. Please check the file format or contents.")
        return(NULL)
      }
      
      gr <- gr[gr$type == "gene"]
      GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
      gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
      gr <- GenomicRanges::trim(gr)
      
      # 提取基因名称
      gname <- if (!is.null(gr$gene_name)) gr$gene_name else if (!is.null(gr$Name)) gr$Name else gr$gene_id
      mcols(gr)$gene <- as.character(gname)
      
      return(gr)
    } else {
      message("File is not in .gtf.gz format. Please provide a valid .gtf.gz file.")
      return(NULL)
    }
  }
  
  # 在 server 中检查基因注释，确保使用 .gtf.gz 文件
  observeEvent(input$build, {
    build_disp <- .norm_build_disp(input$build)
    base <- file.path(
      FILES_DIR,
      paste0("gencode.v", if (build_disp == "37") "19" else "38", ".annotation.gtf")
    )
    
    f <- if (file.exists(paste0(base, ".gz"))) paste0(base, ".gz") else base
    
    if (!file.exists(f)) {
      genes_gr_rv(NULL)
      showNotification("未找到基因注释文件", type = "error")
      return()
    }
    
    withProgress(message = "Loading Gene DB...", {
      gene_gr <- .load_gene_gr(f)
      if (is.null(gene_gr)) {
        genes_gr_rv(NULL)
        showNotification("基因注释加载失败，请检查文件格式或路径", type = "error")
      } else {
        genes_gr_rv(gene_gr)
      }
    })
  }, ignoreInit = FALSE)
  
  .autodetect_mapping <- function(cols){
    cu <- toupper(trimws(sub("^#","", cols)))
    pick <- function(cands){
      hit <- match(toupper(cands), cu)
      hit <- hit[!is.na(hit)]
      if(length(hit)) cols[ hit[1] ] else NA_character_
    }
    list(
      SNP = pick(c("SNP","MARKERNAME","ID","RSID","VARIANT","MARKER","VARIANT_ID")),
      CHR = pick(c("CHR","CHROM","#CHROM","CHROMOSOME")),
      POS = pick(c("POS","BP","POSITION","BASE_PAIR_LOCATION","BP_HG19","POS37")),
      P   = pick(c("P","PVAL","P_VALUE","PVALUE","NEG_LOG_PVALUE","-LOG10P","LOGP","MLP"))
    )
  }
  .read_gwas_with_mapping <- function(path, map){
    stopifnot(is.list(map), all(c("SNP","CHR","POS","P") %in% names(map)))
    dt <- data.table::fread(
      file = path, header = TRUE, data.table = TRUE, fill = TRUE,
      na.strings = c("NA","NaN","","Inf","-Inf",".")
    )
    stopifnot(is.data.frame(dt), nrow(dt) > 0)
    
    # 原始列名 -> 不动；下面只通过 map$* 定位真实列
    nm <- names(dt)
    req_cols <- unlist(map)
    if(any(!req_cols %in% nm))
      stop("映射指向的列在文件中不存在：", paste(setdiff(req_cols, nm), collapse=", "))
    
    # 取出四列，并做必要转换
    snp <- dt[[ map$SNP ]]
    chr <- dt[[ map$CHR ]]
    pos <- dt[[ map$POS ]]
    pv  <- dt[[ map$P   ]]
    
    # CHR/POS 正规化
    chr <- gsub("(?i)^chr", "", as.character(chr))
    chr[tolower(chr)=="x"] <- "23"
    chr[tolower(chr)=="y"] <- "24"
    chr[tolower(chr) %in% c("mt","m")] <- "25"
    suppressWarnings(chr <- as.integer(chr))
    suppressWarnings(pos <- as.integer(pos))
    
    # SNP：只清空空白；不强行改名
    snp <- as.character(snp)
    
    # ---------- 改这里：P 值鲁棒识别 ----------
    pv_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(pv))))
    name_hint  <- grepl("NEG|LOG|MLP|LOG10", map$P, ignore.case = TRUE)
    finite     <- is.finite(pv_num)
    maxv       <- suppressWarnings(max(pv_num[finite], na.rm = TRUE))
    q99        <- suppressWarnings(stats::quantile(pv_num[finite], 0.99, na.rm = TRUE))
    frac_gt1   <- mean(pv_num[finite] > 1, na.rm = TRUE)
    
    if (is.finite(maxv) && maxv <= 1 + 1e-12) {
      # 明确是概率 P（0~1）
      P <- pv_num
    } else if (name_hint || (is.finite(q99) && q99 > 3) || frac_gt1 > 0.10) {
      # 明显像 -log10P
      P <- 10^(-pv_num)
    } else {
      # 默认按 P 处理
      P <- pv_num
    }
    
    # 清洗 + 下限裁切（与你原逻辑一致）
    P[!is.finite(P) | P <= 0] <- NA_real_
    P[P >= 1] <- 1 - 1e-16
    P <- pmax(P, 1e-50)
    # ---------- 到此为止 ----------
    
    out <- data.table::as.data.table(dt)     # 保留原始其他列
    out[, `:=`(SNP = snp, CHR = chr, POS = pos, P = P)]
    out <- out[is.finite(CHR) & is.finite(POS)]
    # 若缺 CHR/POS，尝试从 SNP 解析
    if (!nrow(out) || any(!is.finite(out$CHR)) || any(!is.finite(out$POS))) {
      s <- as.character(snp)
      if (any(grepl(":", s))) {
        sp <- data.table::tstrsplit(gsub("(?i)^CHR","", s), "[:_]", keep = 1:2)
        out[, CHR := suppressWarnings(as.integer(sp[[1]]))]
        out[, POS := suppressWarnings(as.integer(sp[[2]]))]
      }
    }
    out <- out[is.finite(CHR) & is.finite(POS)]
    out
  }
  h2_total_rv <- reactiveVal(NULL)   # list(h2=..., se=...)
  
  applied_once <- reactiveVal(FALSE)
  observeEvent(input$apply_block, {
    applied_once(TRUE)
  }, ignoreInit = TRUE)
  observeEvent(list(input$gwas_file, input$gwas_builtin), {
    path <- if (!is.null(input$gwas_file)) input$gwas_file$datapath else if(nzchar(input$gwas_builtin)) input$gwas_builtin else NULL
    req(path)
    
    applied_once(FALSE)
    gwas_data_rv(NULL) # 这里是变量赋值，不加括号
    
    # 获取列名
    hdr <- names(data.table::fread(path, nrows = 0, header = TRUE))
    
    # 更新下拉框
    updateSelectInput(session, "col_snp", choices = hdr, selected = .pick(hdr, c("SNP", "rsid", "MarkerName")))
    updateSelectInput(session, "col_chr", choices = hdr, selected = .pick(hdr, c("CHR", "CHROM", "Chromosome")))
    updateSelectInput(session, "col_pos", choices = hdr, selected = .pick(hdr, c("POS", "BP", "Position")))
    updateSelectInput(session, "col_p",   choices = hdr, selected = .pick(hdr, c("P", "PVAL", "p-value")))
  }, ignoreInit = TRUE)
  observeEvent(input$apply_selection, {
    path <- if (!is.null(input$gwas_file)) input$gwas_file$datapath else input$gwas_builtin
    req(path)
    
    # 获取用户确认的选择
    user_map <- list(
      snp = input$col_snp, chr = input$col_chr,
      pos = input$col_pos, p = input$col_p
    )
    
    withProgress(message = '正在解析 GWAS 数据...', value = 0.5, {
      dat <- .read_cached_tsv(path, user_map)
      gwas_data_rv(dat)
    })
    showNotification("列名映射已应用！数据加载完成。", type = "message")
  })
  
  catalog_r <- shiny::reactiveFileReader(
    5000, session,
    file.path(FILES_DIR, "gwas_catalog.tsv"),
    function(p){
      cat <- data.table::fread(p, sep="\t", header=TRUE)
      stopifnot(all(c("info","file") %in% names(cat)))
      cat[, file := normalizePath(file.path(APP_ROOT, file), winslash="/", mustWork=FALSE)]
      cat
    }
  )
  observe({
    cat <- catalog_r()
    lbl <- cat$info
    choices <- stats::setNames(cat$file, lbl)
    
    # 默认选中 BMI（按 info 字段匹配，不区分大小写）
    i <- grep("bmi", cat$info, ignore.case = TRUE)
    sel <- if (length(i)) cat$file[i[1]] else cat$file[1]
    
    updateSelectInput(session, "gwas_builtin", choices = choices, selected = sel)
  })
  
  observeEvent(input$gwas_builtin, {
    req(is.character(input$gwas_builtin), nzchar(input$gwas_builtin))
    gw  <- .resolve_gwas(input$gwas_file, input$gwas_builtin)
    bld <- .norm_build_disp(isolate(input$build))
    pop <- .norm_pop(isolate(input$ancestry))
    if (is.data.frame(loc <- .load_loci_from_cache(gw$name, bld, pop)))  loci_rv(loc)
    if (is.data.frame(blk <- .load_blocks_from_cache(gw$name, bld, pop))) blocks_rv(blk)
  }, once = TRUE, ignoreInit = FALSE)
  observeEvent(list(input$gwas_file, input$gwas_builtin), {
    gw <- .resolve_gwas(input$gwas_file, input$gwas_builtin)
    # 只读表头
    hdr <- tryCatch(names(data.table::fread(gw$path, nrows = 0L, header = TRUE)), error = function(e) NULL)
    if (is.null(hdr) || !length(hdr)) {
      updateSelectInput(session, "map_snp", choices = character(0))
      updateSelectInput(session, "map_chr", choices = character(0))
      updateSelectInput(session, "map_pos", choices = character(0))
      updateSelectInput(session, "map_p",   choices = character(0))
      return(invisible(NULL))
    }
    det <- .autodetect_mapping(hdr)
    updateSelectInput(session, "map_snp", choices = hdr, selected = det$SNP)
    updateSelectInput(session, "map_chr", choices = hdr, selected = det$CHR)
    updateSelectInput(session, "map_pos", choices = hdr, selected = det$POS)
    updateSelectInput(session, "map_p",   choices = hdr, selected = det$P)
  }, ignoreInit = FALSE)
  message("[init] WD = ", getwd())
  message("[init] CACHE_ROOT = ", CACHE_ROOT)
  
  # --------- 工具函数：统一用绝对路径 ----------
  .img_cache_for <- function(gwas_name){
    if (is.null(gwas_name) || !nzchar(gwas_name)) gwas_name <- "unknown"
    dir <- file.path(IMG_CACHE_ROOT, gwas_name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    list(
      dir = dir,
      man = file.path(dir, "manhattan.png"),
      qq  = file.path(dir, "qq.png"),
      locuszoom = file.path(dir, "locuszoom.png"),
      blockzoom = file.path(dir, "blockzoom.png")
    )
  }
  .file_ok <- function(p){
    is.character(p) && length(p)==1 && nzchar(p) &&
      file.exists(p) && isTRUE(try(file.info(p)$size > 0, silent = TRUE))
  }

  .cache_paths_abs <- function(gwas_name, build_disp, pop){
    cache_dir  <- file.path(CACHE_ROOT, gwas_name)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    list(
      loci   = file.path(cache_dir, sprintf("significant loci_%s_%s.tsv",   build_disp, pop)),
      blocks = file.path(cache_dir, sprintf("significant blocks_%s_%s.tsv", build_disp, pop)),
      blocks_from_loci = file.path(cache_dir, sprintf("blocks_from_loci_%s_%s.tsv", build_disp, pop))
    )
  }
  .norm_build_disp <- function(x){
    x <- tolower(as.character(x))
    if (x %in% c("grch37","hg19","37")) return("37")
    if (x %in% c("grch38","hg38","38")) return("38")
    gsub("[^0-9a-z]+","", x)
  }
  .norm_build_hdl <- function(disp){ if (disp=="37") "GRCh37" else if (disp=="38") "GRCh38" else disp }
  .norm_pop <- function(x){
    x <- toupper(as.character(x)); key <- gsub("[^A-Z]","", x)
    map <- c(EUR="EUR", EUROPEAN="EUR", EAS="EAS", EASTASIAN="EAS",
             AFR="AFR", AFRICAN="AFR", SAS="SAS", SOUTHASIAN="SAS",
             AMR="AMR", AMERICAN="AMR", MENA="MENA", OCE="OCE")
    if (!is.null(map[[key]])) map[[key]] else gsub("[^A-Z0-9]","", x)
  }
  .resolve_gwas <- function(input_file, builtin) {
    .safe_name <- function(x) {
      x <- gsub("[^A-Za-z0-9._-]+", "_", x)
      x <- gsub("_+", "_", x)
      x <- sub("^_+", "", x)
      x <- sub("_+$", "", x)
      if (!nzchar(x)) "unknown" else x
    }
    
    # 1) 优先使用用户上传的文件
    is_up <- !is.null(input_file) && is.list(input_file) &&
      "datapath" %in% names(input_file) &&
      length(input_file$datapath) >= 1 &&
      nzchar(as.character(input_file$datapath[[1]]))
    
    if (is_up) {
      p  <- normalizePath(as.character(input_file$datapath[[1]]),
                          winslash = "/", mustWork = TRUE)
      nm <- if (!is.null(input_file$name) && length(input_file$name) >= 1)
        as.character(input_file$name[[1]]) else basename(p)
      
      return(list(
        path = p,
        name = .safe_name(tools::file_path_sans_ext(basename(nm)))
      ))
    }
    
    # 2) 否则使用内置 GWAS（来自 builtin 选择）
    req(builtin)
    
    # 强制拿到一个纯字符路径
    f_chr <- as.character(builtin)[1]
    
    # normalizePath 失败时退回原字符
    f <- tryCatch(
      normalizePath(f_chr, winslash = "/", mustWork = FALSE),
      error = function(e) f_chr
    )
    
    # 不再用 validate/need，直接用普通错误（保证是纯字符）
    if (!file.exists(f)) {
      stop(sprintf("GWAS 文件不存在：%s", f), call. = FALSE)
    }
    
    list(
      path = f,
      name = .safe_name(tools::file_path_sans_ext(basename(f)))
    )
  }
  
  loci_rv   <- reactiveVal(NULL)
  blocks_rv <- reactiveVal(NULL)
  .load_loci_from_cache <- function(gwas_name, build_disp, pop){
    paths <- .cache_paths_abs(gwas_name, build_disp, pop)
    if (!file.exists(paths$loci)) return(NULL)
    
    hdr <- tryCatch(
      names(data.table::fread(paths$loci, nrows = 0L, showProgress = FALSE)),
      error = function(e) character(0)
    )
    if (!length(hdr)) return(NULL)
    
    # 不区分大小写匹配：返回“实际表头名”
    pick <- function(...) {
      alts <- unlist(list(...))
      if (!length(alts)) return(NA_character_)
      hdrU <- toupper(trimws(hdr))
      hit  <- match(toupper(trimws(alts)), hdrU)
      hit  <- hit[!is.na(hit)]
      if (length(hit)) hdr[ hit[1] ] else NA_character_
    }
    
    # 需要的列（给出别名备选；返回“实际表头名”）
    cols <- c(
      pick("loc"),
      pick("lead snp","lead_snp","snp","SNP"),
      pick("chr","CHR"),
      pick("pos","POS"),
      pick("p","P","pval","p_value","pvalue","NEG_LOG_PVALUE"),
      pick("h2","H2"),
      pick("nearest gene","NEAREST GENE","nearest_gene"),
      pick("blk","BLK"),
      pick("blk_h5","BLK_H5","blk h5"),
      pick("se","SE")
    )
    sel <- unique(cols[!is.na(cols)])
    if (!length(sel)) return(NULL)
    
    # 仅对“实际存在”的列声明类型（名字必须和 sel 完全一致）
    cc <- c()
    if (!is.na(pick("chr","CHR"))) cc[pick("chr","CHR")] <- "integer"
    if (!is.na(pick("p","P","pval","p_value","pvalue","NEG_LOG_PVALUE"))) cc[pick("p","P","pval","p_value","pvalue","NEG_LOG_PVALUE")] <- "numeric"
    if (!is.na(pick("h2","H2"))) cc[pick("h2","H2")] <- "numeric"
    
    dt <- .read_cached_tsv(paths$loci, sep = "\t", select = sel, colClasses = cc)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    dt[]
  }
  .load_blocks_from_cache <- function(gwas_name, build_disp, anc){
    paths <- .cache_paths_abs(gwas_name, build_disp, anc)
    if (!file.exists(paths$blocks)) return(NULL)
    
    hdr <- tryCatch(
      names(data.table::fread(paths$blocks, nrows = 0L, showProgress = FALSE)),
      error = function(e) character(0)
    )
    if (!length(hdr)) return(NULL)
    
    pick <- function(...) {
      alts <- unlist(list(...))
      if (!length(alts)) return(NA_character_)
      hdrU <- toupper(trimws(hdr))
      hit  <- match(toupper(trimws(alts)), hdrU)
      hit  <- hit[!is.na(hit)]
      if (length(hit)) hdr[hit[1]] else NA_character_
    }
    
    lead_col    <- pick("lead snp","lead_snp","snp","SNP")
    nearest_col <- pick("nearest gene","NEAREST GENE","nearest_gene")
    
    cols <- c(
      pick("blk","BLK"),
      lead_col,
      pick("chr","CHR"),
      pick("pos","POS"),
      pick("p","P","pval","p_value","pvalue","NEG_LOG_PVALUE"),
      pick("h2","H2"),
      nearest_col,                      # ← 新增
      pick("blk_h5","BLK_H5","blk h5"),
      pick("se","SE")
    )
    sel <- unique(cols[!is.na(cols)])
    if (!length(sel)) return(NULL)
    
    cc <- c()
    if (!is.na(pick("blk","BLK"))) cc[pick("blk","BLK")] <- "integer"
    if (!is.na(lead_col))          cc[lead_col]          <- "character"
    if (!is.na(pick("chr","CHR"))) cc[pick("chr","CHR")] <- "integer"
    if (!is.na(pick("pos","POS"))) cc[pick("pos","POS")] <- "character"
    if (!is.na(pick("p","P","pval","p_value","pvalue","NEG_LOG_PVALUE")))
      cc[pick("p","P","pval","p_value","pvalue","NEG_LOG_PVALUE")] <- "numeric"
    if (!is.na(pick("h2","H2")))   cc[pick("h2","H2")]   <- "numeric"
    if (!is.na(pick("blk_h5","BLK_H5","blk h5")))
      cc[pick("blk_h5","BLK_H5","blk h5")] <- "integer"
    
    dt <- .read_cached_tsv(paths$blocks, sep = "\t", select = sel, colClasses = cc)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    
    if (!"lead snp" %in% names(dt) && !is.na(lead_col)) {
      data.table::setnames(dt, lead_col, "lead snp", skip_absent = TRUE)
    }
    if ("se" %in% names(dt)) dt[, se := NULL]
    
    dt[]
  }
  
  observeEvent(list(input$gwas_file, input$gwas_builtin, input$build, input$ancestry), {
    gw  <- .resolve_gwas(input$gwas_file, input$gwas_builtin)
    bld <- .norm_build_disp(input$build)
    pop <- .norm_pop(input$ancestry)
    
    loc <- .load_loci_from_cache(gw$name, bld, pop)
    if (is.data.frame(loc)) loci_rv(loc)
    
    blk <- .load_blocks_from_cache(gw$name, bld, pop)
    if (is.data.frame(blk)) blocks_rv(blk)
  }, ignoreInit = FALSE)
  
  # 定义ldblk目录
  .norm_token <- function(x) toupper(gsub("[^A-Za-z0-9]+","",x))
  .find_chr_h5_in_dir <- function(dir){
    out <- rep(NA_character_, 22)
    for(i in 1:22){
      pat <- sprintf("chr[[:space:]]*%d.*\\.(hdf5|h5)$", i)
      hit <- list.files(dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
      if(length(hit)) out[i] <- hit[1]
    }
    names(out) <- paste0("chr",1:22); out
  }
  .top_peak <- function(dt_std){
    p <- suppressWarnings(as.numeric(dt_std$P))
    p[!is.finite(p)] <- Inf
    i <- which.min(p)
    list(chr = as.integer(dt_std$CHR[i]), pos = as.numeric(dt_std$POS[i]))
  }
  .center <- reactiveVal(NULL)
  resolve_ld_h5_files <- function(base_dir, build, anc){
    build_t <- .norm_token(build)  # “37/38/GRCh37/GRCh38”等都规整
    anc_t   <- .norm_token(anc)    # “EUR/EAS/AFR/...” 规整
    cand_dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
    # 同时匹配 build & ancestry 优先
    hit <- cand_dirs[ grepl(build_t,.norm_token(cand_dirs),fixed=TRUE) &
                        grepl(anc_t,  .norm_token(cand_dirs),fixed=TRUE) ]
    if(!length(hit)){
      hit <- cand_dirs[ grepl(build_t,.norm_token(cand_dirs),fixed=TRUE) |
                          grepl(anc_t,  .norm_token(cand_dirs),fixed=TRUE) ]
    }
    best <- NULL; score <- -1
    for (d in unique(as.character(hit))) {
      chrvec <- .find_chr_h5_in_dir(d)
      sc <- sum(!is.na(chrvec))
      if(sc > score){ best <- chrvec; score <- sc }
      if(sc == 22) break
    }
    if(score <= 0) best <- .find_chr_h5_in_dir(base_dir)
    best
  }
  
  # 1) GWAS 源（上传优先，否则选择 data/ 下文件）
  gw_current <- reactive({
    req(input$gwas_builtin)  # 确保 catalog 观察器已把默认 BMI 填进去
    .resolve_gwas(input$gwas_file, input$gwas_builtin)
  })
  gwas_raw <- reactive({
    gw <- gw_current()
    
    # 四个下拉至少选中 3 个以上才启用映射；否则回退到原有自动识别
    sel <- list(
      SNP = input$map_snp,
      CHR = input$map_chr,
      POS = input$map_pos,
      P   = input$map_p
    )
    have <- vapply(sel, function(x) is.character(x) && length(x)==1 && nzchar(x), logical(1))
    if (sum(have) >= 3 && all(have)) {
      # 尝试映射读取；失败则回退
      dt <- tryCatch(.read_gwas_with_mapping(gw$path, sel),
                     error = function(e) NULL)
      if (is.data.frame(dt) && nrow(dt) > 0) return(dt)
    }
    
    # 根据当前 build 选择 37/38 坐标
    bld <- .norm_build_disp(input$build)   # 返回 "37" 或 "38"
    
    # 回退到自动识别的 read_gwas，同时把 build 传进去
    read_gwas(gw$path, build_disp = bld)
  })
  
  # 2) Manhattan：点了 apply 才生成/更新
  output$manhattan <- renderImage({
    gw <- gw_current();  # 不再依赖 gw_chosen
    paths <- .img_cache_for(gw$name)
    
    # 1) 有缓存就直接读
    if (.file_ok(paths$man)) {
      return(list(src = paths$man, contentType = "image/png"))
    }
    
    # 2) 无缓存就生成一次并写入缓存
    dt <- gwas_raw()[, .(SNP, CHR, POS, P)]
    dt <- dt[CHR %in% 1:25 & is.finite(POS)]
    dt <- dt[is.finite(P)]                             # 丢掉非有限
    p_raw <- suppressWarnings(as.numeric(dt$P))
    p_raw[!is.finite(p_raw)] <- NA
    p_plot <- scale_p(p_raw, min_p = 1e-324, max_log10 = 50)
    p_plot[p_plot >= 1] <- 1 - 1e-16
    d_for_cm <- data.frame(SNP=dt$SNP, CHR=dt$CHR, POS=dt$POS, P=p_plot)
    
    .save_png_atomic(paths$man, width = 1100, height = 1100, res = 150, draw = function(){
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mai = c(0,0,0,0), omi = c(0,0,0,0), plt = c(0,1,0,1),
          xaxs = "i", yaxs = "i", xpd = NA, bg = "white")
      CMplot::CMplot(
        d_for_cm, type="p", plot.type="c", LOG10=TRUE,
        threshold=c(5e-8,1e-5), threshold.lty=c(1,2),
        threshold.col=c("red","grey40"),
        main = NULL,
        file.output=FALSE, verbose=FALSE
      )
    })
    list(src = paths$man, contentType = "image/png")
  }, deleteFile = FALSE)
  
  # 3) QQ：点了 apply 才生成/更新
  output$qq_plot <- renderImage({
    gw <- gw_current()
    paths <- .img_cache_for(gw$name)
    
    # 1) 有缓存就直接读
    if (.file_ok(paths$qq)) {
      return(list(src = paths$qq, contentType = "image/png"))
    }
    
    # 2) QQ plot 数据
    p_all <- suppressWarnings(as.numeric(gwas_raw()[["P"]]))
    p_raw <- p_all[is.finite(p_all)]
    p_raw[p_raw <= 0] <- 1e-300
    p_raw[p_raw >= 1] <- 1 - 1e-16
    
    lambda_all <- median(stats::qchisq(1 - p_raw, df = 1), na.rm = TRUE) / 0.456
    
    p_plot <- scale_p(p_raw, min_p = 1e-324, max_log10 = 50)
    p_plot <- sort(p_plot)
    
    n   <- length(p_plot)
    obs <- -log10(p_plot)
    expx <- -log10((1:n)/(n+1))
    df  <- data.frame(exp = expx, obs = obs)
    
    k   <- 4000L
    idx <- unique(as.integer(seq(1, n, length.out = k)))
    ci <- data.frame(
      exp = -log10((idx)/(n+1)),
      lo  = -log10(qbeta(0.975, idx, n - idx + 1)),
      hi  = -log10(qbeta(0.025, idx, n - idx + 1))
    )
    max_lim <- ceiling(max(df$exp, df$obs, na.rm = TRUE))
    
    # —— h²：只从 significant blocks_%s_%s.tsv 读取并求和 ——
    build_disp  <- .norm_build_disp(input$build)
    anc         <- .norm_pop(input$ancestry)
    cache_paths <- .cache_paths_abs(gw$name, build_disp, anc)
    
    h2_total <- NA_real_
    bt <- tryCatch(
      data.table::fread(cache_paths$blocks, sep = "\t", header = TRUE, data.table = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(bt) && "h2" %in% names(bt)) {
      h2_total <- sum(as.numeric(bt$h2), na.rm = TRUE)
    }
    
    # —— GWAS 整体 p：Fisher 合并检验（all SNPs）——
    p_for_global <- p_raw
    p_for_global[p_for_global <= 0] <- 1e-300
    p_for_global[p_for_global >= 1] <- 1 - 1e-16
    p_for_global <- p_for_global[is.finite(p_for_global)]
    
    p_global <- NA_real_
    if (length(p_for_global) >= 2) {
      chisq_global <- -2 * sum(log(p_for_global))
      df_global    <- 2 * length(p_for_global)
      p_global     <- stats::pchisq(chisq_global, df = df_global, lower.tail = FALSE)
    }
    
    p_show <- if (is.finite(p_global)) {
      if (p_global == 0) "P<1e-300"
      else if (p_global < 1e-6) "P<1e-6"
      else sprintf("P=%.2f", p_global)
    } else {
      "P=NA"
    }
    
    h2_label <- if (is.finite(h2_total)) {
      sprintf("\nh\u00B2=%.3f\n%s\n", h2_total, p_show)
    } else {
      NA_character_
    }
    
    # —— 画图 ——
    .save_png_atomic(paths$qq, width = 1200, height = 1200, res = 144, draw = function(){
      use_rast <- requireNamespace("ggrastr", quietly = TRUE)
      pt_layer <- if (use_rast) {
        ggrastr::geom_point_rast(data = df, ggplot2::aes(exp, obs),
                                 size = 0.25, alpha = 0.7)
      } else {
        ggplot2::geom_point(data = df, ggplot2::aes(exp, obs),
                            size = 0.25, alpha = 0.7)
      }
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = ci,
          ggplot2::aes(x = exp, ymin = lo, ymax = hi),
          fill = "grey92"
        ) +
        ggplot2::geom_abline(
          slope = 1, intercept = 0,
          linetype = "dashed", linewidth = 0.6, color = "grey50"
        ) +
        pt_layer +
        ggplot2::coord_fixed(ratio = 1, clip = "off") +
        ggplot2::scale_x_continuous(
          limits = c(0, max_lim),
          breaks = scales::pretty_breaks(7),
          expand = ggplot2::expansion(mult = c(0.01, 0.02))
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, max_lim),
          breaks = scales::pretty_breaks(7),
          expand = ggplot2::expansion(mult = c(0.01, 0.02))
        ) +
        ggplot2::labs(
          x = "Expected -log10(P)",
          y = "Observed -log10(P)",
          title = "QQ plot",
          subtitle = sprintf("\u03BB = %.3f, n = %s",
                             lambda_all, scales::number(n, big.mark = " "))
        ) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          plot.title       = ggplot2::element_text(face = "bold", hjust = .5),
          plot.subtitle    = ggplot2::element_text(color = "#6B7280", hjust = .5),
          axis.title       = ggplot2::element_text(face = "bold"),
          plot.margin      = grid::unit(c(6, 6, 6, 6), "pt")
        )
      
      if (is.character(h2_label) && nzchar(h2_label)) {
        p <- p + ggplot2::annotate(
          "label",
          # 坐标稍微往里收一点，给巨大的边框留出空间
          x = max_lim * 0.98,  
          y = max_lim * 0.02,  
          label = h2_label,
          hjust = 1, 
          vjust = 0, 
          size = 10,       
          
          # --- 关键：将 padding 设为 2 到 3 左右，配合上面的空格 ---
          label.padding = grid::unit(1.5, "lines"), 
          
          label.r = grid::unit(0.3, "lines"), # 圆角
          fill = "white",
          color = "black",
          linewidth = 0.8   # 增加边框线条粗细，使其视觉感更强
        )
      }
      print(p)
    })
    
    list(src = paths$qq, contentType = "image/png")
  }, deleteFile = FALSE)
  
  # 4) significant loci
  observeEvent(input$apply_block, {
    gw         <- .resolve_gwas(input$gwas_file, input$gwas_builtin)
    build_disp <- .norm_build_disp(input$build)
    pop        <- .norm_pop(input$ancestry)
    paths      <- .cache_paths_abs(gw$name, build_disp, pop)
    build_hdl <- .norm_build_hdl(.norm_build_disp(input$build))
    anc       <- .norm_pop(input$ancestry)
    h5map     <- resolve_ld_h5_files(LD_BASE_DIR, build_hdl, anc)
    snplist_all <- .get_union_snplist_from_h5map(h5map)
   
    # 命中缓存直接返回
    if (.is_nonempty_tsv(paths$loci)) {
      loci_rv(data.table::fread(paths$loci, sep="\t", header=TRUE, data.table=TRUE))
      return(invisible(NULL))
    }
    
    # —— 统一用 read_gwas() 的清洗结果 —— 
    thr <- .to_num(input$p_threshold, 5e-8)
    if (!is.finite(thr)) thr <- 5e-8
    dt  <- data.table::as.data.table(gwas_raw())   # 已含标准列 SNP/CHR/POS/P，且保留原始列
    
    # —— 列名匹配：大小写无关，返回“实际列名” —— 
    nm <- names(dt)
    mapU <- stats::setNames(nm, toupper(nm))
    pickU <- function(alts){
      hits <- mapU[toupper(alts)]
      hits <- hits[!is.na(hits)]
      if (length(hits)) hits[1] else NA_character_
    }
    
    # —— 识别效应列（BETA/EFFECT；若没有则用 OR 转 log(OR)）——
    b_col <- pickU(c("BETA","B","EFFECT","EFFECT_SIZE","LOG_OR","LOG(OR)","LOGODDS","BETA_SNP"))
    if (is.na(b_col)) {
      or_col <- pickU(c("OR","ODDSRATIO","ODDS_RATIO"))
      if (!is.na(or_col)) {
        dt[, BETA := log(suppressWarnings(as.numeric(get(or_col))))]
        b_col <- "BETA"
      }
    }
    
    # —— 识别频率列（支持 alt_allele_freq 等别名）——
    freq_col <- pickU(c(
      "EAF","MAF","AF","ALT_AF","ALT_ALLELE_FREQ","ALT_ALLELE_FRQ",
      "A1FREQ","A1F","ALLELE1_FRQ","FRQ","FREQ","FREQ1.HAPMAP","FREQ1"
    ))
    
    # —— 计算每个 SNP 的 h2 = 2 * MAF * (1-MAF) * BETA^2 —— 
    if (!is.na(b_col) && !is.na(freq_col)) {
      dt[, BETA := suppressWarnings(as.numeric(get(b_col)))]
      dt[, EAF  := suppressWarnings(as.numeric(get(freq_col)))]
      
      # 若频率像是百分比（中位数 > 1 且 ≤ 100），转成 0–1
      med_eaf <- suppressWarnings(stats::median(dt[['EAF']], na.rm = TRUE))
      if (is.finite(med_eaf) && med_eaf > 1 && med_eaf <= 100) dt[, EAF := EAF/100]
      
      # 夹取到 [0,1]，并生成 MAF
      dt[, EAF := pmin(pmax(EAF, 0), 1)]
      dt[, MAF := pmin(EAF, 1 - EAF)]
      
      dt[, h2  := 2 * MAF * (1 - MAF) * (BETA^2)]
    } else {
      dt[, h2 := NA_real_]
    }
    
    # 计算前做诊断日志（只 message，不影响 UI）
    msg_counts <- function(tag, dt, thr){
      p <- suppressWarnings(as.numeric(dt$P))
      n_all <- nrow(dt)
      n_fp  <- sum(is.finite(p))
      n_sig <- sum(is.finite(p) & p <= thr)
      message(sprintf("[loci/%s] rows=%s, finiteP=%s, P<=%g=%s",
                      tag, n_all, n_fp, thr, n_sig))
    }
    thr <- .to_num(input$p_threshold, 5e-8); if (!is.finite(thr)) thr <- 5e-8
    msg_counts("raw", dt, thr)
    
    dt <- dt[CHR %in% 1:25 & is.finite(P) & is.finite(POS)]
    msg_counts("filtered", dt, thr)
    
    if (length(snplist_all)) {
      dt <- dt[SNP %in% snplist_all]
      msg_counts("after_snplist", dt, thr)
    }

    lead <- collapse_to_independent_loci(dt, ldb = NULL, pthr = thr, window_kb = 500)
    data.table::setorder(lead, CHR, POS)
    
    res <- if (nrow(lead)) lead[, .(
      loc = .I,
      `lead snp` = SNP,
      chr = as.integer(CHR),
      pos = sprintf("%.2f Mb", as.numeric(POS)/1e6),
      p   = as.numeric(P),
      h2  = as.numeric(H2)   
    )] else data.table::data.table(
      loc=integer(), snp=character(), chr=integer(), pos=character(), p=numeric(), h2=numeric()
    )
    res[, p := scale_p(p, min_p = 1e-324, max_log10 = 50)]
    res[p >= 1, p := 1 - 1e-16]
    
    # ---- 追加最近基因列（放在 h2 后面）----
    # 先把 "xx.xx Mb" 转为 bp
    pos_bp <- suppressWarnings(as.numeric(sub(" .*", "", res$pos))) * 1e6
    
    # 染色体：把 23/24/25 映射回 X/Y/MT，且都是字符
    chr_chr <- as.character(res$chr)
    chr_chr[chr_chr %in% c("23","x","X")] <- "X"
    chr_chr[chr_chr %in% c("24","y","Y")] <- "Y"
    chr_chr[chr_chr %in% c("25","mt","MT","M")] <- "MT"
    
    # 取注释
    gref <- genes_gr_rv()
    nearest_gene <- rep(NA_character_, nrow(res))
    nearest_dist <- rep(NA_real_,       nrow(res))
    
    if (!is.null(gref) && nrow(res) > 0) {
      # 构造点位 GRanges
      suppressWarnings({
        pts <- GenomicRanges::GRanges(
          seqnames = chr_chr,
          ranges   = IRanges::IRanges(start = as.integer(pos_bp), width = 1)
        )
      })
      # 仅保留在注释中存在的染色体，避免 “不在序列里” 的警告
      pts <- GenomeInfoDb::keepSeqlevels(
        pts, intersect(GenomeInfoDb::seqlevels(pts), GenomeInfoDb::seqlevels(gref)),
        pruning.mode = "coarse"
      )
      
      if (length(pts) > 0) {
        hits <- GenomicRanges::distanceToNearest(pts, gref, ignore.strand = TRUE)
        if (length(hits) > 0) {
          qh <- S4Vectors::queryHits(hits)
          sh <- S4Vectors::subjectHits(hits)
          nearest_gene[qh] <- mcols(gref)$gene[sh]
          nearest_dist[qh] <- S4Vectors::mcols(hits)$distance
        }
      }
    } else {
      showNotification("Gene annotation not ready. Nearest gene will be empty.", type="warning", duration=6)
    }
    
    res[, `nearest gene` := nearest_gene]
    data.table::setcolorder(
      res,
      c("loc","lead snp","chr","pos","p","h2","nearest gene")
    )
    
    readr::write_tsv(res, paths$loci)
    loci_rv(res)
    
  }, ignoreInit = TRUE)
  output$top_table <- DT::renderDT({
    dat <- loci_rv()
    req(!is.null(dat) && nrow(dat) > 0, cancelOutput = TRUE)
    
    DT::datatable(
      dat,
      rownames = FALSE, escape = FALSE,
      class = "compact stripe hover match-rows",
      selection = "single",
      options = list(
        paging = TRUE, pageLength = 10, lengthChange = FALSE,
        searching = FALSE, info = FALSE, dom = "tp", autoWidth = FALSE
      ),
      width = "100%"
    ) |>
      DT::formatSignif(columns = "p", digits = 3) |>
      DT::formatRound(columns = "h2", digits = 4) |>
      DT::formatStyle(
        "lead snp",
        color = "#2563EB",     # 蓝色
        fontWeight = "bold"    # 可选：更醒目
      )
  })
  
  # 5) significant blocks
  observeEvent(input$apply_block, {
    
    # =========================
    # 0) 直接使用已生成的 loci 表
    # =========================
    loci <- data.table::as.data.table(loci_rv())
    req(!is.null(loci) && nrow(loci) > 0)
    
    # =========================
    # 1) loci -> loci4（用于 BED 映射）
    # =========================
    loci4 <- loci[, .(
      SNP = as.character(`lead snp`),
      CHR = suppressWarnings(as.integer(chr)),
      POS = suppressWarnings(as.integer(round(as.numeric(gsub(" Mb", "", pos)) * 1e6))),
      P   = suppressWarnings(as.numeric(p)),
      H2  = suppressWarnings(as.numeric(h2))
    )]
    loci4 <- loci4[is.finite(CHR) & is.finite(POS) & is.finite(P)]
    loci4[!is.finite(H2), H2 := 0]   # NA 的 h2 当 0
    req(nrow(loci4) > 0)
    
    # =========================
    # 2) loci → BED blocks 映射，并在每 blk 取最小 P，同时 h2 求和
    # =========================
    bed_file <- LD_BLOCK_BED
    bed <- .read_ld_blocks_bed(bed_file)
    bed <- data.table::as.data.table(bed)[, .(
      CHR   = suppressWarnings(as.integer(CHR)),
      START = suppressWarnings(as.integer(START)),
      END   = suppressWarnings(as.integer(END)),
      BLK   = suppressWarnings(as.integer(BLK))
    )]
    bed <- bed[is.finite(CHR) & is.finite(START) & is.finite(END) & is.finite(BLK)]
    
    m <- bed[loci4,
             on = .(CHR, START <= POS, END >= POS),
             nomatch = 0L,
             .(blk        = BLK,
               `lead snp` = i.SNP,
               chr        = i.CHR,
               pos_bp     = i.POS,
               p          = i.P,
               h2         = i.H2)
    ]

    if (nrow(m) == 0) {
      blk <- data.table::data.table(
        blk = integer(),
        `lead snp` = character(),
        chr = integer(),
        pos = character(),
        p = numeric(),
        h2 = numeric()
      )
    } else {
      
      # lead：每 blk 取 p 最小
      lead <- m[, .SD[which.min(p)], by = .(blk)]
      lead[, h2 := NULL]  # 避免后面和 h2sum 重名
      
      # h2：每 blk 求和
      h2sum <- m[, .(h2 = sum(h2, na.rm = TRUE)), by = .(blk)]
      
      # 先合并 lead + h2
      blk <- merge(lead, h2sum, by = "blk", all.x = TRUE)
      
      # 再从 bed 里把区间按 blk 合并回来（最稳）
      bed_rng <- bed[, .(blk = BLK,
                         start_bp = START,
                         end_bp   = END)]
      blk <- merge(blk, bed_rng, by = "blk", all.x = TRUE)
      
      # 用 bed 区间生成 pos
      blk[, pos := sprintf("%.2f–%.2f Mb", start_bp / 1e6, end_bp / 1e6)]
      
      # 清理
      blk[, c("start_bp","end_bp","pos_bp") := NULL]
      
      data.table::setorder(blk, chr, blk)
    }
    
    # ---- 把 loci 的 nearest gene 抄到 blk（按 lead snp 对齐）----
    gene_map <- loci[, .(`lead snp` = as.character(`lead snp`),
                         `nearest gene` = as.character(`nearest gene`))]
    
    blk <- merge(blk, gene_map, by = "lead snp", all.x = TRUE)

    data.table::setorder(blk, blk)
    data.table::setcolorder(blk, c("blk","lead snp","chr","pos","p","h2","nearest gene"))
    blocks_rv(blk)
    
    # =========================
    # 3) 写缓存（保持你原行为）
    # =========================
    paths <- .cache_paths_abs(
      .resolve_gwas(input$gwas_file, input$gwas_builtin)$name,
      .norm_build_disp(input$build),
      .norm_pop(input$ancestry)
    )
    .atomic_write_tsv(loci, paths$loci)
    .atomic_write_tsv(blk,  paths$blocks)
    
    showNotification("Blocks computed from existing loci (min-P per block, h2 summed within block).",
                     type = "message")
    
  }, ignoreInit = TRUE)
  output$signal_block_table <- DT::renderDT({
    dat <- blocks_rv()
    req(!is.null(dat) && nrow(dat) > 0, cancelOutput = TRUE)
    
    DT::datatable(
      dat,
      rownames = FALSE, escape = FALSE,
      class = "compact stripe hover match-rows",
      selection = "single",
      options = list(
        paging = TRUE, pageLength = 10, lengthChange = FALSE,
        searching = FALSE, info = FALSE, dom = "tp", autoWidth = FALSE
      ),
      width = "100%"
    ) |>
      DT::formatSignif(columns = "p", digits = 3) |>
      DT::formatRound(columns = "h2", digits = 4)
  })
  
  # 6）locuszoom
  .pick_default_lead_from_loci <- function(loci_dt) {
    loci_dt <- data.table::as.data.table(loci_dt)
    if (!("lead snp" %in% names(loci_dt))) stop("loci_dt 缺少列: lead snp")
    if (!("p" %in% names(loci_dt))) stop("loci_dt 缺少列: p")
    
    loci_dt[, p_num := suppressWarnings(as.numeric(p))]
    loci_dt <- loci_dt[is.finite(p_num)]
    if (nrow(loci_dt) == 0) stop("loci 表里没有可用 p")
    
    as.character(loci_dt[which.min(p_num), `lead snp`])
  }
  observeEvent(input$top_table_rows_selected, {
    req(gwas_raw())
    loci_dt <- data.table::as.data.table(loci_rv())
    sel <- input$top_table_rows_selected
    req(!is.null(loci_dt), nrow(loci_dt) > 0)
    if (is.null(sel) || length(sel) < 1) return(invisible(NULL))
    
    i <- suppressWarnings(as.integer(sel[1]))
    if (!is.finite(i) || i < 1 || i > nrow(loci_dt)) return(invisible(NULL))
    
    rowi <- loci_dt[i, ]
    
    lead0 <- as.character(rowi[["lead snp"]])
    chr0  <- suppressWarnings(as.integer(rowi[["chr"]]))
    
    pos0 <- rowi[["pos"]]
    if (is.character(pos0)) pos0 <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pos0)) * 1e6)
    pos0 <- as.integer(pos0)
    
    .center(list(chr = chr0, pos = pos0, title = "LocusZoom", dt = gwas_raw(), `lead snp` = lead0))
  }, ignoreInit = TRUE)
  output$locuszoom_plot <- renderImage({
    # --- 1. 显式声明依赖：确保点击表格时，这个块会重新执行 ---
    sel <- input$top_table_rows_selected
    cen <- .center() # .center 内部包含了点击后更新的 lead snp
    req(cen, cen[["lead snp"]])
    
    gw <- gw_current()
    # 【重要】修改缓存路径生成逻辑：文件名必须带上 lead snp，否则永远只读第一张缓存
    cache_dir <- file.path("cache", gw$name)
    if(!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    
    # 动态路径：例如 cache/my_gwas/locuszoom_rs12345.png
    specific_path <- file.path(cache_dir, paste0("locuszoom_", cen[["lead snp"]], ".png"))
    
    # --- 2. 检查这个特定 SNP 的缓存是否存在 ---
    if (file.exists(specific_path)) {
      return(list(src = specific_path, contentType = "image/png", width = "100%"))
    }
    
    # --- 3. 如果不存在，则绘图并保存 ---
    lead <- as.character(cen[["lead snp"]])
    anc <- .norm_pop(input$ancestry)
    h5_dir <- file.path(LD_BASE_DIR, paste0("ldblk_ukbb_", tolower(anc)))
    snpinfo_file <- file.path(h5_dir, "snpinfo_ukbb_hm3")
    
    tryCatch({
      .save_png_atomic(specific_path, width = 1200, height = 1200, res = 150, draw = function() {
        generate_locuszoom_by_rsid(
          rsid = lead,
          win_bp = 200e3,
          h5_dir = h5_dir,
          snpinfo_file = snpinfo_file,
          gwas_dt = cen$dt,
          gwas_snp_col = "SNP",
          gwas_p_col = "P",
          ens_db = EnsDb.Hsapiens.v75
        )
      })
      
      list(src = specific_path, contentType = "image/png", width = "100%")
    }, error = function(e) {
      stop(e)
    })
    
  }, deleteFile = FALSE)
  outputOptions(output, "locuszoom_plot", suspendWhenHidden = FALSE)
  observeEvent(list(gwas_raw(), loci_rv()), {
    dt_std <- gwas_raw()
    loci_dt <- loci_rv()
    req(!is.null(dt_std), nrow(dt_std) > 0, !is.null(loci_dt), nrow(loci_dt) > 0)
    
    lead0 <- .pick_default_lead_from_loci(loci_dt)
    
    loci_dt2 <- data.table::as.data.table(loci_dt)
    row0 <- loci_dt2[`lead snp` == lead0][1]
    req(nrow(row0) > 0)
    
    chr0 <- suppressWarnings(as.integer(row0$chr))
    pos0 <- row0$pos
    if (is.character(pos0)) pos0 <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pos0)) * 1e6)
    pos0 <- as.integer(pos0)
    
    .center(list(chr = chr0, pos = pos0, title = "LocusZoom", dt = dt_std, `lead snp` = lead0))
  }, ignoreInit = FALSE)
  
  # 7) blockzoom
  .H5_RSID2BLK_CACHE <- new.env(parent = emptyenv())
  .get_rsid2blk_map <- function(h5file){
    key <- normalizePath(h5file, winslash = "/", mustWork = FALSE)
    if (exists(key, envir = .H5_RSID2BLK_CACHE, inherits = FALSE)) {
      return(get(key, envir = .H5_RSID2BLK_CACHE, inherits = FALSE))
    }
    
    grps <- .get_h5_grps(h5file)  # must exist in your server
    if (!length(grps)) stop("HDF5 内未发现可用 block groups: ", basename(h5file))
    
    mp <- new.env(parent = emptyenv(), hash = TRUE)
    
    for (g in grps) {
      snp <- .h5_read_safe(h5file, paste0(g, "/snplist"))  # must exist in your server
      if (is.null(snp)) next
      
      ids <- NULL
      if (is.data.frame(snp)) {
        nms <- tolower(names(snp))
        idcol <- intersect(c("snp","rsid","id"), nms)[1]
        if (!is.na(idcol) && nzchar(idcol)) {
          ids <- as.character(snp[[which(nms == idcol)[1]]])
        }
      } else if (is.matrix(snp) && ncol(snp) >= 1) {
        ids <- as.character(snp[, 1])
      } else if (is.character(snp)) {
        ids <- as.character(snp)
      }
      if (is.null(ids) || !length(ids)) next
      
      # Parse block_id from group name (handles "/blk_123" or ".../123")
      gid <- suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", gsub(".*/", "", g))))
      if (!is.finite(gid)) gid <- suppressWarnings(as.integer(gsub(".*/", "", g)))
      if (!is.finite(gid)) next
      
      for (rs in ids) {
        if (!is.na(rs) && nzchar(rs)) assign(toupper(rs), gid, envir = mp)
      }
    }
    
    assign(key, mp, envir = .H5_RSID2BLK_CACHE)
    mp
  }
  .find_block_id_from_h5_cached <- function(chr_i, rsid, h5map){
    chr_i <- suppressWarnings(as.integer(chr_i))
    if (!is.finite(chr_i)) stop("chr 无效，无法在 HDF5 中定位 block")
    
    rsid0 <- toupper(trimws(as.character(rsid)[1]))
    if (!nzchar(rsid0)) stop("rsid 无效，无法在 HDF5 中定位 block")
    
    # h5map is expected to be a named vector like chr1..chr22
    h5 <- unname(h5map[paste0("chr", chr_i)])
    if (is.na(h5) || !file.exists(h5)) stop("找不到 chr", chr_i, " 对应的 HDF5 文件")
    
    mp <- .get_rsid2blk_map(h5)
    if (exists(rsid0, envir = mp, inherits = FALSE)) {
      return(get(rsid0, envir = mp, inherits = FALSE))
    }
    
    stop("在 HDF5 的 snplist 中找不到 rsid：", rsid0)
  }
  .draw_blockzoom <- function(row_dt, gwas_dt, build, ancestry) {
    
    chr0 <- suppressWarnings(as.integer(row_dt$chr[1]))
    lead <- as.character(row_dt[[ intersect(names(row_dt), c("lead_snp","lead snp","leadsnp"))[1] ]][1])
    
    if (!is.finite(chr0) || is.na(lead) || !nzchar(lead)) {
      stop("block 表缺少有效的 chr 或 lead snp（列名可能是 lead snp / lead_snp）")
    }

    # build/ancestry -> HDF5 map (server logic)
    build_hdl <- .norm_build_hdl(.norm_build_disp(build))  # must exist
    anc       <- .norm_pop(ancestry)                       # must exist
    
    h5map    <- resolve_ld_h5_files(LD_BASE_DIR, build_hdl, anc)  # must exist
    valid_h5 <- h5map[!is.na(h5map) & file.exists(h5map)]
    if (length(valid_h5) < 1) stop("未找到对应 LD HDF5")
    
    # blockzoomf dependencies (global assignment as you requested)
    ld_dir <<- dirname(valid_h5[1])
    bn       <- basename(valid_h5[1])
    hapref <<- sub("^ldblk_([^_]+)_chr.*$", "\\1", bn)
    
    # locate center block by lead_snp (fast cached)
    block_id <- .find_block_id_from_h5_cached(chr0, lead, h5map)
    
    # triplet logic is inside get_LD / block_plot (blockzoomf logic)
    LD <- get_LD(block_id = block_id, gwas_dt = gwas_dt)
    if (!is.list(LD) || is.null(LD$gwas) || !is.data.frame(LD$gwas) || nrow(LD$gwas) < 1) {
      stop("get_LD 无可画数据")
    }
    
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(cex = 0.7)
    block_plot(gwas = LD$gwas, block_id = block_id, LD = LD, gene_df = NULL)
  }
  output$blockzoom_plot <- renderImage({
    # --- 1. 显式依赖：点击表格行 ---
    dat <- blocks_rv()
    sel <- input$signal_block_table_rows_selected
    req(dat)
    
    # 确定当前选中的行索引
    idx <- if (is.null(sel) || length(sel) < 1) 1L else as.integer(sel[1])
    if (idx > nrow(dat)) idx <- 1L
    row_dt <- dat[idx, , drop = FALSE]
    
    # 获取 lead snp 用于命名缓存
    lead_col <- intersect(names(row_dt), c("lead_snp","lead snp","leadsnp"))[1]
    lead_name <- as.character(row_dt[[lead_col]][1])
    
    gw <- gw_current()
    cache_dir <- file.path("cache", gw$name)
    if(!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    
    # 动态路径：例如 cache/my_gwas/blockzoom_rs999.png
    specific_path <- file.path(cache_dir, paste0("blockzoom_", lead_name, ".png"))
    
    # --- 2. 检查缓存 ---
    if (file.exists(specific_path)) {
      return(list(src = specific_path, contentType = "image/png", width = "100%"))
    }
    
    # --- 3. 绘图 ---
    gwas_dt <- gwas_raw()
    req(gwas_dt)
    
    tryCatch({
      .save_png_atomic(specific_path, width = 1200, height = 1200, res = 150, draw = function() {
        .draw_blockzoom(
          row_dt   = row_dt,
          gwas_dt  = gwas_dt,
          build    = input$build,
          ancestry = input$ancestry
        )
      })
      list(src = specific_path, contentType = "image/png", width = "100%")
    }, error = function(e) {
      stop(e)
    })
    
  }, deleteFile = FALSE)
  outputOptions(output, "blockzoom_plot", suspendWhenHidden = FALSE)
  output$blockzoom_title <- renderText({
    
    dat <- blocks_rv()
    sel <- input$signal_block_table_rows_selected
    req(!is.null(dat), nrow(dat) > 0)
    
    if (is.null(sel) || length(sel) < 1) sel <- 1L
    sel <- suppressWarnings(as.integer(sel[1]))
    if (!is.finite(sel) || sel < 1L || sel > nrow(dat)) sel <- 1L
    
    row_dt <- dat[sel, , drop = FALSE]
    
    # ---- chr ----
    chr_now <- if ("chr" %in% names(row_dt))
      suppressWarnings(as.integer(row_dt$chr[1]))
    else NA_integer_
    
    if (!is.finite(chr_now)) return("BlockZoom")
    
    # ---- lead SNP（兼容列名）----
    lead_col <- intersect(names(row_dt), c("lead_snp","lead snp","leadsnp"))[1]
    if (is.na(lead_col)) {
      return(sprintf("BlockZoom · CHR%d", chr_now))
    }
    
    lead <- as.character(row_dt[[lead_col]][1])
    if (is.na(lead) || !nzchar(lead)) {
      return(sprintf("BlockZoom·CHR%d", chr_now))
    }
    
    # ---- 用与绘图完全一致的逻辑找 block_id ----
    build_hdl <- .norm_build_hdl(.norm_build_disp(input$build))
    anc       <- .norm_pop(input$ancestry)
    h5map     <- resolve_ld_h5_files(LD_BASE_DIR, build_hdl, anc)
    
    block_id <- tryCatch(
      .find_block_id_from_h5_cached(chr_now, lead, h5map),
      error = function(e) NA_integer_
    )
    
    if (!is.finite(block_id)) {
      return(sprintf("BlockZoom·CHR%d", chr_now))
    }
    
    # ---- neighbors（triplet 逻辑）----
    neighbors <- c(block_id - 1L, block_id, block_id + 1L)
    
    sprintf(
      "BlockZoom·CHR%d:blk%d",
      chr_now,
      block_id
    )

  })

}