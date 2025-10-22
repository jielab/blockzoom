options(shiny.maxRequestSize = 2 * 1024^3)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(shiny, data.table, readr, DT, ggplot2, CMplot, rhdf5,
               ggrepel, scales, ggrastr, rtracklayer, GenomicRanges, GenomeInfoDb)
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
try(rhdf5::h5closeAll(), silent = TRUE)
# 限制底层并行线程，避免 native 崩溃
Sys.setenv(OMP_NUM_THREADS = "1",MKL_NUM_THREADS = "1",OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",RCPP_PARALLEL_NUM_THREADS = "1",NUMEXPR_NUM_THREADS = "1")

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
rank_blocks_topN_HDL <- function(
    h5file, sumstats_path, N_eff = 1e5, top_n = Inf, save_csv = NULL, map_file = NULL,
    harmonize = TRUE, max_eig = 12000L, ridge = 1e-4, min_overlap = 50L,
    prefilter = FALSE  # ← 新增：默认不过滤任何块
){
  on.exit(rhdf5::h5closeAll(), add = TRUE)
  
  ss <- data.table::as.data.table(read_gwas(sumstats_path))
  ss[, Z := .ensure_Z(ss)]
  snp_col <- .pick(names(ss), c("SNP","ID","RSID","MARKER","MARKERNAME","VARIANT","VARIANT_ID"))
  if (is.na(snp_col)) stop("sumstats 缺少 SNP/VARIANT 列")
  ss[, SNP_CLEAN := .clean_snp(get(snp_col))]
  data.table::setkey(ss, SNP_CLEAN)
  
  a1_col <- .pick(names(ss), c("A1","EA","ALLELE1"))
  a2_col <- .pick(names(ss), c("A2","NEA","ALLELE2"))
  
  have_chrpos <- all(c("CHR","POS") %in% names(ss))
  if (!have_chrpos && !is.null(map_file)){
    map <- data.table::fread(map_file)
    if (!"SNP" %in% names(map)) stop("map 需包含列 SNP")
    chr <- .pick(names(map), c("CHR","Chr","chr","CHROM","Chromosome","chrom"))
    pos <- .pick(names(map), c("BP","bp","POS","Pos","Bp"))
    if (is.na(chr) || is.na(pos)) stop("map 需提供 CHR 与 BP/POS")
    data.table::setnames(map, c("SNP", chr, pos), c("SNP","CHR","POS"))
    map[, SNP_CLEAN := .clean_snp(SNP)]
    ss <- merge(ss, data.table::unique(map[, .(SNP_CLEAN, CHR, POS)]),
                by = "SNP_CLEAN", all.x = TRUE, sort = FALSE)
    ss <- .norm_chr_pos(ss)
    have_chrpos <- all(c("CHR","POS") %in% names(ss))
  }
  
  grps <- .get_h5_grps(h5file)
  if (!length(grps)) stop("HDF5 未找到包含 snplist 与 ldblk 的分组")
  
  .quick_check_block <- function(h5file, g, ss, snp_col, min_overlap = 50L, pthr_pre = 1e-4){
    sp  <- .h5_read_safe(h5file, file.path(g,"snplist"))
    ref <- if (is.character(sp)) data.table::data.table(rsid=sp) else data.table::as.data.table(sp)
    rn  <- .pick(names(ref), c("SNP","RSID","ID","rsid","id"))
    if (is.na(rn)) return(FALSE)
    data.table::setnames(ref, rn, "RSID")
    ref[, SNP_CLEAN := .clean_snp(RSID)]
    
    idx <- match(ref$SNP_CLEAN, ss$SNP_CLEAN)
    keep <- which(!is.na(idx))
    if (length(keep) < min_overlap) return(FALSE)
    
    pvec <- ss$P[idx[keep]]
    if (!any(is.finite(pvec))) return(FALSE)
    min(pvec, na.rm=TRUE) <= pthr_pre
  }
  
  mle_block <- function(g){
    sp  <- .h5_read_safe(h5file, file.path(g,"snplist"))
    ref <- if (is.character(sp)) data.table::data.table(rsid=sp) else data.table::as.data.table(sp)
    
    nm0 <- names(ref)
    if (!is.null(nm0)) {
      nmU <- toupper(trimws(sub("^\ufeff","", nm0)))
      data.table::setnames(ref, nm0, nmU, skip_absent = TRUE)
    }
    rn  <- .pick(names(ref), c("SNP","RSID","ID","rsid","id"))
    if (is.na(rn)) return(NULL)
    data.table::setnames(ref, rn, "RSID")
    ref[, SNP_CLEAN := .clean_snp(RSID)]
    
    have_ref_chrpos <- all(c("CHR","POS") %in% names(ref))
    if (!have_ref_chrpos) {
      cp <- .parse_chrpos_from_str(ref$RSID)
      if (any(is.finite(cp$CHR) & is.finite(cp$POS))) {
        ref[, `:=`(CHR = cp$CHR, POS = cp$POS)]
        have_ref_chrpos <- TRUE
      }
    }
    
    a1r <- .pick(names(ref), c("A1","EA","ALLELE1","a1","ea"))
    a2r <- .pick(names(ref), c("A2","NEA","ALLELE2","a2","nea"))
    if (!is.na(a1r) && !is.na(a2r)) {
      data.table::setnames(ref, c(a1r,a2r), c("A1","A2"), skip_absent = TRUE)
    } else {
      ref[,`:=`(A1=NA_character_,A2=NA_character_)]
    }
    
    chr_r <- .pick(names(ref), c("CHR","#CHROM","CHROM","CHROMOSOME"))
    pos_r <- .pick(names(ref), c("POS","BP","POSITION","BASE_PAIR_LOCATION"))
    if (!is.na(chr_r) && !is.na(pos_r)) {
      data.table::setnames(ref, c(chr_r,pos_r), c("CHR","POS"), skip_absent = TRUE)
      ref <- .norm_chr_pos(ref)
    }
    
    ss <- data.table::as.data.table(ss)
    idx  <- match(ref$SNP_CLEAN, ss$SNP_CLEAN)
    
    keep <- which(!is.na(idx))
    if (length(keep) < min_overlap) {
      have_ref_chrpos <- all(c("CHR","POS") %in% names(ref))
      have_ss_chrpos  <- all(c("CHR","POS") %in% names(ss))
      if (have_ref_chrpos && have_ss_chrpos) {
        ref[, `:=`(CHR = suppressWarnings(as.integer(gsub("^chr","", as.character(CHR), ignore.case=TRUE))),
                   POS = suppressWarnings(as.integer(POS)))]
        ss[,  `:=`(CHR = suppressWarnings(as.integer(gsub("^chr","", as.character(CHR),  ignore.case=TRUE))),
                   POS = suppressWarnings(as.integer(POS)))]
        ref[, KEY_CP := paste0(CHR, ":", POS)]
        ss[,  KEY_CP := paste0(CHR, ":", POS)]
        idx  <- match(ref$KEY_CP, ss$KEY_CP)
        keep <- which(!is.na(idx))
      }
    }
    if (length(keep) < min_overlap) return(NULL)
    
    R <- as.matrix(.h5_read_safe(h5file, file.path(g,"ldblk"), index = list(keep, keep)))
    if (nrow(R) != ncol(R)) return(NULL)
    idx_ss <- idx[keep]
    z <- ss$Z[idx_ss]
    
    ii <- which(is.finite(z))
    if (length(ii) < min_overlap) return(NULL)
    z <- z[ii]; R <- R[ii, ii, drop=FALSE]; idx_ss <- idx_ss[ii]
    R[!is.finite(R)] <- 0
    
    if (isTRUE(harmonize) && !is.na(a1_col) && !is.na(a2_col) && all(!is.na(ref$A1[keep[ii]]))) {
      a1 <- toupper(as.character(ss[[a1_col]][idx_ss])); a2 <- toupper(as.character(ss[[a2_col]][idx_ss]))
      r1 <- toupper(as.character(ref$A1[keep[ii]]));     r2 <- toupper(as.character(ref$A2[keep[ii]]))
      flip <- (a1==r2 & a2==r1)
      ambi <- ((a1 %in% c("A","T") & a2 %in% c("A","T")) | (a1 %in% c("C","G") & a2 %in% c("C","G")))
      jj <- which(flip & !ambi); if (length(jj)) z[jj] <- -z[jj]
    }
    
    m <- length(z)
    if (m > max_eig) {
      ord <- order(abs(z), decreasing = TRUE)
      sel <- sort(ord[1:max_eig])
      R <- R[sel, sel, drop = FALSE]
      z <- z[sel]
      idx_ss <- idx_ss[sel]
      m <- length(z)
    }
    
    diag(R) <- pmin(1, pmax(0, diag(R)))
    R <- R + diag(ridge, m)
    z <- pmax(pmin(z, 12), -12)
    
    eg  <- try(eigen(R, symmetric=TRUE), silent=TRUE)
    if (inherits(eg, "try-error")) return(NULL)
    lam <- pmax(0, as.numeric(eg$values))
    keep_eig <- which(lam > 1e-8)
    if (length(keep_eig) < 10L) return(NULL)
    U   <- eg$vectors[, keep_eig, drop=FALSE]
    lam <- lam[keep_eig]
    y   <- as.numeric(crossprod(U, z))
    
    nll <- function(s){
      d <- 1 + s * lam
      if (any(d <= 0) || any(!is.finite(d))) return(Inf)
      0.5 * (sum(log(d)) + sum((y*y)/d))
    }
    
    num <- sum(z * as.numeric(R %*% z)) - sum(diag(R))
    den <- sum(R * R); if (!is.finite(den) || den <= 0) den <- 1e-12
    s_mom   <- max(0, num/den)
    lam_max <- max(lam)
    s_upper <- min(1e4, if (lam_max > 0) 1e3/lam_max else 1e4)
    opt <- suppressWarnings(try(optimize(nll, interval = c(0, s_upper)), silent = TRUE))
    if (inherits(opt, "try-error")) return(NULL)
    s_hat <- max(0, min(opt$minimum, s_upper))
    s_hat <- max(s_hat, s_mom, 0)
    
    h  <- max(1e-6, 1e-3 * (1 + s_hat))
    n2 <- (nll(s_hat + h) - 2*nll(s_hat) + nll(s_hat - h)) / (h*h)
    se <- if (is.finite(n2) && n2 > 0 && is.finite(N_eff) && N_eff > 0) sqrt(1/n2)/N_eff else NA_real_
    h2 <- s_hat / N_eff
    h2 <- pmin(pmax(h2, 0), 1)
    
    chr <- start <- end <- NA_real_
    if (have_chrpos){
      cv <- suppressWarnings(as.numeric(ss$CHR[idx_ss])); pv <- suppressWarnings(as.numeric(ss$POS[idx_ss]))
      if (any(is.finite(cv))) chr <- as.integer(names(sort(table(cv[is.finite(cv)]), decr=TRUE))[1])
      if (any(is.finite(pv))) { r <- range(pv[is.finite(pv)], na.rm=TRUE); start <- r[1]; end <- r[2] }
    }
    pvec <- 2*pnorm(-abs(z)); imin <- which.min(pvec)
    lead_snp <- ss[[snp_col]][ idx_ss[imin] ]
    p_min <- suppressWarnings(min(pvec, na.rm = TRUE))
    
    data.table::data.table(
      blk=sub("^/","", g), chr=chr, start=start, end=end,
      lead_snp=lead_snp, p_min=p_min, h2=h2, se=se
    )
  }
  
  res <- data.table::rbindlist(lapply(grps, function(g){
    if (isTRUE(prefilter)) {
      ok <- .quick_check_block(h5file, g, ss, snp_col,
                               min_overlap = min_overlap, pthr_pre = 1e-4)
      if (!ok) return(NULL)  # 只有 prefilter=TRUE 才筛
    }
    tryCatch(mle_block(g), error=function(e) NULL)
  }), fill=TRUE)
  
  if (!nrow(res)) {
    out <- data.table::data.table(
      blk=character(), chr=integer(), start=integer(), end=integer(),
      lead_snp=character(), p_min=numeric(), h2=numeric(), se=numeric()
    )
    if (!is.null(save_csv)) data.table::fwrite(out, save_csv)
    return(out[])
  }
  
  data.table::setorder(res, -h2)
  if (is.infinite(top_n) || is.na(top_n)) top_n <- .Machine$integer.max
  out <- res[1:min(top_n, .N), .(blk, chr, start, end, lead_snp, p_min, h2, se)]
  if (!is.null(save_csv)) data.table::fwrite(out, save_csv)
  out[]
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
read_gwas <- function(path){
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
  pos_col <- pick1(c("POS","BP","POSITION","BASE_PAIR_LOCATION","BP_HG19","POS37"))
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
  P <- pmax(P, 1e-50)   # ← 关键
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
.plot_locus <- function(dt_std, chr0, pos0, win_bp=200e3, title_prefix="LocusZoom"){
  win <- dt_std[CHR==chr0 & POS>=pos0-win_bp & POS<=pos0+win_bp]
  validate(need(nrow(win)>0, "该窗口内无数据"))
  
  # 用原始 P 计算可视化用 p_plot（经过 scale_p 压缩）
  p_raw <- suppressWarnings(as.numeric(win$P))
  p_raw[!is.finite(p_raw) | p_raw <= 0] <- NA
  p_plot <- scale_p(p_raw, min_p = 1e-324, max_log10 = 50)
  
  win[, mlog10P_plot := -log10(p_plot)]
  win[, POS_Mb := POS/1e6]
  
  sig_thr <- 1e-5
  if (nrow(win[P > sig_thr]) > 2000) {
    win <- rbind(win[P <= sig_thr], win[P > sig_thr][sample(.N, 2000)])
  }
  lead_i <- which.min(abs(win$POS - pos0))
  
  use_rast <- requireNamespace("ggrastr", quietly=TRUE)
  gp <- ggplot2::ggplot(win, ggplot2::aes(POS_Mb, mlog10P_plot)) +
    { if (use_rast) ggrastr::geom_point_rast(ggplot2::aes(color = mlog10P_plot), alpha=.6, size=1)
      else ggplot2::geom_point(ggplot2::aes(color = mlog10P_plot), alpha=.6, size=1) } +
    ggplot2::scale_color_gradient(low="blue", high="red") +
    ggplot2::labs(
      title=sprintf("%s · CHR%s:%.2f Mb ± %skb", title_prefix, chr0, pos0/1e6, win_bp/1000),
      x=sprintf("Genomic position on chr%s (Mb)", chr0), y=expression(-log[10](P))
    ) +
    ggplot2::geom_hline(yintercept=-log10(5e-8), linetype="dashed", linewidth=.3) +
    ggplot2::geom_hline(yintercept=-log10(1e-5), linetype="dotted", linewidth=.3) +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), legend.position="none")
  
  if (length(lead_i) == 1) {
    lead <- win[lead_i]
    gp <- gp +
      ggplot2::geom_point(data=lead, ggplot2::aes(POS_Mb, mlog10P_plot),
                          shape=21, stroke=.7, size=2.8, fill="white", color="#D55E00") +
      ggrepel::geom_label_repel(
        data=lead, ggplot2::aes(POS_Mb, mlog10P_plot, label=SNP), size=2.2,
        label.size=.08, label.r=grid::unit(0.05,"lines"), box.padding=.25, point.padding=.15,
        max.overlaps = 30, seed=1
      )
  }
  gp
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
  dt <- tryCatch(
    data.table::fread(p, sep="\t", header=TRUE, nThread=1),
    error=function(e) NULL
  )
  is.data.frame(dt) && nrow(dt) > 0
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

.parse_blk_pos_bp <- function(pos_col){
  s <- gsub("[[:space:]]", "", as.character(pos_col))
  s <- gsub("Mb", "", s, ignore.case = TRUE)
  sp <- data.table::tstrsplit(s, "[-–]", fixed = FALSE)
  start_mb <- suppressWarnings(as.numeric(sp[[1]]))
  end_mb   <- suppressWarnings(as.numeric(sp[[2]]))
  start_bp <- as.numeric(start_mb) * 1e6
  end_bp   <- as.numeric(end_mb)   * 1e6
  data.table::data.table(start_bp = start_bp, end_bp = end_bp,
                         center_bp = (start_bp + end_bp)/2)
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
  
  keep <- logical(nrow(x))
  for (i in seq_len(nrow(x))) {
    if (!keep[i] && !is.na(keep[i])) {
      keep[i] <- TRUE
      idx_same <- which(x$chr == x$chr[i] & !keep & !is.na(keep))
      if (length(idx_same)) {
        dist <- abs(x$center_bp[idx_same] - x$center_bp[i])
        keep[idx_same[dist < window_kb*1000]] <- NA
      }
    }
  }
  # —— 需要的列：blk、SNP、chr、pos、p、h2、blk_h5（可选）——
  base_cols <- c("blk", snp_col, "chr", "pos", "p", "h2", "blk_h5")
  base_cols <- intersect(base_cols, names(x))
  out <- x[keep %in% TRUE, ..base_cols]
  data.table::setnames(out, snp_col, "lead snp", skip_absent = TRUE)
  
  # ① h2 阈
  if (is.finite(min_h2)) out <- out[is.finite(h2) & h2 >= min_h2]
  # ② 每条染色体上限
  if (is.finite(top_per_chr)) { out[, ord := seq_len(.N), by = chr]; out <- out[ord <= top_per_chr][, ord := NULL] }
  # ③ 全局上限
  if (is.finite(top_total) && nrow(out) > top_total) out <- out[order(p, -h2)][seq_len(top_total)]
  
  # —— 最终固定按 blk（数字）升序（若有 chr 就先 chr 再 blk）——
  if ("blk" %in% names(out)) suppressWarnings(out[, blk := as.integer(blk)])
  if (all(c("chr","blk") %in% names(out))) {
    data.table::setorder(out, chr, blk)
  } else if ("blk" %in% names(out)) {
    data.table::setorder(out, blk)
  }
  
  out[]
}

.chr_cache_paths_abs <- function(gwas_name, build_disp, anc, chr_i){
  cache_dir <- file.path(CACHE_ROOT, gwas_name)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- sprintf("%s_%s_chr%02d", build_disp, anc, as.integer(chr_i))
  list(tsv = file.path(cache_dir, paste0(stem, ".tsv")))
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
  #    注意：CHR/POS 本来上游已处理，但这里再保险一次
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
  
  data.table::setorder(DT, CHR, POS, P)
  out <- data.table::data.table(SNP=character(), CHR=integer(),
                                POS=integer(), P=double(), H2=double())
  
  # ⑤ 按染色体逐段聚类
  for (c in unique(DT$CHR)) {
    x <- DT[CHR == c]
    if (!nrow(x)) next
    
    # 选择当前簇的领头（最左且最显著的那行）
    cur_lead <- x[1]
    # 注意：这里若 POS 是 NA 会出错，所以上面已确保 finite(POS)
    cur_left  <- cur_lead$POS - kbp
    cur_right <- cur_lead$POS + kbp
    
    for (i in 2:nrow(x)) {
      # 用 isTRUE 防止 NA 进入 if 条件
      if (isTRUE(x$POS[i] > cur_right)) {
        out <- rbind(out, cur_lead)
        cur_lead <- x[i]
        cur_left  <- cur_lead$POS - kbp
        cur_right <- cur_lead$POS + kbp
      } else {
        # 仍在窗口内：更新更显著的 lead，拓展右边界
        if (isTRUE(x$P[i] < cur_lead$P)) cur_lead <- x[i]
        # max() 对 NA 不友好，这里保证入参有限
        r2 <- x$POS[i] + kbp
        if (!is.finite(r2)) r2 <- cur_right
        cur_right <- max(cur_right, r2, na.rm = TRUE)
      }
    }
    out <- rbind(out, cur_lead)
  }
  
  unique(out, by = c("CHR","POS","SNP"))[]
}
collapse_to_independent_loci <- function(dt, ldb = NULL, pthr = 5e-8, window_kb = 500) {
  collapse_to_loci(dt, p_thr = pthr, window_kb = window_kb)
}

# ----------------------------【Server 主体】------------------------------------
server <- function(input, output, session){
  
  GTF_PATH <- file.path(FILES_DIR, "gencode.v19.annotation.gtf.gz")  # 按需改路径
  genes_gr_rv <- reactiveVal(NULL)
  .load_gene_gr <- function(gtf_path){
    stopifnot(file.exists(gtf_path))
    gr <- rtracklayer::import(gtf_path)
    gr <- gr[gr$type == "gene"]
    
    # ★ 关键：改命名风格而不是手动改 seqnames 向量
    GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"  # chr1→1, chrM→MT
    
    # 仅保留标准染色体
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
    gr <- GenomicRanges::trim(gr)
    
    # gene 名
    gname <- gr$gene_name
    gname[is.na(gname) | !nzchar(gname)] <- NA
    if (all(is.na(gname)) && !is.null(gr$Name)) gname <- gr$Name
    if (all(is.na(gname))) gname <- gr$gene_id
    mcols(gr)$gene <- as.character(gname)
    gr
  }
  observeEvent(TRUE, {
    if (file.exists(GTF_PATH)) {
      genes_gr_rv(.load_gene_gr(GTF_PATH))
    } else {
      showNotification(sprintf("GTF not found: %s", GTF_PATH), type="warning", duration=8)
    }
  }, once = TRUE)
  
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
    applied_once(FALSE)
  }, ignoreInit = TRUE)
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
      qq  = file.path(dir, "qq.png")
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
      blocks = file.path(cache_dir, sprintf("significant blocks_%s_%s.tsv", build_disp, pop))
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
  .resolve_gwas <- function(input_file, builtin){
    .safe_name <- function(x){
      x <- gsub("[^A-Za-z0-9._-]+", "_", x)
      x <- gsub("_+", "_", x)
      x <- sub("^_+", "", x)
      x <- sub("_+$", "", x)
      if (!nzchar(x)) "unknown" else x
    }
    
    # 1) 用户上传优先
    is_up <- !is.null(input_file) && is.list(input_file) &&
      "datapath" %in% names(input_file) &&
      length(input_file$datapath) >= 1 &&
      nzchar(as.character(input_file$datapath[[1]]))
    if (is_up) {
      p  <- normalizePath(as.character(input_file$datapath[[1]]), winslash="/", mustWork=TRUE)
      nm <- if (!is.null(input_file$name) && length(input_file$name)>=1)
        as.character(input_file$name[[1]]) else basename(p)
      return(list(path=p, name=.safe_name(tools::file_path_sans_ext(basename(nm)))))
    }
    
    # 2) 内置（来自 catalog）
    req(builtin)
    f <- normalizePath(as.character(builtin), winslash="/", mustWork=FALSE)
    validate(need(file.exists(f), paste0("GWAS 文件不存在：", f)))
    list(path=f, name=.safe_name(tools::file_path_sans_ext(basename(f))))
  }
  
  loci_rv   <- reactiveVal(NULL)
  blocks_rv <- reactiveVal(NULL)
  .load_loci_from_cache <- function(gwas_name, build_disp, pop){
    paths <- .cache_paths_abs(gwas_name, build_disp, pop)
    if (!file.exists(paths$loci)) return(NULL)
    dt <- tryCatch(data.table::fread(paths$loci, sep="\t", header=TRUE, data.table=TRUE),
                   error=function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) dt else NULL
  }
  .load_blocks_from_cache <- function(gwas_name, build_disp, anc){
    paths <- .cache_paths_abs(gwas_name, build_disp, anc)
    if (!file.exists(paths$blocks)) return(NULL)
    dt <- tryCatch(data.table::fread(paths$blocks, sep="\t", header=TRUE, data.table=TRUE),
                   error=function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) {
      # 统一成 "lead snp"
      if (!"lead snp" %in% names(dt)) {
        for (cand in c("lead_snp", "snp", "SNP")) {
          if (cand %in% names(dt)) {
            data.table::setnames(dt, cand, "lead snp", skip_absent = TRUE)
            break
          }
        }
      }
      # ← 就加这句：彻底移除 se
      if ("se" %in% names(dt)) dt[, se := NULL]
      dt[]
    } else NULL
  }
  observeEvent(list(input$gwas_file, input$gwas_builtin, input$build, input$ancestry, input$p_threshold), {
    gw   <- .resolve_gwas(input$gwas_file, input$gwas_builtin)
    bld  <- .norm_build_disp(input$build)
    pop  <- .norm_pop(input$ancestry)
    
    loci <- .load_loci_from_cache(gw$name, bld, pop)
    if (!is.null(loci)) loci_rv(loci)
    
    blk  <- .load_blocks_from_cache(gw$name, bld, pop)
    if (!is.null(blk))  blocks_rv(blk)
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
  .top_peak <- function(dt_std){
    p <- suppressWarnings(as.numeric(dt_std$P))
    p[!is.finite(p)] <- Inf
    i <- which.min(p)
    list(chr = as.integer(dt_std$CHR[i]), pos = as.numeric(dt_std$POS[i]))
  }
  .center <- reactiveVal(NULL)
  
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
    # 回退到你原有的 read_gwas（自带一套别名自动识别）
    read_gwas(gw$path)
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
    
    # 2) 无缓存就生成一次并写入缓存
    p_all <- suppressWarnings(as.numeric(gwas_raw()[["P"]]))
    p_all <- suppressWarnings(as.numeric(gwas_raw()[["P"]]))
    p_raw <- p_all[is.finite(p_all)]
    p_raw[p_raw <= 0] <- 1e-300            # 数值稳定性，供lambda/CI
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
    
    # —— 从 blocks 汇总总体 h2（优先内存，其次 JSON，再 ALL 文件，最后代表集）——
    build_disp  <- .norm_build_disp(input$build)
    anc         <- .norm_pop(input$ancestry)
    cache_paths <- .cache_paths_abs(gw$name, build_disp, anc)
    
    h2_total <- NA_real_
    se_total <- NA_real_
    
    rv <- h2_total_rv()
    if (is.list(rv) && is.finite(rv$h2)) {
      h2_total <- rv$h2
      se_total <- rv$se
    } else {
      # 先尝试读 JSON
      json_path <- file.path(dirname(cache_paths$blocks), "H2_TOTAL.json")
      if (file.exists(json_path)) {
        jj <- tryCatch(jsonlite::read_json(json_path, simplifyVector = TRUE),
                       error = function(e) NULL)
        if (is.list(jj) && is.finite(jj$h2)) {
          h2_total <- jj$h2
          se_total <- if (is.null(jj$se)) NA_real_ else jj$se
        }
      }
      # 还不行再退到 ALL 或代表集
      if (!is.finite(h2_total)) {
        all_path <- file.path(dirname(cache_paths$blocks),
                              sprintf("ALL_blocks_%s_%s.tsv", build_disp, anc))
        src <- if (file.exists(all_path)) all_path else cache_paths$blocks
        suppressWarnings({
          bt <- tryCatch(data.table::fread(src, sep = "\t", header = TRUE, data.table = FALSE),
                         error = function(e) NULL)
          if (is.data.frame(bt) && "h2" %in% names(bt)) {
            h2_total <- suppressWarnings(sum(as.numeric(bt$h2), na.rm = TRUE))
            if ("se" %in% names(bt)) {
              se_total <- suppressWarnings(sqrt(sum((as.numeric(bt$se))^2, na.rm = TRUE)))
            }
          }
        })
      }
    }
    
    # 计算总体现象学指标
    z_total <- if (is.finite(h2_total) && is.finite(se_total) && se_total > 0) h2_total / se_total else NA_real_
    p_total <- if (is.finite(z_total)) 2 * pnorm(-abs(z_total)) else NA_real_
    
    p_show <- if (is.finite(p_total)) {
      if (p_total < 1e-6) "p<1e-6" else sprintf("p=%.2f", p_total)
    } else "p=NA"
    
    h2_label <- if (is.finite(h2_total)) sprintf("h\u00B2=%.3f\n%s", h2_total, p_show) else NA_character_
    
    .save_png_atomic(paths$qq, width = 1200, height = 1200, res = 144, draw = function(){
      use_rast <- requireNamespace("ggrastr", quietly = TRUE)
      pt_layer <- if (use_rast) {
        ggrastr::geom_point_rast(data = df, ggplot2::aes(exp, obs), size = 0.25, alpha = 0.7)
      } else {
        ggplot2::geom_point(data = df, ggplot2::aes(exp, obs), size = 0.25, alpha = 0.7)
      }
      p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = exp, ymin = lo, ymax = hi), fill = "grey92") +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey50") +
        pt_layer +
        ggplot2::coord_fixed(ratio = 1, clip = "off") +
        ggplot2::scale_x_continuous(limits = c(0, max_lim),
                                    breaks = scales::pretty_breaks(7),
                                    expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
        ggplot2::scale_y_continuous(limits = c(0, max_lim),
                                    breaks = scales::pretty_breaks(7),
                                    expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
        ggplot2::labs(
          x = "Expected -log10(P)", y = "Observed -log10(P)", title = "QQ plot",
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
          x = max_lim * 0.98, y = max_lim * 0.05,
          label = h2_label,
          hjust = 1, vjust = 0,
          size = 8, 
          label.r = grid::unit(0.15, "lines"),
          fill = "white", color = "black"
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
    lead <- collapse_to_independent_loci(dt, ldb = NULL, pthr = thr, window_kb = 500)
    
    res <- if (nrow(lead)) lead[, .(
      loc = .I,
      `lead snp` = SNP,
      chr = as.integer(CHR),
      pos = sprintf("%.2f Mb", as.numeric(POS)/1e6),
      p   = as.numeric(P),
      h2  = as.numeric(H2)   # ← 用 lead SNP 的 H2
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
      DT::formatSignif(columns = "p", digits = 3) |>      # p 依旧用科学计数
      DT::formatRound(columns = "h2", digits = 4)         # h2 固定小数点 4 位
  })
  
  # 5) significant blocks
  observeEvent(input$apply_block, {
    gw         <- .resolve_gwas(input$gwas_file, input$gwas_builtin)
    build_disp <- .norm_build_disp(input$build)
    build_hdl  <- .norm_build_hdl(build_disp)
    anc        <- .norm_pop(input$ancestry)
    paths      <- .cache_paths_abs(gw$name, build_disp, anc)
    
    # 若最终成品已存在，直接读出（并按 chr, blk 排序）
    if (.is_nonempty_tsv(paths$blocks)) {
      dt <- data.table::fread(paths$blocks, sep="\t", header=TRUE, data.table=TRUE)
      if ("blk" %in% names(dt)) {
        suppressWarnings(dt[, blk := as.integer(blk)])
        if ("chr" %in% names(dt)) {
          data.table::setorder(dt, chr, blk)
        } else {
          data.table::setorder(dt, blk)
        }
      }
      blocks_rv(dt[])
      return(invisible(NULL))
    }
    
    
    # 准备 H5 映射
    h5map <- resolve_ld_h5_files(LD_BASE_DIR, build_hdl, anc)
    validate(need(!all(is.na(h5map)), sprintf("HDF5 not found: please check LD directory, build=%s, ancestry=%s", build_hdl, anc)))
    
    # 估计 Neff
    Neff <- tryCatch(.infer_neff_from_dt(gwas_raw()), error=function(e) NA_real_)
    if (!is.finite(Neff) || Neff <= 0) Neff <- 1e5
    
    # 逐 chr 计算（跳过已有 .ok 的）
    for (i in 1:22) {
      chr_h5 <- h5map[[i]]
      if (is.na(chr_h5) || !file.exists(chr_h5)) {
        message(sprintf("[Missing HDF5] %s/%s chr%d", build_hdl, anc, i))
        next
      }
      pc <- .chr_cache_paths_abs(gw$name, build_disp, anc, i)
      
      if (.file_nonempty(pc$tsv)) {
        message(sprintf("[skip] CHR%02d already completed (resume supported)", i))
        next
      }
      
      message(sprintf("[run ] CHR%02d start...", i))
      chr_res <- tryCatch(
        rank_blocks_topN_HDL(
          h5file        = chr_h5,
          sumstats_path = gw$path,
          N_eff         = Neff,
          top_n         = .Machine$integer.max,  # 取全部，再在外面筛
          harmonize     = TRUE,
          prefilter     = FALSE                  # 如需预筛可改 TRUE
        ),
        error = function(e) {
          message(sprintf("[HDL error] %s/%s chr%d: %s",
                          build_hdl, anc, i, conditionMessage(e)))
          NULL
        }
      )
      
      
      if (is.null(chr_res) || !nrow(chr_res)) {
        message(sprintf("[warn] CHR%02d no usable block; skipping write", i))
        next
      }
      
      
      data.table::setorder(chr_res, chr, start, end)
      chr_res[, blk_h5 := as.integer(sub("^/?blk[_ ]?(\\d+)$", "\\1", blk))]
      out_chr <- chr_res[, .(
        blk_h5,                                   # <— HDF5 组编号保存在 blk_h5
        `lead snp` = lead_snp,
        chr,
        pos = sprintf("%.2f–%.2f Mb", as.numeric(start)/1e6, as.numeric(end)/1e6),
        p   = as.numeric(p_min),
        h2  = as.numeric(h2),
        se  = as.numeric(se)
      )]
      
      
      .atomic_write_tsv(out_chr, pc$tsv)
      message(sprintf("[done] CHR%02d wrote %s", i, basename(pc$tsv)))
    }
    
    # 合并所有已有 chrXX.tsv（无论本轮是否新算）
    cache_dir <- dirname(paths$blocks)
    patt <- sprintf("^%s_%s_chr[0-9]{2}\\.tsv$", build_disp, anc)
    chr_files <- list.files(cache_dir, pattern = patt, full.names = TRUE)
    if (!length(chr_files)) {
      showNotification("No per-chromosome results to merge: please check HDF5 and sumstats alignment.", type = "error", duration = 10)
      blocks_rv(data.table::data.table(blk_h5=integer(), snp=character(), chr=integer(), pos=character(), p=numeric(), h2=numeric()))
      return(invisible(NULL))
    }
    
    parts <- lapply(chr_files, function(f){
      tryCatch(data.table::fread(f, sep = "\t", header = TRUE, data.table = TRUE),
               error = function(e) NULL)
    })
    parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, parts)
    if (!length(parts)) {
      showNotification("All per-chromosome results are empty; cannot merge.", type="error", duration=10)
      blocks_rv(data.table::data.table(blk_h5=integer(), snp=character(), chr=integer(), pos=character(), p=numeric(), h2=numeric()))
      return(invisible(NULL))
    }
    
    merged <- data.table::rbindlist(parts, fill = TRUE, use.names = TRUE)
    # 解析 merged$pos = "a–b Mb" 为 bp 区间
    bp <- .parse_blk_pos_bp(merged$pos)
    merged[, `:=`(START = as.integer(round(bp$start_bp)),
                  END   = as.integer(round(bp$end_bp)))]
    
    # 读 bed 并建 GRanges
    bed <- .read_ld_blocks_bed(LD_BLOCK_BED)
    
    # 用 GenomicRanges 做区间重叠
    suppressWarnings({
      gr_m <- GenomicRanges::GRanges(
        seqnames = as.integer(merged$chr),
        ranges   = IRanges::IRanges(start = merged$START, end = merged$END)
      )
      gr_b <- GenomicRanges::GRanges(
        seqnames = as.integer(bed$CHR),
        ranges   = IRanges::IRanges(start = bed$START, end = bed$END)
      )
    })
    hits <- GenomicRanges::findOverlaps(gr_m, gr_b, ignore.strand = TRUE)
    
    # 如果有多对多，选“重叠长度最大”的那个 bed 块
    pick_best <- function(iq, sq){
      # 返回每个 query (merged 行) 对应的单个 subject (bed 行)
      dt <- data.table::data.table(q=iq, s=sq)
      dt[, ovl := as.integer(pmin(merged$END[q], bed$END[s]) - pmax(merged$START[q], bed$START[s]) + 1L)]
      dt <- dt[ovl > 0]
      if (!nrow(dt)) return(integer())
      dt[order(q, -ovl)][, .SD[1], by = q]$s
    }
    best_s <- pick_best(S4Vectors::queryHits(hits), S4Vectors::subjectHits(hits))
    
    # 先全部设 NA，然后填入 bed 的块编号
    merged[, blk := NA_character_]
    if (length(best_s)) {
      qset <- unique(S4Vectors::queryHits(hits))
      sset <- best_s
      merged$blk[qset] <- bed$BLK[sset]
    }
    
    # 排个序：chr, blk（注意 blk 可能是字符；如果是纯数字字符串可以转成整数）
    suppressWarnings({
      blk_num <- as.integer(merged$blk)
      if (all(is.finite(blk_num[!is.na(blk_num)]))) {
        merged[, blk := blk_num]
      }
    })
    if ("blk" %in% names(merged)) data.table::setorder(merged, chr, blk)
    
    # 排序与去重（同一块可能多次出现，按 p/h2 选更优）
    merged[, ord := data.table::frank(list(chr, p, -h2), ties.method="first")]
    data.table::setorder(merged, chr, p, -h2)
    merged <- unique(merged[, -"ord"], by = c("blk","lead snp","chr","pos"))
    
    # ---- 汇总全基因组 HDL block 的 h2/SE，并写入 JSON（QQ 右下角会读）----
    h2_all <- suppressWarnings(sum(as.numeric(merged$h2), na.rm = TRUE))
    se_all <- suppressWarnings(sqrt(sum((as.numeric(merged$se))^2, na.rm = TRUE)))
    if (!is.finite(h2_all)) h2_all <- NA_real_
    if (!is.finite(se_all) || se_all <= 0) se_all <- NA_real_
    
    z_all <- if (is.finite(h2_all) && is.finite(se_all)) h2_all / se_all else NA_real_
    p_all <- if (is.finite(z_all)) 2 * pnorm(-abs(z_all)) else NA_real_
    
    h2_total_rv(list(h2 = h2_all, se = se_all, z = z_all, p = p_all))
    hg_path <- file.path(dirname(paths$blocks), "H2_TOTAL.json")
    jsonlite::write_json(list(h2 = h2_all, se = se_all, z = z_all, p = p_all),
                         path = hg_path, auto_unbox = TRUE)
    
    
    # —— 代表性 block 筛选 —— 
    thr <- .to_num(input$p_threshold, 5e-8); if (!is.finite(thr)) thr <- 5e-8
    rep_blocks <- collapse_blocks_window(merged, pthr = thr, window_kb = 500)
    # —— 代表集按 blk 升序（先 chr 后 blk）——
    if ("blk" %in% names(rep_blocks)) {
      data.table::setorder(rep_blocks, chr, blk)
    }
    full_path <- file.path(dirname(paths$blocks), sprintf("ALL_blocks_%s_%s.tsv", build_disp, anc))
    .atomic_write_tsv(merged, full_path)
    
    gref <- genes_gr_rv()
    if (!is.null(gref) && nrow(rep_blocks) > 0) {
      # 1) 解析 "a–b Mb" 得到区间中心（bp）
      bp <- .parse_blk_pos_bp(rep_blocks$pos)   # 已在前面定义过
      center_bp <- suppressWarnings(as.integer(bp$center_bp))
      
      # 2) 染色体数字转字符（X/Y/MT）
      chr_chr <- as.character(rep_blocks$chr)
      chr_chr[chr_chr %in% c("23","x","X")] <- "X"
      chr_chr[chr_chr %in% c("24","y","Y")] <- "Y"
      chr_chr[chr_chr %in% c("25","mt","MT","M")] <- "MT"
      
      # 3) 构造点位并限制到注释里存在的染色体，避免 seqlevels 警告
      suppressWarnings({
        pts <- GenomicRanges::GRanges(
          seqnames = chr_chr,
          ranges   = IRanges::IRanges(start = center_bp, width = 1)
        )
      })
      pts <- GenomeInfoDb::keepSeqlevels(
        pts, intersect(GenomeInfoDb::seqlevels(pts), GenomeInfoDb::seqlevels(gref)),
        pruning.mode = "coarse"
      )
      
      nearest_gene <- rep(NA_character_, nrow(rep_blocks))
      if (length(pts) > 0) {
        hits <- GenomicRanges::distanceToNearest(pts, gref, ignore.strand = TRUE)
        if (length(hits) > 0) {
          qh <- S4Vectors::queryHits(hits)
          sh <- S4Vectors::subjectHits(hits)
          nearest_gene[qh] <- S4Vectors::mcols(gref)$gene[sh]
          # 如需距离可加：dist_kb <- round(S4Vectors::mcols(hits)$distance/1000, 1)
        }
      }
      
      rep_blocks[, `nearest gene` := nearest_gene]
      
      # 放在 h2 后（兼容是否含 se 列）
      if ("se" %in% names(rep_blocks)) rep_blocks[, se := NULL]
      want_cols <- c("blk","lead snp","chr","pos","p","h2","nearest gene")
      data.table::setcolorder(rep_blocks, intersect(want_cols, names(rep_blocks)))
    }  
    # —— 写盘前锁定排序 —— 
    if ("blk" %in% names(rep_blocks)) {
      suppressWarnings(rep_blocks[, blk := as.integer(blk)])
      if ("chr" %in% names(rep_blocks)) {
        data.table::setorder(rep_blocks, chr, blk)
      } else {
        data.table::setorder(rep_blocks, blk)
      }
    }
    
    # 渲染代表集
    .atomic_write_tsv(rep_blocks, paths$blocks)
    blocks_rv(rep_blocks)
    
    showNotification(sprintf("Blocks: %d → %d (window=%dkb)", 
                             nrow(merged), nrow(rep_blocks), 500),
                     type="message", duration=6)
    
  }, ignoreInit = TRUE)
  output$signal_block_table <- DT::renderDT({
    dat <- blocks_rv()
    req(!is.null(dat) && nrow(dat) > 0, cancelOutput = TRUE)
    
    # 先统一排序：只按 blk（若有 chr 就 chr, blk）
    if ("blk" %in% names(dat)) suppressWarnings(dat[, blk := as.integer(blk)])
    if (all(c("chr","blk") %in% names(dat))) {
      data.table::setorder(dat, chr, blk)
    } else if ("blk" %in% names(dat)) {
      data.table::setorder(dat, blk)
    }
    
    # 若存在 blk_h5，则隐藏该列
    cols <- names(dat)
    hide_idx <- which(cols == "blk_h5") - 1L  # DT 列索引从 0 开始
    
    DT::datatable(
      dat,
      rownames = FALSE, escape = FALSE,
      class = "compact stripe hover match-rows",
      selection = "single",
      options = list(
        paging = TRUE, pageLength = 10, lengthChange = FALSE,
        searching = FALSE, info = FALSE, dom = "tp", autoWidth = FALSE,
        order = list(list(0, 'asc')),
        columnDefs = c(
          if (length(hide_idx) == 1) list(list(visible = FALSE, targets = hide_idx))
        )
      ),
      width = "100%"
    ) |>
      DT::formatSignif(columns = "p", digits = 3) |>
      DT::formatRound(columns = "h2", digits = 4)
  })
  
  # 6）locuszoom
  observeEvent(input$apply_block, {
    dt_std <- gwas_raw()
    pk <- .top_peak(dt_std)
    .center(c(pk, list(title="LocusZoom", dt=dt_std)))
    # 清空选中，让“默认最强”生效
    try(DT::selectRows(DT::dataTableProxy("top_table"), NULL), silent=TRUE)
  }, ignoreInit=TRUE)
  observeEvent(input$top_table_rows_selected, {
    idx <- input$top_table_rows_selected
    if (length(idx)==1) {
      dat <- loci_rv(); req(!is.null(dat) && nrow(dat)>=idx)
      chr0 <- as.integer(dat$chr[idx])
      pos0 <- suppressWarnings(as.numeric(sub(" .*", "", dat$pos[idx]))) * 1e6
      dt_std <- gwas_raw()
      .center(list(chr=chr0, pos=pos0, title="LocusZoom", dt=dt_std))
    }
  }, ignoreInit=TRUE)
  output$locuszoom_plot <- renderPlot({
    cen <- .center(); req(!is.null(cen$dt), is.finite(cen$chr), is.finite(cen$pos))
    .plot_locus(cen$dt, cen$chr, cen$pos, win_bp=200e3, title_prefix=cen$title)
  }, height=480, res=120)
  outputOptions(output, "locuszoom_plot", suspendWhenHidden=FALSE)
  observeEvent(input$gwas_builtin, {
    dt_std <- gwas_raw()
    pk <- .top_peak(dt_std)
    .center(list(chr = pk$chr, pos = pk$pos, title = "LocusZoom", dt = dt_std))
  }, ignoreInit = FALSE)
  
  # 7) blockzoom
  MAX_SNP_PER_BLOCK <- 100
  output$blockzoom_title <- renderText("BlockZoom")
  myplot <- function(
    gwas,
    ld,
    gap = 0.2,
    r = 0.10,
    h1 = 1.0,
    h2 = 1.5,
    col = c("darkgreen", "yellow", "red"),
    show_legend = TRUE
  ){
    if(!is.list(ld)) stop("blocks of ld must be provided via list!")
    n_intv <- length(ld)
    if(n_intv < 2) stop("at leat 2 blocks should be provided!")
    if(length(unique(gwas[, 1])) != nrow(gwas)) stop("duplicate names of markers detected in gwas data!")
    
    # 清理位置列
    suppressWarnings(non_dgt_pos <- is.na(as.numeric(gwas[, 3])))
    if(sum(non_dgt_pos)) gwas <- gwas[!non_dgt_pos, , drop = FALSE]
    suppressWarnings(numeric.chr <- as.numeric(gwas[, 2]))
    suppressWarnings(max.chr <- max(numeric.chr, na.rm = TRUE))
    if(is.infinite(max.chr)) max.chr <- 0
    suppressWarnings(map.xy.index <- which(!numeric.chr %in% c(0:max.chr)))
    if(length(map.xy.index) != 0){
      chr.xy <- unique(gwas[map.xy.index, 2])
      for(i in 1:length(chr.xy)){
        gwas[gwas[, 2] == chr.xy[i], 2] <- max.chr + i
      }
    }
    gwas <- gwas[order(as.numeric(gwas[, 2]), as.numeric(gwas[, 3])), ]
    
    # 与各 block 对齐（按 LD 列名的顺序抽 GWAS 行，保证一一同序）
    gwas_list <- vector("list", n_intv)
    for (i in seq_len(n_intv)) {
      ld_id <- colnames(ld[[i]])
      if (length(unique(ld_id)) != length(ld_id)) {
        stop(sprintf("第 %d 个 block 的 LD 列名存在重复", i))
      }
      idx <- match(ld_id, gwas[, 1])
      if (anyNA(idx)) {
        miss <- ld_id[is.na(idx)]
        stop(sprintf("第 %d 个 block：GWAS 缺少 %d 个 LD 标记，例如：%s",
                     i, length(miss), paste(head(miss, 3), collapse = ", ")))
      }
      gwas_list[[i]] <- gwas[idx, , drop = FALSE]
    }
    
    col_fun <- colorRamp(col)
    line_cross <- function(l1, l2){
      x1 <- l1[1]; y1 <- l1[2]; x2 <- l1[3]; y2 <- l1[4]
      x3 <- l2[1]; y3 <- l2[2]; x4 <- l2[3]; y4 <- l2[4]
      denom <- (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
      px <- ((x1*y2 - y1*x2) * (x3 - x4) - (x1 - x2) * (x3*y4 - y3*x4)) / denom
      py <- ((x1*y2 - y1*x2) * (y3 - y4) - (y1 - y2) * (x3*y4 - y3*x4)) / denom
      c(px, py)
    }
    
    .clamp_p <- function(v){
      v <- suppressWarnings(as.numeric(v))
      v <- scale_p(v, min_p = 1e-324, max_log10 = 50)
      v[v >= 1] <- 1 - 1e-16
      v
    }
    
    min_pval <- min(sapply(gwas_list, function(x) {
      px <- .clamp_p(x[, ncol(x)])
      min(px, na.rm = TRUE)
    }))
    
    min_pval_pos <- sapply(gwas_list, function(x){
      px <- .clamp_p(x[, ncol(x)])
      which.min(px)
    })
    
    gwas_pval_log_scale <- lapply(gwas_list, function(x){
      px <- .clamp_p(x[, ncol(x)])
      -log10(px) * h2 / -log10(min_pval)
    })
    
    
    base_ang <- 1 / n_intv
    pts_x <- vector("list", n_intv)
    pts_y <- vector("list", n_intv)
    pts_cross_x <- vector("list", n_intv)
    pts_cross_y <- vector("list", n_intv)
    
    for(i in 1:n_intv){
      r_x <- r * sin(2*pi*((i - 0.5) * base_ang))
      r_y <- r * cos(2*pi*((i - 0.5) * base_ang))
      
      base_pts <- seq(0, h1, length = nrow(gwas_list[[i]]) + 1)
      # 两边的径向边
      x1 <- r_x + sin(2*pi*((i - 1) * base_ang)) * base_pts
      y1 <- r_y + cos(2*pi*((i - 1) * base_ang)) * base_pts
      x2 <- r_x + sin(2*pi*(i * base_ang)) * base_pts
      y2 <- r_y + cos(2*pi*(i * base_ang)) * base_pts
      
      pts_x[[i]] <- list(); pts_y[[i]] <- list()
      pts_x[[i]][[1]] <- x1; pts_y[[i]][[1]] <- y1
      pts_x[[i]][[2]] <- x2; pts_y[[i]][[2]] <- y2
      
      # 圆弧上的等距点
      base_ang_i <- seq(0, base_ang, length = nrow(gwas_list[[i]]) + 1)
      x_arc <- r_x + sin(2*pi*(base_ang_i + (i - 1) * base_ang)) * h1
      y_arc <- r_y + cos(2*pi*(base_ang_i + (i - 1) * base_ang)) * h1
      pts_x[[i]][[3]] <- x_arc
      pts_y[[i]][[3]] <- y_arc
      
      # 外圈（高度由 -log10 P 决定）
      x_out <- r_x + sin(2*pi*(base_ang_i[-1] + (i - 1) * base_ang - base_ang_i[2] / 2)) *
        (h1 + gap + gwas_pval_log_scale[[i]])
      y_out <- r_y + cos(2*pi*(base_ang_i[-1] + (i - 1) * base_ang - base_ang_i[2] / 2)) *
        (h1 + gap + gwas_pval_log_scale[[i]])
      pts_x[[i]][[4]] <- x_out
      pts_y[[i]][[4]] <- y_out
      
      # 预计算多边形交点
      pts_cross_x[[i]] <- vector("list", length(base_ang_i))
      pts_cross_y[[i]] <- vector("list", length(base_ang_i))
      for(j in 1:(length(base_ang_i) - 2)){
        if(j == 1){
          pts_cross_x[[i]][[j]] <- rev(x1)
          pts_cross_y[[i]][[j]] <- rev(y1)
        } else {
          cx <- x_arc[j]; cy <- y_arc[j]
          l1 <- c(x_arc[j], y_arc[j], x2[j], y2[j])
          for(k in j:(length(base_ang_i) - 2)){
            l2 <- c(x_arc[k + 1], y_arc[k + 1], x1[length(base_ang_i) - k], y1[length(base_ang_i) - k])
            pxy <- line_cross(l1, l2)
            cx <- c(cx, pxy[1]); cy <- c(cy, pxy[2])
          }
          cx <- c(cx, x2[j]); cy <- c(cy, y2[j])
          pts_cross_x[[i]][[j]] <- cx
          pts_cross_y[[i]][[j]] <- cy
        }
      }
      pts_cross_x[[i]][[length(base_ang_i) - 1]] <- c(x_arc[length(base_ang_i) - 1], x2[length(base_ang_i) - 1])
      pts_cross_y[[i]][[length(base_ang_i) - 1]] <- c(y_arc[length(base_ang_i) - 1], y2[length(base_ang_i) - 1])
      pts_cross_x[[i]][[length(base_ang_i)]] <- x_arc[length(base_ang_i)]
      pts_cross_y[[i]][[length(base_ang_i)]] <- y_arc[length(base_ang_i)]
    }
    
    # —— 零边距防报错 —— 
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mai = rep(0,4), omi = rep(0,4), xaxs = "i", yaxs = "i", xpd = NA)
    
    plot.new()
    xr <- c(-r - h1 - h2 - gap, r + h1 + h2 + gap)
    yr <- xr
    plot.window(xlim = xr, ylim = yr, asp = 1)
    
    # 多边形+散点
    for(i in 1:n_intv){
      index <- match(gwas_list[[i]][, 1], colnames(ld[[i]]))
      # 多边形填充（LD 色带）
      for(j in 1:(nrow(gwas_list[[i]]) - 1)){
        for(k in 1:(length(pts_cross_x[[i]][[j]]) - 2)){
          polygon(
            c(pts_cross_x[[i]][[j]][k + 1],
              pts_cross_x[[i]][[j + 1]][k],
              pts_cross_x[[i]][[j + 1]][k + 1],
              pts_cross_x[[i]][[j]][k + 2]),
            c(pts_cross_y[[i]][[j]][k + 1],
              pts_cross_y[[i]][[j + 1]][k],
              pts_cross_y[[i]][[j + 1]][k + 1],
              pts_cross_y[[i]][[j]][k + 2]),
            col = rgb(col_fun(abs(ld[[i]][index[j], index[k + 1]])) / 255),
            border = "gray", lwd = 0.1
          )
        }
      }
      # 外圈散点（用该块最显著位点那一列着色）
      points(pts_x[[i]][[4]], pts_y[[i]][[4]], pch = 16,
             col = rgb(col_fun(abs(ld[[i]][, index[min_pval_pos[i]]])) / 255),
             cex = 0.6)
    }
    
    # 小型色标（不使用 image）
    if (isTRUE(show_legend)) {
      w <- diff(xr); h <- diff(yr)
      # 调整这两行，控制色标的宽度
      x0 <- xr[1] + 0.86*w   # 起点更靠右
      x1 <- xr[1] + 0.92*w   # 终点保持
      y0 <- yr[1] + 0.62*h
      y1 <- yr[1] + 0.90*h
      
      rect(x0, y0, x1, y1, border = "grey40", lwd = 0.6)
      nstrip <- 120
      yy <- seq(y0, y1, length.out = nstrip + 1)
      cc <- colorRampPalette(col)(nstrip)
      for (ii in seq_len(nstrip)) rect(x0, yy[ii], x1, yy[ii+1], col = cc[ii], border = NA)
      
      ticky <- pretty(c(0,1), n = 5); ticky <- ticky[ticky >= 0 & ticky <= 1]
      ytick <- y0 + (y1 - y0) * ticky
      segments(x1, ytick, x1 + 0.01*w, ytick, col = "grey25", lwd = 0.6)
      text(x1 + 0.018*w, ytick, labels = format(ticky), adj = c(0, 0.5), cex = 0.7)
    }
    
  }
  
  .block_h5_for_chr <- function(chr_i){
    build_hdl <- .norm_build_hdl(.norm_build_disp(input$build))
    anc       <- .norm_pop(input$ancestry)
    h5map     <- resolve_ld_h5_files(LD_BASE_DIR, build_hdl, anc)
    req(!all(is.na(h5map)), "未找到 LD HDF5 映射")
    h5map[[ sprintf("chr%d", as.integer(chr_i)) ]]
  }
  .read_block_aligned <- function(h5file, blk_group, gwas_dt, topn = MAX_SNP_PER_BLOCK){
    on.exit(rhdf5::h5closeAll(), add = TRUE)
    
    # 读 snplist
    sp  <- rhdf5::h5read(h5file, file.path(blk_group, "snplist"))
    ref <- if (is.character(sp)) data.table::data.table(rsid = sp) else data.table::as.data.table(sp)
    
    rn  <- .pick(names(ref), c("SNP","RSID","ID","rsid","id"))
    validate(need(!is.na(rn), "snplist 缺少 SNP 列"))
    data.table::setnames(ref, rn, "rsid", skip_absent = TRUE)
    ref[, rsid := as.character(rsid)]
    ref[, SNP_CLEAN := .clean_snp(rsid)]
    
    # GWAS 数据
    dt <- data.table::as.data.table(gwas_dt)
    sc <- .pick(names(dt), c("SNP","ID","RSID","MARKER","MARKERNAME","VARIANT","VARIANT_ID"))
    if (is.na(sc)) stop("GWAS 缺少 SNP/ID 列以对齐 LD")
    dt[, SNP_CLEAN := .clean_snp(get(sc))]
    data.table::setkey(dt, SNP_CLEAN)
    
    # ① 先用 rsID/标准化ID 对齐
    idx  <- match(ref$SNP_CLEAN, dt$SNP_CLEAN)
    keep <- which(!is.na(idx))
    
    # ② 若太少，回退到 chr:pos 对齐（snplist 里解析）
    if (length(keep) < 2) {
      cp <- .parse_chrpos_from_str(ref$rsid)  # 确保上面已定义该函数
      have_dt_cp <- all(c("CHR","POS") %in% names(dt))
      if (have_dt_cp && any(is.finite(cp$CHR) & is.finite(cp$POS))) {
        ref[, `:=`(CHR = cp$CHR, POS = cp$POS)]
        ref[, KEY_CP := paste0(CHR, ":", POS)]
        dt[,  KEY_CP := paste0(as.integer(CHR), ":", as.integer(POS))]
        idx  <- match(ref$KEY_CP, dt$KEY_CP)
        keep <- which(!is.na(idx))
      }
    }
    validate(need(length(keep) >= 2, paste0("块 ", blk_group, " 与 GWAS 交集过小")))
    
    # 取 LD 子矩阵
    R <- as.matrix(rhdf5::h5read(h5file, file.path(blk_group, "ldblk")))
    validate(need(nrow(R) == ncol(R), "LD 不是方阵"))
    R <- R[keep, keep, drop = FALSE]
    
    # 统一命名：让 GWAS$SNP 与 LD 列名一致
    rs <- as.character(ref$SNP_CLEAN[keep])
    colnames(R) <- rownames(R) <- rs
    
    g_sub <- dt[idx[keep], .(CHR, POS, P)]
    g_sub[, `:=`(
      CHR = suppressWarnings(as.integer(CHR)),
      POS = suppressWarnings(as.numeric(POS)),
      P   = suppressWarnings(as.numeric(P))
    )]
    g_sub <- g_sub[is.finite(CHR) & is.finite(POS) & is.finite(P)]
    g_sub[, SNP := rs]
    data.table::setcolorder(g_sub, c("SNP","CHR","POS","P"))
    
    # 只保留最显著的前 topn，并同步裁剪 LD
    g_sub <- g_sub[is.finite(P)]
    data.table::setorder(g_sub, P)
    if (is.finite(topn) && topn > 0 && nrow(g_sub) > topn) g_sub <- g_sub[1:topn]
    sel <- match(g_sub$SNP, rs)
    R   <- R[sel, sel, drop = FALSE]
    
    validate(need(nrow(g_sub) >= 2, paste0("块 ", blk_group, " 剩余 SNP < 2")))
    list(gwas = g_sub, ld = R)
  }
  .h5_block_ids_for_chr <- function(h5file){
    on.exit(rhdf5::h5closeAll(), add = TRUE)
    info <- tryCatch(rhdf5::h5ls(h5file, all = TRUE), error = function(e) NULL)
    if (is.null(info) || !nrow(info)) return(integer())
    g <- as.character(info$group)
    # 只保留形如 "/blk_38" 的分组，抽出数字 38
    ids <- suppressWarnings(as.integer(sub("^/blk[_ ]?(\\d+)$", "\\1", g)))
    ids <- ids[is.finite(ids)]
    sort(unique(ids))
  }
  .prepare_triplet <- function(row_dt){
    validate(need(!is.null(row_dt) && nrow(row_dt)==1, "未找到所选 block"))
    chr_i <- as.integer(row_dt$chr[1])
    
    # —— 取行里的块编号：优先 HDF5 组编号 blk_h5，其次退回 bed 编号 blk —— 
    blk_raw <- if ("blk_h5" %in% names(row_dt)) row_dt$blk_h5[1]
    else if ("blk" %in% names(row_dt)) row_dt$blk[1]
    else NA
    if (is.numeric(blk_raw) || grepl("^[0-9]+$", as.character(blk_raw))) {
      blk_id <- as.integer(blk_raw)
    } else {
      blk_id <- suppressWarnings(as.integer(sub("^/?blk[_ ]?", "", as.character(blk_raw))))
    }
    validate(need(is.finite(blk_id), "无法解析所选 block 编号"))
    
    # —— 用 H5 的实际组编号来找邻居 —— 
    h5 <- .block_h5_for_chr(chr_i); validate(need(file.exists(h5), paste0("HDF5 不存在：", h5)))
    all_ids <- .h5_block_ids_for_chr(h5)
    validate(need(length(all_ids) > 0, sprintf("CHR%s 在 H5 中未找到任何块", chr_i)))
    validate(need(blk_id %in% all_ids, sprintf("H5 中不存在 /blk_%d（编号体系可能不同）", blk_id)))
    
    n <- length(all_ids)
    center_pos <- match(blk_id, all_ids)
    
    # 先取中心±1，并做边界裁剪
    lo <- max(1, center_pos - 1)
    hi <- min(n, center_pos + 1)
    cand_pos <- lo:hi
    
    # 若不足 3 个，再向两侧补充；用条件判断避免 seq 的方向冲突
    if (length(cand_pos) < 3) {
      extra_left  <- if (lo > 1) seq(lo - 1, 1, by = -1) else integer(0)
      extra_right <- if (hi < n) seq(hi + 1, n, by = 1)  else integer(0)
      cand_pos <- unique(c(cand_pos, extra_left, extra_right))
    }
    
    neigh_ids <- all_ids[cand_pos][1:min(3, n)]
    validate(need(length(neigh_ids) >= 2, "相邻块不足以绘三连块"))
    
    dt_gwas <- gwas_raw()
    
    parts <- lapply(neigh_ids, function(b){
      gname <- sprintf("/blk_%d", as.integer(b))
      tryCatch(.read_block_aligned(h5, gname, dt_gwas), error = function(e) NULL)
    })
    # 过滤掉异常/过小块，至少要有 2 个 SNP 才能画
    parts <- Filter(function(x) is.list(x) && is.data.frame(x$gwas) && nrow(x$gwas) >= 2, parts)
    validate(need(length(parts) >= 2, "可视化需要至少两个有效相邻块"))
    
    g_merge <- data.table::rbindlist(lapply(parts, `[[`, "gwas"), use.names = TRUE, fill = TRUE)
    g_merge <- unique(g_merge, by = "SNP")
    ld_list <- lapply(parts, `[[`, "ld")
    
    list(
      gwas = g_merge,
      ld_list = ld_list,
      info  = data.table::data.table(chr = chr_i, center = blk_id, triplet = list(neigh_ids))
    )
  }
  .draw_blockzoom <- function(row_dt){
    trip <- .prepare_triplet(row_dt)
    gdf  <- as.data.frame(trip$gwas)
    validate(need(ncol(gdf) >= 4, "GWAS 合并失败：列不足"))
    
    info <- trip$info
    nb   <- paste(info$triplet[[1]], collapse = ",")  # H5 邻居，仅用于内部参考
    
    # —— 新：标题用 BED 的 blk —— 
    bed_blk <- if ("blk" %in% names(row_dt)) as.character(row_dt$blk[1]) else NA
    ttl  <- if (is.na(bed_blk)) {
      sprintf("BlockZoom · CHR%s : block ? (HDF5 center %s; neighbors %s)",
              info$chr[1], info$center[1], nb)
    } else {
      sprintf("BlockZoom · CHR%s : blk %s (neighbors by HDF5: %s)",
              info$chr[1], bed_blk, nb)
    }
    output$blockzoom_title <- shiny::renderText(ttl)
    
    myplot(gdf, trip$ld_list, h1 = 2, r = 0.06)
  }
  observeEvent(blocks_rv(), {
    dat <- blocks_rv()
    req(!is.null(dat), nrow(dat) > 0)
    data.table::setorder(dat, chr, blk) 
    
    output$blockzoom_plot <- renderPlot({
      .draw_blockzoom(dat[1, ])
    }, res = 144)
    
    try(DT::selectRows(DT::dataTableProxy("signal_block_table"), 1), silent = TRUE)
  }, ignoreInit = FALSE)
  observeEvent(input$signal_block_table_rows_selected, {
    sel <- input$signal_block_table_rows_selected
    dat <- blocks_rv()
    if (is.null(dat) || !length(sel) || sel[1] < 1 || sel[1] > nrow(dat)) return(invisible(NULL))
    row_dt <- dat[sel[1], ]
    output$blockzoom_plot <- renderPlot({
      .draw_blockzoom(row_dt)
    }, res = 144)
  }, ignoreInit = TRUE)
  
}