library(rhdf5)
library(data.table)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)

generate_locuszoom_by_rsid <- function(
    rsid,
    win_bp = 200e3,
    h5_dir,
    snpinfo_file,
    gwas_dt = NULL,             # ✅新增：优先用内存 GWAS（Shiny用）
    gwas_file = NULL,           # 可选：仍可支持文件模式
    gwas_snp_col = "SNP",
    gwas_p_col   = "p",
    ens_db = EnsDb.Hsapiens.v75
) {
  
  
  rsid <- toupper(trimws(as.character(rsid)))
  
  ## 1) snpinfo 定位 rsid -> chr0/pos0
  snpinfo <- fread(snpinfo_file)
  setnames(snpinfo, c("CHR","SNP","BP"), c("CHR","SNP","POS"), skip_absent = TRUE)
  
  if (!all(c("CHR","SNP","POS") %in% names(snpinfo))) {
    stop("snpinfo 缺少必要列：CHR/SNP/POS（或 BP）")
  }
  
  snpinfo[, CHR := suppressWarnings(as.integer(CHR))]
  snpinfo[, POS := suppressWarnings(as.integer(POS))]
  snpinfo[, SNP := toupper(trimws(as.character(SNP)))]
  
  hit <- snpinfo[SNP == rsid]
  
  if (nrow(hit) == 0) stop("snpinfo 中找不到该 rsid：", rsid)
  if (nrow(hit) > 1) {
    # 极少数情况下同一 rsid 可能重复，取第一个并提示
    message("警告：snpinfo 中 rsid 重复，取第一条记录")
    hit <- hit[1]
  }
  
  chr0 <- hit$CHR[1]
  pos0 <- hit$POS[1]
  
  message(sprintf("从 snpinfo 定位：%s -> chr%d:%d", rsid, chr0, pos0))
  
  ## 2) 窗口 SNP
  left  <- max(1L, pos0 - as.integer(win_bp))
  right <- pos0 + as.integer(win_bp)
  
  win <- snpinfo[CHR == chr0 & POS >= left & POS <= right, .(CHR, POS, SNP)]
  if (nrow(win) == 0) stop("窗口内没有 SNP（不太可能）")
  
  ## 3) 找 lead block
  h5file <- file.path(h5_dir, paste0("ldblk_ukbb_chr", chr0, ".hdf5"))
  if (!file.exists(h5file)) stop("找不到 hdf5：", h5file)
  
  blk_ids <- h5ls(h5file, recursive = FALSE)$name
  blk_ids <- blk_ids[grepl("^blk_", blk_ids)]
  if (length(blk_ids) == 0) stop("hdf5 顶层未找到 blk_* group")
  
  lead_blk <- NA_character_
  for (b in blk_ids) {
    snps <- toupper(as.character(h5read(h5file, paste0(b, "/snplist"))))
    if (rsid %in% snps) { lead_blk <- b; break }
  }
  if (is.na(lead_blk)) stop("lead SNP 不在任何 block：", rsid)
  
  message("Lead SNP 在 ", lead_blk)
  
  ## 4) 计算 r2（lead vs block）
  snps_blk <- toupper(as.character(h5read(h5file, paste0(lead_blk, "/snplist"))))
  R <- h5read(h5file, paste0(lead_blk, "/ldblk"))
  
  idx <- which(snps_blk == rsid)
  if (length(idx) != 1) stop("在 lead block 的 snplist 中找不到或不唯一：", rsid)
  
  r2dt <- data.table(SNP = snps_blk, r2 = as.numeric(R[idx, ]^2))
  win  <- merge(win, r2dt, by = "SNP", all.x = TRUE)
  
  ## 5) P 值：优先用 gwas_dt，否则读 gwas_file，否则 P=1
  win[, P := 1]
  
  gwas <- NULL
  if (!is.null(gwas_dt) && nrow(as.data.table(gwas_dt)) > 0) {
    gwas <- as.data.table(gwas_dt)
  } else if (!is.null(gwas_file)) {
    gwas <- fread(gwas_file)
  }
  
  if (!is.null(gwas)) {
    if (!(gwas_snp_col %in% names(gwas))) stop("GWAS 缺少 SNP 列：", gwas_snp_col)
    if (!(gwas_p_col   %in% names(gwas))) stop("GWAS 缺少 P 列：", gwas_p_col)
    
    gwas[, SNP := toupper(trimws(as.character(get(gwas_snp_col))))]
    gwas[, P_gwas := suppressWarnings(as.numeric(get(gwas_p_col)))]
    gwas <- gwas[is.finite(P_gwas) & P_gwas > 0, .(SNP, P_gwas)]
    
    setkey(win, SNP)
    setkey(gwas, SNP)
    win <- gwas[win]
    win[, P := fifelse(is.finite(P_gwas), P_gwas, 1)]
    win[, P_gwas := NULL]
  }
  
  
  message(sprintf(
    "窗口SNP=%d | lead块内可上色(r2非NA)=%d | GWAS匹配到P=%d",
    nrow(win),
    sum(!is.na(win$r2)),
    sum(win$P < 1)
  ))
  
  ## 6) locuszoom
  loc <- locus(
    data = win,
    seqname = as.character(chr0),
    xrange  = c(left, right),
    ens_db  = ens_db,
    chrom   = "CHR",
    pos     = "POS",
    p       = "P",
    LD      = "r2",
    index_snp = rsid
  )
  
  locus_plot(loc)
}

# generate_locuszoom_by_rsid(
#   rsid = "rs571312",
#   win_bp = 200000,
#   h5_dir = "D:/hello/data/ldblk/ldblk_ukbb_eur",
#   snpinfo_file = "D:/hello/data/ldblk/ldblk_ukbb_eur/snpinfo_ukbb_hm3",
#   gwas_file = "D:/hello/data/gwas/bmi2015.gz"   # 可选：想要非一条线就传
# )
