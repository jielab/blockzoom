library(data.table)
library(rhdf5)

ld_path_for_chr = function(chr) {paste0(ld_dir, "/ldblk_", hapref, "_chr", chr, ".hdf5")}

.pick <- function(nms, alts){
  z <- intersect(alts, nms)
  if (length(z)) z[1] else NA_character_
}

.clean_snp <- function(x){
  x <- toupper(trimws(as.character(x)))
  is_rsid <- grepl("^RS\\d+$", x)
  y <- x
  y[!is_rsid] <- sub("^CHR", "", y[!is_rsid])
  y[!is_rsid] <- sub("^([0-9]{1,2}|X|Y|MT|M)[:|_/ ]+([0-9]+).*$",
                     "\\1:\\2", y[!is_rsid], perl = TRUE)
  y
}

.parse_chrpos_from_str <- function(x){
  s <- gsub(",", "", as.character(x))
  m <- regexec("(?i)^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[[:punct:]_ ]+([0-9]+)", s, perl = TRUE)
  mt <- regmatches(s, m)
  chr <- vapply(mt, function(v) if (length(v) >= 3) v[2] else NA_character_, "")
  pos <- vapply(mt, function(v) if (length(v) >= 3) v[3] else NA_character_, "")
  chr[toupper(chr) == "X"]  <- "23"
  chr[toupper(chr) == "Y"]  <- "24"
  chr[toupper(chr) %in% c("MT","M")] <- "25"
  list(
    CHR = suppressWarnings(as.integer(chr)),
    POS = suppressWarnings(as.integer(pos))
  )
}

## ---------------------------------------------------------
## 3. GWAS reader (SNP / CHR / POS / P)
## ---------------------------------------------------------
gwas_raw <- function(f){
  dt <- fread(
    file = f,
    header = TRUE,
    data.table = TRUE,
    na.strings = c("NA","NaN","","Inf","-Inf",".")
  )
  stopifnot(nrow(dt) > 0)

  nm0 <- names(dt)
  nmU <- toupper(trimws(sub("^\ufeff","", nm0)))
  setnames(dt, nm0, nmU, skip_absent = TRUE)

  pick1 <- function(alts){
    hit <- intersect(alts, names(dt))
    if (length(hit)) hit[1] else NA_character_
  }

  snp_col <- pick1(c("SNP","MARKER","MARKERNAME","ID","RSID","VARIANT"))
  chr_col <- pick1(c("CHR","CHROM","#CHROM","CHROMOSOME"))
  pos_col <- pick1(c("POS.38","POS38","POS","POS.37","POS37","BP",
                     "POSITION","BASE_PAIR_LOCATION"))
  p_col   <- pick1(c("P","PVAL","P_VALUE","PVALUE","NEG_LOG_PVALUE"))

  if (any(is.na(c(snp_col, chr_col, pos_col, p_col))))
    stop("GWAS 文件必须包含 SNP / CHR / POS / P（或等价列名）")

  CHR <- gsub("(?i)^CHR", "", as.character(dt[[chr_col]]))
  CHR[tolower(CHR) == "x"]  <- "23"
  CHR[tolower(CHR) == "y"]  <- "24"
  CHR[tolower(CHR) %in% c("mt","m")] <- "25"

  if (toupper(p_col) == "NEG_LOG_PVALUE") {
    P <- 10^(-suppressWarnings(as.numeric(dt[[p_col]])))
  } else {
    P <- suppressWarnings(as.numeric(dt[[p_col]]))
  }
  P[!is.finite(P) | P <= 0] <- NA_real_
  P[P >= 1] <- 1-1e-16
  P <- pmax(P, 1e-300)

  out <- data.table(
    SNP = as.character(dt[[snp_col]]),
    CHR = suppressWarnings(as.integer(CHR)),
    POS = suppressWarnings(as.numeric(dt[[pos_col]])),
    P   = P
  )
  out[is.finite(CHR) & is.finite(POS) & is.finite(P)]
}

## ---------------------------------------------------------
## 4. HDF5 block helpers
## ---------------------------------------------------------
MAX_SNP_PER_BLOCK <- 100

.read_block_aligned <- function(h5file, blk_group, gwas_dt, topn = MAX_SNP_PER_BLOCK){
  on.exit(rhdf5::h5closeAll(), add = TRUE)

  sp  <- rhdf5::h5read(h5file, file.path(blk_group, "snplist"))
  ref <- if (is.character(sp)) data.table(rsid = sp) else as.data.table(sp)

  rn  <- .pick(names(ref), c("SNP","RSID","ID","rsid","id"))
  if (is.na(rn)) stop("snplist 缺少 SNP 列")
  setnames(ref, rn, "rsid", skip_absent = TRUE)
  ref[, rsid := as.character(rsid)]
  ref[, SNP_CLEAN := .clean_snp(rsid)]

  dt <- as.data.table(gwas_dt)
  sc <- .pick(names(dt), c("SNP","ID","RSID","MARKER","MARKERNAME","VARIANT","VARIANT_ID"))
  if (is.na(sc)) stop("GWAS 缺少 SNP/ID 列对齐 LD")
  dt[, SNP_CLEAN := .clean_snp(get(sc))]
  setkey(dt, SNP_CLEAN)

  idx  <- match(ref$SNP_CLEAN, dt$SNP_CLEAN)
  keep <- which(!is.na(idx))

  if (length(keep) < 2){
    cp <- .parse_chrpos_from_str(ref$rsid)
    have_dt_cp <- all(c("CHR","POS") %in% names(dt))
    if (have_dt_cp && any(is.finite(cp$CHR) & is.finite(cp$POS))){
      ref[, `:=`(CHR = cp$CHR, POS = cp$POS)]
      ref[, KEY_CP := paste0(CHR, ":", POS)]
      dt[,  KEY_CP := paste0(as.integer(CHR), ":", as.integer(POS))]
      idx  <- match(ref$KEY_CP, dt$KEY_CP)
      keep <- which(!is.na(idx))
    }
  }
  if (length(keep) < 2) stop(blk_group, " 与 GWAS 交集过小")

  R <- as.matrix(rhdf5::h5read(h5file, file.path(blk_group, "ldblk")))
  if (nrow(R) != ncol(R)) stop("LD 不是方阵")
  R <- R[keep, keep, drop = FALSE]

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
  setcolorder(g_sub, c("SNP","CHR","POS","P"))
  g_sub <- g_sub[is.finite(P)]
  setorder(g_sub, P)

  if (is.finite(topn) && topn > 0 && nrow(g_sub) > topn)
    g_sub <- g_sub[1:topn]

  sel <- match(g_sub$SNP, rs)
  R   <- R[sel, sel, drop = FALSE]
  if (nrow(g_sub) < 2) stop(blk_group, " 剩余 SNP < 2")

  list(gwas = g_sub, ld = R)
}

.h5_block_ids_for_chr <- function(h5file){
  on.exit(rhdf5::h5closeAll(), add = TRUE)
  info <- tryCatch(rhdf5::h5ls(h5file, all = TRUE), error = function(e) NULL)
  if (is.null(info) || !nrow(info)) return(integer())
  ids <- suppressWarnings(
    as.integer(sub("^blk_","",
                   info$name[info$otype=="H5I_GROUP" & grepl("^blk_", info$name)]))
  )
  ids[is.finite(ids)]
}

## ---- prepare 3 blocks: center, left=center-1, right=center+1 ----
.prepare_triplet <- function(chr_i, blk_id, gwas_dt){
  chr_i  <- as.integer(chr_i)
  blk_id <- as.integer(blk_id)
  if (!is.finite(chr_i) || !is.finite(blk_id))
    stop("chr_i / blk_id 必须是整数")

  h5 <- ld_path_for_chr(chr_i)
  if (!file.exists(h5)) stop("HDF5 不存在：", h5)
  all_ids <- .h5_block_ids_for_chr(h5)
  if (!length(all_ids)) stop("CHR", chr_i, " 中无 block")
  if (!(blk_id %in% all_ids))
    stop("在 CHR", chr_i, " 中找不到中心块 blk_", blk_id)

  left_id  <- blk_id - 1L
  right_id <- blk_id + 1L
  if (!(left_id  %in% all_ids) ||
      !(right_id %in% all_ids))
    stop("目前要求左右相邻块都存在：", blk_id-1, ", ", blk_id+1)

  id_vec <- c(blk_id, left_id, right_id)   # 顺时针顺序：中心、左、右

  dt_gwas <- as.data.table(gwas_dt)
  parts   <- vector("list", 3L)
  for (i in 1:3){
    gname <- sprintf("/blk_%d", id_vec[i])
    parts[[i]] <- .read_block_aligned(h5, gname, dt_gwas)
  }

  g_merge <- rbindlist(lapply(parts, `[[`, "gwas"),
                       use.names = TRUE, fill = TRUE)
  g_merge <- unique(g_merge, by = "SNP")

  list(
    gwas    = g_merge,
    ld_list = lapply(parts, `[[`, "ld"),
    info    = data.table(
      chr     = chr_i,
      center  = blk_id,
      triplet = list(id_vec)      # c(center, left, right)
    )
  )
}

get_LD <- function(block_id, gwas_dt){
  block_id <- as.integer(block_id)
  if (!is.finite(block_id)) stop("block_id 必须是整数")

  found_chr <- NA_integer_
  for (chr_i in 1:26){
    h5file <- ld_path_for_chr(chr_i)
    if (!file.exists(h5file)) next
    ids <- tryCatch(.h5_block_ids_for_chr(h5file), error = function(e) integer())
    if (block_id %in% ids){ found_chr <- chr_i; break }
  }
  if (is.na(found_chr))
    stop("在任何 H5 中都没找到 block_id = ", block_id)

  .prepare_triplet(found_chr, block_id, gwas_dt)
}

## ---------------------------------------------------------
## 5. Gene list reader (optional)
## ---------------------------------------------------------
.read_glist <- function(f){
  if (is.null(f) || is.na(f)) return(NULL)
  gl <- fread(f, header = FALSE)
  if (ncol(gl) < 4) stop("glist 至少需要 4 列: CHR START END GENE")
  setnames(gl, 1:4, c("CHR","START","END","GENE"))
  gl[, CHR := suppressWarnings(as.integer(CHR))]
  gl
}

## ---------------------------------------------------------
## 6. Plot core: myplot()
## ---------------------------------------------------------
myplot <- function(gwas, ld,
                   gap = 0.2, r = 0.10, h1 = 1.0, h2 = 1.5,
                   col = c("darkgreen", "yellow", "red"),
                   block_ids = NULL,
                   gene_df = NULL){

  if (!is.list(ld)) stop("ld 必须是 list")
  n_intv <- length(ld)
  if (n_intv != 3) stop("当前版本假定正好 3 个块（中心 + 左 + 右）")

  g <- data.table::as.data.table(gwas)
  if (!all(c("SNP","CHR","POS","P") %in% names(g)))
    stop("gwas 需要列: SNP / CHR / POS / P")
  if (length(unique(g$SNP)) != nrow(g))
    stop("GWAS SNP 名有重复")
  data.table::setorder(g, CHR, POS)

  ## ---- GWAS per block (顺序: center, left, right) ----
  g_list <- vector("list", n_intv)
  for (i in seq_len(n_intv)) {
    ids <- colnames(ld[[i]])
    idx <- match(ids, g$SNP)
    if (anyNA(idx))
      stop("第 ", i, " 个 block 有 SNP 不在 GWAS 中")
    g_list[[i]] <- as.data.frame(g[idx, .(SNP, CHR, POS, P)])
  }

  ## ---- block length -> angle (clockwise) ----
  pos_range <- lapply(g_list, function(x){
    pos <- as.numeric(x$POS)
    c(min(pos, na.rm = TRUE), max(pos, na.rm = TRUE))
  })
  blk_start <- sapply(pos_range, `[`, 1)
  blk_end   <- sapply(pos_range, `[`, 2)
  blk_len   <- pmax(blk_end - blk_start, 1)

  total_len   <- sum(blk_len)
  span_vec    <- 2*pi*blk_len/total_len
  center_span <- span_vec[1]
  rem_span    <- 2*pi - center_span
  side_sum    <- blk_len[2] + blk_len[3]
  right_span  <- rem_span * blk_len[3] / side_sum
  left_span   <- rem_span * blk_len[2] / side_sum

  alpha_start  <- alpha_end <- alpha_center <- numeric(3)
  alpha_center[1] <- 0
  alpha_start[1]  <- -center_span/2
  alpha_end[1]    <-  center_span/2
  alpha_start[3]  <- alpha_end[1]
  alpha_end[3]    <- alpha_start[3] + right_span
  alpha_center[3] <- (alpha_start[3] + alpha_end[3]) / 2
  alpha_end[2]    <- alpha_start[1]
  alpha_start[2]  <- alpha_end[2] - left_span
  alpha_center[2] <- (alpha_start[2] + alpha_end[2]) / 2

  ## ---- -log10(P) scaling ----
  p_min <- max(0, p_min)
  if (p_max <= p_min) stop("p_max 必须大于 p_min")

  logp_list <- lapply(g_list, function(x){
    p  <- as.numeric(x$P)
    lp <- -log10(p)
    lp[!is.finite(lp)] <- NA_real_
    lp <- pmin(lp, p_max)
    lp[lp < p_min] <- NA_real_
    lp
  })
  raw_max   <- max(unlist(logp_list), na.rm = TRUE)
  if (!is.finite(raw_max)) raw_max <- p_max
  scale_max <- ceiling(raw_max)
  col_fun   <- colorRamp(col)
  cols_all  <- colorRampPalette(c("darkgreen","yellow","red"))(100)

  sign_ang <- -1                      # clockwise
  r0       <- h1 + gap
  Rmax     <- r0 + h2                 # 最外层浅蓝圈半径

  line_cross <- function(l1,l2){
    x1<-l1[1];y1<-l1[2];x2<-l1[3];y2<-l1[4]
    x3<-l2[1];y3<-l2[2];x4<-l2[3];y4<-l2[4]
    d<-(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
    c(((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/d,
      ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/d)
  }

  boundary_label <- function(i, angle, value, rad_offset){
    th    <- sign_ang * angle
    r_lab <- h1 + rad_offset
    x0 <- sin(th)*r_lab
    y0 <- cos(th)*r_lab

    tx <- cos(th)
    ty <- -sin(th)
    interior_sign <- if (alpha_center[i] > angle) 1 else -1

    off_tan <- 0.10 * interior_sign   # 离辐射线更远一点
    x <- x0 + off_tan*tx
    y <- y0 + off_tan*ty

    text(x, y,
         labels = sprintf("%.1f", value),
         cex = 0.7, col = "black", font = 2)
  }

  ## ---- plot frame ----
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mai = rep(0,4), omi = rep(0,4), xaxs = "i", yaxs = "i", xpd = NA)
  plot.new()
  xr <- c(-r - h1 - h2 - gap, r + h1 + h2 + gap); yr <- xr
  plot.window(xlim = xr, ylim = yr, asp = 1)

  ## ---- 3 darkblue 辐射线 ----
  rad_angles <- c(0,  2*pi/3, -2*pi/3)
  ticks <- pretty(c(p_min, scale_max), n = 4)
  ticks <- ticks[ticks >= p_min & ticks <= scale_max]

  for (a in rad_angles){
    th <- sign_ang * a
    segments(0, 0, sin(th)*Rmax, cos(th)*Rmax,
             col = "darkblue", lwd = 1.4)
    for (tk in ticks){
      rr <- r0 + (tk/scale_max)*(h2 - 0.05)  # margin inside blue circle
      points(sin(th)*rr, cos(th)*rr,
             pch = 16, col = "darkblue", cex = 0.45)
      text  (sin(th)*rr, cos(th)*rr,
             labels = tk, pos = 4, cex = 0.75, col = "black")
    }
  }

  ## ---- each block (LD + p dots + ID) ----
  for (i in 1:3){
    ang_c <- sign_ang * alpha_center[i]
    ang_l <- sign_ang * alpha_start[i]
    ang_r <- sign_ang * alpha_end[i]
    cx <- r * sin(ang_c); cy <- r * cos(ang_c)

    gg <- g_list[[i]]
    n  <- nrow(gg)
    base_pts <- seq(0, h1, length.out = n+1)

    x1 <- cx + sin(ang_l)*base_pts
    y1 <- cy + cos(ang_l)*base_pts
    x2 <- cx + sin(ang_r)*base_pts
    y2 <- cy + cos(ang_r)*base_pts

    th_seq <- seq(alpha_start[i], alpha_end[i], length.out = n+1)
    th_arc <- sign_ang*th_seq
    xa <- cx + sin(th_arc)*h1
    ya <- cy + cos(th_arc)*h1

    px <- vector("list", n+1); py <- vector("list", n+1)
    for (j in 1:(n-1)){
      if (j == 1){
        px[[j]] <- rev(x1); py[[j]] <- rev(y1)
      } else {
        xx <- xa[j]; yy <- ya[j]
        l1 <- c(xa[j], ya[j], x2[j], y2[j])
        for (k in j:(n-1)){
          l2 <- c(xa[k+1], ya[k+1], x1[n+2-k], y1[n+2-k])
          p  <- line_cross(l1, l2)
          xx <- c(xx, p[1]); yy <- c(yy, p[2])
        }
        xx <- c(xx, x2[j]); yy <- c(yy, y2[j])
        px[[j]] <- xx; py[[j]] <- yy
      }
    }
    px[[n]]   <- c(xa[n], x2[n])
    py[[n]]   <- c(ya[n], y2[n])
    px[[n+1]] <- xa[n+1]; py[[n+1]] <- ya[n+1]

    idx <- match(gg$SNP, colnames(ld[[i]]))
    for (j in 1:(n-1)){
      for (k in 1:(length(px[[j]])-2)){
        polygon(
          c(px[[j]][k+1], px[[j+1]][k], px[[j+1]][k+1], px[[j]][k+2]),
          c(py[[j]][k+1], py[[j+1]][k], py[[j+1]][k+1], py[[j]][k+2]),
          col    = rgb(col_fun(abs(ld[[i]][idx[j], idx[k+1]]))/255),
          border = "gray", lwd = 0.1
        )
      }
    }

    ## solar dots: keep inside blue circle
    lp  <- logp_list[[i]]
    pos <- gg$POS
    len <- max(pos) - min(pos); if (!is.finite(len) || len <= 0) len <- 1
    frac  <- (pos - min(pos))/len
    theta <- alpha_start[i] + frac*(alpha_end[i] - alpha_start[i])
    thp   <- sign_ang*theta
    rr    <- r0 + (lp/scale_max)*(h2 - 0.05)

    xp <- cx + sin(thp)*rr
    yp <- cy + cos(thp)*rr
    ic <- cut(lp, 100, labels = FALSE)
    colp <- cols_all[ic]
    ok <- which(!is.na(colp) & is.finite(xp) & is.finite(yp))
    if (length(ok))
      points(xp[ok], yp[ok], pch = 16, col = colp[ok], cex = 0.5)

    ## block ID at centre of each sector
    if (!is.null(block_ids) && length(block_ids) >= i){
      text(cx + sin(ang_c)*h1*0.35,
           cy + cos(ang_c)*h1*0.35,
           labels = block_ids[i],
           cex = 1.2, font = 2)
    }
  }

  ## ---- Mb labels on block edges ----
  for (i in 1:3){
    boundary_label(i, alpha_start[i], blk_start[i]/1e6, rad_offset = 0.08)
    boundary_label(i, alpha_end[i],   blk_end[i]/1e6,   rad_offset = 0.20)
  }

  ## ---- purple 5e-8 circle ----
  thr   <- min(-log10(5e-8), scale_max)
  r_thr <- r0 + (thr/scale_max)*(h2 - 0.05)
  tt    <- seq(0, 2*pi, length.out = 361)
  lines(sin(tt)*r_thr, cos(tt)*r_thr, col = "purple", lwd = 1.4)

  ## ---- lightblue outer circle (max -log10(P)) ----
  Rmax <- r0 + h2
  lines(sin(tt)*Rmax, cos(tt)*Rmax, col = "lightblue", lwd = 1.2)

## ---- gene arcs (keep only 1 no-overlap curve per bin) ----
if (!is.null(gene_df) && nrow(gene_df) > 0) {
  
  region_start <- min(blk_start)
  region_end   <- max(blk_end)

  # filter genes that at least partially enter this block-triplet
  genes <- gene_df[START <= region_end & END >= region_start]
  if (nrow(genes) > 0) {

    midpos <- pmax(genes$START, region_start)
    midpos <- (midpos + pmin(genes$END, region_end)) / 2

    blk_idx   <- integer(nrow(genes))
    theta_mid <- rep(NA_real_, nrow(genes))

    # 计算角度
    for (b in 1:3) {
      in_b <- which(midpos >= blk_start[b] & midpos <= blk_end[b])
      if (!length(in_b)) next
      fracm <- (midpos[in_b] - blk_start[b])/(blk_end[b] - blk_start[b])
      theta_mid[in_b] <- alpha_start[b] + fracm*(alpha_end[b] - alpha_start[b])
      blk_idx[in_b]   <- b
    }

    keep <- which(is.finite(theta_mid) & blk_idx > 0)
    if (length(keep) > 0) {
      genes     <- genes[keep]
      theta_mid <- theta_mid[keep]
      blk_idx   <- blk_idx[keep]
      midpos    <- midpos[keep]

      ##---------------------------------------
      ## COLLAPSE OVERLAPPING GENES → one per angle-bin
      ##---------------------------------------
      bin_width <- 2*pi/240
      bins <- floor((theta_mid + pi) / bin_width)

      gene_len <- genes$END - genes$START
      dt_bins  <- data.table(idx = seq_along(bins), bin = bins, len = gene_len)
      dt_sel   <- dt_bins[, .SD[which.max(len)], by = bin]
      sel      <- sort(dt_sel$idx)

      if (length(sel) > 0) {
        genes     <- genes[sel]
        theta_mid <- theta_mid[sel]
        blk_idx   <- blk_idx[sel]
        midpos    <- midpos[sel]

        r_gene <- (r0 + h2) - 0.25   # place below outer circle

        ## draw each gene arc
        for (ii in seq_along(theta_mid)) {
          gs <- max(genes$START[ii], region_start)
          ge <- min(genes$END[ii],   region_end)
          if (ge <= gs) next

          for (b in 1:3) {
            s <- max(gs, blk_start[b])
            e <- min(ge, blk_end[b])
            if (e <= s) next

            frac1 <- (s - blk_start[b])/(blk_end[b] - blk_start[b])
            frac2 <- (e - blk_start[b])/(blk_end[b] - blk_start[b])
            th1   <- alpha_start[b] + frac1*(alpha_end[b] - alpha_start[b])
            th2   <- alpha_start[b] + frac2*(alpha_end[b] - alpha_start[b])

            tseg <- seq(th1, th2, length.out = 40)
            xs <- sin(sign_ang*tseg)*r_gene
            ys <- cos(sign_ang*tseg)*r_gene
            lines(xs, ys, col = "grey25", lwd = 1.0)
          }

          thm <- theta_mid[ii]
          xm  <- sin(sign_ang*thm)*(r_gene - 0.05)
          ym  <- cos(sign_ang*thm)*(r_gene - 0.05)
          text(xm, ym, labels = genes$GENE[ii],
               cex = 0.7, col = "black", srt = (-sign_ang * thm * 180/pi) - 90)
        }
      }
    }
  }
}


  ## ---- LD colour legend bar on top ----
  usr <- par("usr"); xr <- usr[1:2]; yr <- usr[3:4]
  w <- diff(xr); h <- diff(yr)
  x0 <- xr[1] + 0.30*w
  x1 <- xr[1] + 0.70*w
  y1 <- yr[2] - 0.02*h
  y0 <- y1 - 0.04*h
  rect(x0, y0, x1, y1, border = "grey40", lwd = 0.6)
  nstrip <- 180
  xx <- seq(x0, x1, length.out = nstrip+1)
  cc <- colorRampPalette(col)(nstrip)
  for (ii in seq_len(nstrip)){
    rect(xx[ii], y0, xx[ii+1], y1, col = cc[ii], border = NA)
  }
}

## ---------------------------------------------------------
## 7. block_plot & blockzoom
## ---------------------------------------------------------
block_plot <- function(gwas, block_id, LD, gene_df = NULL){

  if (is.null(LD$gwas) || is.null(LD$ld_list) || is.null(LD$info))
    stop("LD 对象结构不完整")

  gdf <- as.data.frame(LD$gwas)

  all_ids   <- LD$info$triplet[[1]]
  center_id <- as.integer(block_id)
  left_id   <- center_id - 1L
  right_id  <- center_id + 1L

  target_ids <- c(center_id, left_id, right_id)    # 顺时钟：中、左、右
  idx_ord    <- match(target_ids, all_ids)
  if (any(is.na(idx_ord)))
    stop("LD$info$triplet 中缺少中心或相邻块：", paste(target_ids, collapse = ","))

  ld_ord   <- LD$ld_list[idx_ord]
  blk_ids  <- target_ids
  chr_now  <- LD$info$chr[1]

  gene_sub <- NULL
  if (!is.null(gene_df))
    gene_sub <- gene_df[CHR == chr_now]

  myplot(gdf, ld_ord, h1 = 2, r = 0.06, gap = 0.2, block_ids = blk_ids, gene_df = gene_sub)
  invisible(NULL)
}


blockzoom <- function(gwas_file, block_id){

  if (missing(gwas_file) || missing(block_id)) stop("需要指定 gwas_file 和 block_id")

  gwas_dt <- gwas_raw(gwas_file)
  LD      <- get_LD(block_id = block_id, gwas_dt = gwas_dt)
  gene_df <- .read_glist(glist)

  block_plot(gwas_dt, block_id = block_id, LD, gene_df = gene_df)

  invisible(list(gwas = gwas_dt, LD = LD, genes = gene_df))
}
