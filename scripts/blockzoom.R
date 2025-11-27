
get_LD <- function("ukbb", "race", "grch", block) {

}

block_plot <- function (gwas, block, LD)  { 
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



block_id = 37
grch = 38
race = "EUR"
hapref = "ukbb"

gwas <- read.table("bmi.gz")
LD <- get_LD(hapref, race, grch, block_id)
block_plot(gwas, block, LD)