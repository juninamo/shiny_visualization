# =============================================================================
# Seurat Viewer Shiny App
# Seurat objectのRDSファイルを読み込み、遺伝子発現を可視化する
# RDS fileがある場所にこのapp.Rをおいて、R内で
#> library(shiny);
#> runApp(".")
# と動かしてください。
# =============================================================================

library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(DT)

# --- アプリのディレクトリ内のRDSファイルを検出 ---
app_dir <- dirname(sys.frame(1)$ofile %||% ".")
if (app_dir == ".") app_dir <- getwd()
rds_files <- list.files(app_dir, pattern = "\\.rds$", ignore.case = TRUE)

# --- 外部データベースURL生成ヘルパー ---
ncbi_gene_url <- function(gene_name) {
  paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", gene_name)
}
immunexut_url <- function(gene_name) {
  paste0("https://www.immunexut.org/eqtlGenes?gene_symbol=", gene_name)
}
ampra2_url <- function(gene_name) {
  paste0("https://immunogenomics.io/ampra2/app/?ds=fibroblast&gene=",gene_name,"&groupby=none")
}

# --- 組成プロット用のカラーパレット ---
manual_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFD92F", "#A65628", "#F781BF",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#E5C494", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A",
  "#33A02C", "#FB9A99", "#FDBF6F", "#CAB2D6", "#6A3D9A", "#B15928", "#8DD3C7", "#BEBADA",
  "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F"
)

# --- 系統 (lineage) ごとのグラデーション配色 ---
# クラスター名末尾の系統名 (例: "0_B_Plasma" -> "B_Plasma") を取り出し、
# 系統内で数字プレフィックス順に明→暗のグラデーションを割り当てる。
lineage_color_ramp <- function(category, n) {
  if (n == 0) return(character(0))
  ramps <- list(
    Epithelial          = c("#FFEBED99", "#B71B1B99"),
    Stromal_Endothelial = c("#E0E0E099", "#42424299"),
    TNK_ILC             = c("#ACD9AE99", "#358A3999"),
    B_Plasma            = c("#C8B9E599", "#311A9299"),
    Myeloid             = c("#FFECB399", "#E6510099")
  )
  if (!is.null(ramps[[category]])) {
    return(grDevices::colorRampPalette(ramps[[category]], alpha = TRUE)(n))
  }
  rep("grey60", n)
}

# 末尾の系統名を取り出す（"0_B_Plasma" -> "B_Plasma"）。
# regmatches はマッチしない要素を取り除いて長さが変わるため、ここでは
# マッチしない要素は "" にして必ず入力と同じ長さのベクトルを返す。
extract_lineage <- function(x) {
  x <- as.character(x)
  m <- regexpr("[A-Za-z_]+$", x)
  out <- rep("", length(x))
  hit <- m != -1
  if (any(hit)) out[hit] <- regmatches(x, m)
  sub("^_+", "", out)
}

generate_cluster_colors <- function(levels) {
  raw <- as.character(levels)
  clean <- raw[raw != "NA" & !is.na(raw)]
  if (length(clean) == 0) {
    out <- setNames(character(0), character(0))
  } else {
    num_prefix <- suppressWarnings(as.numeric(sub("^(\\d+).*", "\\1", clean)))
    # 末尾の系統名（長さを保つ）
    coarse <- extract_lineage(clean)
    df <- data.frame(cluster = clean, num_prefix = num_prefix,
                     coarse = coarse, color = NA_character_,
                     stringsAsFactors = FALSE)
    df <- df[order(df$coarse, df$num_prefix), ]
    for (ct in unique(df$coarse)) {
      idx <- df$coarse == ct
      df$color[idx] <- lineage_color_ramp(ct, sum(idx))
    }
    out <- setNames(df$color, df$cluster)
  }
  # NA 用の色を追加
  if (any(raw == "NA" | is.na(raw))) {
    out <- c(out, setNames("grey80", "NA"))
  }
  out
}

# 既知系統リスト
known_lineages <- c("Epithelial", "Stromal_Endothelial", "TNK_ILC",
                    "B_Plasma", "Myeloid")

# クラスター変数が「数字_系統」形式で既知系統に該当するか判定する。
# 半分以上のレベルが既知系統に一致すれば名前付きカラーベクトルを、
# そうでなければ NULL を返す（その場合はデフォルト配色にフォールバック）。
lineage_colors_or_null <- function(levels) {
  levs <- as.character(unique(levels))
  levs <- levs[!is.na(levs) & levs != "NA"]
  if (length(levs) == 0) return(NULL)
  coarse <- extract_lineage(levs)
  if (length(coarse) == 0 || mean(coarse %in% known_lineages) < 0.5) return(NULL)
  generate_cluster_colors(levels)
}

# クラスターレベルの並び順を決める。
# 「数字_系統」形式（既知系統が半数以上）の場合は、まず系統ごと
# （known_lineages の順 → 未知系統はアルファベット順）に、その中で先頭の
# 数字の昇順に並べる。形式が合わなければ自然順（数字を考慮したソート）。
cluster_level_order <- function(levels) {
  levs <- as.character(unique(levels))
  na_lev <- levs[is.na(levs) | levs == "NA"]
  levs <- levs[!is.na(levs) & levs != "NA"]
  if (length(levs) == 0) return(c(levs, na_lev))
  coarse <- extract_lineage(levs)
  num_prefix <- suppressWarnings(as.numeric(sub("^(\\d+).*", "\\1", levs)))
  if (length(coarse) > 0 && mean(coarse %in% known_lineages) >= 0.5) {
    # 系統順（既知系統を優先）→ 系統名 → 先頭数字
    lin_rank <- match(coarse, known_lineages)
    lin_rank[is.na(lin_rank)] <- length(known_lineages) + 1L
    ord <- order(lin_rank, coarse, num_prefix, levs)
  } else if (all(!is.na(num_prefix))) {
    ord <- order(num_prefix, levs)   # 数字主体のクラスター名は数値順
  } else {
    ord <- order(levs)
  }
  c(levs[ord], na_lev)
}

# =============================================================================
# よく使うマーカー遺伝子セット（Heatmap / Dot plot 用に事前登録）
# =============================================================================
# 細胞種ごとのマーカー（グループ付き: Dot plot のファセットに使用）
celltype_markers <- list(
  "Adipocyte"                  = c("ADIPOQ", "PLIN1", "PLIN4"),
  "B cell"                     = c("CD19", "CD79A", "MS4A1"),
  "Dendritic cell"             = c("WDFY4", "CLEC9A", "CLEC10A", "CD1C"),
  "Endothelial"                = c("CDH5", "FLT1", "PECAM1", "VWF"),
  "Lining fibroblast"          = c("CD55", "PRG4"),
  "Lymphatic endothelial cell" = c("CCL21", "EGFL7", "FLT4"),
  "Macrophage"                 = c("C1QA", "CD68", "CTSB", "MARCO"),
  "Mast cell"                  = c("HDC", "KIT"),
  "Mural cell"                 = c("NOTCH3"),
  "Neutrophil"                 = c("MPO", "TMC8"),
  "pDC"                        = c("IRF7", "TCF4"),
  "Plasma cell"                = c("MZB1", "PRDM1", "XBP1"),
  "Proliferating cell"         = c("MKI67", "TOP2A"),
  "Sublining fibroblast"       = c("COL6A1", "CXCL12", "PDPN", "THY1"),
  "T cell"                     = c("CD3E", "CD4", "CD8A", "IL7R", "NCAM1", "NKG7"),
  "Epithelial"                 = c("EPCAM", "LYZ", "PIP", "MUC7", "MUC5B", "CEACAM6",
                                   "BPIFB2", "TFCP2L1", "KRT19", "KRT5", "KRT14",
                                   "TP63", "CAV1", "ACTA2")
)

# 系統別の詳細マーカー（フラットなベクトル。重複は unique() で除去）
t_markers <- unique(c(
  "SELL","CCR7","TCF7","CD69","IL7R","GIMAP5","LEF1","CD3E","CD4","CD8A",
  "CXCR3","CXCR4","CXCR5","CXCR6","CX3CR1","CCR2","CCR4","CCR6","CCL5","CCL4",
  "CXCL13","FOXP3","IL2RA","PDCD1","CTLA4","ICOS","CD40LG","TIGIT","LAG3","MAF",
  "RORC","GATA3","TBX21","HLA-DRA","HLA-DRB1","B3GAT1","CD27","CD38","ITGAE",
  "KLRB1","SLC4A10","IFI44L","MX1","IFIT3","ZNF683","XCL1","GZMB","GZMK","GZMA",
  "GZMH","GNLY","PRF1","NKG7","IFNG","TRDV1","TRDV2","TRGV9","AQP3","MKI67",
  "MALAT1","CCR5","IL21","NR3C1","PCNA","ZBTB16","TRDC"
))

b_markers <- unique(c(
  "CD19","MS4A1","CD79A","CR2","FCER2","TCL1A","CD27","IGHD","IGHG1","IGHA1",
  "SDC1","RPN2","XBP1","PRDX4","MZB1","CD38","ITGAX","TBX21","ZEB2","EMP3",
  "FGR","HSPA1B","TNFRSF13B","S100A10","S100A4","NR4A1","NFKBIA","NFKBIZ","BCL6",
  "SLAMF7","CD24","MX1","AICDA","CD86","CD1C","CD22","CD5","BCL2","CXCR5","FAS",
  "TLR9","FCRL4","FCRL5","CD83","CD79B","CD40","MME","PAX5","ISG15","MKI67",
  "PRDM1","HLA-DRB1","ITGAM","IGHG3"
))

m_markers <- unique(c(
  "CD14","FCGR3A","CD163","CFLAR","LYZ","CLEC1B","S100A9","S100A8","ACTB","FTL",
  "FCGR3B","MERTK","FOLR2","CXCL10","STMN1","MKI67","WDFY4","CLEC10A","LAMP3",
  "CD1C","CLEC9A","THBD","FBP1","MRC1","LPL","IL1RN","CCR2","CD300E","VEGFA",
  "CXCL9","TREM2","KIT","GATA2","CSF1R","PDGFA","PDGFB","PDGFC","PDGFD","AREG",
  "IL1B","C3AR1","CXCR6","LYVE1","SELENOP","HBEGF","SPP1","APOE","C1QA","C1QB",
  "C1QC","STAT1","FCN1","PLCG2","FCER1A","AXL","MIR155","KLF4","ATF3","IFNAR1",
  "ALDOA","FABP5","CD36","NR4A3","IL6","CSF1","FZD1","LRP5","LRP6","NOTCH1",
  "HES1","CD68"
))

n_markers <- unique(c(
  "NCAM1","FCER1G","KLRF1","NKG7","GNLY","PRF1","CX3CR1","GZMB","GZMK","IFNG",
  "PCNA","MKI67","IL7R","RORC","KLRB1","GATA3"
))

tissue_markers <- unique(c(
  "THY1","FAP","PRG4","CD34","DKK3","NOTCH3","NOTCH4","CXCL12","POSTN","COL4A1",
  "KDR","PLVAP","CLDN5","EFNB2","NRP2","SELE","LIFR","PDPN","MYOT","MYH2","CCL21",
  "FLT4","PROX1","CLU","ACTA2","CDH5","MCAM","ICAM1","VCAM1","VWF","RUNX1",
  "PDGFRA","PDGFRB","EGFR","CCL2","IL6R","TGFBR1","TGFBR2","TGFBR3","CSF1","CSF2",
  "CSF3","IL1R1","IL1R2","CXCL16","C3","IL6","CLIC5","PI16","DPP4","CD74","HLA-DRA",
  "SFRP1","RSPO3","SPARC","COL1A1","COL1A2","COL3A1","LTBP4","DLL4","LYVE1","CTNNB1",
  "DLL1","LGR5","AXIN2","WNT1","WNT2","WNT2B","WNT3","WNT4","TCF4","SOX9"
))

# マーカーセットの登録（UI で選択するパネル）
# - グループ付きリスト: Dot plot のファセット分けに利用
# - フラットベクトル: 単一グループとして扱う
marker_sets <- list(
  "Cell type (grouped)" = celltype_markers,
  "T cell"              = t_markers,
  "B cell"              = b_markers,
  "Myeloid"             = m_markers,
  "NK / ILC"            = n_markers,
  "Stromal / Tissue"    = tissue_markers,
  "All curated"         = unique(c(unlist(celltype_markers, use.names = FALSE),
                                   t_markers, b_markers, m_markers,
                                   n_markers, tissue_markers))
)

# マーカーセットを (feature, group) のデータフレームに展開
marker_set_to_df <- function(set) {
  if (is.list(set)) {
    do.call(rbind, lapply(names(set), function(g) {
      data.frame(feature = set[[g]], group = g, stringsAsFactors = FALSE)
    }))
  } else {
    data.frame(feature = set, group = "markers", stringsAsFactors = FALSE)
  }
}

# Seurat v4/v5 両対応で発現マトリクス (genes x cells) を取得
# Seurat v5 で data レイヤーが複数（未結合の multi-layer assay）の場合、
# GetAssayData(layer="data") はエラーになる。その場合は各 data レイヤーを
# LayerData で取り出して横結合（cbind）する。返す行列の列は obj の細胞順に揃える。
get_expr_matrix <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)
  assay <- SeuratObject::DefaultAssay(obj)
  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "data"),
                    error = function(e) NULL)
  }
  if (is.null(mat)) {
    # v5 複数レイヤー: data.* （無ければ counts.*）レイヤーを結合
    layers <- tryCatch(SeuratObject::Layers(obj, assay = assay),
                       error = function(e) character(0))
    dl <- layers[grepl("^data", layers)]
    if (length(dl) == 0) dl <- layers[grepl("^counts", layers)]
    if (length(dl) == 0) return(NULL)
    parts <- lapply(dl, function(ly) {
      m <- SeuratObject::LayerData(obj, assay = assay, layer = ly)
      m[intersect(genes, rownames(m)), , drop = FALSE]
    })
    mat <- do.call(cbind, parts)
  }
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) return(NULL)
  mat <- mat[genes, , drop = FALSE]
  # 列をオブジェクトの細胞順に揃える（cbind で順序が変わる場合の対策）
  cells <- intersect(colnames(obj), colnames(mat))
  mat[, cells, drop = FALSE]
}

# --- 並列処理ヘルパー ---------------------------------------------------------
# 利用可能なワーカー数（1コア残し、最大8）。環境変数 SHINY_VIZ_THREADS で上書き可。
n_workers <- function() {
  ov <- suppressWarnings(as.integer(Sys.getenv("SHINY_VIZ_THREADS", "")))
  if (!is.na(ov) && ov >= 1) return(ov)
  n <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (is.na(n) || n < 1) n <- 1L
  as.integer(max(1L, min(n - 1L, 8L)))
}
# 各要素を seeds[i] で個別にシードして並列実行（コア数に依らず再現性を担保）。
# fork が使える環境(mac/linux)は parallel::mclapply、それ以外は逐次 lapply。
par_lapply_seeded <- function(X, FUN, seeds) {
  ff <- function(i) { set.seed(seeds[[match(i, X)]]); FUN(i) }
  nw <- n_workers()
  if (nw > 1L && .Platform$OS.type != "windows" &&
      requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(X, ff, mc.cores = nw)
  } else {
    lapply(X, ff)
  }
}

# Seurat オブジェクトから空間座標 (cell, x, y) を取り出す。
# meta.data の x/y 系カラムを優先し、無ければ @images の TissueCoordinates。
spatial_coords <- function(obj) {
  meta <- obj@meta.data
  cand <- list(c("x", "y"), c("x_centroid", "y_centroid"), c("sdimx", "sdimy"),
               c("imagecol", "imagerow"), c("X", "Y"), c("spatial_1", "spatial_2"))
  for (p in cand) {
    if (all(p %in% names(meta)) && is.numeric(meta[[p[1]]]) && is.numeric(meta[[p[2]]])) {
      return(data.frame(cell = rownames(meta), x = meta[[p[1]]], y = meta[[p[2]]],
                        stringsAsFactors = FALSE))
    }
  }
  if (length(obj@images) > 0) {
    xy <- tryCatch(do.call(rbind, lapply(names(obj@images), function(im) {
      tc <- SeuratObject::GetTissueCoordinates(obj@images[[im]])
      if (!("cell" %in% names(tc))) tc$cell <- rownames(tc)
      xc <- if ("x" %in% names(tc)) "x" else names(tc)[1]
      yc <- if ("y" %in% names(tc)) "y" else names(tc)[2]
      data.frame(cell = tc$cell, x = tc[[xc]], y = tc[[yc]], stringsAsFactors = FALSE)
    })), error = function(e) NULL)
    if (!is.null(xy)) return(xy)
  }
  NULL
}

# @images にセグメンテーション(細胞境界ポリゴン)があれば取り出す。
# 戻り値: data.frame(cell, x, y) の頂点列（細胞ごとに複数行）。無ければ NULL。
spatial_segmentation <- function(obj) {
  if (length(obj@images) == 0) return(NULL)
  segs <- tryCatch(do.call(rbind, lapply(names(obj@images), function(im) {
    img <- obj@images[[im]]
    bnames <- tryCatch(SeuratObject::Boundaries(img), error = function(e) character(0))
    if (!("segmentation" %in% bnames)) return(NULL)
    tc <- SeuratObject::GetTissueCoordinates(img, which = "segmentation")
    if (is.null(tc) || !all(c("x", "y") %in% names(tc))) return(NULL)
    if (!("cell" %in% names(tc))) tc$cell <- rownames(tc)
    data.frame(cell = as.character(tc$cell), x = tc$x, y = tc$y, stringsAsFactors = FALSE)
  })), error = function(e) NULL)
  if (is.null(segs) || nrow(segs) == 0) return(NULL)
  segs
}

# =============================================================================
# 翻訳辞書
# =============================================================================
i18n <- list(
  ja = list(
    # サイドバー
    data_load       = "\U0001F4C2 データ読み込み",
    rds_file        = "RDSファイル (複数選択可)",
    load_btn        = "読み込む",
    active_ds       = "アクティブ",
    ref_ds          = "リファレンス (比較用)",
    ref_none        = "（なし）",
    ref_missing     = "リファレンスに存在しません",
    vis_settings    = "\U0001F52C 可視化設定",
    gene_name       = "遺伝子名",
    gene_placeholder = "ファイルを読み込んでください",
    group_var       = "グループ変数",
    ref_group_var   = "グループ変数 (リファレンス)",
    umap_reduction  = "UMAP (reduction)",
    plot_settings   = "\U0001F3A8 プロット設定",
    pt_size         = "点のサイズ",
    umap_label_size = "UMAPラベルの大きさ",
    plot_height     = "プロットの高さ (px)",
    plot_width      = "プロットの幅 (px)",
    ref_plot_settings = "\U0001F3A8 プロット設定 (リファレンス)",
    ref_pt_size     = "点のサイズ (リファレンス)",
    ref_plot_height = "プロットの高さ (リファレンス)",
    ref_plot_width  = "プロットの幅 (リファレンス)",

    # DEG
    deg_settings      = "DEG解析設定",
    deg_group_var     = "グループ変数",
    deg_cat_prefix    = "[カテゴリ] ",
    deg_num_prefix    = "[数値] ",
    deg_percentile    = "Top / Bottom パーセンタイル (%)",
    deg_group1        = "Group 1 (テスト群)",
    deg_group2        = "Group 2 (コントロール群・複数選択可)",
    deg_all_others    = "その他全て (All others)",
    deg_feature_mode  = "解析する遺伝子",
    deg_feat_all      = "全遺伝子",
    deg_feat_set      = "マーカーセット",
    deg_feat_set_custom = "セット + 自分で選択",
    deg_feat_custom   = "自分で選択のみ",
    deg_marker_set    = "マーカーセット",
    deg_custom_genes  = "遺伝子を選択",
    deg_logfc         = "logFC閾値",
    deg_pval          = "p値閾値",
    deg_run           = "DEG解析を実行",
    deg_total_cells   = "全細胞数: %s",
    deg_top_cells     = "Top %d%%: %s 細胞",
    deg_bottom_cells  = "Bottom %d%%: %s 細胞",

    # Composition
    comp_settings        = "組成プロット設定",
    comp_cluster         = "クラスター変数 (塗り分け)",
    comp_x               = "X軸変数",
    comp_x_ref           = "X軸変数 (リファレンス)",
    comp_facets          = "ファセット変数 (任意・複数可)",
    comp_color_scheme    = "配色",
    comp_scheme_lineage  = "系統グラデーション (動的)",
    comp_scheme_manual   = "パレット",
    comp_title           = "クラスター組成",
    comp_yaxis           = "割合",
    comp_no_cat          = "カテゴリ変数が見つかりません。",

    # Heatmap / Dot plot 共通
    ref_settings_label   = "リファレンス設定",
    ref_marker_cluster   = "クラスター変数 (リファレンス)",
    ref_draw_clusters    = "描画するクラスター (リファレンス)",
    marker_settings      = "マーカー設定",
    marker_set           = "マーカーセット",
    marker_cluster       = "クラスター変数",
    draw_clusters        = "描画するクラスター",
    coarse_var           = "絞り込み列 (大分類)",
    coarse_values        = "絞り込む値",
    mk_feature_mode      = "遺伝子の選択",
    mk_feat_set          = "マーカーセット",
    mk_feat_set_custom   = "セット + 自分で選択",
    mk_feat_custom       = "自分で選択のみ",
    mk_custom_genes      = "遺伝子を選択",
    mk_custom_group      = "カスタム",
    mk_set_genes         = "セット内の遺伝子（不要なものは削除可）",
    plot_run             = "描画",
    plot_run_hint        = "設定後「描画」ボタンを押してください",
    hm_scale             = "Z-score (行スケーリング)",
    hm_cluster_rows      = "行をクラスタリング",
    hm_cluster_cols      = "列をクラスタリング",
    hm_combine           = "リファレンスと統合 (似たクラスターが近くに並ぶ)",
    dot_combine          = "リファレンスと統合 (行・列をクラスタリング)",
    corr_title           = "発現パターンの近さ: 各アクティブクラスター → 最も近いリファレンス",
    corr_download        = "対応表をCSVでダウンロード",
    corr_active          = "アクティブクラスター",
    corr_best            = "最も近いリファレンス",
    corr_cor             = "相関",
    corr_second          = "次点リファレンス",
    corr_margin          = "マージン(1位-2位)",
    corr_proposed        = "提案クラスター名",
    corr_conf            = "確信度",
    corr_help            = "各アクティブクラスターの平均発現プロファイル（共通遺伝子をデータセット内でZ-score化）を、各リファレンスクラスターとピアソン相関で比較し、相関が最も高いものを「最も近いリファレンス」としています。マージン = 1位の相関 − 2位の相関（大きいほど対応が明確）。提案クラスター名 = 最も近いリファレンス名（1位と2位が僅差(マージン<0.05)の場合は両方を併記し「?」付き）。確信度は相関の高さとマージンから high/medium/low。",
    subc_settings        = "サブクラスタリング設定",
    subc_init            = "現在のアクティブクラスターから初期化 / リセット",
    subc_target          = "細分化するクラスター (複数選択可)",
    subc_k               = "サブクラスター数 (k)",
    subc_run             = "細分化を実行",
    subc_need_init       = "先に「初期化」ボタンを押してください。",
    subc_no_ref          = "比較するリファレンスを選択してください。",
    subc_running         = "サブクラスタリング中...",
    subc_done            = "✅ 細分化完了: %d 個のサブクラスター",
    subc_placeholder     = "クラスターを選び「細分化を実行」を押してください（相関の低い/曖昧なクラスターを選ぶと有効です）",
    subc_leaf            = "現在の葉クラスター数: %d",
    subc_dl_anno         = "全クラスター注釈(細胞→名前)をCSVで保存",
    subc_help            = "選んだ各クラスターの細胞だけで、発現量の高い遺伝子を使ってサブクラスタリング（kmeans, k個）し、元のクラスター名の入れ子（例 5-TNK_ILC.0, 5-TNK_ILC.1）として命名します。生成したサブクラスターだけを同じリファレンスの遺伝子・クラスターと統合して Heatmap / Dot plot と相関解析を再描画します。さらに細分化したいクラスターを選んで何度でも繰り返せます。",
    hm_no_genes          = "選択したマーカーがデータ内に見つかりません。",
    hm_pkg_missing       = "pheatmap パッケージが必要です: install.packages('pheatmap')",
    dot_scale            = "ドットサイズ",
    dot_facet            = "グループでファセット",
    dot_pct             = "発現割合 (%)",
    dot_avg             = "平均発現 (scaled)",

    # 通知
    notify_loading      = "データを読み込み中...",
    notify_finalizing   = "セレクタを更新中...",
    notify_not_seurat   = "選択されたファイルはSeuratオブジェクトではありません。",
    notify_no_umap      = "注意: UMAPが計算されていません。UMAP表示は利用できません。",
    notify_load_done    = "✅ 読み込み完了: %s 細胞 × %s 遺伝子",
    notify_error        = "エラー: ",
    notify_deg_running  = "DEG解析を実行中...",
    notify_deg_same     = "Group 1とGroup 2は異なるグループを選択してください。",
    notify_deg_no_features = "解析対象の遺伝子が選択されていません。",
    notify_deg_done     = "✅ DEG解析完了: %d Up, %d Down",

    # Enrichr
    enrichr_up        = "Up-DEG → Enrichr",
    enrichr_down      = "Down-DEG → Enrichr",
    enrichr_none      = "有意な遺伝子がありません。",
    enrichr_opening   = "Enrichrにリストを送信中...",
    enrichr_open      = "Enrichrで結果を開く (%s, %d遺伝子)",
    enrichr_err       = "Enrichrへの送信に失敗しました: ",

    # GSEA
    gsea_settings     = "GSEA解析設定 (fgsea × MSigDB)",
    gsea_group_var    = "グループ変数",
    gsea_group1       = "対象グループ (このグループで上昇する順にランキング)",
    gsea_species      = "種",
    gsea_collection   = "MSigDB コレクション",
    gsea_minsize      = "最小サイズ",
    gsea_maxsize      = "最大サイズ",
    gsea_run          = "GSEAを実行",
    gsea_running      = "GSEA実行中... (遺伝子ランキング → fgsea)",
    gsea_done         = "✅ GSEA完了: %d パスウェイ検定, 有意(padj<0.05) %d",
    gsea_no_genes     = "ランキングに使える遺伝子がありません。",
    gsea_no_sets      = "重複する遺伝子セットがありません（パネルが小さい可能性）。",
    gsea_placeholder  = "設定して「GSEAを実行」ボタンを押してください",
    gsea_nes_title    = "上位パスウェイ (NES)",
    gsea_need_deg     = "先に「DEG」タブでDEG解析を実行してください。",
    gsea_uses_deg     = "GSEA は上の比較設定で実行した DEG 結果（avg_log2FC でランキング）を使用します。",

    # Spatial
    spatial_settings  = "空間プロット設定",
    spatial_color     = "色分け変数",
    spatial_sample_col = "サンプル列",
    spatial_sample    = "サンプル",
    spatial_pt        = "点のサイズ",
    spatial_flip      = "Y軸を反転",
    spatial_none      = "このオブジェクトに位置座標(x/y)が見つかりません。",
    spatial_ds        = "空間データ",
    spatial_seg       = "セグメンテーション(細胞境界)を表示",
    spatial_by        = "色分けの基準",
    spatial_by_meta   = "メタデータ変数",
    spatial_by_gene   = "遺伝子発現",
    spatial_gene      = "遺伝子",
    spatial_highlight = "強調・近傍解析の対象クラスター (複数選択可)",
    spatial_hl_alpha  = "強調: 不透明度",
    spatial_hl_size   = "強調: 点サイズ",
    spatial_other_alpha = "その他: 不透明度",
    spatial_other_size  = "その他: 点サイズ",
    spatial_map       = "マップ",
    spatial_nbr       = "近傍距離",
    spatial_nbr_need  = "メタデータ(カテゴリ)で色分けし、対象クラスターを選択してください。",
    spatial_nbr_x     = "距離 (最近接、log scale)",
    spatial_nbr_y     = "対象細胞の累積割合",
    spatial_nbr_title = "%s から各クラスターへの最近接距離 (近いほど左)",
    spatial_nbr_run   = "近傍解析を実行",
    spatial_nbr_help  = "対象クラスターの各細胞から、他クラスターの最も近い細胞までの距離を計算し、その累積分布(ECDF)を描きます。曲線が左にあるほどそのクラスターは対象クラスターの近くに位置します(y=0.5が距離の中央値)。同一サンプル内でのみ距離を計算します。複数サンプル選択時はサンプルごとに計算し、平均±標準偏差を信頼帯(リボン)として表示します。",
    spatial_co        = "Co-occurrence",
    spatial_co_run    = "Co-occurrence を計算",
    spatial_co_x      = "距離",
    spatial_co_y      = "Co-occurrence 比 (P(j|i,d)/P(j))",
    spatial_co_title  = "%s 周辺における各クラスターの共局在 (>1:濃縮, <1:希薄)",
    spatial_co_help   = "squidpy の co_occurrence を参考にした指標です。対象クラスター i の細胞から距離 d の範囲(リング)内に他クラスター j の細胞が見つかる条件付き確率 P(j | i, d) を、全体での j の割合 P(j) で割った比を距離ごとに描きます。1より大きいと i の周囲に j が濃縮、1未満だと希薄を意味します。サンプルごとに計算し平均±標準偏差をリボン表示します。",
    spatial_ne        = "Neighbors enrichment",
    spatial_ne_run    = "Neighbors enrichment を計算",
    spatial_ne_k      = "近傍数 k",
    spatial_ne_perm   = "並べ替え回数",
    spatial_ne_stars  = "有意マーク(*)を表示",
    spatial_ne_prog   = "Neighbors enrichment を計算中",
    spatial_eta       = "%d/%d (残り約 %s)",
    spatial_ne_title  = "Neighbors enrichment z-score (正:隣接, 負:回避)",
    spatial_ne_help   = "squidpy の nhood_enrichment を参考にした指標です。各細胞の k 近傍で空間グラフを作り、クラスター間の隣接エッジ数を、ラベルをランダムに並べ替えた帰無分布と比較して z-score を計算します。正の値はそのクラスター対が予想より隣接、負の値は回避を意味します。z-score 行列をヒートマップで表示します(同一サンプル内のエッジのみ)。 有意なペアにアスタリスク(z→正規近似で両側p値→BH補正; *<0.05, **<0.01, ***<0.001)を表示します。行・列は階層クラスタリングで並べ替えます。",
    spatial_need_run  = "対象を選び「計算」ボタンを押してください。",

    # プレースホルダ
    placeholder_load  = "\U0001F4C2 RDSファイルを選択して「読み込む」ボタンを押してください",
    placeholder_deg   = "上のパネルでグループを設定し「DEG解析を実行」ボタンを押してください",
    no_umap_title     = "⚠️ このSeuratオブジェクトにはUMAP reductionが含まれていません",
    no_umap_msg       = "先にRunUMAP()を実行して保存し直してください。"
  ),

  en = list(
    # Sidebar
    data_load       = "\U0001F4C2 Load Data",
    rds_file        = "RDS File (multi-select)",
    load_btn        = "Load",
    active_ds       = "Active dataset",
    ref_ds          = "Reference (compare)",
    ref_none        = "(none)",
    ref_missing     = "Not present in reference",
    vis_settings    = "\U0001F52C Visualization",
    gene_name       = "Gene",
    gene_placeholder = "Load an RDS file first",
    group_var       = "Group Variable",
    ref_group_var   = "Group Variable (reference)",
    umap_reduction  = "UMAP (reduction)",
    plot_settings   = "\U0001F3A8 Plot Settings",
    pt_size         = "Point Size",
    umap_label_size = "UMAP Label Size",
    plot_height     = "Plot Height (px)",
    plot_width      = "Plot Width (px)",
    ref_plot_settings = "\U0001F3A8 Plot Settings (reference)",
    ref_pt_size     = "Point Size (reference)",
    ref_plot_height = "Plot Height (reference)",
    ref_plot_width  = "Plot Width (reference)",

    # DEG
    deg_settings      = "DEG Analysis Settings",
    deg_group_var     = "Group Variable",
    deg_cat_prefix    = "[Category] ",
    deg_num_prefix    = "[Numeric] ",
    deg_percentile    = "Top / Bottom Percentile (%)",
    deg_group1        = "Group 1 (Test)",
    deg_group2        = "Group 2 (Control, multi-select)",
    deg_all_others    = "All others",
    deg_feature_mode  = "Genes to test",
    deg_feat_all      = "All genes",
    deg_feat_set      = "Marker set",
    deg_feat_set_custom = "Set + custom",
    deg_feat_custom   = "Custom only",
    deg_marker_set    = "Marker Set",
    deg_custom_genes  = "Select genes",
    deg_logfc         = "logFC Threshold",
    deg_pval          = "p-value Threshold",
    deg_run           = "Run DEG Analysis",
    deg_total_cells   = "Total cells: %s",
    deg_top_cells     = "Top %d%%: %s cells",
    deg_bottom_cells  = "Bottom %d%%: %s cells",

    # Composition
    comp_settings        = "Composition Plot Settings",
    comp_cluster         = "Cluster Variable (fill)",
    comp_x               = "X-axis Variable",
    comp_x_ref           = "X-axis Variable (reference)",
    comp_facets          = "Facet Variables (optional)",
    comp_color_scheme    = "Color Scheme",
    comp_scheme_lineage  = "Lineage gradient (dynamic)",
    comp_scheme_manual   = "Manual palette",
    comp_title           = "Cluster Composition",
    comp_yaxis           = "Proportion",
    comp_no_cat          = "No categorical variables found.",

    # Heatmap / Dot plot shared
    ref_settings_label   = "Reference settings",
    ref_marker_cluster   = "Cluster Variable (reference)",
    ref_draw_clusters    = "Clusters to plot (reference)",
    marker_settings      = "Marker Settings",
    marker_set           = "Marker Set",
    marker_cluster       = "Cluster Variable",
    draw_clusters        = "Clusters to plot",
    coarse_var           = "Filter column (coarse)",
    coarse_values        = "Filter values",
    mk_feature_mode      = "Gene selection",
    mk_feat_set          = "Marker set",
    mk_feat_set_custom   = "Set + custom",
    mk_feat_custom       = "Custom only",
    mk_custom_genes      = "Select genes",
    mk_custom_group      = "Custom",
    mk_set_genes         = "Genes in set (remove any you don't need)",
    plot_run             = "Plot",
    plot_run_hint        = "Configure options, then click 'Plot'",
    hm_scale             = "Z-score (row scaling)",
    hm_cluster_rows      = "Cluster rows",
    hm_cluster_cols      = "Cluster columns",
    hm_combine           = "Combine with reference (similar clusters sit together)",
    dot_combine          = "Combine with reference (cluster rows & columns)",
    corr_title           = "Expression similarity: each active cluster → nearest reference",
    corr_download        = "Download correspondence table (CSV)",
    corr_active          = "Active cluster",
    corr_best            = "Nearest reference",
    corr_cor             = "Correlation",
    corr_second          = "2nd reference",
    corr_margin          = "Margin(1st-2nd)",
    corr_proposed        = "Proposed cluster name",
    corr_conf            = "Confidence",
    corr_help            = "Each active cluster's mean-expression profile (over shared genes, z-scored within its dataset) is compared to every reference cluster by Pearson correlation; the highest correlation is the 'nearest reference'. Margin = 1st − 2nd correlation (larger = more confident). Proposed cluster name = the nearest reference label (if 1st and 2nd are within margin <0.05, both are shown with a '?'). Confidence (high/medium/low) is from the correlation strength and margin.",
    subc_settings        = "Sub-clustering settings",
    subc_init            = "Initialize / reset from current active clusters",
    subc_target          = "Clusters to sub-cluster (multi-select)",
    subc_k               = "Number of sub-clusters (k)",
    subc_run             = "Run sub-clustering",
    subc_need_init       = "Click Initialize first.",
    subc_no_ref          = "Select a reference dataset to compare.",
    subc_running         = "Sub-clustering...",
    subc_done            = "✅ Done: %d sub-clusters",
    subc_placeholder     = "Select clusters and click Run sub-clustering (pick the low-correlation / ambiguous ones)",
    subc_leaf            = "current leaf clusters: %d",
    subc_dl_anno         = "Download full annotation (cell→label) CSV",
    subc_help            = "Re-clusters the cells of each selected cluster (kmeans into k, on the most variable genes) and names them nested under the original (e.g. 5-TNK_ILC.0, 5-TNK_ILC.1). The new sub-clusters are then combined with the same reference genes/clusters and re-plotted as Heatmap / Dot plot with the correlation analysis. Pick clusters and repeat as many times as you like.",
    hm_no_genes          = "None of the selected markers were found in the data.",
    hm_pkg_missing       = "The 'pheatmap' package is required: install.packages('pheatmap')",
    dot_scale            = "Dot size",
    dot_facet            = "Facet by group",
    dot_pct             = "Percent expressed (%)",
    dot_avg             = "Avg expression (scaled)",

    # Notifications
    notify_loading      = "Loading data...",
    notify_finalizing   = "Updating selectors...",
    notify_not_seurat   = "The selected file is not a Seurat object.",
    notify_no_umap      = "Note: No UMAP reduction found. UMAP plots are unavailable.",
    notify_load_done    = "✅ Loaded: %s cells × %s genes",
    notify_error        = "Error: ",
    notify_deg_running  = "Running DEG analysis...",
    notify_deg_same     = "Group 1 and Group 2 must be different.",
    notify_deg_no_features = "No genes are selected for analysis.",
    notify_deg_done     = "✅ DEG complete: %d Up, %d Down",

    # Enrichr
    enrichr_up        = "Up-DEGs → Enrichr",
    enrichr_down      = "Down-DEGs → Enrichr",
    enrichr_none      = "No significant genes.",
    enrichr_opening   = "Submitting list to Enrichr...",
    enrichr_open      = "Open results in Enrichr (%s, %d genes)",
    enrichr_err       = "Failed to submit to Enrichr: ",

    # GSEA
    gsea_settings     = "GSEA Settings (fgsea × MSigDB)",
    gsea_group_var    = "Group Variable",
    gsea_group1       = "Target group (rank by up-regulation in this group)",
    gsea_species      = "Species",
    gsea_collection   = "MSigDB Collection",
    gsea_minsize      = "Min size",
    gsea_maxsize      = "Max size",
    gsea_run          = "Run GSEA",
    gsea_running      = "Running GSEA... (ranking genes → fgsea)",
    gsea_done         = "✅ GSEA done: %d pathways tested, %d significant (padj<0.05)",
    gsea_no_genes     = "No genes available for ranking.",
    gsea_no_sets      = "No overlapping gene sets (panel may be small).",
    gsea_placeholder  = "Configure and click 'Run GSEA'",
    gsea_nes_title    = "Top pathways (NES)",
    gsea_need_deg     = "Run a DEG analysis first (on the DEG sub-tab).",
    gsea_uses_deg     = "GSEA uses the DEG result from the comparison above (genes ranked by avg_log2FC).",

    # Spatial
    spatial_settings  = "Spatial plot settings",
    spatial_color     = "Color variable",
    spatial_sample_col = "Sample column",
    spatial_sample    = "Sample",
    spatial_pt        = "Point size",
    spatial_flip      = "Flip Y axis",
    spatial_none      = "No spatial coordinates (x/y) found in this object.",
    spatial_ds        = "Spatial dataset",
    spatial_seg       = "Show segmentation (cell boundaries)",
    spatial_by        = "Color by",
    spatial_by_meta   = "Metadata variable",
    spatial_by_gene   = "Gene expression",
    spatial_gene      = "Gene",
    spatial_highlight = "Target clusters (highlight + neighborhood, multi-select)",
    spatial_hl_alpha  = "Highlight: opacity",
    spatial_hl_size   = "Highlight: point size",
    spatial_other_alpha = "Others: opacity",
    spatial_other_size  = "Others: point size",
    spatial_map       = "Map",
    spatial_nbr       = "Neighborhood",
    spatial_nbr_need  = "Color by a categorical metadata variable and select target clusters.",
    spatial_nbr_x     = "Distance (nearest, log scale)",
    spatial_nbr_y     = "Cumulative fraction of target cells",
    spatial_nbr_title = "Nearest-neighbor distance from %s to each cluster (closer = left)",
    spatial_nbr_run   = "Run neighborhood analysis",
    spatial_nbr_help  = "For each cell of a target cluster, computes the distance to the nearest cell of every other cluster and draws the cumulative distribution (ECDF). A curve further to the left means that cluster sits closer to the target (y=0.5 is the median distance). Distances are computed within the same sample only. With multiple samples it computes per sample and shows mean +/- SD as a confidence band.",
    spatial_co        = "Co-occurrence",
    spatial_co_run    = "Compute co-occurrence",
    spatial_co_x      = "Distance",
    spatial_co_y      = "Co-occurrence ratio (P(j|i,d)/P(j))",
    spatial_co_title  = "Co-occurrence of clusters around %s (>1: enriched, <1: depleted)",
    spatial_co_help   = "Inspired by squidpy's co_occurrence. For target cluster i, the conditional probability P(j | i, d) of finding cluster j within a distance ring d of i is divided by the overall fraction P(j); plotted vs distance. >1 means j is enriched near i, <1 depleted. Computed per sample; mean +/- SD shown as a band.",
    spatial_ne        = "Neighbors enrichment",
    spatial_ne_run    = "Compute neighbors enrichment",
    spatial_ne_k      = "Neighbors k",
    spatial_ne_perm   = "Permutations",
    spatial_ne_stars  = "Show significance asterisks",
    spatial_ne_prog   = "Computing neighbors enrichment",
    spatial_eta       = "%d/%d (about %s left)",
    spatial_ne_title  = "Neighbors-enrichment z-score (positive: adjacent, negative: avoidance)",
    spatial_ne_help   = "Inspired by squidpy's nhood_enrichment. Builds a spatial kNN graph and compares the number of edges between each cluster pair to a permutation null (shuffled labels) as a z-score. Positive = that cluster pair is adjacent more than expected, negative = avoidance. Shown as a z-score matrix heatmap (within-sample edges only). Significant pairs are marked with asterisks (z -> two-sided p via normal approx -> BH; *<0.05, **<0.01, ***<0.001). Rows and columns are hierarchically clustered.",
    spatial_need_run  = "Select targets and click Compute.",

    # Placeholders
    placeholder_load  = "\U0001F4C2 Select an RDS file and click 'Load'",
    placeholder_deg   = "Configure groups above and click 'Run DEG Analysis'",
    no_umap_title     = "⚠️ This Seurat object does not contain a UMAP reduction",
    no_umap_msg       = "Run RunUMAP() first and re-save the object."
  )
)

# =============================================================================
# UI
# =============================================================================
ui <- page_sidebar(
  title = "Seurat Viewer",
  theme = bs_theme(
    version = 5,
    bootswatch = "darkly",
    primary = "#6ea8fe",
    "navbar-bg" = "#1a1d23"
  ),

  # --- 計算中インジケータ（Shinyが計算中=html.shiny-busy のときCSSで表示）---
  # 経過時間をライブ表示（JSで shiny-busy を監視してカウント）。
  tags$head(tags$style(HTML("
    #app-busy { display: none; position: fixed; top: 64px; right: 22px; z-index: 100000;
      background: rgba(33,37,41,0.94); color: #fff; padding: 9px 16px; border-radius: 8px;
      font-size: 14px; box-shadow: 0 3px 10px rgba(0,0,0,0.35); align-items: center; gap: 9px; }
    html.shiny-busy #app-busy { display: inline-flex; }
  "))),
  tags$div(id = "app-busy",
    tags$span(class = "spinner-border spinner-border-sm", role = "status"),
    tags$span("計算中… / Computing…"),
    tags$span(id = "app-busy-secs", style = "opacity:0.85;")),
  tags$script(HTML("
    (function(){
      var t0=null, iv=null;
      function fmt(s){ return s<60 ? s+'s' : Math.floor(s/60)+'m'+(s%60)+'s'; }
      function tick(){ var el=document.getElementById('app-busy-secs');
        if(el && t0!==null){ el.textContent='('+fmt(Math.round((Date.now()-t0)/1000))+')'; } }
      new MutationObserver(function(){
        var busy=document.documentElement.classList.contains('shiny-busy');
        if(busy && iv===null){ t0=Date.now(); tick(); iv=setInterval(tick,500); }
        else if(!busy && iv!==null){ clearInterval(iv); iv=null; t0=null; }
      }).observe(document.documentElement, {attributes:true, attributeFilter:['class']});
    })();
  ")),

  # --- サイドバー ---
  sidebar = sidebar(
    width = 460,

    # ダーク/ブライトモード切替
    input_dark_mode(id = "dark_mode", mode = "dark"),

    # 言語切替
    radioButtons("lang", NULL,
      choices = c("\U0001F1EF\U0001F1F5 日本語" = "ja", "\U0001F1EC\U0001F1E7 English" = "en"),
      selected = "ja", inline = TRUE),

    hr(),

    # 動的サイドバーコンテンツ
    uiOutput("sidebar_content_ui")
  ),

  # --- メインパネル ---
  navset_card_tab(
    id = "main_tabs",

    # Violin / Feature UMAP / Group UMAP \u30921\u3064\u306E\u5927\u30BF\u30D6\u306B\u307E\u3068\u3081\u3001
    # \u907A\u4F1D\u5B50\u30FB\u30B0\u30EB\u30FC\u30D7\u5909\u6570\u306E\u9078\u629E\u3092\u305D\u306E\u4E2D\u3067\u5171\u6709\u3059\u308B
    nav_panel(
      title = "\U0001F3BB Expression",
      value = "expression",
      card_body(
        class = "p-2",
        uiOutput("expr_panel_ui")
      )
    ),

    # DEG と GSEA を1つの大タブにまとめ、同じ比較設定・同じDEG結果を共有する
    nav_panel(
      title = "\U0001F4CA DEG / GSEA",
      value = "deg",
      card_body(
        class = "p-2",
        uiOutput("deg_panel_ui")
      )
    ),

    nav_panel(
      title = "\U0001F9F1 Composition",
      value = "composition",
      card_body(
        class = "p-2",
        uiOutput("comp_panel_ui")
      )
    ),

    # Heatmap と Dot Plot を1つの大タブにまとめ、遺伝子・クラスター選択を
    # 共有しつつ、内部の独立タブで切り替える
    nav_panel(
      title = "\U0001F525 Heatmap / Dot",
      value = "markers",
      card_body(
        class = "p-2",
        uiOutput("marker_panel_ui")
      )
    ),

    # 空間座標プロット（位置情報を含むデータ用）
    nav_panel(
      title = "\U0001F5FA️ Spatial",
      value = "spatial",
      card_body(
        class = "p-2",
        uiOutput("spatial_panel_ui")
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  # --- 翻訳ヘルパー ---
  t <- function(key) {
    lang <- input$lang %||% "ja"
    i18n[[lang]][[key]] %||% key
  }

  # --- テーマ切替をリアクティブに反映 ---
  observe({
    is_dark <- input$dark_mode == "dark"
    if (is_dark) {
      session$setCurrentTheme(
        bs_theme(version = 5, bootswatch = "darkly", primary = "#6ea8fe")
      )
    } else {
      session$setCurrentTheme(
        bs_theme(version = 5, bootswatch = "flatly", primary = "#2c7be5")
      )
    }
  })

  # --- プロットテーマをモードに応じて切替 ---
  plot_theme <- reactive({
    is_dark <- isTRUE(input$dark_mode == "dark")
    if (is_dark) {
      list(
        bg = "#2b3035",
        fg = "#dee2e6",
        fg2 = "#adb5bd",
        accent = "#8ab4f8",
        theme = theme(
          plot.background = element_rect(fill = "#2b3035", color = NA),
          panel.background = element_rect(fill = "#2b3035", color = NA),
          text = element_text(color = "#dee2e6"),
          axis.text = element_text(color = "#adb5bd"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "none",
          plot.title = element_text(size = 16, face = "bold", color = "#8ab4f8")
        ),
        theme_legend = theme(
          plot.background = element_rect(fill = "#2b3035", color = NA),
          panel.background = element_rect(fill = "#2b3035", color = NA),
          text = element_text(color = "#dee2e6"),
          axis.text = element_text(color = "#adb5bd"),
          legend.text = element_text(color = "#dee2e6"),
          plot.title = element_text(size = 16, face = "bold", color = "#8ab4f8")
        )
      )
    } else {
      list(
        bg = "#ffffff",
        fg = "#212529",
        fg2 = "#495057",
        accent = "#2c7be5",
        theme = theme(
          plot.background = element_rect(fill = "#ffffff", color = NA),
          panel.background = element_rect(fill = "#ffffff", color = NA),
          text = element_text(color = "#212529"),
          axis.text = element_text(color = "#495057"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "none",
          plot.title = element_text(size = 16, face = "bold", color = "#2c7be5")
        ),
        theme_legend = theme(
          plot.background = element_rect(fill = "#ffffff", color = NA),
          panel.background = element_rect(fill = "#ffffff", color = NA),
          text = element_text(color = "#212529"),
          axis.text = element_text(color = "#495057"),
          legend.text = element_text(color = "#212529"),
          plot.title = element_text(size = 16, face = "bold", color = "#2c7be5")
        )
      )
    }
  })

  # --- リアクティブ値 ---
  seurat_obj <- reactiveVal(NULL)
  data_loaded <- reactiveVal(FALSE)
  deg_results <- reactiveVal(NULL)
  meta_col_types <- reactiveVal(list(cat = character(0), num = character(0)))
  deg_title_label <- reactiveVal("")
  loaded_objs <- reactiveVal(list())   # 読み込んだデータセット (名前付きリスト)

  # --- リファレンス(比較用)データセット ---
  ref_obj <- reactive({
    objs <- loaded_objs()
    ds <- input$ref_ds
    if (is.null(ds) || identical(ds, "__none__") || is.null(objs[[ds]])) return(NULL)
    objs[[ds]]
  })

  # ==========================================================================
  # 動的サイドバー
  # ==========================================================================
  output$sidebar_content_ui <- renderUI({
    # input$langへの依存で言語切替時に再描画
    lang <- input$lang

    tagList(
      h5(t("data_load"), class = "text-primary mb-2"),

      selectInput(
        "rds_file", t("rds_file"),
        choices = rds_files,
        selected = rds_files[1],
        multiple = TRUE
      ),
      actionButton(
        "load_btn", t("load_btn"),
        class = "btn-primary w-100 mb-3",
        icon = icon("upload")
      ),

      # アクティブ / リファレンス データセット選択（spatial 以外のタブで表示）
      conditionalPanel(
        "input.main_tabs != 'spatial'",
        uiOutput("dataset_selectors_ui")
      ),

      hr(),

      h5(t("plot_settings"), class = "text-primary mb-2"),

      sliderInput("pt_size", t("pt_size"), min = 0, max = 2, value = 0.3, step = 0.1),
      sliderInput("umap_label_size", t("umap_label_size"), min = 2, max = 12, value = 4, step = 0.5),
      sliderInput("plot_height", t("plot_height"), min = 400, max = 3000, value = 600, step = 50),
      sliderInput("plot_width", t("plot_width"), min = 400, max = 4000, value = 800, step = 50),

      # リファレンス用のプロット設定（比較時のみ表示）
      uiOutput("ref_plot_settings_ui"),

      # 空間タブ用: プロットするデータセットの選択
      conditionalPanel(
        "input.main_tabs == 'spatial'",
        hr(),
        h5("\U0001F5FA️ Spatial", class = "text-primary mb-2"),
        uiOutput("spatial_ds_ui")
      )
    )
  })

  # --- 外部データベースリンク ---
  # 遺伝子セレクタ(gene/spatial_gene)はサーバーサイドで選択肢を投入し、
  # 大量の遺伝子でもクライアントが重くならないようにする。
  observe({
    if (!isTRUE(data_loaded())) return()
    input$lang
    if (!identical(input$main_tabs, "expression")) return()
    obj <- seurat_obj(); if (is.null(obj)) return()
    genes <- sort(rownames(obj))
    sel <- isolate(input$gene); if (is.null(sel) || !(sel %in% genes)) sel <- genes[1]
    updateSelectizeInput(session, "gene", choices = genes, selected = sel, server = TRUE)
  })
  observe({
    if (!isTRUE(data_loaded())) return()
    input$lang; input$spatial_by
    if (!identical(input$main_tabs, "spatial")) return()
    obj <- spatial_obj(); if (is.null(obj)) return()
    genes <- sort(rownames(obj))
    sel <- isolate(input$spatial_gene); if (is.null(sel) || !(sel %in% genes)) sel <- genes[1]
    updateSelectizeInput(session, "spatial_gene", choices = genes, selected = sel, server = TRUE)
  })

  output$external_links_ui <- renderUI({
    req(input$gene)
    gene <- input$gene
    tagList(
      div(
        class = "d-grid gap-1 mb-2",
        tags$a(
          href = ncbi_gene_url(gene), target = "_blank",
          class = "btn btn-outline-info btn-sm",
          icon("dna"), paste0(" ", gene, " — NCBI Gene")
        ),
        tags$a(
          href = immunexut_url(gene), target = "_blank",
          class = "btn btn-outline-success btn-sm",
          icon("microscope"), paste0(" ", gene, " — ImmuNexUT")
        ),
        tags$a(
          href = ampra2_url(gene), target = "_blank",
          class = "btn btn-outline-warning btn-sm",
          icon("bone"), paste0(" ", gene, " — AMP RA2")
        )
      )
    )
  })

  # ==========================================================================
  # 言語切替時にすでに読込済みデータの選択肢を維持
  # ==========================================================================
  # 言語切替時、gene/group_var は expr_panel_ui が、deg_group_var は deg_panel_ui が
  # それぞれ再描画して選択を維持するため、ここでの明示更新は不要。

  # ==========================================================================
  # RDSファイル読み込み
  # ==========================================================================
  # --- アクティブにするデータセットで各セレクタ・メタ情報を更新 ---
  # 現在の遺伝子/グループ選択はデータセット切替後も（存在すれば）維持する。
  activate_dataset <- function(obj) {
    seurat_obj(obj)
    data_loaded(TRUE)
    deg_results(NULL)

    # meta.dataの列を分類。gene/group_var は expr_panel_ui が、
    # deg_group_var は deg_panel_ui が選択肢を組み立てるためここでは更新しない。
    meta <- obj@meta.data
    cat_cols <- names(meta)[sapply(meta, function(x) {
      is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 50)
    })]
    num_cols <- names(meta)[sapply(meta, function(x) {
      is.numeric(x) && length(unique(x)) > 50
    })]
    meta_col_types(list(cat = cat_cols, num = num_cols))

    if (is.null(find_umap_reduction(obj))) {
      showNotification(t("notify_no_umap"), type = "warning", duration = 8)
    }
  }

  # --- 複数RDSの読み込み ---
  # withProgress を使うことで、ブロッキングする readRDS の前/最中でも
  # 「読み込み中」インジケータがブラウザに即時表示され、完了まで残り続ける。
  observeEvent(input$load_btn, {
    req(input$rds_file)
    files <- input$rds_file
    n <- length(files)
    withProgress(message = t("notify_loading"), value = 0, {
      tryCatch({
        objs <- list()
        for (i in seq_along(files)) {
          f <- files[i]
          incProgress(0, detail = sprintf("%s (%d/%d)", f, i, n))
          o <- readRDS(file.path(app_dir, f))
          if (inherits(o, "Seurat")) objs[[f]] <- o
          incProgress(1 / n)
        }
        if (length(objs) == 0) {
          showNotification(t("notify_not_seurat"), type = "error", id = "loading")
          return()
        }
        incProgress(0, detail = t("notify_finalizing"))
        loaded_objs(objs)
        activate_dataset(objs[[1]])

        o1 <- objs[[1]]
        showNotification(
          sprintf(t("notify_load_done"),
                  format(ncol(o1), big.mark = ","),
                  format(nrow(o1), big.mark = ",")),
          type = "message", id = "loading", duration = 5
        )
      }, error = function(e) {
        showNotification(paste(t("notify_error"), e$message), type = "error", id = "loading")
      })
    })
  })

  # --- アクティブデータセットの切り替え ---
  observeEvent(input$active_ds, {
    objs <- loaded_objs()
    if (!is.null(input$active_ds) && !is.null(objs[[input$active_ds]])) {
      activate_dataset(objs[[input$active_ds]])
    }
  }, ignoreInit = TRUE)

  # --- アクティブ / リファレンス セレクタ ---
  output$dataset_selectors_ui <- renderUI({
    lang <- input$lang
    objs <- loaded_objs()
    if (length(objs) == 0) return(NULL)
    nms <- names(objs)
    active_sel <- isolate(input$active_ds)
    if (is.null(active_sel) || !(active_sel %in% nms)) active_sel <- nms[1]
    ref_sel <- isolate(input$ref_ds)
    if (is.null(ref_sel) || !(ref_sel %in% c("__none__", nms))) ref_sel <- "__none__"
    tagList(
      selectInput("active_ds", t("active_ds"), choices = nms, selected = active_sel),
      selectInput("ref_ds", t("ref_ds"),
                  choices = c(setNames("__none__", t("ref_none")), setNames(nms, nms)),
                  selected = ref_sel),
      div(class = "mb-3")
    )
  })

  # ==========================================================================
  # DEG解析パネル UI（動的・言語対応）
  # ==========================================================================
  output$deg_panel_ui <- renderUI({
    lang <- input$lang  # 言語切替で再描画

    # DEGグループ変数の選択肢を直接生成
    col_types <- meta_col_types()
    cat_cols <- col_types$cat
    num_cols <- col_types$num
    deg_choices <- NULL
    deg_selected <- NULL
    if (length(cat_cols) > 0 || length(num_cols) > 0) {
      deg_choices <- c(
        setNames(cat_cols, paste0(t("deg_cat_prefix"), cat_cols)),
        setNames(num_cols, paste0(t("deg_num_prefix"), num_cols))
      )
      # 現在の選択を維持
      current <- isolate(input$deg_group_var)
      deg_selected <- if (!is.null(current) && current %in% c(cat_cols, num_cols)) {
        current
      } else if ("seurat_clusters" %in% cat_cols) {
        "seurat_clusters"
      } else if (length(cat_cols) > 0) {
        cat_cols[1]
      } else {
        num_cols[1]
      }
    }

    # 解析対象遺伝子の候補（自分で選ぶ用）
    deg_gene_choices <- if (data_loaded()) sort(rownames(seurat_obj())) else character(0)
    feat_mode_sel <- isolate(input$deg_feature_mode) %||% "all"

    tagList(
      # DEG設定パネル
      div(
        class = "card mb-3",
        div(
          class = "card-body",
          h6(t("deg_settings"), class = "card-title text-primary"),
          fluidRow(
            column(12,
              selectInput("deg_group_var", t("deg_group_var"),
                          choices = deg_choices,
                          selected = deg_selected)
            )
          ),
          # カテゴリカル/数値で動的に切り替わるUI
          uiOutput("deg_group_settings_ui"),

          # 解析する遺伝子: 全 / セット / セット+カスタム / カスタムのみ
          fluidRow(
            column(12,
              radioButtons("deg_feature_mode", t("deg_feature_mode"),
                choices = c(
                  setNames("all",        t("deg_feat_all")),
                  setNames("set",        t("deg_feat_set")),
                  setNames("set_custom", t("deg_feat_set_custom")),
                  setNames("custom",     t("deg_feat_custom"))
                ),
                selected = feat_mode_sel, inline = TRUE)
            )
          ),
          conditionalPanel(
            "input.deg_feature_mode == 'set' || input.deg_feature_mode == 'set_custom'",
            fluidRow(
              column(12,
                selectInput("deg_marker_set", t("deg_marker_set"),
                            choices = names(marker_sets),
                            selected = isolate(input$deg_marker_set) %||% names(marker_sets)[1],
                            multiple = TRUE)
              )
            )
          ),
          conditionalPanel(
            "input.deg_feature_mode == 'custom' || input.deg_feature_mode == 'set_custom'",
            fluidRow(
              column(12,
                selectizeInput("deg_custom_genes", t("deg_custom_genes"),
                               choices = deg_gene_choices,
                               selected = isolate(input$deg_custom_genes),
                               multiple = TRUE,
                               options = list(placeholder = t("gene_placeholder"),
                                              maxOptions = 1000))
              )
            )
          ),

          fluidRow(
            column(4,
              numericInput("deg_logfc", t("deg_logfc"), value = 0.25, min = 0, step = 0.1)
            ),
            column(4,
              numericInput("deg_pval", t("deg_pval"), value = 0.05, min = 0, max = 1, step = 0.01)
            ),
            column(4,
              div(style = "margin-top: 24px;",
                actionButton("run_deg", t("deg_run"),
                             class = "btn-primary w-100",
                             icon = icon("flask"))
              )
            )
          )
        )
      ),

      # 上の比較設定（同じDEG結果）を DEG / GSEA の内部タブで共有
      navset_card_tab(
        id = "deg_subtab",
        nav_panel(title = "\U0001F4CA DEG", value = "deg_res",
                  card_body(class = "p-2", uiOutput("deg_results_ui"))),
        nav_panel(title = "\U0001F9EC GSEA", value = "gsea_res",
                  card_body(class = "p-2", uiOutput("gsea_inner_ui")))
      )
    )
  })

  # --- DEGグループ変数の型判定 ---
  deg_var_is_numeric <- reactive({
    req(input$deg_group_var)
    col_types <- meta_col_types()
    input$deg_group_var %in% col_types$num
  })

  # --- DEGで解析する遺伝子（NULL = 全遺伝子） ---
  deg_features <- reactive({
    mode <- input$deg_feature_mode %||% "all"
    if (mode == "all") return(NULL)
    obj <- seurat_obj()
    req(obj)
    set_genes <- character(0)
    if (mode %in% c("set", "set_custom") && length(input$deg_marker_set) > 0) {
      sets_present <- input$deg_marker_set[input$deg_marker_set %in% names(marker_sets)]
      set_genes <- unique(do.call(c, lapply(sets_present,
        function(s) marker_set_to_df(marker_sets[[s]])$feature)))
    }
    custom <- if (mode %in% c("custom", "set_custom")) input$deg_custom_genes else character(0)
    intersect(unique(c(set_genes, custom)), rownames(obj))
  })

  # --- DEGグループ設定UI（カテゴリ/数値で切替） ---
  output$deg_group_settings_ui <- renderUI({
    req(input$deg_group_var, seurat_obj())
    lang <- input$lang  # 言語依存

    if (deg_var_is_numeric()) {
      # 数値列の場合: パーセンタイルスライダー
      fluidRow(
        column(6,
          sliderInput("deg_percentile", t("deg_percentile"),
                      min = 5, max = 50, value = 25, step = 5,
                      post = "%")
        ),
        column(6,
          div(class = "mt-3",
            uiOutput("deg_numeric_info")
          )
        )
      )
    } else {
      # カテゴリカル列の場合: グループ選択（系統順に並べて選びやすく）
      obj <- seurat_obj()
      groups <- cluster_level_order(obj@meta.data[[input$deg_group_var]])
      group2_choices <- c(
        setNames("__ALL_OTHERS__", t("deg_all_others")),
        setNames(groups, groups)
      )
      fluidRow(
        column(6,
          selectInput("deg_ident1", t("deg_group1"),
                      choices = groups, selected = groups[1])
        ),
        column(6,
          selectInput("deg_ident2", t("deg_group2"),
                      choices = group2_choices,
                      selected = if (length(groups) >= 2) groups[2] else "__ALL_OTHERS__",
                      multiple = TRUE)
        )
      )
    }
  })

  # --- 数値DEGの細胞数表示 ---
  output$deg_numeric_info <- renderUI({
    req(seurat_obj(), input$deg_group_var, deg_var_is_numeric(), input$deg_percentile)
    lang <- input$lang  # 言語依存
    obj <- seurat_obj()
    vals <- obj@meta.data[[input$deg_group_var]]
    n_total <- length(vals)
    pct <- input$deg_percentile / 100
    n_group <- floor(n_total * pct)
    div(
      class = "small",
      p(class = "mb-1", sprintf(t("deg_total_cells"), format(n_total, big.mark = ","))),
      p(class = "mb-1", sprintf(t("deg_top_cells"), input$deg_percentile, format(n_group, big.mark = ","))),
      p(class = "mb-0", sprintf(t("deg_bottom_cells"), input$deg_percentile, format(n_group, big.mark = ",")))
    )
  })

  # --- ヘルパー: 利用可能な UMAP 関連 reduction 名を全て返す ---
  # "umap"（完全一致）を優先し、続けて名前に "umap" を含む reduction を返す。
  # 見つからなければ character(0)。
  list_umap_reductions <- function(obj) {
    if (is.null(obj)) return(character(0))
    red_names <- names(obj@reductions)
    hits <- red_names[grepl("umap", red_names, ignore.case = TRUE)]
    # "umap" 完全一致を先頭に
    exact <- hits[hits == "umap"]
    rest  <- hits[hits != "umap"]
    c(exact, rest)
  }

  # --- ヘルパー: デフォルトで使用する UMAP reduction 名 ---
  # "umap" があればそれを、なければ "umap" を含む最初の reduction を返す。
  find_umap_reduction <- function(obj) {
    hits <- list_umap_reductions(obj)
    if (length(hits) > 0) return(hits[1])
    NULL
  }

  # --- ヘルパー: 利用可能な UMAP reduction 名（リアクティブ） ---
  umap_reductions <- reactive({
    list_umap_reductions(seurat_obj())
  })

  # --- ヘルパー: 使用する UMAP reduction 名（ユーザー選択を優先） ---
  umap_reduction <- reactive({
    choices <- umap_reductions()
    if (length(choices) == 0) return(NULL)
    sel <- input$umap_reduction
    if (!is.null(sel) && sel %in% choices) return(sel)
    choices[1]
  })

  # --- UMAP reduction セレクタ（複数ある場合のみ表示） ---
  output$umap_reduction_ui <- renderUI({
    choices <- umap_reductions()
    if (length(choices) < 2) return(NULL)
    selectInput(
      "umap_reduction", t("umap_reduction"),
      choices = choices,
      selected = isolate(umap_reduction())
    )
  })

  # --- ヘルパー: UMAPの有無を確認 ---
  has_umap <- reactive({
    length(umap_reductions()) > 0
  })

  # --- プレースホルダUI ---
  placeholder_ui <- function() {
    div(
      class = "text-center text-muted py-5",
      h4(t("placeholder_load"))
    )
  }

  no_umap_ui <- function() {
    div(
      class = "text-center text-warning py-5",
      h4(t("no_umap_title")),
      p(t("no_umap_msg"))
    )
  }

  # ==========================================================================
  # 比較表示（アクティブ vs リファレンス）ヘルパー
  # ==========================================================================
  active_name <- reactive({ input$active_ds %||% "" })
  ref_name    <- reactive({ if (is.null(ref_obj())) NULL else input$ref_ds })

  # 欠損時のプレースホルダ ggplot
  empty_panel <- function(msg) {
    ggplot() +
      annotate("text", x = 0, y = 0, label = msg, size = 5, color = "grey50") +
      theme_void()
  }

  # --- サイズ・点サイズのアクセサ（アクティブ / リファレンス別） ---
  act_pt <- reactive({ input$pt_size %||% 0.3 })
  ref_pt <- reactive({ input$ref_pt_size %||% (input$pt_size %||% 0.3) })
  act_h  <- reactive({ paste0(input$plot_height %||% 600, "px") })
  act_w  <- reactive({ paste0(input$plot_width %||% 800, "px") })
  ref_h  <- reactive({ paste0(input$ref_plot_height %||% (input$plot_height %||% 600), "px") })
  ref_w  <- reactive({ paste0(input$ref_plot_width %||% (input$plot_width %||% 800), "px") })

  # リファレンス用プロット設定スライダー（比較時のみ表示）
  output$ref_plot_settings_ui <- renderUI({
    lang <- input$lang
    if (is.null(ref_obj())) return(NULL)
    tagList(
      hr(),
      h6(t("ref_plot_settings"), class = "text-primary mb-2"),
      sliderInput("ref_pt_size", t("ref_pt_size"), min = 0, max = 2,
                  value = isolate(input$ref_pt_size) %||% (isolate(input$pt_size) %||% 0.3), step = 0.1),
      sliderInput("ref_plot_height", t("ref_plot_height"), min = 400, max = 3000,
                  value = isolate(input$ref_plot_height) %||% (isolate(input$plot_height) %||% 600), step = 50),
      sliderInput("ref_plot_width", t("ref_plot_width"), min = 400, max = 4000,
                  value = isolate(input$ref_plot_width) %||% (isolate(input$plot_width) %||% 800), step = 50)
    )
  })

  # --- アクティブ/リファレンスの2プロットを左右に並べるUI（個別サイズ） ---
  # have_ref = FALSE なら active 単独。各パネルにデータセット名ラベル。
  side_by_side_ui <- function(active_out, ref_out, have_ref, plotly = FALSE) {
    out_fn <- if (plotly) {
      function(id, height, width) plotly::plotlyOutput(id, height = height, width = width)
    } else {
      function(id, height, width) plotOutput(id, height = height, width = width)
    }
    if (!have_ref) {
      return(out_fn(active_out, act_h(), act_w()))
    }
    div(
      style = "display:flex; gap:14px; overflow-x:auto; align-items:flex-start;",
      div(
        h6(active_name(), class = "text-center text-muted mb-1"),
        out_fn(active_out, act_h(), act_w())
      ),
      div(
        h6(ref_name() %||% "", class = "text-center text-muted mb-1"),
        out_fn(ref_out, ref_h(), ref_w())
      )
    )
  }

  # 2つのggplotを左右に並べる（refがNULLなら単独）。各パネルに名前を付与。
  combine_gg <- function(p_active, p_ref, name_active = NULL, name_ref = NULL) {
    if (is.null(p_ref)) return(p_active)
    if (!is.null(name_active)) p_active <- p_active + ggtitle(name_active)
    if (!is.null(name_ref))    p_ref    <- p_ref + ggtitle(name_ref)
    patchwork::wrap_plots(p_active, p_ref, ncol = 2)
  }

  # --- リファレンスのカテゴリ列・グループ変数（アクティブと別に選択可） ---
  ref_cat_cols <- reactive({
    rb <- ref_obj()
    if (is.null(rb)) return(character(0))
    meta <- rb@meta.data
    names(meta)[sapply(meta, function(x) {
      is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 50)
    })]
  })

  # リファレンスで使うグループ変数（未指定ならアクティブと同名→seurat_clusters→先頭）
  ref_group_var <- reactive({
    cols <- ref_cat_cols()
    if (length(cols) == 0) return(NULL)
    sel <- input$ref_group_var
    if (!is.null(sel) && sel %in% cols) return(sel)
    if (!is.null(input$group_var) && input$group_var %in% cols) return(input$group_var)
    if ("seurat_clusters" %in% cols) return("seurat_clusters")
    cols[1]
  })

  # Composition のリファレンス X軸（複数選択可。未指定ならアクティブ→先頭）
  ref_comp_x <- reactive({
    cols <- ref_cat_cols()
    if (length(cols) == 0) return(NULL)
    sel <- input$comp_x_ref
    if (!is.null(sel) && length(sel) > 0 && all(sel %in% cols)) return(sel)
    ax <- input$comp_x
    if (!is.null(ax) && all(ax %in% cols)) return(ax)
    cols[1]
  })

  output$comp_x_ref_ui <- renderUI({
    lang <- input$lang
    if (is.null(ref_obj())) return(NULL)
    cols <- ref_cat_cols()
    if (length(cols) == 0) return(NULL)
    sel <- isolate(input$comp_x_ref)
    if (is.null(sel) || length(sel) == 0 || !all(sel %in% cols)) {
      ax <- isolate(input$comp_x)
      sel <- if (!is.null(ax) && all(ax %in% cols)) ax else cols[1]
    }
    selectInput("comp_x_ref", t("comp_x_ref"), choices = cols, selected = sel,
                multiple = TRUE)
  })

  output$ref_group_var_ui <- renderUI({
    lang <- input$lang
    cols <- ref_cat_cols()
    if (length(cols) == 0) return(NULL)
    sel <- isolate(input$ref_group_var)
    if (is.null(sel) || !(sel %in% cols)) {
      ag <- isolate(input$group_var)
      sel <- if (!is.null(ag) && ag %in% cols) ag
             else if ("seurat_clusters" %in% cols) "seurat_clusters" else cols[1]
    }
    selectInput("ref_group_var", t("ref_group_var"), choices = cols, selected = sel)
  })

  # ==========================================================================
  # Expression 大タブ（Violin / Feature UMAP / Group UMAP の共通設定 + 内部タブ）
  # ==========================================================================
  output$expr_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    ct <- meta_col_types()
    cat_cols <- ct$cat
    if (length(cat_cols) == 0) {
      return(div(class = "text-center text-muted py-4", h5(t("comp_no_cat"))))
    }
    genes <- sort(rownames(seurat_obj()))
    cur_gene <- isolate(input$gene)
    gsel <- if (!is.null(cur_gene) && cur_gene %in% genes) cur_gene else genes[1]
    cur_group <- isolate(input$group_var)
    gv_sel <- if (!is.null(cur_group) && cur_group %in% cat_cols) {
      cur_group
    } else if ("seurat_clusters" %in% cat_cols) "seurat_clusters" else cat_cols[1]

    tagList(
      div(class = "card mb-3", div(class = "card-body",
        fluidRow(
          column(6,
            # 遺伝子リストはサーバーサイドで投入（多数の選択肢でも高速）
            selectizeInput("gene", t("gene_name"),
                           choices = if (!is.null(gsel)) stats::setNames(gsel, gsel) else NULL,
                           selected = gsel,
                           options = list(placeholder = t("gene_placeholder"),
                                          maxOptions = 1000)),
            uiOutput("external_links_ui")
          ),
          column(6,
            selectInput("group_var", t("group_var"), choices = cat_cols, selected = gv_sel),
            uiOutput("ref_group_var_ui"),
            uiOutput("umap_reduction_ui")
          )
        )
      )),
      navset_card_tab(
        id = "expr_subtab",
        nav_panel(title = "\U0001F3BB Violin", value = "violin",
                  card_body(class = "p-2", uiOutput("violin_ui"))),
        nav_panel(title = "\U0001F5FA️ Feature UMAP", value = "feature_umap",
                  card_body(class = "p-2", uiOutput("feature_umap_ui"))),
        nav_panel(title = "\U0001F3F7️ Group UMAP", value = "group_umap",
                  card_body(class = "p-2", uiOutput("group_umap_ui")))
      )
    )
  })

  # ==========================================================================
  # Violin Plot
  # ==========================================================================
  output$violin_ui <- renderUI({
    if (!data_loaded()) return(placeholder_ui())
    side_by_side_ui("violin_plot", "violin_plot_ref", !is.null(ref_obj()))
  })

  build_violin <- function(obj, group_var, pt_size) {
    pt <- plot_theme()
    if (!(input$gene %in% rownames(obj))) return(empty_panel(t("ref_missing")))
    if (is.null(group_var) || !(group_var %in% names(obj@meta.data))) return(empty_panel(t("ref_missing")))
    # クラスターを系統順 factor に並べ替え（軸・色順を統一）
    Idents(obj) <- factor(
      as.character(obj@meta.data[[group_var]]),
      levels = cluster_level_order(obj@meta.data[[group_var]])
    )
    p <- VlnPlot(obj, features = input$gene, pt.size = pt_size) + pt$theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    lin_cols <- lineage_colors_or_null(levels(Idents(obj)))
    if (!is.null(lin_cols)) p <- p + scale_fill_manual(values = lin_cols)
    p
  }

  output$violin_plot <- renderPlot({
    req(seurat_obj(), input$gene, input$group_var)
    build_violin(seurat_obj(), input$group_var, act_pt())
  }, bg = "transparent")

  output$violin_plot_ref <- renderPlot({
    req(ref_obj(), input$gene)
    build_violin(ref_obj(), ref_group_var(), ref_pt())
  }, bg = "transparent")

  # ==========================================================================
  # Feature UMAP
  # ==========================================================================
  output$feature_umap_ui <- renderUI({
    if (!data_loaded()) return(placeholder_ui())
    if (!has_umap()) return(no_umap_ui())
    side_by_side_ui("feature_umap_plot", "feature_umap_plot_ref", !is.null(ref_obj()))
  })

  build_feature <- function(obj, reduction, pt_size) {
    pt <- plot_theme()
    if (is.null(reduction)) return(empty_panel(t("no_umap_title")))
    if (!(input$gene %in% rownames(obj))) return(empty_panel(t("ref_missing")))
    FeaturePlot(obj, features = input$gene, reduction = reduction,
                pt.size = pt_size) + pt$theme_legend +
      theme(legend.position = "bottom")
  }

  output$feature_umap_plot <- renderPlot({
    req(seurat_obj(), input$gene, has_umap())
    build_feature(seurat_obj(), umap_reduction(), act_pt())
  }, bg = "transparent")

  output$feature_umap_plot_ref <- renderPlot({
    req(ref_obj(), input$gene)
    build_feature(ref_obj(), find_umap_reduction(ref_obj()), ref_pt())
  }, bg = "transparent")

  # ==========================================================================
  # Group UMAP
  # ==========================================================================
  output$group_umap_ui <- renderUI({
    if (!data_loaded()) return(placeholder_ui())
    if (!has_umap()) return(no_umap_ui())
    # plotly があればインタラクティブ（点にホバーでクラスター名を表示）
    if (requireNamespace("plotly", quietly = TRUE)) {
      side_by_side_ui("group_umap_plotly", "group_umap_plotly_ref",
                      !is.null(ref_obj()), plotly = TRUE)
    } else {
      side_by_side_ui("group_umap_plot", "group_umap_plot_ref", !is.null(ref_obj()))
    }
  })

  build_group_umap <- function(obj, group_var, reduction, pt_size, label_size = 4) {
    pt <- plot_theme()
    if (is.null(reduction)) return(empty_panel(t("no_umap_title")))
    if (is.null(group_var) || !(group_var %in% names(obj@meta.data))) return(empty_panel(t("ref_missing")))
    # クラスターを系統順 factor に並べ替え（凡例・色順を統一）
    obj@meta.data[[group_var]] <- factor(
      as.character(obj@meta.data[[group_var]]),
      levels = cluster_level_order(obj@meta.data[[group_var]])
    )
    n_lev <- nlevels(obj@meta.data[[group_var]])
    p <- DimPlot(obj, reduction = reduction, group.by = group_var,
                 label = FALSE, pt.size = pt_size) + pt$theme_legend +
      # クラスター数が多いと凡例が巨大になり描画領域が足りずエラーになるため、
      # 多い場合は凡例を非表示（プロット上のラベルで識別できる）
      theme(legend.position = if (n_lev > 30) "none" else "bottom")
    lin_cols <- lineage_colors_or_null(obj@meta.data[[group_var]])
    if (!is.null(lin_cols)) p <- p + scale_color_manual(values = lin_cols)
    # ラベルは自前で repel 描画（max.overlaps=Inf で全ラベル表示・重なり回避）。
    # Seurat の repel=TRUE はクラスター数が多いと構造化エラーを出すため使わない。
    emb <- tryCatch(Embeddings(obj, reduction = reduction)[, 1:2, drop = FALSE],
                    error = function(e) NULL)
    if (!is.null(emb)) {
      lab_df <- data.frame(Ux = emb[, 1], Uy = emb[, 2],
                           Cl = as.character(obj@meta.data[[group_var]]),
                           stringsAsFactors = FALSE)
      lab_df <- lab_df[!is.na(lab_df$Cl), , drop = FALSE]
      if (nrow(lab_df) > 0) {
        cent <- aggregate(cbind(Ux, Uy) ~ Cl, data = lab_df, FUN = median)
        p <- p + ggrepel::geom_text_repel(
          data = cent, aes(x = Ux, y = Uy, label = Cl),
          size = label_size, max.overlaps = Inf, inherit.aes = FALSE,
          color = pt$fg, bg.color = pt$bg, bg.r = 0.12, seed = 1)
      }
    }
    p
  }

  output$group_umap_plot <- renderPlot({
    req(seurat_obj(), input$group_var, has_umap())
    build_group_umap(seurat_obj(), input$group_var, umap_reduction(), act_pt(),
                     input$umap_label_size %||% 4)
  }, bg = "transparent")

  output$group_umap_plot_ref <- renderPlot({
    req(ref_obj(), has_umap())
    build_group_umap(ref_obj(), ref_group_var(), find_umap_reduction(ref_obj()), ref_pt(),
                     input$umap_label_size %||% 4)
  }, bg = "transparent")

  # --- インタラクティブ版（点にホバーでクラスター名を表示） ---
  # ggplotly は10万点規模で非常に遅いため、plotly(scattergl)で直接描く。
  if (requireNamespace("plotly", quietly = TRUE)) {
    plotly_msg <- function(msg) {
      plotly::layout(plotly::plotly_empty(type = "scatter", mode = "markers"),
        annotations = list(text = msg, showarrow = FALSE,
                           font = list(size = 16, color = "grey50")))
    }
    build_group_umap_plotly <- function(obj, group_var, reduction, pt_size) {
      pt <- plot_theme()
      if (is.null(reduction)) return(plotly_msg(t("no_umap_title")))
      if (is.null(group_var) || !(group_var %in% names(obj@meta.data))) {
        return(plotly_msg(t("ref_missing")))
      }
      emb <- tryCatch(Embeddings(obj, reduction = reduction)[, 1:2, drop = FALSE],
                      error = function(e) NULL)
      if (is.null(emb)) return(plotly_msg(t("no_umap_title")))
      levs <- cluster_level_order(obj@meta.data[[group_var]])
      df <- data.frame(X = emb[, 1], Y = emb[, 2],
                       Cl = factor(as.character(obj@meta.data[[group_var]]), levels = levs))
      df <- df[!is.na(df$Cl), , drop = FALSE]
      n_lev <- nlevels(droplevels(df$Cl))
      axn <- colnames(emb)
      lin_cols <- lineage_colors_or_null(obj@meta.data[[group_var]])
      pal <- if (!is.null(lin_cols)) unname(lin_cols[levs]) else scales::hue_pal()(length(levs))
      msize <- max(2, (pt_size %||% 0.3) * 6)
      add_layout <- function(p, showleg) {
        # シンプルに: x=0/y=0 のゼロライン・グリッドを消す
        ax <- function(ttl) list(title = ttl, zeroline = FALSE, showgrid = FALSE,
                                 showline = FALSE, ticks = "")
        plotly::layout(p,
          xaxis = ax(axn[1]), yaxis = ax(axn[2]),
          paper_bgcolor = pt$bg, plot_bgcolor = pt$bg,
          font = list(color = pt$fg), showlegend = showleg,
          legend = list(font = list(color = pt$fg)))
      }
      if (n_lev <= 30) {
        # クラスターごとにトレース（凡例あり）
        p <- plotly::plot_ly(df, x = ~X, y = ~Y, color = ~Cl, colors = pal,
                             type = "scattergl", mode = "markers",
                             marker = list(size = msize, opacity = 0.85),
                             text = ~Cl, hoverinfo = "text")
        add_layout(p, TRUE)
      } else {
        # 多クラスターは単一トレース（高速・凡例なし、ホバーで識別）
        cmap <- setNames(pal, levs)
        df$.col <- cmap[as.character(df$Cl)]
        p <- plotly::plot_ly(df, x = ~X, y = ~Y, type = "scattergl", mode = "markers",
                             marker = list(color = ~I(.col), size = msize, opacity = 0.85),
                             text = ~Cl, hoverinfo = "text")
        add_layout(p, FALSE)
      }
    }
    output$group_umap_plotly <- plotly::renderPlotly({
      req(seurat_obj(), input$group_var, has_umap())
      build_group_umap_plotly(seurat_obj(), input$group_var, umap_reduction(), act_pt())
    })
    output$group_umap_plotly_ref <- plotly::renderPlotly({
      req(ref_obj(), has_umap())
      build_group_umap_plotly(ref_obj(), ref_group_var(),
                              find_umap_reduction(ref_obj()), ref_pt())
    })
  }

  # ==========================================================================
  # Composition (組成プロット)
  # ==========================================================================
  output$comp_panel_ui <- renderUI({
    lang <- input$lang  # 言語切替で再描画
    if (!data_loaded()) return(placeholder_ui())

    col_types <- meta_col_types()
    cat_cols <- col_types$cat
    if (length(cat_cols) == 0) {
      return(div(class = "text-center text-muted py-4", h5(t("comp_no_cat"))))
    }

    # 現在の選択を維持しつつデフォルトを決定
    cur_cluster <- isolate(input$comp_cluster)
    cluster_sel <- if (!is.null(cur_cluster) && cur_cluster %in% cat_cols) {
      cur_cluster
    } else if ("seurat_clusters" %in% cat_cols) {
      "seurat_clusters"
    } else {
      cat_cols[1]
    }
    cur_x <- isolate(input$comp_x)
    x_sel <- if (!is.null(cur_x) && all(cur_x %in% cat_cols) && length(cur_x) > 0) cur_x else cat_cols[1]
    cur_facets <- isolate(input$comp_facets)
    facet_sel <- cur_facets[cur_facets %in% cat_cols]
    cur_scheme <- isolate(input$comp_color_scheme)
    scheme_sel <- if (!is.null(cur_scheme)) cur_scheme else "manual"

    tagList(
      div(
        class = "card mb-3",
        div(
          class = "card-body",
          h6(t("comp_settings"), class = "card-title text-primary"),
          fluidRow(
            column(4,
              selectInput("comp_cluster", t("comp_cluster"),
                          choices = cat_cols, selected = cluster_sel)
            ),
            column(4,
              selectInput("comp_x", t("comp_x"),
                          choices = cat_cols, selected = x_sel, multiple = TRUE)
            ),
            column(4,
              selectInput("comp_facets", t("comp_facets"),
                          choices = cat_cols, selected = facet_sel,
                          multiple = TRUE)
            )
          ),
          # 大分類で絞り込み → 描画するクラスター（選択範囲内で割合を再計算）
          coarse_block("comp", cat_cols),
          uiOutput("comp_clusters_ui"),
          # リファレンス用のクラスター選択（比較時のみ）
          ref_cluster_controls_ui("comp"),
          uiOutput("comp_x_ref_ui"),
          fluidRow(
            column(12,
              radioButtons("comp_color_scheme", t("comp_color_scheme"),
                           choices = c(setNames("manual",  t("comp_scheme_manual")),
                                       setNames("lineage", t("comp_scheme_lineage"))),
                           selected = scheme_sel, inline = TRUE)
            )
          ),
          plot_run_btn("comp_run")
        )
      ),
      uiOutput("comp_plot_ui")
    )
  })

  output$comp_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    if ((input$comp_run %||% 0) == 0) return(run_hint_ui())
    # plotly があればインタラクティブ（バーにホバーでクラスター表示）
    if (requireNamespace("plotly", quietly = TRUE)) {
      if (is.null(ref_obj())) {
        plotly::plotlyOutput("comp_plotly", height = act_h(), width = act_w())
      } else {
        # リファレンス比較は2つのplotlyを左右に並べる（個別サイズ）
        div(
          style = "display:flex; gap:14px; overflow-x:auto; align-items:flex-start;",
          div(h6(active_name(), class = "text-center text-muted mb-1"),
              plotly::plotlyOutput("comp_plotly", height = act_h(), width = act_w())),
          div(h6(ref_name() %||% "", class = "text-center text-muted mb-1"),
              plotly::plotlyOutput("comp_plotly_ref", height = ref_h(), width = ref_w()))
        )
      }
    } else {
      plotOutput("comp_plot",
                 height = paste0(input$plot_height, "px"),
                 width = paste0(input$plot_width, "px"))
    }
  })

  # 組成プロットを1データセットから構築
  # fill_var/clusters_sel を引数化（active と reference で別々のクラスター選択に対応）
  build_comp <- function(obj, interactive = FALSE,
                         fill_var = input$comp_cluster,
                         clusters_sel = input$comp_clusters_sel,
                         x_var = input$comp_x) {
    pt <- plot_theme()
    meta <- obj@meta.data

    if (is.null(fill_var) || is.null(x_var) || !all(c(fill_var, x_var) %in% names(meta))) {
      return(empty_panel(t("ref_missing")))
    }
    facet_vars <- input$comp_facets
    facet_vars <- facet_vars[facet_vars %in% names(meta)]
    x_cols <- x_var                              # X軸（複数選択可）
    x_label <- paste(x_cols, collapse = "_")
    key_cols <- unique(c(fill_var, x_cols, facet_vars))

    # --- 割合を計算 (X×facet 内で proportion を算出) ---
    df <- meta[key_cols]
    for (cc in key_cols) df[[cc]] <- as.character(df[[cc]])
    # 複数変数のX軸は paste して1列にまとめる
    x_use <- if (length(x_cols) == 1) x_cols else "__x_combined"
    if (length(x_cols) > 1) df[[x_use]] <- do.call(paste, c(df[x_cols], sep = "_"))
    group_vars <- unique(c(x_use, facet_vars))
    # 描画クラスターに限定（選択範囲内で割合を再計算 → 例: T細胞内での割合）
    sel_clusters <- selected_clusters_v(clusters_sel, obj, fill_var)
    df <- df[df[[fill_var]] %in% sel_clusters, , drop = FALSE]
    if (nrow(df) == 0) return(empty_panel(t("ref_missing")))
    agg_cols <- unique(c(fill_var, x_use, facet_vars))
    counts <- aggregate(list(n = rep(1L, nrow(df))), by = df[agg_cols], FUN = sum)
    grp_key <- do.call(paste, c(counts[group_vars], sep = "\r"))
    totals <- tapply(counts$n, grp_key, sum)
    counts$proportion <- counts$n / as.numeric(totals[grp_key])

    # --- クラスター(塗り分け)とX軸を系統順 factor に ---
    levs <- cluster_level_order(unique(counts[[fill_var]]))
    # plotly はスタックの先頭(=factorの最初)を下に積むため、interactive 時は
    # factor を反転し凡例を反転して戻す → バーの上端と凡例の先頭を一致させる
    fill_levels <- if (interactive) rev(levs) else levs
    counts[[fill_var]] <- factor(counts[[fill_var]], levels = fill_levels)
    counts[[x_use]] <- factor(counts[[x_use]],
                              levels = cluster_level_order(unique(counts[[x_use]])))

    # --- 配色 ---
    if (identical(input$comp_color_scheme, "manual")) {
      pal <- manual_colors
      if (length(levs) > length(pal)) {
        pal <- grDevices::colorRampPalette(manual_colors)(length(levs))
      }
      fill_colors <- setNames(pal[seq_along(levs)], levs)
    } else {
      fill_colors <- generate_cluster_colors(levs)
    }
    # 念のため未割当のレベルを補完
    missing <- setdiff(levs, names(fill_colors))
    if (length(missing) > 0) fill_colors[missing] <- "grey60"

    # ホバー用テキスト（クラスター名・X値・割合）
    counts$.hover <- paste0(
      fill_var, ": ", as.character(counts[[fill_var]]),
      "\n", x_label, ": ", as.character(counts[[x_use]]),
      "\n", t("comp_yaxis"), ": ", sprintf("%.1f%%", counts$proportion * 100),
      " (n=", counts$n, ")"
    )

    bar_aes <- if (interactive) aes(text = .hover) else aes()
    p <- ggplot(counts, aes(x = .data[[x_use]], y = proportion,
                            fill = .data[[fill_var]])) +
      geom_bar(mapping = bar_aes, stat = "identity") +
      scale_fill_manual(values = fill_colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = t("comp_title"), x = x_label, y = t("comp_yaxis"), fill = fill_var) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        text = element_text(color = pt$fg),
        axis.text = element_text(color = pt$fg2),
        legend.text = element_text(color = pt$fg),
        legend.title = element_text(color = pt$fg),
        strip.text = element_text(color = pt$fg),
        strip.background = element_rect(fill = pt$bg, color = pt$fg2),
        plot.title = element_text(size = 16, face = "bold", color = pt$accent)
      )

    # --- ファセット ---
    # interactive(plotly)時は facet_grid（ggh4x::facet_nested は ggplotly 非対応）
    if (length(facet_vars) > 0) {
      fct <- stats::as.formula(paste("~", paste(facet_vars, collapse = " + ")))
      if (!interactive && requireNamespace("ggh4x", quietly = TRUE)) {
        p <- p + ggh4x::facet_nested(fct, scales = "free_x",
                                     space = "free_x", nest_line = TRUE)
      } else {
        p <- p + facet_grid(fct, scales = "free_x", space = "free_x")
      }
    }

    p
  }

  # 「描画」ボタンを押したときだけ計算（アクティブ + リファレンス）
  comp_plot_obj <- eventReactive(input$comp_run, {
    req(seurat_obj(), input$comp_cluster, input$comp_x)
    rb <- ref_obj()
    pa <- build_comp(seurat_obj())
    pr <- if (!is.null(rb)) {
      build_comp(rb, fill_var = ref_tab_cluster("comp"),
                 clusters_sel = input$comp_clusters_sel_ref,
                 x_var = ref_comp_x())
    } else NULL
    combine_gg(pa, pr, active_name(), ref_name())
  }, ignoreInit = TRUE)

  # 静的版（plotly が無い場合のフォールバック）
  output$comp_plot <- renderPlot({ comp_plot_obj() }, bg = "transparent")

  # インタラクティブ版（バーにホバーでクラスター名・割合を表示）
  # リファレンス比較時は2つの独立した plotly を左右に並べる（subplot は facet と相性が悪い）
  if (requireNamespace("plotly", quietly = TRUE)) {
    comp_gg_titled <- function(obj, title = NULL,
                               fill_var = input$comp_cluster,
                               clusters_sel = input$comp_clusters_sel,
                               x_var = input$comp_x) {
      p <- build_comp(obj, interactive = TRUE, fill_var = fill_var,
                      clusters_sel = clusters_sel, x_var = x_var)
      if (!is.null(title)) p <- p + ggtitle(title)
      gp <- plotly::ggplotly(p, tooltip = "text")
      # スタックは factor 反転で上端=系統先頭。凡例も系統先頭が上に来るよう反転。
      plotly::layout(gp, legend = list(traceorder = "reversed"))
    }
    comp_plotly_obj <- eventReactive(input$comp_run, {
      req(seurat_obj(), input$comp_cluster, input$comp_x)
      title <- if (is.null(ref_obj())) NULL else active_name()
      comp_gg_titled(seurat_obj(), title)
    }, ignoreInit = TRUE)
    comp_plotly_ref_obj <- eventReactive(input$comp_run, {
      rb <- ref_obj()
      if (is.null(rb)) return(NULL)
      comp_gg_titled(rb, ref_name(),
                     fill_var = ref_tab_cluster("comp"),
                     clusters_sel = input$comp_clusters_sel_ref,
                     x_var = ref_comp_x())
    }, ignoreInit = TRUE)

    output$comp_plotly     <- plotly::renderPlotly({ comp_plotly_obj() })
    output$comp_plotly_ref <- plotly::renderPlotly({ req(comp_plotly_ref_obj()) })
  }

  # ==========================================================================
  # マーカー設定UI（Heatmap / Dot plot 共通）
  # ==========================================================================
  # prefix で input ID を分け、tab ごとに状態を保持する
  # feature_mode = TRUE のとき「セット / セット+カスタム / 自分で選択のみ」を選べる
  marker_settings_card <- function(prefix, extra = NULL, feature_mode = FALSE) {
    col_types <- meta_col_types()
    cat_cols <- col_types$cat

    cur_clu <- isolate(input[[paste0(prefix, "_cluster")]])
    clu_sel <- if (!is.null(cur_clu) && cur_clu %in% cat_cols) {
      cur_clu
    } else if ("seurat_clusters" %in% cat_cols) {
      "seurat_clusters"
    } else {
      cat_cols[1]
    }
    cur_set <- isolate(input[[paste0(prefix, "_set")]])
    set_sel <- if (!is.null(cur_set) && length(cur_set) > 0 && all(cur_set %in% names(marker_sets))) {
      cur_set
    } else {
      names(marker_sets)[1]
    }

    id_mode   <- paste0(prefix, "_feature_mode")
    id_set    <- paste0(prefix, "_set")
    id_custom <- paste0(prefix, "_custom_genes")
    gene_choices <- if (data_loaded()) sort(rownames(seurat_obj())) else character(0)

    if (feature_mode) {
      mode_sel <- isolate(input[[id_mode]]) %||% "set"
      gene_block <- tagList(
        fluidRow(
          column(12,
            radioButtons(id_mode, t("mk_feature_mode"),
              choices = c(
                setNames("set",        t("mk_feat_set")),
                setNames("set_custom", t("mk_feat_set_custom")),
                setNames("custom",     t("mk_feat_custom"))
              ),
              selected = mode_sel, inline = TRUE)
          )
        ),
        conditionalPanel(
          sprintf("input.%s == 'set' || input.%s == 'set_custom'", id_mode, id_mode),
          fluidRow(column(12,
            selectInput(id_set, t("marker_set"),
                        choices = names(marker_sets), selected = set_sel, multiple = TRUE)
          )),
          # セット内の遺伝子（不要なものを削除可能）
          uiOutput(paste0(prefix, "_set_genes_ui"))
        ),
        conditionalPanel(
          sprintf("input.%s == 'custom' || input.%s == 'set_custom'", id_mode, id_mode),
          fluidRow(column(12,
            selectizeInput(id_custom, t("mk_custom_genes"),
                           choices = gene_choices,
                           selected = isolate(input[[id_custom]]),
                           multiple = TRUE,
                           options = list(placeholder = t("gene_placeholder"),
                                          maxOptions = 1000))
          ))
        )
      )
      top_row <- fluidRow(column(12,
        selectInput(paste0(prefix, "_cluster"), t("marker_cluster"),
                    choices = cat_cols, selected = clu_sel)
      ))
    } else {
      gene_block <- NULL
      top_row <- fluidRow(
        column(6,
          selectInput(paste0(prefix, "_cluster"), t("marker_cluster"),
                      choices = cat_cols, selected = clu_sel)
        ),
        column(6,
          selectInput(id_set, t("marker_set"),
                      choices = names(marker_sets), selected = set_sel, multiple = TRUE)
        )
      )
    }

    div(
      class = "card mb-3",
      div(
        class = "card-body",
        h6(t("marker_settings"), class = "card-title text-primary"),
        top_row,
        gene_block,
        # 大分類で絞り込み → 描画するクラスターの選択（2段階）
        coarse_block(prefix, cat_cols),
        uiOutput(paste0(prefix, "_clusters_ui")),
        # リファレンス用のクラスター選択（比較時のみ）
        ref_cluster_controls_ui(prefix),
        extra
      )
    )
  }

  # --- マーカー描画用の (feature, group) を解決（セット / セット+カスタム / カスタムのみ） ---
  resolve_marker_set_df <- function(prefix) {
    mode <- input[[paste0(prefix, "_feature_mode")]] %||% "set"
    setn <- input[[paste0(prefix, "_set")]]   # 複数セット可
    parts <- list()
    if (mode %in% c("set", "set_custom") && length(setn) > 0) {
      sets_present <- setn[setn %in% names(marker_sets)]
      if (length(sets_present) > 0) {
        sdf <- do.call(rbind, lapply(sets_present, function(s) marker_set_to_df(marker_sets[[s]])))
        # ユーザーがセット内から削除した遺伝子を反映（kept が NULL のときは全て）
        kept <- input[[paste0(prefix, "_set_genes")]]
        if (!is.null(kept)) sdf <- sdf[sdf$feature %in% kept, , drop = FALSE]
        parts[[length(parts) + 1]] <- sdf
      }
    }
    if (mode %in% c("custom", "set_custom")) {
      cg <- input[[paste0(prefix, "_custom_genes")]]
      if (length(cg) > 0) {
        parts[[length(parts) + 1]] <- data.frame(
          feature = cg, group = t("mk_custom_group"), stringsAsFactors = FALSE)
      }
    }
    if (length(parts) == 0) {
      return(data.frame(feature = character(0), group = character(0),
                        stringsAsFactors = FALSE))
    }
    df <- do.call(rbind, parts)
    df[!duplicated(df$feature), , drop = FALSE]   # セットの groupを優先
  }

  # --- マーカー描画用の遺伝子のみ ---
  resolve_marker_genes <- function(prefix) {
    unique(resolve_marker_set_df(prefix)$feature)
  }

  # --- セット内の遺伝子の編集UI（選んだセットの遺伝子を削除可能） ---
  # セットを変えると新しいセットの全遺伝子で再初期化される。
  set_genes_select_ui <- function(prefix) {
    setn <- input[[paste0(prefix, "_set")]]
    setn <- setn[setn %in% names(marker_sets)]
    if (length(setn) == 0) return(NULL)
    all_genes <- unique(do.call(c, lapply(setn, function(s) marker_set_to_df(marker_sets[[s]])$feature)))
    id_sel <- paste0(prefix, "_set_genes")
    prev <- isolate(input[[id_sel]])
    sel <- if (!is.null(prev) && length(prev) > 0 && all(prev %in% all_genes)) prev else all_genes
    selectizeInput(id_sel, t("mk_set_genes"), choices = all_genes, selected = sel,
                   multiple = TRUE,
                   options = list(plugins = list("remove_button"), maxOptions = 2000))
  }
  output$mk_set_genes_ui  <- renderUI({ set_genes_select_ui("mk") })

  # --- 描画クラスター選択UI（クラスター変数に追従、デフォルトは全選択） ---
  # 大分類(coarse)列で絞ったあとの対象クラスター（ネスト2段階選択の1段目）
  clusters_in_coarse <- function(obj, cluster_var, coarse_var, coarse_sel) {
    levs <- cluster_level_order(obj@meta.data[[cluster_var]])
    if (is.null(coarse_var) || identical(coarse_var, "__none__") ||
        !(coarse_var %in% names(obj@meta.data)) ||
        is.null(coarse_sel) || length(coarse_sel) == 0) return(levs)
    cl <- as.character(obj@meta.data[[cluster_var]])
    co <- as.character(obj@meta.data[[coarse_var]])
    keep <- unique(cl[co %in% coarse_sel])
    levs[levs %in% keep]
  }

  # coarse 列の選択 + 値選択UI（suffix="_ref" で reference 用）
  coarse_block <- function(prefix, cat_cols, suffix = "") {
    id_var <- paste0(prefix, "_coarse_var", suffix)
    tagList(
      selectInput(id_var, t("coarse_var"),
                  choices = c(setNames("__none__", t("ref_none")), cat_cols),
                  selected = isolate(input[[id_var]]) %||% "__none__"),
      uiOutput(paste0(prefix, "_coarse", suffix, "_ui"))
    )
  }
  # coarse 値の選択UI（coarse_var に追従）
  coarse_values_ui <- function(prefix, obj, suffix = "") {
    if (is.null(obj)) return(NULL)
    cvar <- input[[paste0(prefix, "_coarse_var", suffix)]]
    if (is.null(cvar) || identical(cvar, "__none__") || !(cvar %in% names(obj@meta.data))) return(NULL)
    levs <- cluster_level_order(obj@meta.data[[cvar]])
    id_sel <- paste0(prefix, "_coarse_sel", suffix)
    prev <- isolate(input[[id_sel]])
    sel <- if (!is.null(prev) && length(prev) > 0 && all(prev %in% levs)) prev else levs
    selectizeInput(id_sel, t("coarse_values"), choices = levs, selected = sel,
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  }

  cluster_select_ui <- function(prefix) {
    obj <- seurat_obj()
    if (is.null(obj)) return(NULL)
    cvar <- input[[paste0(prefix, "_cluster")]]
    if (is.null(cvar) || !(cvar %in% names(obj@meta.data))) return(NULL)
    levs <- clusters_in_coarse(obj, cvar,
                               input[[paste0(prefix, "_coarse_var")]],
                               input[[paste0(prefix, "_coarse_sel")]])
    prev <- isolate(input[[paste0(prefix, "_clusters_sel")]])
    sel <- if (!is.null(prev) && length(prev) > 0 && all(prev %in% levs)) prev else levs
    selectizeInput(
      paste0(prefix, "_clusters_sel"), t("draw_clusters"),
      choices = levs, selected = sel, multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  }

  # --- 選択中の描画クラスター（未指定なら全クラスター） ---
  selected_clusters <- function(prefix, obj, cluster_var) {
    selected_clusters_v(input[[paste0(prefix, "_clusters_sel")]], obj, cluster_var)
  }
  # sel ベクトルを直接受け取る版（リファレンス用）
  selected_clusters_v <- function(sel, obj, cluster_var) {
    levs <- cluster_level_order(obj@meta.data[[cluster_var]])
    if (is.null(sel) || length(sel) == 0) return(levs)
    levs[levs %in% sel]   # 並び順は cluster_level_order を維持
  }

  # --- リファレンス用のクラスター変数（タブごと） ---
  ref_tab_cluster <- function(prefix) {
    cols <- ref_cat_cols()
    if (length(cols) == 0) return(NULL)
    sel <- input[[paste0(prefix, "_cluster_ref")]]
    if (!is.null(sel) && sel %in% cols) return(sel)
    act <- input[[paste0(prefix, "_cluster")]]
    if (!is.null(act) && act %in% cols) return(act)
    if ("seurat_clusters" %in% cols) return("seurat_clusters")
    cols[1]
  }

  # --- リファレンス用のクラスター変数＋描画クラスター選択UI（比較時のみ） ---
  ref_cluster_controls_ui <- function(prefix) {
    rb <- ref_obj()
    if (is.null(rb)) return(NULL)
    cols <- ref_cat_cols()
    if (length(cols) == 0) return(NULL)
    cur <- isolate(input[[paste0(prefix, "_cluster_ref")]])
    sel <- if (!is.null(cur) && cur %in% cols) {
      cur
    } else {
      act <- isolate(input[[paste0(prefix, "_cluster")]])
      if (!is.null(act) && act %in% cols) act
      else if ("seurat_clusters" %in% cols) "seurat_clusters" else cols[1]
    }
    tagList(
      hr(),
      h6(t("ref_settings_label"), class = "text-muted mb-2"),
      selectInput(paste0(prefix, "_cluster_ref"), t("ref_marker_cluster"),
                  choices = cols, selected = sel),
      coarse_block(prefix, cols, suffix = "_ref"),
      uiOutput(paste0(prefix, "_clusters_ref_ui"))
    )
  }
  # リファレンスの描画クラスター サブ選択UI（ref_tab_cluster + coarse に追従）
  ref_cluster_select_ui <- function(prefix) {
    rb <- ref_obj()
    if (is.null(rb)) return(NULL)
    cvar <- ref_tab_cluster(prefix)
    if (is.null(cvar) || !(cvar %in% names(rb@meta.data))) return(NULL)
    levs <- clusters_in_coarse(rb, cvar,
                               input[[paste0(prefix, "_coarse_var_ref")]],
                               input[[paste0(prefix, "_coarse_sel_ref")]])
    prev <- isolate(input[[paste0(prefix, "_clusters_sel_ref")]])
    sel <- if (!is.null(prev) && length(prev) > 0 && all(prev %in% levs)) prev else levs
    selectizeInput(
      paste0(prefix, "_clusters_sel_ref"), t("ref_draw_clusters"),
      choices = levs, selected = sel, multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  }

  # --- マーカーセットの平均発現マトリクス (gene x cluster) を計算 ---
  # genes はセット順を維持。clusters で描画クラスターを限定・並び替え。
  marker_avg_matrix <- function(obj, cluster_var, genes, clusters = NULL) {
    mat <- get_expr_matrix(obj, genes)        # genes x cells（存在する遺伝子のみ）
    if (is.null(mat)) return(NULL)
    # 列(細胞)順に合わせてラベルを取得（位置ではなく細胞名で対応付け）
    labels <- as.character(obj@meta.data[colnames(mat), cluster_var])
    keep <- !is.na(labels)
    mat <- mat[, keep, drop = FALSE]
    labels <- labels[keep]
    # 描画クラスターの順序（指定があればその順、なければ系統順）
    cl_levels <- if (!is.null(clusters)) clusters else cluster_level_order(unique(labels))
    cl_levels <- cl_levels[cl_levels %in% labels]
    if (length(cl_levels) == 0) return(NULL)
    avg <- sapply(cl_levels, function(cl) {
      Matrix::rowMeans(mat[, labels == cl, drop = FALSE])
    })
    if (is.null(dim(avg))) avg <- matrix(avg, nrow = nrow(mat),
                                         dimnames = list(rownames(mat), cl_levels))
    avg  # genes x clusters
  }

  # 描画クラスター選択UI（mk = Heatmap/Dot 共通、comp = Composition）
  output$mk_clusters_ui   <- renderUI({ cluster_select_ui("mk") })
  output$comp_clusters_ui <- renderUI({ cluster_select_ui("comp") })
  # リファレンス側のサブ選択
  output$mk_clusters_ref_ui   <- renderUI({ ref_cluster_select_ui("mk") })
  output$comp_clusters_ref_ui <- renderUI({ ref_cluster_select_ui("comp") })
  # coarse(大分類)の値選択UI（active / reference）
  output$mk_coarse_ui   <- renderUI({ coarse_values_ui("mk", seurat_obj()) })
  output$comp_coarse_ui <- renderUI({ coarse_values_ui("comp", seurat_obj()) })
  output$mk_coarse_ref_ui   <- renderUI({ coarse_values_ui("mk", ref_obj(), "_ref") })
  output$comp_coarse_ref_ui <- renderUI({ coarse_values_ui("comp", ref_obj(), "_ref") })

  # 描画実行ボタン（自動描画せず、押したときだけ描画）
  plot_run_btn <- function(id) {
    div(class = "d-grid mt-3",
        actionButton(id, t("plot_run"), class = "btn-primary",
                     icon = icon("play")))
  }
  # 未実行時に表示するヒント
  run_hint_ui <- function() {
    div(class = "text-center text-muted py-5", h5(t("plot_run_hint")))
  }

  # ==========================================================================
  # Heatmap
  # ==========================================================================
  # Heatmap / Dot を1つの大タブにまとめる。上部に共通のマーカー・クラスター
  # 設定（遺伝子・クラスター選択を両者で共有）、その下に Heatmap / Dot の内部タブ。
  output$marker_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    col_types <- meta_col_types()
    if (length(col_types$cat) == 0) {
      return(div(class = "text-center text-muted py-4", h5(t("comp_no_cat"))))
    }
    tagList(
      marker_settings_card("mk", NULL, feature_mode = TRUE),
      navset_card_tab(
        id = "marker_subtab",
        nav_panel(
          title = "\U0001F525 Heatmap", value = "heatmap",
          card_body(class = "p-2", uiOutput("heatmap_panel_ui"))
        ),
        nav_panel(
          title = "\U0001F535 Dot Plot", value = "dotplot",
          card_body(class = "p-2", uiOutput("dotplot_panel_ui"))
        ),
        nav_panel(
          title = "\U0001FAB5 Sub-cluster", value = "subcluster",
          card_body(class = "p-2", uiOutput("subc_panel_ui"))
        )
      )
    )
  })

  output$heatmap_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    col_types <- meta_col_types()
    if (length(col_types$cat) == 0) {
      return(div(class = "text-center text-muted py-4", h5(t("comp_no_cat"))))
    }

    tagList(
      div(class = "card mb-3", div(class = "card-body",
        fluidRow(
          column(4, checkboxInput("hm_scale", t("hm_scale"),
                                  value = isolate(input$hm_scale) %||% TRUE)),
          column(4, checkboxInput("hm_cluster_rows", t("hm_cluster_rows"),
                                  value = isolate(input$hm_cluster_rows) %||% TRUE)),
          column(4, checkboxInput("hm_cluster_cols", t("hm_cluster_cols"),
                                  value = isolate(input$hm_cluster_cols) %||% TRUE))
        ),
        # リファレンスがある時のみ: 1枚に統合してクラスター対応を見やすく
        if (!is.null(ref_obj())) {
          checkboxInput("hm_combine", t("hm_combine"),
                        value = isolate(input$hm_combine) %||% FALSE)
        },
        plot_run_btn("hm_run")
      )),
      uiOutput("heatmap_plot_ui")
    )
  })

  output$heatmap_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    if ((input$hm_run %||% 0) == 0) return(run_hint_ui())
    # 統合モードは1枚を全幅で表示、それ以外は左右に並べる
    if (isTRUE(input$hm_combine) && !is.null(ref_obj())) {
      plotOutput("heatmap_plot", height = act_h(), width = act_w())
    } else {
      side_by_side_ui("heatmap_plot", "heatmap_plot_ref", !is.null(ref_obj()))
    }
  })

  # 「描画」ボタンを押したときだけ計算（自動描画しない）
  # 1データセットの平均発現マトリクス（Z-score・クリップ済み）を構築
  build_hm_mat <- function(obj, cluster_var = input$mk_cluster,
                           clusters_sel = input$mk_clusters_sel) {
    genes <- resolve_marker_genes("mk")
    if (length(genes) == 0) return(NULL)
    if (is.null(cluster_var) || !(cluster_var %in% names(obj@meta.data))) return(NULL)
    clusters <- selected_clusters_v(clusters_sel, obj, cluster_var)
    avg <- marker_avg_matrix(obj, cluster_var, genes, clusters)
    if (is.null(avg) || nrow(avg) == 0) return(NULL)
    mat <- avg
    # Z-score はクラスター(列)が2つ以上ないと計算できない（1列は全てNaNになる）
    if (isTRUE(input$hm_scale) && ncol(mat) >= 2) {
      # NOTE: server scope の t() は翻訳ヘルパーなので転置は base::t() を使う
      mat <- base::t(scale(base::t(mat)))   # 遺伝子(行)ごとに Z-score
      mat <- mat[stats::complete.cases(mat), , drop = FALSE]
      mat[mat > 2] <- 2
      mat[mat < -2] <- -2
    }
    if (nrow(mat) == 0) return(NULL)
    mat
  }

  # 純 ggplot(geom_tile) でヒートマップを構築。
  # pheatmap は grid に直接描画してしまい（RStudio のプロットペインに出る等）
  # Shiny のデバイス制御と相性が悪いため、ggplot で描く。
  # mat: genes x clusters（Z-score 済み）。行=クラスター, 列=遺伝子 で表示。
  # row_source: クラスター名 -> データセット名 の名前付きベクトル。指定すると
  # 左端にデータセット別のアノテーションカラー列を追加する（統合ヒートマップ用）。
  build_hm_ggplot <- function(mat, main, cluster_rows, cluster_cols, row_source = NULL) {
    pt <- plot_theme()
    m <- base::t(mat)                      # clusters x genes
    row_order <- rownames(m)               # clusters
    col_order <- colnames(m)               # genes
    if (isTRUE(cluster_rows) && nrow(m) > 2) {
      row_order <- rownames(m)[stats::hclust(stats::dist(m))$order]
    }
    if (isTRUE(cluster_cols) && ncol(m) > 2) {
      col_order <- colnames(m)[stats::hclust(stats::dist(base::t(m)))$order]
    }

    has_ann <- !is.null(row_source) && requireNamespace("ggnewscale", quietly = TRUE)
    ann_lab <- "▮"                    # アノテーション列のラベル（細い縦棒記号）
    gene_levels <- if (has_ann) c(ann_lab, col_order) else col_order

    df <- data.frame(
      cluster = factor(rep(rownames(m), times = ncol(m)), levels = row_order),
      gene    = factor(rep(colnames(m), each = nrow(m)),  levels = gene_levels),
      value   = as.vector(m),
      stringsAsFactors = FALSE
    )
    lim <- max(abs(df$value), na.rm = TRUE)
    p <- ggplot() +
      geom_tile(data = df, aes(x = gene, y = cluster, fill = value),
                color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#0072B5FF", mid = "white", high = "#BC3C29FF",
                           midpoint = 0, limits = c(-lim, lim), name = NULL)

    if (has_ann) {
      srcs <- unique(unname(row_source))
      src_pal <- setNames(c("#4C78A8", "#F58518", "#54A24B", "#E45756")[seq_along(srcs)], srcs)
      ann <- data.frame(
        cluster = factor(rownames(m), levels = row_order),
        gene    = factor(ann_lab, levels = gene_levels),
        src     = factor(unname(row_source[rownames(m)]), levels = srcs),
        stringsAsFactors = FALSE
      )
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_tile(data = ann, aes(x = gene, y = cluster, fill = src),
                  color = "white", linewidth = 0.3) +
        scale_fill_manual(values = src_pal, name = "dataset", drop = FALSE) +
        # x の並び順を明示（指定しないと2層目の "▮" が右端に回るため左端に固定）
        scale_x_discrete(limits = gene_levels)
    }

    p +
      labs(title = main, x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        text = element_text(color = pt$fg),
        axis.text = element_text(color = pt$fg2),
        legend.text = element_text(color = pt$fg),
        legend.title = element_text(color = pt$fg),
        plot.title = element_text(size = 14, face = "bold", color = pt$accent)
      )
  }

  # アクティブとリファレンスのクラスター平均を共通遺伝子で1枚に統合。
  # 各データセット内で遺伝子ごとに Z-score 化（プラットフォーム差を吸収）してから結合。
  combine_hm <- function(aa, ar, na, nr, scale_within) {
    if (is.null(aa) || is.null(ar)) return(NULL)
    common <- intersect(rownames(aa), rownames(ar))
    if (length(common) < 2) return(NULL)
    aa <- aa[common, , drop = FALSE]; ar <- ar[common, , drop = FALSE]
    if (scale_within) {
      sc <- function(m) {
        if (ncol(m) >= 2) { z <- base::t(scale(base::t(m))); z[!is.finite(z)] <- 0; z } else m
      }
      aa <- sc(aa); ar <- sc(ar)
    }
    colnames(aa) <- paste0(na, " | ", colnames(aa))
    colnames(ar) <- paste0(nr, " | ", colnames(ar))
    mat <- cbind(aa, ar)
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    if (scale_within) { mat[mat > 2] <- 2; mat[mat < -2] <- -2 }
    mat
  }

  heatmap_spec <- eventReactive(input$hm_run, {
    req(seurat_obj(), input$mk_cluster)
    rb <- ref_obj()
    combine <- isTRUE(input$hm_combine) && !is.null(rb)
    if (combine) {
      genes <- resolve_marker_genes("mk")
      aa <- marker_avg_matrix(seurat_obj(), input$mk_cluster, genes,
                              selected_clusters("mk", seurat_obj(), input$mk_cluster))
      rcv <- ref_tab_cluster("mk")
      ar <- if (!is.null(rcv)) {
        marker_avg_matrix(rb, rcv, genes,
                          selected_clusters_v(input$mk_clusters_sel_ref, rb, rcv))
      } else NULL
      combined <- combine_hm(aa, ar, active_name(), ref_name() %||% "ref",
                             isTRUE(input$hm_scale))
      validate(need(!is.null(combined) && nrow(combined) > 0 && ncol(combined) >= 2,
                    t("hm_no_genes")))
      # 各列(クラスター)のデータセット名を "名前 | クラスター" 接頭辞から取得
      src <- sub(" \\| .*$", "", colnames(combined))
      list(combine = TRUE, combined = combined,
           source = setNames(src, colnames(combined)),
           cluster_cols = isTRUE(input$hm_cluster_cols))
    } else {
      ma <- build_hm_mat(seurat_obj())
      validate(need(!is.null(ma) && nrow(ma) > 0, t("hm_no_genes")))
      list(combine = FALSE, active = ma,
           ref = if (!is.null(rb)) {
             build_hm_mat(rb, cluster_var = ref_tab_cluster("mk"),
                          clusters_sel = input$mk_clusters_sel_ref)
           } else NULL,
           ref_present = !is.null(rb),
           names = c(active_name(), ref_name() %||% ""),
           cluster_rows = isTRUE(input$hm_cluster_rows),
           cluster_cols = isTRUE(input$hm_cluster_cols))
    }
  }, ignoreInit = TRUE)

  output$heatmap_plot <- renderPlot({
    spec <- heatmap_spec()
    if (isTRUE(spec$combine)) {
      # 統合モード: 行(クラスター)を必ず階層クラスタリングして似たもの同士を隣接。
      # 左端にデータセット別アノテーションカラーを表示。
      build_hm_ggplot(spec$combined, "", TRUE, spec$cluster_cols, row_source = spec$source)
    } else {
      build_hm_ggplot(spec$active, spec$names[1], spec$cluster_rows, spec$cluster_cols)
    }
  }, bg = "transparent")

  output$heatmap_plot_ref <- renderPlot({
    spec <- heatmap_spec()
    req(isFALSE(spec$combine), spec$ref_present)
    if (is.null(spec$ref)) {
      empty_panel(t("ref_missing"))
    } else {
      build_hm_ggplot(spec$ref, spec$names[2], spec$cluster_rows, spec$cluster_cols)
    }
  }, bg = "transparent")

  # ==========================================================================
  # Dot plot
  # ==========================================================================
  output$dotplot_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    col_types <- meta_col_types()
    if (length(col_types$cat) == 0) {
      return(div(class = "text-center text-muted py-4", h5(t("comp_no_cat"))))
    }

    tagList(
      div(class = "card mb-3", div(class = "card-body",
        fluidRow(
          column(6, sliderInput("dot_scale", t("dot_scale"), min = 2, max = 12,
                                value = isolate(input$dot_scale) %||% 6, step = 1)),
          column(6, checkboxInput("dot_facet", t("dot_facet"),
                                  value = isolate(input$dot_facet) %||% TRUE))
        ),
        if (!is.null(ref_obj())) {
          checkboxInput("dot_combine", t("dot_combine"),
                        value = isolate(input$dot_combine) %||% FALSE)
        },
        plot_run_btn("dot_run")
      )),
      uiOutput("dotplot_plot_ui")
    )
  })

  output$dotplot_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    if ((input$dot_run %||% 0) == 0) return(run_hint_ui())
    if (isTRUE(input$dot_combine) && !is.null(ref_obj())) {
      # 統合モード: 1枚の統合ドットプロット + 対応表 + CSVダウンロード
      tagList(
        plotOutput("dot_combined_plot", height = act_h(), width = act_w()),
        div(class = "mt-3",
          h6(class = "text-primary", t("corr_title"), " ",
             bslib::tooltip(
               tags$span(icon("circle-question"), style = "cursor: help;"),
               t("corr_help"), placement = "right")),
          downloadButton("dot_corr_dl", t("corr_download"), class = "btn-outline-success btn-sm mb-2"),
          DTOutput("dot_corr_table")
        )
      )
    } else if (requireNamespace("plotly", quietly = TRUE)) {
      side_by_side_ui("dot_plotly", "dot_plotly_ref", !is.null(ref_obj()), plotly = TRUE)
    } else {
      side_by_side_ui("dotplot_plot", "dotplot_plot_ref", !is.null(ref_obj()))
    }
  })

  # 「描画」ボタンを押したときだけ計算
  build_dot <- function(obj, cluster_var = input$mk_cluster,
                        clusters_sel = input$mk_clusters_sel,
                        interactive = FALSE) {
    pt <- plot_theme()
    if (is.null(cluster_var) || !(cluster_var %in% names(obj@meta.data))) return(empty_panel(t("ref_missing")))

    set_df <- resolve_marker_set_df("mk")
    genes_all <- unique(set_df$feature)
    validate(need(length(genes_all) > 0, t("hm_no_genes")))
    mat <- get_expr_matrix(obj, genes_all)     # genes x cells
    if (is.null(mat) || nrow(mat) == 0) return(empty_panel(t("ref_missing")))

    # 列(細胞)順に合わせてラベルを取得（位置ではなく細胞名で対応付け）
    labels <- as.character(obj@meta.data[colnames(mat), cluster_var])
    keep <- !is.na(labels)
    mat <- mat[, keep, drop = FALSE]
    labels <- labels[keep]
    # 描画クラスターを選択・系統順に並べ、対象細胞だけに限定
    cl_levels <- selected_clusters_v(clusters_sel, obj, cluster_var)
    cl_levels <- cl_levels[cl_levels %in% labels]
    if (length(cl_levels) == 0) return(empty_panel(t("ref_missing")))
    sel_idx <- labels %in% cl_levels
    mat <- mat[, sel_idx, drop = FALSE]
    labels <- labels[sel_idx]

    # クラスターごとに発現割合(pct.exp)と平均発現(avg.exp)を算出
    rows <- lapply(cl_levels, function(cl) {
      idx <- labels == cl
      sub <- mat[, idx, drop = FALSE]
      data.frame(
        feature = rownames(mat),
        id = cl,
        pct.exp = Matrix::rowSums(sub > 0) / sum(idx) * 100,
        avg.exp = Matrix::rowMeans(sub),
        stringsAsFactors = FALSE
      )
    })
    dot <- do.call(rbind, rows)

    # 遺伝子ごとに avg.exp を Z-score 化 (Seurat DotPlot と同様)
    dot$avg.exp.scaled <- ave(dot$avg.exp, dot$feature, FUN = function(x) {
      s <- as.vector(scale(x))
      s[is.na(s)] <- 0
      s
    })

    # グループ情報を付与し、セット順で x 軸を並べる
    grp_map <- set_df[!duplicated(set_df$feature), ]
    dot$group <- grp_map$group[match(dot$feature, grp_map$feature)]
    gene_order <- set_df$feature[set_df$feature %in% rownames(mat)]
    gene_order <- gene_order[!duplicated(gene_order)]
    dot$feature <- factor(dot$feature, levels = gene_order)
    dot$id <- factor(dot$id, levels = cl_levels)
    dot$group <- factor(dot$group, levels = unique(grp_map$group))

    # ホバー用テキスト（遺伝子・クラスター・発現割合・平均発現）
    dot$.hover <- paste0(
      "gene: ", as.character(dot$feature),
      "\n", cluster_var, ": ", as.character(dot$id),
      "\n", t("dot_pct"), ": ", sprintf("%.1f%%", dot$pct.exp),
      "\n", t("dot_avg"), ": ", sprintf("%.2f", dot$avg.exp.scaled)
    )
    pt_aes <- if (interactive) {
      aes(size = pct.exp, color = avg.exp.scaled, text = .hover)
    } else {
      aes(size = pct.exp, color = avg.exp.scaled)
    }

    p <- ggplot(dot, aes(x = feature, y = id)) +
      geom_point(mapping = pt_aes) +
      scale_radius(range = c(0, input$dot_scale), limits = c(0, 100)) +
      scale_color_gradient2(midpoint = 0, low = "#3C5488FF",
                            mid = "grey90", high = "#DC0000FF", space = "Lab") +
      labs(x = NULL, y = cluster_var,
           size = t("dot_pct"), color = t("dot_avg")) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        panel.grid = element_line(color = paste0(pt$fg2, "30")),
        text = element_text(color = pt$fg),
        axis.text = element_text(color = pt$fg2),
        legend.text = element_text(color = pt$fg),
        legend.title = element_text(color = pt$fg),
        strip.text = element_text(color = pt$fg, face = "bold"),
        panel.spacing = unit(0.6, "lines")
      )

    # グループでファセット（グループが複数あり、ファセット有効時のみ）
    if (isTRUE(input$dot_facet) && nlevels(dot$group) > 1) {
      p <- p + facet_grid(cols = vars(group), scales = "free_x",
                          space = "free_x") +
        theme(strip.background = element_blank())
    }

    p
  }

  # 「描画」ボタンを押したときだけ計算（アクティブ / リファレンス別出力）
  dot_active_obj <- eventReactive(input$dot_run, {
    req(seurat_obj(), input$mk_cluster)
    build_dot(seurat_obj())
  }, ignoreInit = TRUE)
  dot_ref_obj <- eventReactive(input$dot_run, {
    rb <- ref_obj()
    if (is.null(rb)) return(NULL)
    build_dot(rb, cluster_var = ref_tab_cluster("mk"),
              clusters_sel = input$mk_clusters_sel_ref)
  }, ignoreInit = TRUE)

  output$dotplot_plot     <- renderPlot({ dot_active_obj() }, bg = "transparent")
  output$dotplot_plot_ref <- renderPlot({ req(dot_ref_obj()) }, bg = "transparent")

  # --- 統合ドット用: クラスター×遺伝子の平均発現と発現割合 ---
  dot_avg_pct <- function(obj, cluster_var, clusters, genes) {
    mat <- get_expr_matrix(obj, genes)
    if (is.null(mat) || nrow(mat) == 0) return(NULL)
    labels <- as.character(obj@meta.data[colnames(mat), cluster_var])
    keep <- !is.na(labels) & labels %in% clusters
    mat <- mat[, keep, drop = FALSE]; labels <- labels[keep]
    cls <- clusters[clusters %in% labels]
    if (length(cls) == 0) return(NULL)
    avg <- vapply(cls, function(cl) Matrix::rowMeans(mat[, labels == cl, drop = FALSE]),
                  numeric(nrow(mat)))
    pct <- vapply(cls, function(cl) Matrix::rowSums(mat[, labels == cl, drop = FALSE] > 0) /
                                    sum(labels == cl) * 100, numeric(nrow(mat)))
    dimnames(avg) <- list(rownames(mat), cls); dimnames(pct) <- list(rownames(mat), cls)
    list(avg = avg, pct = pct)
  }

  # --- A(avg,pct) と R(avg,pct) から 統合ドット用データ + 対応表 を計算（共通） ---
  combine_dot_corr <- function(A, R, na, nr) {
    if (is.null(A) || is.null(R)) return(NULL)
    shared <- intersect(rownames(A$avg), rownames(R$avg))
    if (length(shared) < 2) return(NULL)
    zsc <- function(m) { z <- base::t(scale(base::t(m))); z[!is.finite(z)] <- 0; z }
    zA <- zsc(A$avg[shared, , drop = FALSE]); zR <- zsc(R$avg[shared, , drop = FALSE])
    cmat <- suppressWarnings(stats::cor(zA, zR)); cmat[!is.finite(cmat)] <- 0
    corr <- do.call(rbind, lapply(rownames(cmat), function(ac) {
      v <- cmat[ac, ]; ord <- order(v, decreasing = TRUE)
      best_name <- colnames(cmat)[ord[1]]; c1 <- v[ord[1]]
      c2 <- if (length(ord) >= 2) v[ord[2]] else NA
      sec <- if (length(ord) >= 2) colnames(cmat)[ord[2]] else NA
      mg  <- if (!is.na(c2)) c1 - c2 else NA
      conf <- if (c1 >= 0.5 && (is.na(mg) || mg >= 0.05)) "high"
              else if (c1 >= 0.3) "medium" else "low"
      proposed <- if (!is.na(mg) && mg < 0.05 && !is.na(sec)) paste0(best_name, " / ", sec, " ?") else best_name
      data.frame(active = ac, proposed_name = proposed, confidence = conf,
                 best = best_name, cor = round(c1, 3), second = sec,
                 cor2 = if (!is.na(c2)) round(c2, 3) else NA,
                 margin = if (!is.na(mg)) round(mg, 3) else NA, stringsAsFactors = FALSE)
    }))
    corr <- corr[order(-corr$cor), ]
    ids <- c(paste0(na, " | ", colnames(zA)), paste0(nr, " | ", colnames(zR)))
    src <- c(rep(na, ncol(zA)), rep(nr, ncol(zR)))
    avg_comb <- cbind(zA, zR); pct_comb <- cbind(A$pct[shared, , drop = FALSE], R$pct[shared, , drop = FALSE])
    colnames(avg_comb) <- ids; colnames(pct_comb) <- ids
    gene_ord <- shared[stats::hclust(stats::dist(avg_comb))$order]
    clu_ord  <- ids[stats::hclust(stats::dist(base::t(avg_comb)))$order]
    long <- data.frame(
      feature = factor(rep(shared, times = ncol(avg_comb)), levels = gene_ord),
      id      = factor(rep(ids, each = length(shared)), levels = clu_ord),
      avg     = as.vector(avg_comb), pct = as.vector(pct_comb), stringsAsFactors = FALSE)
    list(long = long, corr = corr, src = setNames(src, ids), clu_ord = clu_ord)
  }

  # 統合ドット(左端にデータセットアノテーション色)を ggplot で描く（共通）
  combined_dot_ggplot <- function(cc, dot_scale) {
    pt <- plot_theme()
    ann_lab <- "▮"
    gene_levels <- c(ann_lab, levels(cc$long$feature))
    dots <- cc$long; dots$feature <- factor(as.character(dots$feature), levels = gene_levels)
    ann <- data.frame(
      feature = factor(ann_lab, levels = gene_levels),
      id      = factor(names(cc$src), levels = cc$clu_ord),
      src     = factor(unname(cc$src), levels = unique(unname(cc$src))), stringsAsFactors = FALSE)
    src_pal <- setNames(c("#4C78A8", "#F58518", "#54A24B", "#E45756")[seq_along(levels(ann$src))],
                        levels(ann$src))
    ggplot() +
      geom_tile(data = ann, aes(x = feature, y = id, fill = src), width = 0.95, height = 0.95) +
      scale_fill_manual(values = src_pal, name = "dataset", drop = FALSE) +
      geom_point(data = dots, aes(x = feature, y = id, size = pct, color = avg)) +
      scale_x_discrete(limits = gene_levels) +   # "▮" を左端に固定
      scale_radius(range = c(0, dot_scale), limits = c(0, 100)) +
      scale_color_gradient2(midpoint = 0, low = "#3C5488FF", mid = "grey90",
                            high = "#DC0000FF", space = "Lab") +
      labs(x = NULL, y = NULL, size = t("dot_pct"), color = t("dot_avg")) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid = element_line(color = paste0(pt$fg2, "30")),
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        text = element_text(color = pt$fg), axis.text = element_text(color = pt$fg2),
        legend.text = element_text(color = pt$fg), legend.title = element_text(color = pt$fg))
  }

  # 対応表 DT（共通）
  corr_datatable <- function(corr) {
    datatable(corr, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Blfrtip"),
              colnames = c(t("corr_active"), t("corr_proposed"), t("corr_conf"),
                           t("corr_best"), t("corr_cor"),
                           t("corr_second"), paste0(t("corr_cor"), "2"), t("corr_margin")))
  }

  # 統合ドット + 対応表（描画ボタン押下時に計算）
  dot_combined <- eventReactive(input$dot_run, {
    req(seurat_obj(), input$mk_cluster)
    rb <- ref_obj(); if (is.null(rb)) return(NULL)
    genes <- resolve_marker_genes("mk")
    A <- dot_avg_pct(seurat_obj(), input$mk_cluster,
                     selected_clusters("mk", seurat_obj(), input$mk_cluster), genes)
    rcv <- ref_tab_cluster("mk")
    R <- if (!is.null(rcv)) {
      dot_avg_pct(rb, rcv, selected_clusters_v(input$mk_clusters_sel_ref, rb, rcv), genes)
    } else NULL
    combine_dot_corr(A, R, active_name(), ref_name() %||% "ref")
  }, ignoreInit = TRUE)

  output$dot_combined_plot <- renderPlot({
    cc <- dot_combined(); req(cc)
    combined_dot_ggplot(cc, input$dot_scale %||% 6)
  }, bg = "transparent")

  output$dot_corr_table <- renderDT({ cc <- dot_combined(); req(cc); corr_datatable(cc$corr) })

  output$dot_corr_dl <- downloadHandler(
    filename = function() "active_vs_reference_correspondence.csv",
    content = function(file) {
      cc <- dot_combined()
      if (!is.null(cc)) utils::write.csv(cc$corr, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # Sub-cluster（反復サブクラスタリング）
  # ==========================================================================
  work_labels  <- reactiveVal(NULL)   # 細胞名 -> 現在の作業ラベル
  subc_parents <- reactiveVal(NULL)   # 直近に細分化して生成したサブクラスター名

  # ラベルベースの平均発現/発現割合（作業ラベルから集計）
  avg_pct_labels <- function(obj, labels_vec, clusters, genes) {
    mat <- get_expr_matrix(obj, genes)
    if (is.null(mat) || nrow(mat) == 0) return(NULL)
    lab <- as.character(labels_vec[colnames(mat)])
    keep <- !is.na(lab) & lab %in% clusters
    mat <- mat[, keep, drop = FALSE]; lab <- lab[keep]
    cls <- clusters[clusters %in% lab]
    if (!length(cls)) return(NULL)
    avg <- vapply(cls, function(cl) Matrix::rowMeans(mat[, lab == cl, drop = FALSE]), numeric(nrow(mat)))
    pct <- vapply(cls, function(cl) Matrix::rowSums(mat[, lab == cl, drop = FALSE] > 0) /
                                    sum(lab == cl) * 100, numeric(nrow(mat)))
    dimnames(avg) <- list(rownames(mat), cls); dimnames(pct) <- list(rownames(mat), cls)
    list(avg = avg, pct = pct)
  }

  # 指定細胞を k 個にサブクラスタリング（変動の大きい遺伝子 + kmeans）
  quick_subcluster <- function(obj, cells, k) {
    mat <- get_expr_matrix(obj, rownames(obj))
    if (is.null(mat)) return(NULL)
    cells <- intersect(cells, colnames(mat))
    if (length(cells) <= k) return(NULL)
    mat <- mat[, cells, drop = FALSE]
    n <- ncol(mat); rmean <- Matrix::rowMeans(mat)
    rv <- (Matrix::rowSums(mat^2) - n * rmean^2) / (n - 1)
    rv[is.na(rv)] <- 0
    topg <- names(sort(rv, decreasing = TRUE))[seq_len(min(1000, sum(rv > 0)))]
    if (length(topg) < 2) return(NULL)
    m <- as.matrix(mat[topg, , drop = FALSE])
    m <- base::t(scale(base::t(m))); m[!is.finite(m)] <- 0
    set.seed(1)
    km <- stats::kmeans(base::t(m), centers = k, nstart = 5, iter.max = 50)
    setNames(km$cluster, cells)
  }

  output$subc_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    if (is.null(ref_obj())) {
      return(div(class = "text-center text-warning py-4", h5(t("subc_no_ref"))))
    }
    leaf_n <- if (!is.null(work_labels())) length(unique(work_labels())) else 0
    tagList(
      div(class = "alert alert-secondary py-2 small mb-2", icon("circle-info"), " ", t("subc_help")),
      div(class = "card mb-3", div(class = "card-body",
        h6(t("subc_settings"), class = "card-title text-primary"),
        actionButton("subc_init", t("subc_init"), class = "btn-outline-secondary btn-sm mb-2",
                     icon = icon("rotate-left")),
        if (leaf_n > 0) div(class = "small text-muted mb-2", sprintf(t("subc_leaf"), leaf_n)),
        uiOutput("subc_target_ui"),
        fluidRow(
          column(6, numericInput("subc_k", t("subc_k"), value = isolate(input$subc_k) %||% 2,
                                 min = 2, max = 12, step = 1)),
          column(6, div(style = "margin-top: 24px;",
            actionButton("subc_run", t("subc_run"), class = "btn-primary w-100",
                         icon = icon("sitemap"))))
        ),
        downloadButton("subc_anno_dl", t("subc_dl_anno"), class = "btn-outline-success btn-sm mt-2")
      )),
      uiOutput("subc_results_ui")
    )
  })

  output$subc_target_ui <- renderUI({
    wl <- work_labels()
    if (is.null(wl)) return(div(class = "small text-warning", t("subc_need_init")))
    leaves <- sort(unique(wl))
    selectizeInput("subc_target", t("subc_target"), choices = leaves,
                   selected = isolate(input$subc_target), multiple = TRUE,
                   options = list(plugins = list("remove_button")))
  })

  observeEvent(input$subc_init, {
    req(seurat_obj(), input$mk_cluster)
    obj <- seurat_obj()
    if (!(input$mk_cluster %in% names(obj@meta.data))) return()
    lab <- as.character(obj@meta.data[[input$mk_cluster]])
    names(lab) <- rownames(obj@meta.data)
    # アクティブの選択クラスターに限定（未選択なら全て）
    sel <- selected_clusters("mk", obj, input$mk_cluster)
    lab <- lab[lab %in% sel]
    work_labels(lab)
    subc_parents(NULL)
  })

  observeEvent(input$subc_run, {
    req(seurat_obj())
    wl <- work_labels()
    if (is.null(wl)) { showNotification(t("subc_need_init"), type = "warning"); return() }
    targets <- input$subc_target; k <- input$subc_k %||% 2
    req(length(targets) > 0, k >= 2)
    showNotification(t("subc_running"), id = "subc", type = "message", duration = NULL)
    tryCatch({
      obj <- seurat_obj(); children <- character(0)
      for (tg in targets) {
        cells <- names(wl)[wl == tg]
        sub <- quick_subcluster(obj, cells, k)
        if (is.null(sub)) next
        childlab <- paste0(tg, ".", sub - 1)        # 入れ子命名: tg.0, tg.1, ...
        wl[names(sub)] <- childlab
        children <- c(children, unique(childlab))
      }
      work_labels(wl); subc_parents(sort(unique(children)))
      showNotification(sprintf(t("subc_done"), length(unique(children))),
                       id = "subc", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste(t("notify_error"), conditionMessage(e)), type = "error", id = "subc")
    })
  })

  # 生成したサブクラスターを同じリファレンスと統合して再計算
  subc_view <- reactive({
    parents <- subc_parents(); wl <- work_labels()
    if (is.null(parents) || length(parents) == 0 || is.null(wl)) return(NULL)
    rb <- ref_obj(); if (is.null(rb)) return(NULL)
    genes <- resolve_marker_genes("mk")
    A <- avg_pct_labels(seurat_obj(), wl, parents, genes)
    rcv <- ref_tab_cluster("mk")
    R <- if (!is.null(rcv)) {
      dot_avg_pct(rb, rcv, selected_clusters_v(input$mk_clusters_sel_ref, rb, rcv), genes)
    } else NULL
    cc <- combine_dot_corr(A, R, active_name(), ref_name() %||% "ref")
    if (is.null(cc)) return(NULL)
    hm <- combine_hm(A$avg, R$avg, active_name(), ref_name() %||% "ref", TRUE)
    list(cc = cc,
         hm = hm,
         hm_src = if (!is.null(hm)) setNames(sub(" \\| .*$", "", colnames(hm)), colnames(hm)) else NULL)
  })

  output$subc_results_ui <- renderUI({
    if (is.null(subc_view())) {
      return(div(class = "text-center text-muted py-4", h5(t("subc_placeholder"))))
    }
    tagList(
      h6("\U0001F525 Heatmap", class = "text-primary mt-2"),
      plotOutput("subc_heatmap", height = act_h(), width = act_w()),
      h6("\U0001F535 Dot plot", class = "text-primary mt-3"),
      plotOutput("subc_dot", height = act_h(), width = act_w()),
      div(class = "mt-3",
        h6(class = "text-primary", t("corr_title"), " ",
           bslib::tooltip(tags$span(icon("circle-question"), style = "cursor: help;"),
                          t("corr_help"), placement = "right")),
        downloadButton("subc_corr_dl", t("corr_download"), class = "btn-outline-success btn-sm mb-2"),
        DTOutput("subc_corr_table"))
    )
  })

  output$subc_heatmap <- renderPlot({
    v <- subc_view(); req(v, !is.null(v$hm))
    build_hm_ggplot(v$hm, "", TRUE, TRUE, row_source = v$hm_src)
  }, bg = "transparent")
  output$subc_dot <- renderPlot({
    v <- subc_view(); req(v)
    combined_dot_ggplot(v$cc, input$subc_dot_scale %||% 6)
  }, bg = "transparent")
  output$subc_corr_table <- renderDT({ v <- subc_view(); req(v); corr_datatable(v$cc$corr) })
  output$subc_corr_dl <- downloadHandler(
    filename = function() "subcluster_vs_reference_correspondence.csv",
    content = function(file) {
      v <- subc_view(); if (!is.null(v)) utils::write.csv(v$cc$corr, file, row.names = FALSE)
    }
  )
  output$subc_anno_dl <- downloadHandler(
    filename = function() "subcluster_annotation.csv",
    content = function(file) {
      wl <- work_labels()
      if (!is.null(wl)) utils::write.csv(
        data.frame(cell = names(wl), label = unname(wl), stringsAsFactors = FALSE),
        file, row.names = FALSE)
    }
  )

  # インタラクティブ版（ドットにホバーで発現割合などを表示）
  if (requireNamespace("plotly", quietly = TRUE)) {
    # 指定 % に対応する plotly マーカーサイズ（同じ size スケールを再現して取得）
    dot_size_px <- function(breaks, dot_scale) {
      d <- data.frame(x = seq_along(breaks), y = 1, v = breaks)
      g <- plotly::plotly_build(plotly::ggplotly(
        ggplot(d, aes(x, y)) + geom_point(aes(size = v)) +
          scale_radius(range = c(0, dot_scale), limits = c(0, 100))))
      g$x$data[[1]]$marker$size
    }
    # ggplotly は facet の space="free_x" を無視し全パネルを等幅にするため、
    # 各パネル幅をグループの遺伝子数に比例させ、ドット間隔をグループ間で揃える。
    adjust_facet_widths <- function(gp, p) {
      if (is.null(p$data) || is.null(p$data$group)) return(gp)
      xaxes <- grep("^xaxis", names(gp$x$layout), value = TRUE)
      if (length(xaxes) < 2) return(gp)
      gl <- levels(droplevels(factor(p$data$group)))
      counts <- vapply(gl, function(g) {
        length(unique(as.character(p$data$feature[as.character(p$data$group) == g])))
      }, integer(1))
      counts <- counts[counts > 0]
      if (length(counts) != length(xaxes)) return(gp)
      # 現在のドメイン開始位置で左→右順にソート
      starts_cur <- vapply(xaxes, function(nm) gp$x$layout[[nm]]$domain[1], numeric(1))
      xord <- xaxes[order(starts_cur)]
      n <- length(xord); gap <- 0.02
      w <- counts / sum(counts) * (1 - gap * (n - 1))
      new_starts <- c(0, utils::head(cumsum(w + gap), -1))
      doms <- Map(function(s, wi) c(s, s + wi), new_starts, w)
      centers <- vapply(doms, mean, numeric(1))
      for (i in seq_len(n)) gp$x$layout[[xord[i]]]$domain <- doms[[i]]
      # strip ラベル(annotation)を新パネル中心へ移動
      ann <- gp$x$layout$annotations
      if (!is.null(ann)) {
        for (j in seq_along(ann)) {
          idx <- match(ann[[j]]$text, names(counts))
          if (!is.na(idx) && identical(ann[[j]]$xref, "paper")) {
            gp$x$layout$annotations[[j]]$x <- centers[idx]
          }
        }
      }
      gp
    }
    # ggplotly はサイズ凡例を落とすため、ダミーの凡例マーカーを追加する
    dot_to_plotly <- function(p, dot_scale) {
      gp <- plotly::plotly_build(plotly::ggplotly(p, tooltip = "text"))
      gp <- adjust_facet_widths(gp, p)
      breaks <- c(25, 50, 75, 100)
      sizes <- tryCatch(dot_size_px(breaks, dot_scale), error = function(e) NULL)
      if (!is.null(sizes)) {
        for (i in seq_along(breaks)) {
          gp$x$data[[length(gp$x$data) + 1]] <- list(
            x = list(NA), y = list(NA), type = "scatter", mode = "markers",
            marker = list(size = sizes[i], color = "grey50", line = list(width = 0)),
            name = paste0(breaks[i], "%"), legendgroup = "sizeleg",
            showlegend = TRUE, hoverinfo = "skip")
        }
      }
      gp
    }
    dot_plotly_active <- eventReactive(input$dot_run, {
      req(seurat_obj(), input$mk_cluster)
      dot_to_plotly(build_dot(seurat_obj(), interactive = TRUE), input$dot_scale)
    }, ignoreInit = TRUE)
    dot_plotly_refobj <- eventReactive(input$dot_run, {
      rb <- ref_obj()
      if (is.null(rb)) return(NULL)
      dot_to_plotly(
        build_dot(rb, cluster_var = ref_tab_cluster("mk"),
                  clusters_sel = input$mk_clusters_sel_ref, interactive = TRUE),
        input$dot_scale)
    }, ignoreInit = TRUE)
    output$dot_plotly     <- plotly::renderPlotly({ dot_plotly_active() })
    output$dot_plotly_ref <- plotly::renderPlotly({ req(dot_plotly_refobj()) })
  }

  # ==========================================================================
  # DEG解析
  # ==========================================================================
  observeEvent(input$run_deg, {
    req(seurat_obj(), input$deg_group_var)

    showNotification(t("notify_deg_running"), type = "message", id = "deg_run")

    tryCatch({
      obj <- seurat_obj()

      # 解析対象遺伝子（NULL = 全遺伝子）。限定時は min.pct=0 で選択遺伝子を確実に表示
      features <- deg_features()
      if (!is.null(features) && length(features) == 0) {
        showNotification(t("notify_deg_no_features"), type = "error", id = "deg_run")
        return()
      }
      deg_min_pct <- if (is.null(features)) 0.1 else 0

      if (deg_var_is_numeric()) {
        # --- 数値列: Top X% vs Bottom X% ---
        req(input$deg_percentile)
        vals <- obj@meta.data[[input$deg_group_var]]
        pct <- input$deg_percentile / 100
        top_threshold <- quantile(vals, 1 - pct, na.rm = TRUE)
        bottom_threshold <- quantile(vals, pct, na.rm = TRUE)

        label_col <- paste0("__deg_group_", input$deg_group_var)
        obj@meta.data[[label_col]] <- ifelse(
          vals >= top_threshold, paste0("Top ", input$deg_percentile, "%"),
          ifelse(vals <= bottom_threshold, paste0("Bottom ", input$deg_percentile, "%"),
                 "Middle")
        )
        Idents(obj) <- obj@meta.data[[label_col]]

        ident1_label <- paste0("Top ", input$deg_percentile, "%")
        ident2_label <- paste0("Bottom ", input$deg_percentile, "%")
        deg_title_label(paste0(input$deg_group_var, ": ", ident1_label, " vs ", ident2_label))

        markers <- FindMarkers(obj,
                               ident.1 = ident1_label,
                               ident.2 = ident2_label,
                               features = features,
                               logfc.threshold = 0,
                               min.pct = deg_min_pct)
      } else {
        # --- カテゴリカル列 ---
        req(input$deg_ident1, input$deg_ident2)

        # Group 2 (コントロール群) は複数選択可。
        # 「その他全て」が含まれる、または未選択の場合は ident.2 = NULL（残り全て）。
        ident2_sel <- input$deg_ident2
        use_all_others <- "__ALL_OTHERS__" %in% ident2_sel || length(ident2_sel) == 0
        ident2_groups <- setdiff(ident2_sel, "__ALL_OTHERS__")

        if (!use_all_others) {
          # ident.1 がコントロール群に含まれていないか確認
          if (input$deg_ident1 %in% ident2_groups) {
            showNotification(t("notify_deg_same"), type = "error")
            return()
          }
        }

        Idents(obj) <- obj@meta.data[[input$deg_group_var]]
        ident2_val <- if (use_all_others) NULL else ident2_groups
        ident2_display <- if (use_all_others) "All others" else paste(ident2_groups, collapse = ", ")
        deg_title_label(paste0(input$deg_ident1, " vs ", ident2_display))

        markers <- FindMarkers(obj,
                               ident.1 = input$deg_ident1,
                               ident.2 = ident2_val,
                               features = features,
                               logfc.threshold = 0,
                               min.pct = deg_min_pct)
      }

      markers$gene <- rownames(markers)
      markers$neg_log10_pval <- -log10(markers$p_val_adj + 1e-300)
      markers$significance <- ifelse(
        abs(markers$avg_log2FC) >= input$deg_logfc & markers$p_val_adj < input$deg_pval,
        ifelse(markers$avg_log2FC > 0, "Up", "Down"),
        "NS"
      )

      deg_results(markers)

      n_up <- sum(markers$significance == "Up")
      n_down <- sum(markers$significance == "Down")
      showNotification(
        sprintf(t("notify_deg_done"),  n_up, n_down),
        type = "message", id = "deg_run", duration = 5
      )

    }, error = function(e) {
      showNotification(paste(t("notify_error"), e$message), type = "error", id = "deg_run")
    })
  })

  # --- DEG結果UI ---
  output$deg_results_ui <- renderUI({
    res <- deg_results()
    if (is.null(res)) {
      return(div(
        class = "text-center text-muted py-4",
        h5(t("placeholder_deg"))
      ))
    }

    # plotly があればインタラクティブ（ホバーで遺伝子名表示）、無ければ静的
    volcano_out <- if (requireNamespace("plotly", quietly = TRUE)) {
      plotly::plotlyOutput("volcano_plotly", height = "500px")
    } else {
      plotOutput("volcano_plot", height = "500px",
                 width = paste0(input$plot_width, "px"))
    }

    tagList(
      # Enrichr へ送るボタン（有意な Up / Down DEG）
      div(class = "mb-2 d-flex gap-2 align-items-center flex-wrap",
        actionButton("enrichr_up_btn", t("enrichr_up"),
                     class = "btn-outline-danger btn-sm", icon = icon("up-long")),
        actionButton("enrichr_down_btn", t("enrichr_down"),
                     class = "btn-outline-primary btn-sm", icon = icon("down-long")),
        uiOutput("enrichr_link_ui", inline = TRUE)
      ),
      # Volcano Plot
      div(class = "mb-3", volcano_out),
      # DEGテーブル
      div(
        DTOutput("deg_table")
      )
    )
  })

  # --- Enrichr 連携（有意な Up/Down DEG を Enrichr に送信しリンクを表示） ---
  enrichr_link <- reactiveVal(NULL)
  submit_enrichr <- function(direction) {
    res <- deg_results()
    if (is.null(res)) return()
    genes <- res$gene[res$significance == direction]
    if (length(genes) == 0) {
      showNotification(t("enrichr_none"), type = "warning")
      return()
    }
    showNotification(t("enrichr_opening"), id = "enrichr", type = "message")
    tryCatch({
      r <- httr::POST(
        "https://maayanlab.cloud/Enrichr/addList",
        body = list(list = paste(genes, collapse = "\n"),
                    description = paste0("shiny_", direction)),
        encode = "multipart")
      httr::stop_for_status(r)
      # Enrichr は content-type が text/html だが本文は JSON なので明示的にパース
      sid <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))$shortId
      if (is.null(sid) || !nzchar(sid)) stop("no shortId")
      url <- paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=", sid)
      enrichr_link(list(url = url, dir = direction, n = length(genes)))
      removeNotification("enrichr")
    }, error = function(e) {
      showNotification(paste0(t("enrichr_err"), conditionMessage(e)),
                       type = "error", id = "enrichr")
    })
  }
  observeEvent(input$enrichr_up_btn,   { submit_enrichr("Up") })
  observeEvent(input$enrichr_down_btn, { submit_enrichr("Down") })
  output$enrichr_link_ui <- renderUI({
    l <- enrichr_link()
    if (is.null(l)) return(NULL)
    tags$a(href = l$url, target = "_blank",
           class = "btn btn-success btn-sm",
           icon("up-right-from-square"),
           sprintf(t("enrichr_open"), l$dir, l$n))
  })

  # --- Volcano Plot ビルダー（静的/インタラクティブ共通） ---
  build_volcano <- function(interactive = FALSE) {
    res <- deg_results()
    pt <- plot_theme()
    colors <- c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#95a5a6")

    # ラベル用: top 10遺伝子
    top_genes <- res[res$significance != "NS", ]
    top_genes <- top_genes[order(top_genes$p_val_adj), ]
    top_genes <- head(top_genes, 10)

    # インタラクティブ時はホバー用テキスト（遺伝子名＋統計量）
    gp_aes <- if (interactive) {
      aes(text = paste0(gene,
                        "\navg_log2FC: ", round(avg_log2FC, 3),
                        "\n-log10(padj): ", round(neg_log10_pval, 2),
                        "\nstatus: ", significance))
    } else {
      aes()
    }

    p <- ggplot(res, aes(x = avg_log2FC, y = neg_log10_pval, color = significance)) +
      geom_point(mapping = gp_aes, alpha = 0.6, size = 1.5) +
      scale_color_manual(values = colors) +
      geom_vline(xintercept = c(-input$deg_logfc, input$deg_logfc),
                 linetype = "dashed", color = pt$fg2, linewidth = 0.5) +
      geom_hline(yintercept = -log10(input$deg_pval),
                 linetype = "dashed", color = pt$fg2, linewidth = 0.5) +
      labs(
        title = paste0("Volcano Plot: ", deg_title_label()),
        x = "avg_log2FC",
        y = "-log10(p_val_adj)",
        color = ""
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        panel.grid.major = element_line(color = paste0(pt$fg2, "40")),
        panel.grid.minor = element_blank(),
        text = element_text(color = pt$fg),
        axis.text = element_text(color = pt$fg2),
        legend.text = element_text(color = pt$fg),
        plot.title = element_text(size = 16, face = "bold", color = pt$accent)
      )

    # 静的のみ: トップ遺伝子ラベル（plotly はホバーで表示するため不要）
    if (!interactive && nrow(top_genes) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        color = pt$fg,
        size = 3.5,
        max.overlaps = 20,
        segment.color = pt$fg2
      )
    }
    p
  }

  # 静的版（plotly が無い場合のフォールバック）
  output$volcano_plot <- renderPlot({
    req(deg_results())
    build_volcano(interactive = FALSE)
  }, bg = "transparent")

  # インタラクティブ版（ドットにホバーで遺伝子名を表示）
  if (requireNamespace("plotly", quietly = TRUE)) {
    output$volcano_plotly <- plotly::renderPlotly({
      req(deg_results())
      p <- build_volcano(interactive = TRUE)
      gp <- plotly::ggplotly(p, tooltip = "text")
      bg <- plot_theme()$bg
      plotly::layout(gp, paper_bgcolor = bg, plot_bgcolor = bg)
    })
  }

  # --- DEGテーブル ---
  output$deg_table <- renderDT({
    req(deg_results())
    res <- deg_results()

    # 有意な遺伝子のみフィルタ（全体も閲覧可能にするが初期表示は有意なもの）
    display_cols <- c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "significance")
    df <- res[, display_cols, drop = FALSE]

    # 外部データベースリンクを追加
    df$Links <- paste0(
      '<a href="', sapply(df$gene, ncbi_gene_url),
      '" target="_blank" class="btn btn-outline-info btn-sm btn-xs me-1" title="NCBI Gene">',
      'NCBI</a>',
      '<a href="', sapply(df$gene, immunexut_url),
      '" target="_blank" class="btn btn-outline-success btn-sm btn-xs me-1" title="ImmuNexUT">',
      'ImmuNexUT</a>',
      '<a href="', sapply(df$gene, ampra2_url),
      '" target="_blank" class="btn btn-outline-warning btn-sm btn-xs" title="AMP RA2">',
      'AMP RA2</a>'
    )

    # 数値の丸め（p_val_adjは数値のまま渡し、JSで表示時に指数表記にする）
    df$avg_log2FC <- round(df$avg_log2FC, 4)
    df$pct.1 <- round(df$pct.1, 3)
    df$pct.2 <- round(df$pct.2, 3)

    datatable(
      df,
      escape = FALSE,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 20,
        order = list(list(4, "asc")),
        dom = "Blfrtip",
        scrollX = TRUE,
        columnDefs = list(
          list(
            targets = 4,
            render = DT::JS(
              "function(data, type, row, meta) {",
              "  if (type === 'display' || type === 'filter') {",
              "    if (data === null) return '';",
              "    return Number(data).toExponential(2);",
              "  }",
              "  return data;",
              "}"
            )
          )
        )
      ),
      colnames = c("Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "Status", "External DB")
    )
  })

  # ==========================================================================
  # GSEA (fgsea × MSigDB)
  # ==========================================================================
  gsea_results <- reactiveVal(NULL)

  # MSigDB コレクション選択肢（H, C2:..., C5:... など）
  gsea_collections <- reactive({
    if (!requireNamespace("msigdbr", quietly = TRUE)) return(NULL)
    colls <- tryCatch(msigdbr::msigdbr_collections(), error = function(e) NULL)
    if (is.null(colls)) return(NULL)
    labs <- ifelse(colls$gs_subcat == "", colls$gs_cat,
                   paste0(colls$gs_cat, ": ", colls$gs_subcat))
    setNames(paste(colls$gs_cat, colls$gs_subcat, sep = "|"), labs)
  })

  # GSEA は DEG タブと同じ比較設定で得た deg_results() をランキングに使う
  output$gsea_inner_ui <- renderUI({
    lang <- input$lang
    if (!requireNamespace("fgsea", quietly = TRUE) ||
        !requireNamespace("msigdbr", quietly = TRUE)) {
      return(div(class = "text-center text-warning py-4",
        h5("fgsea / msigdbr が必要です: install.packages('msigdbr'); BiocManager::install('fgsea')")))
    }
    colls <- gsea_collections()
    coll_sel <- isolate(input$gsea_collection)
    if (is.null(coll_sel) || !(coll_sel %in% colls)) coll_sel <- "H|"

    tagList(
      div(class = "alert alert-secondary py-2 small mb-2", icon("circle-info"),
          " ", t("gsea_uses_deg")),
      div(class = "card mb-3", div(class = "card-body",
        h6(t("gsea_settings"), class = "card-title text-primary"),
        fluidRow(
          column(4, selectInput("gsea_species", t("gsea_species"),
                                choices = c("Homo sapiens", "Mus musculus"),
                                selected = isolate(input$gsea_species) %||% "Homo sapiens")),
          column(8, selectInput("gsea_collection", t("gsea_collection"),
                                choices = colls, selected = coll_sel))
        ),
        fluidRow(
          column(4, numericInput("gsea_minsize", t("gsea_minsize"),
                                value = isolate(input$gsea_minsize) %||% 5, min = 1, step = 1)),
          column(4, numericInput("gsea_maxsize", t("gsea_maxsize"),
                                value = isolate(input$gsea_maxsize) %||% 500, min = 10, step = 10)),
          column(4, div(style = "margin-top: 24px;",
            actionButton("gsea_run", t("gsea_run"), class = "btn-primary w-100",
                         icon = icon("dna"))))
        )
      )),
      uiOutput("gsea_results_ui")
    )
  })

  observeEvent(input$gsea_run, {
    res <- deg_results()
    if (is.null(res)) {
      showNotification(t("gsea_need_deg"), type = "warning", id = "gsea"); return()
    }
    showNotification(t("gsea_running"), id = "gsea", type = "message", duration = NULL)
    tryCatch({
      # DEG 結果の avg_log2FC で遺伝子をランキング
      ranks <- stats::setNames(res$avg_log2FC, res$gene)
      ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & !duplicated(names(ranks))]
      if (length(ranks) < 5) {
        showNotification(t("gsea_no_genes"), type = "error", id = "gsea"); return()
      }
      ranks <- sort(ranks, decreasing = TRUE)
      sel <- strsplit(input$gsea_collection, "\\|")[[1]]
      catg <- sel[1]; subcat <- if (length(sel) >= 2) sel[2] else ""
      msig <- msigdbr::msigdbr(species = input$gsea_species, category = catg,
                               subcategory = if (nzchar(subcat)) subcat else NULL)
      pathways <- split(msig$gene_symbol, msig$gs_name)
      set.seed(1)   # 再現性（fgsea の内部乱数を固定）
      fg <- fgsea::fgsea(pathways = pathways, stats = ranks,
                         minSize = max(1, input$gsea_minsize %||% 5),
                         maxSize = max(10, input$gsea_maxsize %||% 500),
                         nproc = n_workers())   # マルチスレッド（利用可能時）
      if (is.null(fg) || nrow(fg) == 0) {
        showNotification(t("gsea_no_sets"), type = "warning", id = "gsea"); return()
      }
      fg <- fg[order(fg$padj), ]
      gsea_results(fg)
      nsig <- sum(fg$padj < 0.05, na.rm = TRUE)
      showNotification(sprintf(t("gsea_done"), nrow(fg), nsig),
                       id = "gsea", type = "message", duration = 6)
    }, error = function(e) {
      showNotification(paste(t("notify_error"), conditionMessage(e)),
                       type = "error", id = "gsea")
    })
  })

  output$gsea_results_ui <- renderUI({
    if (is.null(gsea_results())) {
      return(div(class = "text-center text-muted py-4", h5(t("gsea_placeholder"))))
    }
    tagList(
      div(class = "mb-3", plotOutput("gsea_plot", height = "460px")),
      div(DTOutput("gsea_table"))
    )
  })

  output$gsea_plot <- renderPlot({
    fg <- gsea_results(); req(fg)
    pt <- plot_theme()
    top <- utils::head(fg[order(fg$padj), ], 20)
    top <- top[order(top$NES), ]
    top$pathway <- factor(as.character(top$pathway), levels = as.character(top$pathway))
    ggplot(top, aes(x = NES, y = pathway, fill = NES > 0)) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"), guide = "none") +
      labs(title = t("gsea_nes_title"), x = "NES", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        panel.grid.major.y = element_blank(),
        text = element_text(color = pt$fg),
        axis.text = element_text(color = pt$fg2),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 15, face = "bold", color = pt$accent)
      )
  }, bg = "transparent")

  output$gsea_table <- renderDT({
    fg <- gsea_results(); req(fg)
    df <- data.frame(
      pathway = as.character(fg$pathway),
      NES = round(fg$NES, 3),
      pval = fg$pval,
      padj = fg$padj,
      size = fg$size,
      leadingEdge = vapply(fg$leadingEdge,
                           function(x) paste(utils::head(x, 30), collapse = ", "),
                           character(1)),
      stringsAsFactors = FALSE
    )
    df$Link <- paste0(
      '<a href="https://www.gsea-msigdb.org/gsea/msigdb/cards/', df$pathway,
      '.html" target="_blank" class="btn btn-outline-info btn-sm btn-xs">MSigDB</a>')
    datatable(
      df, escape = FALSE, rownames = FALSE, filter = "top",
      options = list(
        pageLength = 20, order = list(list(3, "asc")), dom = "Blfrtip", scrollX = TRUE,
        columnDefs = list(list(
          targets = c(2, 3),
          render = DT::JS(
            "function(data, type, row, meta) {",
            "  if (type === 'display' || type === 'filter') {",
            "    if (data === null) return '';",
            "    return Number(data).toExponential(2);",
            "  }",
            "  return data;",
            "}")))),
      colnames = c("Pathway", "NES", "pval", "padj", "Size", "Leading edge", "DB")
    )
  })

  # ==========================================================================
  # Spatial（空間座標プロット）
  # ==========================================================================
  # プロット対象データセット（サイドバーで選択。デフォルトはアクティブ）
  spatial_obj <- reactive({
    objs <- loaded_objs()
    ds <- input$spatial_ds
    if (!is.null(ds) && !is.null(objs[[ds]])) return(objs[[ds]])
    seurat_obj()
  })

  output$spatial_ds_ui <- renderUI({
    lang <- input$lang
    objs <- loaded_objs()
    if (length(objs) == 0) return(NULL)
    nms <- names(objs)
    sel <- isolate(input$spatial_ds)
    if (is.null(sel) || !(sel %in% nms)) sel <- input$active_ds %||% nms[1]
    selectInput("spatial_ds", t("spatial_ds"), choices = nms, selected = sel)
  })

  # 選択データの空間座標
  spatial_xy <- reactive({
    obj <- spatial_obj(); req(obj)
    spatial_coords(obj)
  })
  # 選択データのセグメンテーション（あれば）
  spatial_seg <- reactive({
    obj <- spatial_obj(); req(obj)
    spatial_segmentation(obj)
  })

  output$spatial_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    xy <- spatial_xy()
    if (is.null(xy)) {
      return(div(class = "text-center text-warning py-4", h5(t("spatial_none"))))
    }
    obj <- spatial_obj()
    meta <- obj@meta.data
    cat_cols <- names(meta)[sapply(meta, function(x) {
      is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 50)
    })]
    color_choices <- names(meta)   # カテゴリ/数値どちらでも色分け可
    col_sel <- isolate(input$spatial_color)
    if (is.null(col_sel) || !(col_sel %in% color_choices)) {
      col_sel <- if ("seurat_clusters" %in% color_choices) "seurat_clusters"
                 else if (length(cat_cols) > 0) cat_cols[1] else color_choices[1]
    }
    # サンプル列の候補（複数サンプルを含む空間データを1サンプルに絞るため）
    samp_cands <- intersect(c("sample_id", "sample_id_full", "orig.ident", "PatientID",
                              "donor", "slide", "sample", "library", "fov", "Sample_Type"),
                            cat_cols)
    samp_cands <- c("(なし)" = "__none__", setNames(samp_cands, samp_cands))
    samp_col_sel <- isolate(input$spatial_sample_col)
    if (is.null(samp_col_sel) || !(samp_col_sel %in% samp_cands)) {
      samp_col_sel <- if (length(samp_cands) > 1) samp_cands[[2]] else "__none__"
    }
    genes <- sort(rownames(obj))
    gene_sel <- isolate(input$spatial_gene)
    if (is.null(gene_sel) || !(gene_sel %in% genes)) gene_sel <- genes[1]
    by_sel <- isolate(input$spatial_by) %||% "meta"
    tagList(
      div(class = "card mb-3", div(class = "card-body",
        h6(t("spatial_settings"), class = "card-title text-primary"),
        radioButtons("spatial_by", t("spatial_by"),
          choices = c(setNames("meta", t("spatial_by_meta")),
                      setNames("gene", t("spatial_by_gene"))),
          selected = by_sel, inline = TRUE),
        fluidRow(
          column(6,
            conditionalPanel("input.spatial_by == 'meta'",
              selectInput("spatial_color", t("spatial_color"),
                          choices = color_choices, selected = col_sel)),
            conditionalPanel("input.spatial_by == 'gene'",
              selectizeInput("spatial_gene", t("spatial_gene"),
                             choices = if (!is.null(gene_sel)) stats::setNames(gene_sel, gene_sel) else NULL,
                             selected = gene_sel,
                             options = list(placeholder = t("gene_placeholder"),
                                            maxOptions = 1000)))
          ),
          column(6, selectInput("spatial_sample_col", t("spatial_sample_col"),
                                choices = samp_cands, selected = samp_col_sel))
        ),
        uiOutput("spatial_sample_ui"),
        fluidRow(
          column(6, sliderInput("spatial_pt", t("spatial_pt"),
                                min = 0.1, max = 4, value = isolate(input$spatial_pt) %||% 0.8, step = 0.1)),
          column(6, div(style = "margin-top: 30px;",
            checkboxInput("spatial_flip", t("spatial_flip"),
                          value = isolate(input$spatial_flip) %||% TRUE)))
        ),
        # セグメンテーションがあれば表示切替（デフォルトON）
        if (!is.null(spatial_seg())) {
          checkboxInput("spatial_use_seg", t("spatial_seg"),
                        value = isolate(input$spatial_use_seg) %||% TRUE)
        },
        # 強調・近傍解析の対象クラスター（メタデータ・カテゴリ色分け時のみ）
        conditionalPanel("input.spatial_by == 'meta'",
          uiOutput("spatial_highlight_ui"),
          # 強調クラスター/その他クラスターの 不透明度・点サイズ 個別指定
          fluidRow(
            column(3, sliderInput("spatial_hl_alpha", t("spatial_hl_alpha"),
                                  min = 0.1, max = 1, value = isolate(input$spatial_hl_alpha) %||% 1, step = 0.05)),
            column(3, sliderInput("spatial_hl_size", t("spatial_hl_size"),
                                  min = 0.1, max = 4, value = isolate(input$spatial_hl_size) %||% 1.2, step = 0.1)),
            column(3, sliderInput("spatial_other_alpha", t("spatial_other_alpha"),
                                  min = 0, max = 1, value = isolate(input$spatial_other_alpha) %||% 0.3, step = 0.05)),
            column(3, sliderInput("spatial_other_size", t("spatial_other_size"),
                                  min = 0.1, max = 4, value = isolate(input$spatial_other_size) %||% 0.6, step = 0.1))
          )),
        # 近傍解析(Neighborhood/Co-occurrence/Enrichment)用のサンプル選択
        # （3タブで共有。Map以外のサブタブ表示時のみ）
        conditionalPanel("input.spatial_subtab && input.spatial_subtab != 'map'",
          uiOutput("spatial_nbr_sample_ui"))
      )),
      navset_card_tab(
        id = "spatial_subtab",
        nav_panel(title = "\U0001F5FA️ Map", value = "map",
                  card_body(class = "p-2", uiOutput("spatial_plot_ui"))),
        nav_panel(title = "\U0001F4CF Neighborhood", value = "nbr",
                  card_body(class = "p-2", uiOutput("spatial_nbr_ui"))),
        nav_panel(title = "\U0001F517 Co-occurrence", value = "co",
                  card_body(class = "p-2", uiOutput("spatial_co_ui"))),
        nav_panel(title = "\U0001F525 Enrichment", value = "ne",
                  card_body(class = "p-2", uiOutput("spatial_ne_ui")))
      )
    )
  })

  # 対象クラスター選択UI（現在の色分け変数のレベル＝カテゴリ時のみ）
  output$spatial_highlight_ui <- renderUI({
    obj <- spatial_obj(); req(obj)
    cv <- input$spatial_color
    if (identical(input$spatial_by, "gene") || is.null(cv) || !(cv %in% names(obj@meta.data))) return(NULL)
    v <- obj@meta.data[[cv]]
    if (!(is.factor(v) || is.character(v) || (is.numeric(v) && length(unique(v)) <= 50))) return(NULL)
    levs <- cluster_level_order(v)
    prev <- isolate(input$spatial_highlight)
    sel <- if (!is.null(prev)) prev[prev %in% levs] else NULL
    selectizeInput("spatial_highlight", t("spatial_highlight"), choices = levs, selected = sel,
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  })

  # サンプル選択UI（サンプル列に追従。複数選択可）
  output$spatial_sample_ui <- renderUI({
    obj <- spatial_obj(); req(obj)
    sc <- input$spatial_sample_col
    if (is.null(sc) || identical(sc, "__none__") || !(sc %in% names(obj@meta.data))) return(NULL)
    vals <- sort(unique(as.character(obj@meta.data[[sc]])))
    prev <- isolate(input$spatial_sample)
    sel <- if (!is.null(prev) && all(prev %in% vals) && length(prev) > 0) prev else vals[1]
    selectizeInput("spatial_sample", t("spatial_sample"), choices = vals, selected = sel,
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  })

  # プロット出力は常にDOMに置く（動的挿入による htmlwidget 初期化失敗を回避）。
  # 計算は「描画」ボタン押下時のみ（observeEvent → reactiveVal）。
  output$spatial_plot_ui <- renderUI({
    tagList(
      div(class = "mb-2",
        actionButton("spatial_map_run", t("plot_run"), class = "btn-primary btn-sm",
                     icon = icon("play")),
        span(class = "text-muted small ms-2", t("plot_run_hint"))),
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput("spatial_plotly", height = act_h(), width = act_w())
      } else {
        plotOutput("spatial_plot", height = act_h(), width = act_w())
      }
    )
  })

  # 「描画」ボタン押下時に、その時点の設定をスナップショットして保存。
  spatial_map_spec_v <- reactiveVal(NULL)
  observeEvent(input$spatial_map_run, {
    spatial_map_spec_v(list(
      pr = spatial_prep(),
      pt_sz = input$spatial_pt %||% 0.8,
      flip = isTRUE(input$spatial_flip),
      colorvar = spatial_colorvar(),
      highlight = input$spatial_highlight,
      hl_alpha = input$spatial_hl_alpha %||% 1,
      hl_size = input$spatial_hl_size %||% 1.2,
      other_alpha = input$spatial_other_alpha %||% 0.3,
      other_size = input$spatial_other_size %||% 0.6))
  }, ignoreInit = TRUE)
  spatial_map_spec <- reactive({ req(spatial_map_spec_v()) })

  # 色分けに使う名前（メタデータ変数名 or 遺伝子名）
  spatial_colorvar <- reactive({
    if (identical(input$spatial_by, "gene")) input$spatial_gene else input$spatial_color
  })

  # 描画用データの準備（点/ポリゴン・plotly/ggplot で共通）
  spatial_prep <- reactive({
    obj <- spatial_obj(); req(obj)
    xy <- spatial_xy(); req(xy)
    meta <- obj@meta.data
    by_gene <- identical(input$spatial_by, "gene")
    if (by_gene) req(input$spatial_gene) else req(input$spatial_color)
    df <- xy[xy$cell %in% rownames(meta), , drop = FALSE]
    if (by_gene) {
      # 選択遺伝子の発現量（data レイヤー）で色分け
      g <- input$spatial_gene
      em <- get_expr_matrix(obj, g)
      validate(need(!is.null(em), t("spatial_none")))
      ev <- as.numeric(em[1, ])
      names(ev) <- colnames(em)
      df <- df[df$cell %in% names(ev), , drop = FALSE]
      df$col <- ev[df$cell]
    } else {
      df$col <- meta[df$cell, input$spatial_color]
    }
    sc <- input$spatial_sample_col
    facet <- FALSE
    if (!is.null(sc) && !identical(sc, "__none__") && sc %in% names(meta)) {
      df$sample <- as.character(meta[df$cell, sc])
      if (!is.null(input$spatial_sample) && length(input$spatial_sample) > 0) {
        df <- df[df$sample %in% input$spatial_sample, , drop = FALSE]
      }
      facet <- length(unique(df$sample)) > 1
    }
    validate(need(nrow(df) > 0, t("spatial_none")))
    # 遺伝子発現は常に連続値、メタデータは型から判定
    is_cat <- !by_gene && (is.factor(df$col) || is.character(df$col) ||
              (is.numeric(df$col) && length(unique(df$col)) <= 50))
    if (is_cat) df$col <- factor(as.character(df$col), levels = cluster_level_order(df$col))
    seg <- spatial_seg()
    use_seg <- isTRUE(input$spatial_use_seg) && !is.null(seg)
    sg <- if (use_seg) seg[seg$cell %in% df$cell, , drop = FALSE] else NULL
    list(df = df, sg = sg, use_seg = use_seg, is_cat = is_cat, facet = facet)
  })

  # 色(hex/名前)に不透明度を付与して rgba 文字列にする
  hex_rgba <- function(col, a) {
    rgb <- grDevices::col2rgb(col)
    sprintf("rgba(%d,%d,%d,%.3f)", rgb[1, ], rgb[2, ], rgb[3, ], a)
  }

  # --- インタラクティブ版（点/ポリゴンにホバーでラベル表示）---
  if (requireNamespace("plotly", quietly = TRUE)) {
    # 1サンプル分の plotly を構築。hl= 強調クラスター集合, hl_alpha/size・other_alpha/size で
    # 透明度と点サイズを強調/非強調ごとに指定。
    spatial_one_plotly <- function(d, sg, use_seg, is_cat, flip, colorvar, disc_pal, showlegend,
                                   hl, hl_alpha, hl_size, other_alpha, other_size, base_size) {
      yf <- if (flip) function(v) -v else function(v) v
      d$hover <- paste0(colorvar, ": ", as.character(d$col), "<br>cell: ", d$cell)
      has_hl <- is_cat && !is.null(hl) && length(hl) > 0
      d$ishl <- if (has_hl) as.character(d$col) %in% hl else TRUE
      fig <- plotly::plot_ly()
      if (use_seg && !is.null(sg) && nrow(sg) > 0 && is_cat) {
        # カテゴリごとに1トレース（ポリゴンをNAで区切り fill='toself'）
        for (lev in levels(droplevels(d$col))) {
          cl_cells <- d$cell[d$col == lev]
          s <- sg[sg$cell %in% cl_cells, , drop = FALSE]
          if (!nrow(s)) next
          a <- if (!has_hl) 1 else if (lev %in% hl) hl_alpha else other_alpha
          if (a <= 0) next
          xs <- unlist(lapply(split(s$x, factor(s$cell, unique(s$cell))), function(z) c(z, NA)))
          ys <- unlist(lapply(split(s$y, factor(s$cell, unique(s$cell))), function(z) c(yf(z), NA)))
          fig <- plotly::add_trace(fig, x = xs, y = ys, type = "scatter", mode = "lines",
            fill = "toself", fillcolor = hex_rgba(unname(disc_pal[lev]), a),
            line = list(width = 0.4, color = "rgba(255,255,255,0.5)"),
            name = lev, legendgroup = lev, showlegend = showlegend, hoverinfo = "skip")
        }
        fig <- plotly::add_trace(fig, data = d, x = ~x, y = yf(d$y), type = "scattergl",
          mode = "markers", marker = list(size = 3, color = "rgba(0,0,0,0)"),
          text = ~hover, hoverinfo = "text", showlegend = FALSE)
      } else if (is_cat) {
        add_group <- function(dd, sz, a, sl) {
          if (nrow(dd) == 0) return(fig)
          dd$mc <- hex_rgba(unname(disc_pal[as.character(dd$col)]), a)
          plotly::add_trace(fig, data = dd, x = ~x, y = yf(dd$y), type = "scattergl",
            mode = "markers", marker = list(size = sz * 5, color = ~mc),
            text = ~hover, hoverinfo = "text", showlegend = FALSE)
        }
        if (has_hl) {
          # 非強調を下(背景)、強調を上(前面)に描く
          fig <- add_group(d[!d$ishl, , drop = FALSE], other_size, other_alpha, FALSE)
          fig <- add_group(d[d$ishl, , drop = FALSE], hl_size, hl_alpha, FALSE)
          # 凡例用の色付きダミー（強調クラスターのみ）
          if (showlegend) for (lev in hl[hl %in% levels(d$col)]) {
            fig <- plotly::add_trace(fig, x = c(NA), y = c(NA), type = "scattergl", mode = "markers",
              marker = list(size = 8, color = unname(disc_pal[lev])), name = lev,
              showlegend = TRUE, hoverinfo = "skip")
          }
        } else {
          fig <- add_group(d, base_size, 1, showlegend)
          if (showlegend) for (lev in levels(droplevels(d$col))) {
            fig <- plotly::add_trace(fig, x = c(NA), y = c(NA), type = "scattergl", mode = "markers",
              marker = list(size = 8, color = unname(disc_pal[lev])), name = lev,
              showlegend = TRUE, hoverinfo = "skip")
          }
        }
      } else {
        fig <- plotly::add_trace(fig, data = d, x = ~x, y = yf(d$y), type = "scattergl",
          mode = "markers",
          marker = list(size = base_size * 5, color = ~col, colorscale = "Viridis",
                        showscale = showlegend, colorbar = list(title = colorvar)),
          text = ~hover, hoverinfo = "text", showlegend = FALSE)
      }
      plotly::layout(fig,
        xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
        yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE,
                     scaleanchor = "x", scaleratio = 1))
    }

    output$spatial_plotly <- plotly::renderPlotly({
      spec <- spatial_map_spec(); pr <- spec$pr; pt <- plot_theme()
      df <- pr$df; flip <- spec$flip; colorvar <- spec$colorvar
      hl <- spec$highlight
      disc_pal <- NULL
      if (pr$is_cat) {
        lin <- lineage_colors_or_null(df$col)
        levs <- levels(df$col)
        disc_pal <- if (!is.null(lin)) lin[levs] else
          setNames(scales::hue_pal()(length(levs)), levs)
        # 強調指定時、その他クラスターは灰色に
        if (!is.null(hl) && length(hl) > 0) {
          disc_pal[!(names(disc_pal) %in% hl)] <- "#9E9E9E"
        }
      }
      showleg <- !(pr$is_cat && nlevels(df$col) > 30)
      samples <- if (pr$facet) sort(unique(df$sample)) else NA
      build1 <- function(s, first) {
        d <- if (is.na(s)) df else df[df$sample == s, , drop = FALSE]
        sg <- if (!is.null(pr$sg)) pr$sg[pr$sg$cell %in% d$cell, , drop = FALSE] else NULL
        f <- spatial_one_plotly(d, sg, pr$use_seg, pr$is_cat, flip, colorvar, disc_pal,
          showlegend = (first && showleg), hl = hl, hl_alpha = spec$hl_alpha,
          hl_size = spec$hl_size, other_alpha = spec$other_alpha, other_size = spec$other_size,
          base_size = spec$pt_sz)
        if (!is.na(s)) f <- plotly::layout(f, annotations = list(text = s, x = 0.5, y = 1.0,
          xref = "paper", yref = "paper", showarrow = FALSE, font = list(color = pt$fg)))
        f
      }
      fig <- if (pr$facet) {
        # サンプル数が多いときは正方形に近いグリッドに配置して見やすく
        figs <- lapply(seq_along(samples), function(i) build1(samples[i], i == 1))
        nr <- max(1, round(sqrt(length(samples))))
        plotly::subplot(figs, nrows = nr, margin = 0.015, titleX = FALSE, titleY = FALSE)
      } else build1(NA, TRUE)
      plotly::layout(fig, paper_bgcolor = pt$bg, plot_bgcolor = pt$bg,
                     font = list(color = pt$fg), legend = list(font = list(color = pt$fg)))
    })
  }

  output$spatial_plot <- renderPlot({
    spec <- spatial_map_spec(); pt <- plot_theme(); pr <- spec$pr
    df <- pr$df; sg <- pr$sg; use_seg <- pr$use_seg; is_cat <- pr$is_cat; facet <- pr$facet
    pt_sz <- spec$pt_sz
    if (use_seg && !is.null(sg)) {
      cmap <- setNames(df$col, df$cell); sg$col <- cmap[sg$cell]
      if (!is.null(df$sample)) sg$sample <- setNames(df$sample, df$cell)[sg$cell]
    }

    eff_seg <- use_seg && !is.null(sg)
    p <- ggplot()
    if (eff_seg) {
      p <- p +
        geom_polygon(data = sg, aes(x = x, y = y, group = cell, fill = col),
                     color = pt$bg, linewidth = 0.1) +
        labs(fill = spatial_colorvar())
    } else {
      p <- p +
        geom_point(data = df, aes(x = x, y = y, color = col),
                   size = pt_sz, shape = 16) +
        labs(color = spatial_colorvar())
    }
    p <- p + coord_fixed() + labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        plot.background = element_rect(fill = pt$bg, color = NA),
        panel.background = element_rect(fill = pt$bg, color = NA),
        text = element_text(color = pt$fg),
        legend.text = element_text(color = pt$fg),
        legend.title = element_text(color = pt$fg),
        strip.text = element_text(color = pt$fg)
      )
    if (isTRUE(input$spatial_flip)) p <- p + scale_y_reverse()

    # 配色（塗り=fill / 点=color を共通の関数で適用）
    apply_disc <- function(p, aes_fn_name) {
      lin <- lineage_colors_or_null(df$col)
      if (!is.null(lin)) {
        p <- p + (if (aes_fn_name == "fill") scale_fill_manual(values = lin)
                  else scale_color_manual(values = lin))
      }
      if (nlevels(df$col) > 30) p <- p + theme(legend.position = "none")
      else if (aes_fn_name == "color")
        p <- p + guides(color = guide_legend(override.aes = list(size = 3)))
      p
    }
    if (is_cat) {
      p <- apply_disc(p, if (eff_seg) "fill" else "color")
    } else {
      p <- p + (if (eff_seg) scale_fill_viridis_c() else scale_color_viridis_c())
    }
    if (facet) p <- p + facet_wrap(~ sample)
    p
  }, bg = "transparent")

  # ==========================================================================
  # Spatial: 近傍距離（対象クラスター → 他クラスターへの最近接距離 ECDF）
  # ==========================================================================
  output$spatial_nbr_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    tagList(
      div(class = "d-flex align-items-center gap-2 mb-2",
        actionButton("spatial_nbr_run", t("spatial_nbr_run"),
                     class = "btn-primary btn-sm", icon = icon("ruler")),
        bslib::tooltip(tags$span(icon("circle-question"), style = "cursor: help;"),
                       t("spatial_nbr_help"), placement = "right")),
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput("spatial_nbr_plot", height = act_h(), width = act_w())
      } else plotOutput("spatial_nbr_plot_static", height = act_h(), width = act_w())
    )
  })

  # 近傍解析専用のサンプル選択UI（サンプル列に追従。複数選択可）
  output$spatial_nbr_sample_ui <- renderUI({
    obj <- spatial_obj(); req(obj)
    sc <- input$spatial_sample_col
    if (is.null(sc) || identical(sc, "__none__") || !(sc %in% names(obj@meta.data))) return(NULL)
    vals <- sort(unique(as.character(obj@meta.data[[sc]])))
    prev <- isolate(input$spatial_nbr_sample)
    sel <- if (!is.null(prev) && all(prev %in% vals) && length(prev) > 0) prev else vals[1]
    selectizeInput("spatial_nbr_sample", t("spatial_sample"), choices = vals, selected = sel,
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  })

  # 近傍解析用データ（マップとは独立: 色分け変数=クラスター + 近傍解析用サンプル）
  spatial_nbr_df <- reactive({
    obj <- spatial_obj(); req(obj)
    xy <- spatial_xy(); req(xy)
    meta <- obj@meta.data
    cv <- input$spatial_color
    ok <- !identical(input$spatial_by, "gene") && !is.null(cv) && cv %in% names(meta)
    validate(need(ok, t("spatial_nbr_need")))
    v <- meta[[cv]]
    is_cat <- is.factor(v) || is.character(v) || (is.numeric(v) && length(unique(v)) <= 50)
    validate(need(is_cat, t("spatial_nbr_need")))
    df <- xy[xy$cell %in% rownames(meta), , drop = FALSE]
    df$col <- factor(as.character(meta[df$cell, cv]), levels = cluster_level_order(meta[[cv]]))
    sc <- input$spatial_sample_col
    if (!is.null(sc) && !identical(sc, "__none__") && sc %in% names(meta)) {
      df$sample <- as.character(meta[df$cell, sc])
      sel <- input$spatial_nbr_sample
      if (!is.null(sel) && length(sel) > 0) df <- df[df$sample %in% sel, , drop = FALSE]
    } else df$sample <- "all"
    validate(need(nrow(df) > 0, t("spatial_none")))
    df
  })

  # 対象細胞群から他クラスター細胞群への最近接距離（同一サンプル内）
  # 対象クラスター(A)の各細胞から クラスターB の最近接距離（サンプルごとのベクトル）
  nn_dist_per_sample <- function(df, A, B, cap = 3000) {
    lapply(unique(df$sample), function(smp) {
      a <- df[df$sample == smp & as.character(df$col) == A, c("x", "y"), drop = FALSE]
      b <- df[df$sample == smp & as.character(df$col) == B, c("x", "y"), drop = FALSE]
      if (nrow(a) < 3 || nrow(b) < 1) return(NULL)
      if (nrow(a) > cap) a <- a[sample(nrow(a), cap), , drop = FALSE]
      nn <- tryCatch(FNN::get.knnx(as.matrix(b), as.matrix(a), k = 1)$nn.dist[, 1],
                     error = function(e) NULL)
      if (is.null(nn)) nn <- apply(as.matrix(a), 1, function(p)
        sqrt(min(colSums((base::t(as.matrix(b)) - p)^2))))
      nn <- nn[is.finite(nn) & nn > 0]
      if (length(nn) < 3) NULL else nn
    })
  }

  # 近傍距離 ECDF（距離グリッド上で「最近接 ≤ x」割合、サンプルごと→平均±SD）
  spatial_nbr_obj <- eventReactive(input$spatial_nbr_run, {
    set.seed(1)   # 再現性（細胞の間引きを固定）
    df <- spatial_nbr_df()
    anchors <- input$spatial_highlight
    validate(need(length(anchors) > 0, t("spatial_nbr_need")))
    levs <- levels(df$col)
    recs <- list()
    for (A in anchors) for (B in setdiff(levs, A)) {
      d <- Filter(Negate(is.null), nn_dist_per_sample(df, A, B))
      if (!length(d)) next
      recs[[paste(A, "||", B)]] <- list(A = A, B = B, d = d,
                                        med = stats::median(unlist(d)))
    }
    validate(need(length(recs) > 0, t("spatial_nbr_need")))
    alld <- unlist(lapply(recs, function(r) unlist(r$d)))
    rng <- range(alld[alld > 0]); rng[1] <- max(rng[1], rng[2] / 1e4)
    grid <- exp(seq(log(rng[1]), log(rng[2]), length.out = 60))
    # アンカーごとに中央値が近い順 最大40クラスターに絞る
    byA <- split(recs, vapply(recs, function(r) r$A, ""))
    keepRecs <- unlist(lapply(byA, function(rs) {
      ord <- order(vapply(rs, function(r) r$med, 0))
      rs[ord[seq_len(min(40, length(rs)))]]
    }), recursive = FALSE)
    do.call(rbind, lapply(keepRecs, function(r) {
      fr <- vapply(r$d, function(dd) stats::ecdf(dd)(grid), numeric(length(grid)))
      if (is.null(dim(fr))) fr <- matrix(fr, ncol = 1)
      ym <- rowMeans(fr); ys <- if (ncol(fr) > 1) apply(fr, 1, stats::sd) else rep(0, nrow(fr))
      data.frame(anchor = r$A, cluster = r$B, x = grid, ymean = ym,
                 ylow = pmax(0, ym - ys), yhigh = pmin(1, ym + ys),
                 nsamp = ncol(fr), stringsAsFactors = FALSE)
    }))
  }, ignoreInit = TRUE)

  # 折れ線 + 信頼帯(リボン) の共通描画（co-occurrence でも再利用）
  spatial_curve_ggplot <- function(d, xlab, ylab, title_fmt, logx = TRUE, hline = NULL) {
    pt <- plot_theme()
    d$cluster <- factor(d$cluster, levels = cluster_level_order(unique(d$cluster)))
    d$tip <- paste0(d$cluster, "<br>", xlab, ": ", signif(d$x, 3),
                    "<br>", ylab, ": ", round(d$ymean, 3))
    p <- ggplot(d, aes(x = x, group = cluster))
    if (any(d$nsamp > 1)) p <- p +
      geom_ribbon(aes(ymin = ylow, ymax = yhigh, fill = cluster), alpha = 0.15, color = NA)
    p <- p + geom_line(aes(y = ymean, color = cluster, text = tip), linewidth = 0.6)
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                                             color = pt$fg2, linewidth = 0.4)
    if (logx) p <- p + scale_x_log10()
    p <- p + labs(x = xlab, y = ylab, color = NULL, fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(plot.background = element_rect(fill = pt$bg, color = NA),
            panel.background = element_rect(fill = pt$bg, color = NA),
            panel.grid.minor = element_blank(),
            text = element_text(color = pt$fg), axis.text = element_text(color = pt$fg2),
            legend.text = element_text(color = pt$fg),
            strip.text = element_text(color = pt$fg, face = "bold"))
    lin <- lineage_colors_or_null(levels(d$cluster))
    if (!is.null(lin)) p <- p + scale_color_manual(values = lin) + scale_fill_manual(values = lin)
    if (length(unique(d$anchor)) > 1) p <- p + facet_wrap(~ anchor)
    else p <- p + ggtitle(sprintf(title_fmt, unique(d$anchor)))
    p
  }

  spatial_nbr_ggplot <- reactive({
    spatial_curve_ggplot(spatial_nbr_obj(), t("spatial_nbr_x"), t("spatial_nbr_y"),
                         t("spatial_nbr_title"))
  })
  if (requireNamespace("plotly", quietly = TRUE)) {
    output$spatial_nbr_plot <- plotly::renderPlotly({
      plotly::ggplotly(spatial_nbr_ggplot(), tooltip = "text")
    })
  }
  output$spatial_nbr_plot_static <- renderPlot({ spatial_nbr_ggplot() }, bg = "transparent")

  # ==========================================================================
  # Spatial: Co-occurrence probability（squidpy 参考）
  # ==========================================================================
  output$spatial_co_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    tagList(
      div(class = "d-flex align-items-center gap-2 mb-2",
        actionButton("spatial_co_run", t("spatial_co_run"), class = "btn-primary btn-sm",
                     icon = icon("link")),
        bslib::tooltip(tags$span(icon("circle-question"), style = "cursor: help;"),
                       t("spatial_co_help"), placement = "right")),
      if (requireNamespace("plotly", quietly = TRUE))
        plotly::plotlyOutput("spatial_co_plot", height = act_h(), width = act_w())
      else plotOutput("spatial_co_plot_static", height = act_h(), width = act_w())
    )
  })

  # 1サンプル内: アンカーAから距離リングごとの co-occurrence 比（クラスター別）
  co_one_sample <- function(coords, lab, A, breaks, levs, n_tot, cap = 1500) {
    ai <- which(lab == A)
    if (length(ai) < 3) return(NULL)
    if (length(ai) > cap) ai <- sample(ai, cap)
    Rmax <- max(breaks)
    nn <- tryCatch(RANN::nn2(coords, coords[ai, , drop = FALSE], k = min(nrow(coords), 400),
                             searchtype = "radius", radius = Rmax), error = function(e) NULL)
    if (is.null(nn)) return(NULL)
    fi <- as.vector(nn$nn.idx); fd <- as.vector(nn$nn.dists)
    keep <- fi > 0 & fd > 0
    if (!any(keep)) return(NULL)
    ring <- findInterval(fd[keep], breaks, rightmost.closed = TRUE)
    nbl  <- lab[fi[keep]]
    ok <- ring >= 1 & ring <= (length(breaks) - 1)
    tab <- table(factor(ring[ok], levels = seq_len(length(breaks) - 1)),
                 factor(nbl[ok], levels = levs))
    tot <- rowSums(tab)
    pj <- vapply(levs, function(j) sum(lab == j) / n_tot, 0)   # 全体での割合
    ratio <- sweep(tab, 1, ifelse(tot == 0, NA, tot), "/")     # P(j|A,ring)
    ratio <- sweep(ratio, 2, ifelse(pj == 0, NA, pj), "/")     # / P(j)
    ratio  # ring x cluster
  }

  spatial_co_obj <- eventReactive(input$spatial_co_run, {
    set.seed(1)   # 再現性
    df <- spatial_nbr_df()
    anchors <- input$spatial_highlight
    validate(need(length(anchors) > 0, t("spatial_need_run")))
    levs <- levels(droplevels(df$col))
    # 距離リング: 全細胞の最近接距離中央値を基準に Rmax を決め、線形に分割
    mednn <- stats::median(vapply(unique(df$sample), function(smp) {
      cc <- as.matrix(df[df$sample == smp, c("x", "y")])
      if (nrow(cc) < 5) return(NA_real_)
      stats::median(FNN::get.knn(cc, k = 1)$nn.dist[, 1])
    }, 0), na.rm = TRUE)
    if (!is.finite(mednn) || mednn <= 0) mednn <- 10
    breaks <- seq(0, mednn * 40, length.out = 21)
    mids <- (head(breaks, -1) + breaks[-1]) / 2
    out <- list()
    for (A in anchors) {
      per <- lapply(unique(df$sample), function(smp) {
        sub <- df[df$sample == smp, , drop = FALSE]
        co_one_sample(as.matrix(sub[, c("x", "y")]), as.character(sub$col), A,
                      breaks, levs, nrow(sub))
      })
      per <- Filter(Negate(is.null), per)
      if (!length(per)) next
      for (j in levs[levs != A]) {
        vals <- vapply(per, function(m) as.numeric(m[, j]), numeric(length(mids)))
        if (is.null(dim(vals))) vals <- matrix(vals, ncol = 1)
        ym <- rowMeans(vals, na.rm = TRUE)
        ys <- if (ncol(vals) > 1) apply(vals, 1, stats::sd, na.rm = TRUE) else rep(0, nrow(vals))
        keep <- is.finite(ym)
        if (!any(keep)) next
        out[[paste(A, j)]] <- data.frame(anchor = A, cluster = j, x = mids[keep],
          ymean = ym[keep], ylow = pmax(0, (ym - ys)[keep]), yhigh = (ym + ys)[keep],
          nsamp = ncol(vals), stringsAsFactors = FALSE)
      }
    }
    validate(need(length(out) > 0, t("spatial_need_run")))
    do.call(rbind, out)
  }, ignoreInit = TRUE)

  spatial_co_ggplot <- reactive({
    spatial_curve_ggplot(spatial_co_obj(), t("spatial_co_x"), t("spatial_co_y"),
                         t("spatial_co_title"), logx = FALSE, hline = 1)
  })
  if (requireNamespace("plotly", quietly = TRUE)) {
    output$spatial_co_plot <- plotly::renderPlotly({
      plotly::ggplotly(spatial_co_ggplot(), tooltip = "text")
    })
  }
  output$spatial_co_plot_static <- renderPlot({ spatial_co_ggplot() }, bg = "transparent")

  # ==========================================================================
  # Spatial: Neighbors enrichment（squidpy 参考、置換 z-score）
  # ==========================================================================
  output$spatial_ne_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    tagList(
      fluidRow(
        column(3, numericInput("spatial_ne_k", t("spatial_ne_k"),
                               value = isolate(input$spatial_ne_k) %||% 6, min = 2, max = 30, step = 1)),
        column(3, numericInput("spatial_ne_perm", t("spatial_ne_perm"),
                               value = isolate(input$spatial_ne_perm) %||% 100, min = 20, max = 1000, step = 20)),
        column(3, div(style = "margin-top: 30px;",
          checkboxInput("spatial_ne_stars", t("spatial_ne_stars"),
                        value = isolate(input$spatial_ne_stars) %||% TRUE))),
        column(3, div(style = "margin-top: 24px;",
          div(class = "d-flex align-items-center gap-2",
            actionButton("spatial_ne_run", t("spatial_ne_run"), class = "btn-primary btn-sm",
                         icon = icon("project-diagram")),
            bslib::tooltip(tags$span(icon("circle-question"), style = "cursor: help;"),
                           t("spatial_ne_help"), placement = "right"))))
      ),
      if (requireNamespace("plotly", quietly = TRUE))
        plotly::plotlyOutput("spatial_ne_plot", height = act_h(), width = act_w())
      else plotOutput("spatial_ne_plot_static", height = act_h(), width = act_w())
    )
  })

  spatial_ne_obj <- eventReactive(input$spatial_ne_run, {
    set.seed(1)   # 再現性（間引き + 並べ替え検定を固定）
    df <- spatial_nbr_df()
    validate(need(nlevels(droplevels(df$col)) >= 2, t("spatial_need_run")))
    k <- max(2, input$spatial_ne_k %||% 6); nperm <- max(20, input$spatial_ne_perm %||% 100)
    levs <- levels(droplevels(df$col)); K <- length(levs)
    li <- setNames(seq_len(K), levs)
    e1 <- integer(0); e2 <- integer(0); labv <- integer(0); sampv <- integer(0)
    base <- 0L; sid <- 0L
    for (smp in unique(df$sample)) {
      sub <- df[df$sample == smp, , drop = FALSE]
      cc <- as.matrix(sub[, c("x", "y")]); n <- nrow(cc)
      if (n > 20000) { ix <- sample(n, 20000); cc <- cc[ix, ]; sub <- sub[ix, ]; n <- 20000 }
      if (n <= k) next
      sid <- sid + 1L
      knn <- FNN::get.knn(cc, k = k)$nn.index
      u <- rep(seq_len(n), times = k); v <- as.vector(knn)
      a <- pmin(u, v) + base; b <- pmax(u, v) + base
      e1 <- c(e1, a); e2 <- c(e2, b)
      labv <- c(labv, li[as.character(sub$col)]); sampv <- c(sampv, rep(sid, n))
      base <- base + n
    }
    validate(need(length(e1) > 0, t("spatial_need_run")))
    # 無向ユニークエッジ
    key <- paste(e1, e2); dup <- !duplicated(key); e1 <- e1[dup]; e2 <- e2[dup]
    count_mat <- function(lab) {
      a <- lab[e1]; b <- lab[e2]
      m <- matrix(0, K, K)
      tab <- table(factor(a, levels = seq_len(K)), factor(b, levels = seq_len(K)))
      m <- as.matrix(tab); m + base::t(m)   # 対称化（t は翻訳ヘルパーなので base::t）
    }
    obs <- count_mat(labv)
    sm <- matrix(0, K, K); sm2 <- matrix(0, K, K)
    sample_idx <- split(seq_along(labv), sampv)
    # 各並べ替えに固定シードを割り当て（コア数に依らず再現性を担保）
    perm_seeds <- sample.int(.Machine$integer.max, nperm)
    one_perm <- function(p) {
      lp <- labv
      for (idx in sample_idx) lp[idx] <- sample(labv[idx])   # サンプル内でシャッフル
      count_mat(lp)
    }
    fmt_sec <- function(s) if (s < 60) paste0(round(s), "s") else
      paste0(s %/% 60, "m", round(s %% 60), "s")
    nw <- n_workers()
    # 進捗表示のためチャンク分割（チャンク内は並列、チャンク間で進捗更新）
    nchunk <- min(nperm, max(10L, nw * 4L))
    chunks <- split(seq_len(nperm), cut(seq_len(nperm), nchunk, labels = FALSE))
    msg <- if (nw > 1L) paste0(t("spatial_ne_prog"), " (", nw, " threads)") else t("spatial_ne_prog")
    withProgress(message = msg, value = 0, {
      t0 <- Sys.time(); done <- 0L
      for (ch in chunks) {
        cms <- par_lapply_seeded(ch, one_perm, perm_seeds[ch])
        for (cm in cms) { sm <- sm + cm; sm2 <- sm2 + cm^2 }
        done <- done + length(ch)
        el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        eta <- el / done * (nperm - done)
        setProgress(done / nperm, detail = sprintf(t("spatial_eta"), done, nperm, fmt_sec(eta)))
      }
    })
    mu <- sm / nperm; sdv <- sqrt(pmax(0, sm2 / nperm - mu^2))
    z <- (obs - mu) / ifelse(sdv == 0, NA, sdv)
    dimnames(z) <- list(levs, levs)
    z
  }, ignoreInit = TRUE)

  spatial_ne_ggplot <- reactive({
    z <- spatial_ne_obj(); req(z)
    pt <- plot_theme()
    levs <- rownames(z); K <- length(levs)
    zc <- z; zc[!is.finite(zc)] <- 0
    # 行・列を階層クラスタリングして似たものを近くに（対称行列なので行=列順）
    ord <- if (K > 2) levs[stats::hclust(stats::dist(zc))$order] else levs
    # z-score を正規近似で両側p値に変換 → ユニークなペア(上三角)で BH 補正
    pmat <- 2 * stats::pnorm(-abs(z))
    padj <- matrix(NA_real_, K, K, dimnames = dimnames(z))
    ut <- upper.tri(pmat, diag = TRUE)
    padj[ut] <- stats::p.adjust(pmat[ut], method = "BH")
    padj[lower.tri(padj)] <- base::t(padj)[lower.tri(padj)]   # 対称化（t は翻訳ヘルパー）
    star <- function(p) ifelse(is.na(p), "",
              ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
    df <- data.frame(
      a = factor(rep(levs, times = K), levels = ord),
      b = factor(rep(levs, each = K), levels = ord),
      z = as.vector(z), padj = as.vector(padj), stringsAsFactors = FALSE)
    df$lab <- star(df$padj)
    df$tip <- paste0(df$a, " - ", df$b, "<br>z: ", round(df$z, 2),
                     "<br>padj: ", signif(df$padj, 2), " ", df$lab)
    # 自己ペア(対角)の z は非常に大きく、それに合わせると非対角(=注目したい
    # クラスター間)の色が白付近に潰れる。色の範囲は非対角の値を基準に決め、
    # 範囲外(対角など)は squish で両端の色にクリップして 0付近の解像度を上げる。
    offdiag <- df$z[as.character(df$a) != as.character(df$b)]
    lim <- stats::quantile(abs(offdiag), 0.98, na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- max(abs(df$z), na.rm = TRUE)
    lim <- max(lim, 2)   # 最低でも ±2 は表示
    p <- ggplot(df, aes(x = a, y = b, fill = z, text = tip)) +
      geom_tile(color = "grey80", linewidth = 0.3)   # 白に潰れても枠で見える
    if (isTRUE(input$spatial_ne_stars %||% TRUE))
      p <- p + geom_text(aes(label = lab), size = 3, color = "black", fontface = "bold")
    p +
      scale_fill_gradient2(low = "#3C5488FF", mid = "#F7F7F7", high = "#DC0000FF",
                           midpoint = 0, limits = c(-lim, lim),
                           oob = scales::squish, name = "z") +
      labs(x = NULL, y = NULL, title = t("spatial_ne_title")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 8), panel.grid = element_blank(),
            plot.background = element_rect(fill = pt$bg, color = NA),
            panel.background = element_rect(fill = pt$bg, color = NA),
            text = element_text(color = pt$fg), axis.text = element_text(color = pt$fg2),
            plot.title = element_text(size = 13, face = "bold", color = pt$accent))
  })
  if (requireNamespace("plotly", quietly = TRUE)) {
    output$spatial_ne_plot <- plotly::renderPlotly({
      plotly::ggplotly(spatial_ne_ggplot(), tooltip = "text")
    })
  }
  output$spatial_ne_plot_static <- renderPlot({ spatial_ne_ggplot() }, bg = "transparent")
}

# --- アプリ実行 ---
shinyApp(ui = ui, server = server)
