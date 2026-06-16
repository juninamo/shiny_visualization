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

  # --- サイドバー ---
  sidebar = sidebar(
    width = 320,

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

    nav_panel(
      title = "\U0001F3BB Violin",
      value = "violin",
      card_body(
        class = "p-2",
        uiOutput("violin_ui")
      )
    ),

    nav_panel(
      title = "\U0001F5FA\uFE0F Feature UMAP",
      value = "feature_umap",
      card_body(
        class = "p-2",
        uiOutput("feature_umap_ui")
      )
    ),

    nav_panel(
      title = "\U0001F3F7\uFE0F Group UMAP",
      value = "group_umap",
      card_body(
        class = "p-2",
        uiOutput("group_umap_ui")
      )
    ),

    nav_panel(
      title = "\U0001F4CA DEG",
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

      # アクティブ / リファレンス データセット選択（読み込み後に表示）
      uiOutput("dataset_selectors_ui"),

      hr(),

      h5(t("vis_settings"), class = "text-primary mb-2"),

      selectizeInput(
        "gene", t("gene_name"),
        choices = NULL,
        options = list(
          placeholder = t("gene_placeholder"),
          maxOptions = 50
        )
      ),

      # 外部データベースリンク
      uiOutput("external_links_ui"),

      selectInput(
        "group_var", t("group_var"),
        choices = NULL
      ),

      # リファレンス用グループ変数（比較時のみ表示）
      uiOutput("ref_group_var_ui"),

      # UMAP reduction 選択（umap関連が複数ある場合のみ表示）
      uiOutput("umap_reduction_ui"),

      hr(),

      h5(t("plot_settings"), class = "text-primary mb-2"),

      sliderInput("pt_size", t("pt_size"), min = 0, max = 2, value = 0.3, step = 0.1),
      sliderInput("umap_label_size", t("umap_label_size"), min = 2, max = 12, value = 4, step = 0.5),
      sliderInput("plot_height", t("plot_height"), min = 400, max = 3000, value = 600, step = 50),
      sliderInput("plot_width", t("plot_width"), min = 400, max = 4000, value = 800, step = 50),

      # リファレンス用のプロット設定（比較時のみ表示）
      uiOutput("ref_plot_settings_ui")
    )
  })

  # --- 外部データベースリンク ---
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
  observeEvent(input$lang, {
    if (data_loaded()) {
      obj <- seurat_obj()
      meta <- obj@meta.data
      col_types <- meta_col_types()
      cat_cols <- col_types$cat
      num_cols <- col_types$num

      # 遺伝子リスト: 現在選択中を維持
      current_gene <- input$gene
      genes <- sort(rownames(obj))
      updateSelectizeInput(session, "gene",
                           choices = genes,
                           selected = current_gene,
                           server = TRUE)

      # グループ変数: 現在の選択を維持
      current_group <- input$group_var
      updateSelectInput(session, "group_var",
                        choices = cat_cols,
                        selected = current_group)

      # DEGグループ変数: ラベルを再生成
      current_deg_var <- input$deg_group_var
      deg_choices <- c(
        setNames(cat_cols, paste0(t("deg_cat_prefix"), cat_cols)),
        setNames(num_cols, paste0(t("deg_num_prefix"), num_cols))
      )
      updateSelectInput(session, "deg_group_var",
                        choices = deg_choices,
                        selected = current_deg_var)
    }
  }, ignoreInit = TRUE)

  # ==========================================================================
  # RDSファイル読み込み
  # ==========================================================================
  # --- アクティブにするデータセットで各セレクタ・メタ情報を更新 ---
  # 現在の遺伝子/グループ選択はデータセット切替後も（存在すれば）維持する。
  activate_dataset <- function(obj) {
    seurat_obj(obj)
    data_loaded(TRUE)
    deg_results(NULL)

    # 遺伝子リスト（現在の選択を維持）
    genes <- sort(rownames(obj))
    cur_gene <- isolate(input$gene)
    gsel <- if (!is.null(cur_gene) && cur_gene %in% genes) cur_gene else genes[1]
    updateSelectizeInput(session, "gene", choices = genes,
                         selected = gsel, server = TRUE)

    # meta.dataの列を分類
    meta <- obj@meta.data
    cat_cols <- names(meta)[sapply(meta, function(x) {
      is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 50)
    })]
    num_cols <- names(meta)[sapply(meta, function(x) {
      is.numeric(x) && length(unique(x)) > 50
    })]
    meta_col_types(list(cat = cat_cols, num = num_cols))

    cur_group <- isolate(input$group_var)
    default_group <- if (!is.null(cur_group) && cur_group %in% cat_cols) {
      cur_group
    } else if ("seurat_clusters" %in% cat_cols) {
      "seurat_clusters"
    } else {
      cat_cols[1]
    }
    updateSelectInput(session, "group_var",
                      choices = cat_cols, selected = default_group)

    cur_deg <- isolate(input$deg_group_var)
    deg_choices <- c(
      setNames(cat_cols, paste0(t("deg_cat_prefix"), cat_cols)),
      setNames(num_cols, paste0(t("deg_num_prefix"), num_cols))
    )
    deg_sel <- if (!is.null(cur_deg) && cur_deg %in% c(cat_cols, num_cols)) {
      cur_deg
    } else default_group
    updateSelectInput(session, "deg_group_var",
                      choices = deg_choices, selected = deg_sel)

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
                            selected = isolate(input$deg_marker_set) %||% names(marker_sets)[1])
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

      # 結果: Volcano + テーブル
      uiOutput("deg_results_ui")
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
    if (mode %in% c("set", "set_custom") && !is.null(input$deg_marker_set)) {
      set_genes <- unique(marker_set_to_df(marker_sets[[input$deg_marker_set]])$feature)
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
        plotly::layout(p,
          xaxis = list(title = axn[1]), yaxis = list(title = axn[2]),
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
    set_sel <- if (!is.null(cur_set) && cur_set %in% names(marker_sets)) {
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
                        choices = names(marker_sets), selected = set_sel)
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
                      choices = names(marker_sets), selected = set_sel)
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
    setn <- input[[paste0(prefix, "_set")]]
    parts <- list()
    if (mode %in% c("set", "set_custom") && !is.null(setn)) {
      sdf <- marker_set_to_df(marker_sets[[setn]])
      # ユーザーがセット内から削除した遺伝子を反映（kept が NULL のときは全て）
      kept <- input[[paste0(prefix, "_set_genes")]]
      if (!is.null(kept)) sdf <- sdf[sdf$feature %in% kept, , drop = FALSE]
      parts[[length(parts) + 1]] <- sdf
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
    if (is.null(setn) || !(setn %in% names(marker_sets))) return(NULL)
    all_genes <- unique(marker_set_to_df(marker_sets[[setn]])$feature)
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
        plot_run_btn("hm_run")
      )),
      uiOutput("heatmap_plot_ui")
    )
  })

  output$heatmap_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    if ((input$hm_run %||% 0) == 0) return(run_hint_ui())
    side_by_side_ui("heatmap_plot", "heatmap_plot_ref", !is.null(ref_obj()))
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
  build_hm_ggplot <- function(mat, main, cluster_rows, cluster_cols) {
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
    df <- data.frame(
      cluster = factor(rep(rownames(m), times = ncol(m)), levels = row_order),
      gene    = factor(rep(colnames(m), each = nrow(m)),  levels = col_order),
      value   = as.vector(m),
      stringsAsFactors = FALSE
    )
    lim <- max(abs(df$value), na.rm = TRUE)
    ggplot(df, aes(x = gene, y = cluster, fill = value)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#0072B5FF", mid = "white", high = "#BC3C29FF",
                           midpoint = 0, limits = c(-lim, lim), name = NULL) +
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
        plot.title = element_text(size = 14, face = "bold", color = pt$accent)
      )
  }

  heatmap_spec <- eventReactive(input$hm_run, {
    req(seurat_obj(), input$mk_cluster)
    ma <- build_hm_mat(seurat_obj())
    validate(need(!is.null(ma) && nrow(ma) > 0, t("hm_no_genes")))
    rb <- ref_obj()
    list(active = ma,
         ref = if (!is.null(rb)) {
           build_hm_mat(rb, cluster_var = ref_tab_cluster("mk"),
                        clusters_sel = input$mk_clusters_sel_ref)
         } else NULL,
         ref_present = !is.null(rb),
         names = c(active_name(), ref_name() %||% ""),
         cluster_rows = isTRUE(input$hm_cluster_rows),
         cluster_cols = isTRUE(input$hm_cluster_cols))
  }, ignoreInit = TRUE)

  output$heatmap_plot <- renderPlot({
    spec <- heatmap_spec()
    build_hm_ggplot(spec$active, spec$names[1], spec$cluster_rows, spec$cluster_cols)
  }, bg = "transparent")

  output$heatmap_plot_ref <- renderPlot({
    spec <- heatmap_spec()
    req(spec$ref_present)
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
        plot_run_btn("dot_run")
      )),
      uiOutput("dotplot_plot_ui")
    )
  })

  output$dotplot_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    if ((input$dot_run %||% 0) == 0) return(run_hint_ui())
    if (requireNamespace("plotly", quietly = TRUE)) {
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
      # Volcano Plot
      div(class = "mb-3", volcano_out),
      # DEGテーブル
      div(
        DTOutput("deg_table")
      )
    )
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
}

# --- アプリ実行 ---
shinyApp(ui = ui, server = server)
