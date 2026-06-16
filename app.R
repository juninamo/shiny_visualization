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

generate_cluster_colors <- function(levels) {
  raw <- as.character(levels)
  clean <- raw[raw != "NA" & !is.na(raw)]
  if (length(clean) == 0) {
    out <- setNames(character(0), character(0))
  } else {
    num_prefix <- suppressWarnings(as.numeric(sub("^(\\d+).*", "\\1", clean)))
    coarse_raw <- regmatches(clean, regexpr("[A-Za-z_]+$", clean))
    # 末尾の系統名から先頭のアンダースコアを除去 ("_B_Plasma" -> "B_Plasma")
    coarse <- sub("^_+", "", coarse_raw)
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
  coarse_raw <- regmatches(levs, regexpr("[A-Za-z_]+$", levs))
  coarse <- sub("^_+", "", coarse_raw)
  if (length(coarse) == 0 || mean(coarse %in% known_lineages) < 0.5) return(NULL)
  generate_cluster_colors(levels)
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
get_expr_matrix <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)
  mat <- tryCatch(
    Seurat::GetAssayData(obj, layer = "data"),
    error = function(e) Seurat::GetAssayData(obj, slot = "data")
  )
  mat[genes, , drop = FALSE]
}

# =============================================================================
# 翻訳辞書
# =============================================================================
i18n <- list(
  ja = list(
    # サイドバー
    data_load       = "\U0001F4C2 データ読み込み",
    rds_file        = "RDSファイル",
    load_btn        = "読み込む",
    vis_settings    = "\U0001F52C 可視化設定",
    gene_name       = "遺伝子名",
    gene_placeholder = "ファイルを読み込んでください",
    group_var       = "グループ変数",
    umap_reduction  = "UMAP (reduction)",
    plot_settings   = "\U0001F3A8 プロット設定",
    pt_size         = "点のサイズ",
    plot_height     = "プロットの高さ (px)",
    plot_width      = "プロットの幅 (px)",

    # DEG
    deg_settings      = "DEG解析設定",
    deg_group_var     = "グループ変数",
    deg_cat_prefix    = "[カテゴリ] ",
    deg_num_prefix    = "[数値] ",
    deg_percentile    = "Top / Bottom パーセンタイル (%)",
    deg_group1        = "Group 1 (テスト群)",
    deg_group2        = "Group 2 (コントロール群・複数選択可)",
    deg_all_others    = "その他全て (All others)",
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
    comp_facets          = "ファセット変数 (任意・複数可)",
    comp_color_scheme    = "配色",
    comp_scheme_lineage  = "系統グラデーション (動的)",
    comp_scheme_manual   = "パレット",
    comp_title           = "クラスター組成",
    comp_yaxis           = "割合",
    comp_no_cat          = "カテゴリ変数が見つかりません。",

    # Heatmap / Dot plot 共通
    marker_settings      = "マーカー設定",
    marker_set           = "マーカーセット",
    marker_cluster       = "クラスター変数",
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
    notify_not_seurat   = "選択されたファイルはSeuratオブジェクトではありません。",
    notify_no_umap      = "注意: UMAPが計算されていません。UMAP表示は利用できません。",
    notify_load_done    = "✅ 読み込み完了: %s 細胞 × %s 遺伝子",
    notify_error        = "エラー: ",
    notify_deg_running  = "DEG解析を実行中...",
    notify_deg_same     = "Group 1とGroup 2は異なるグループを選択してください。",
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
    rds_file        = "RDS File",
    load_btn        = "Load",
    vis_settings    = "\U0001F52C Visualization",
    gene_name       = "Gene",
    gene_placeholder = "Load an RDS file first",
    group_var       = "Group Variable",
    umap_reduction  = "UMAP (reduction)",
    plot_settings   = "\U0001F3A8 Plot Settings",
    pt_size         = "Point Size",
    plot_height     = "Plot Height (px)",
    plot_width      = "Plot Width (px)",

    # DEG
    deg_settings      = "DEG Analysis Settings",
    deg_group_var     = "Group Variable",
    deg_cat_prefix    = "[Category] ",
    deg_num_prefix    = "[Numeric] ",
    deg_percentile    = "Top / Bottom Percentile (%)",
    deg_group1        = "Group 1 (Test)",
    deg_group2        = "Group 2 (Control, multi-select)",
    deg_all_others    = "All others",
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
    comp_facets          = "Facet Variables (optional)",
    comp_color_scheme    = "Color Scheme",
    comp_scheme_lineage  = "Lineage gradient (dynamic)",
    comp_scheme_manual   = "Manual palette",
    comp_title           = "Cluster Composition",
    comp_yaxis           = "Proportion",
    comp_no_cat          = "No categorical variables found.",

    # Heatmap / Dot plot shared
    marker_settings      = "Marker Settings",
    marker_set           = "Marker Set",
    marker_cluster       = "Cluster Variable",
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
    notify_not_seurat   = "The selected file is not a Seurat object.",
    notify_no_umap      = "Note: No UMAP reduction found. UMAP plots are unavailable.",
    notify_load_done    = "✅ Loaded: %s cells × %s genes",
    notify_error        = "Error: ",
    notify_deg_running  = "Running DEG analysis...",
    notify_deg_same     = "Group 1 and Group 2 must be different.",
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

    nav_panel(
      title = "\U0001F525 Heatmap",
      value = "heatmap",
      card_body(
        class = "p-2",
        uiOutput("heatmap_panel_ui")
      )
    ),

    nav_panel(
      title = "\U0001F535 Dot Plot",
      value = "dotplot",
      card_body(
        class = "p-2",
        uiOutput("dotplot_panel_ui")
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
        selected = rds_files[1]
      ),
      actionButton(
        "load_btn", t("load_btn"),
        class = "btn-primary w-100 mb-3",
        icon = icon("upload")
      ),

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

      # UMAP reduction 選択（umap関連が複数ある場合のみ表示）
      uiOutput("umap_reduction_ui"),

      hr(),

      h5(t("plot_settings"), class = "text-primary mb-2"),

      sliderInput("pt_size", t("pt_size"), min = 0, max = 2, value = 0.3, step = 0.1),
      sliderInput("plot_height", t("plot_height"), min = 400, max = 1200, value = 600, step = 50),
      sliderInput("plot_width", t("plot_width"), min = 400, max = 1600, value = 800, step = 50)
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
  observeEvent(input$load_btn, {
    req(input$rds_file)

    showNotification(t("notify_loading"), type = "message", id = "loading")

    tryCatch({
      file_path <- file.path(app_dir, input$rds_file)
      obj <- readRDS(file_path)

      # Seuratオブジェクトか確認
      if (!inherits(obj, "Seurat")) {
        showNotification(t("notify_not_seurat"),
                         type = "error", id = "loading")
        return()
      }

      seurat_obj(obj)
      data_loaded(TRUE)
      deg_results(NULL)

      # 遺伝子リストを更新
      genes <- sort(rownames(obj))
      updateSelectizeInput(session, "gene",
                           choices = genes,
                           selected = genes[1],
                           server = TRUE)

      # meta.dataの列を分類
      meta <- obj@meta.data
      cat_cols <- names(meta)[sapply(meta, function(x) {
        is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 50)
      })]
      num_cols <- names(meta)[sapply(meta, function(x) {
        is.numeric(x) && length(unique(x)) > 50
      })]
      meta_col_types(list(cat = cat_cols, num = num_cols))

      # seurat_clustersを優先的にデフォルトにする
      default_group <- if ("seurat_clusters" %in% cat_cols) {
        "seurat_clusters"
      } else {
        cat_cols[1]
      }

      updateSelectInput(session, "group_var",
                        choices = cat_cols,
                        selected = default_group)

      # DEG用のグループ変数: カテゴリカル + 数値を表示
      deg_choices <- c(
        setNames(cat_cols, paste0(t("deg_cat_prefix"), cat_cols)),
        setNames(num_cols, paste0(t("deg_num_prefix"), num_cols))
      )
      updateSelectInput(session, "deg_group_var",
                        choices = deg_choices,
                        selected = default_group)

      # UMAP存在チェック（"umap" がなければ "umap" を含む reduction を探す）
      has_umap <- !is.null(find_umap_reduction(obj))
      if (!has_umap) {
        showNotification(
          t("notify_no_umap"),
          type = "warning", duration = 8
        )
      }

      n_cells <- ncol(obj)
      n_genes <- nrow(obj)
      showNotification(
        sprintf(t("notify_load_done"),
                format(n_cells, big.mark = ","),
                format(n_genes, big.mark = ",")),
        type = "message", id = "loading", duration = 5
      )

    }, error = function(e) {
      showNotification(paste(t("notify_error"), e$message), type = "error", id = "loading")
    })
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
      # カテゴリカル列の場合: グループ選択
      obj <- seurat_obj()
      groups <- sort(unique(as.character(obj@meta.data[[input$deg_group_var]])))
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
  # Violin Plot
  # ==========================================================================
  output$violin_ui <- renderUI({
    if (!data_loaded()) return(placeholder_ui())
    plotOutput("violin_plot",
               height = paste0(input$plot_height, "px"),
               width = paste0(input$plot_width, "px"))
  })

  output$violin_plot <- renderPlot({
    req(seurat_obj(), input$gene, input$group_var)
    obj <- seurat_obj()
    pt <- plot_theme()
    Idents(obj) <- obj@meta.data[[input$group_var]]
    p <- VlnPlot(obj, features = input$gene, pt.size = input$pt_size) + pt$theme
    # 系統(lineage)グラデーション配色（形式が合えば適用、合わなければデフォルト）
    lin_cols <- lineage_colors_or_null(levels(Idents(obj)))
    if (!is.null(lin_cols)) p <- p + scale_fill_manual(values = lin_cols)
    p
  }, bg = "transparent")

  # ==========================================================================
  # Feature UMAP
  # ==========================================================================
  output$feature_umap_ui <- renderUI({
    if (!data_loaded()) return(placeholder_ui())
    if (!has_umap()) return(no_umap_ui())
    plotOutput("feature_umap_plot",
               height = paste0(input$plot_height, "px"),
               width = paste0(input$plot_width, "px"))
  })

  output$feature_umap_plot <- renderPlot({
    req(seurat_obj(), input$gene, has_umap())
    obj <- seurat_obj()
    pt <- plot_theme()
    FeaturePlot(obj, features = input$gene, reduction = umap_reduction(),
                pt.size = input$pt_size) + pt$theme_legend
  }, bg = "transparent")

  # ==========================================================================
  # Group UMAP
  # ==========================================================================
  output$group_umap_ui <- renderUI({
    if (!data_loaded()) return(placeholder_ui())
    if (!has_umap()) return(no_umap_ui())
    plotOutput("group_umap_plot",
               height = paste0(input$plot_height, "px"),
               width = paste0(input$plot_width, "px"))
  })

  output$group_umap_plot <- renderPlot({
    req(seurat_obj(), input$group_var, has_umap())
    obj <- seurat_obj()
    pt <- plot_theme()
    p <- DimPlot(obj, reduction = umap_reduction(), group.by = input$group_var,
                 label = TRUE, pt.size = input$pt_size) + pt$theme_legend
    # 系統(lineage)グラデーション配色（形式が合えば適用、合わなければデフォルト）
    lin_cols <- lineage_colors_or_null(obj@meta.data[[input$group_var]])
    if (!is.null(lin_cols)) p <- p + scale_color_manual(values = lin_cols)
    p
  }, bg = "transparent")

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
    x_sel <- if (!is.null(cur_x) && cur_x %in% cat_cols) cur_x else cat_cols[1]
    cur_facets <- isolate(input$comp_facets)
    facet_sel <- cur_facets[cur_facets %in% cat_cols]
    cur_scheme <- isolate(input$comp_color_scheme)
    scheme_sel <- if (!is.null(cur_scheme)) cur_scheme else "lineage"

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
                          choices = cat_cols, selected = x_sel)
            ),
            column(4,
              selectInput("comp_facets", t("comp_facets"),
                          choices = cat_cols, selected = facet_sel,
                          multiple = TRUE)
            )
          ),
          fluidRow(
            column(12,
              radioButtons("comp_color_scheme", t("comp_color_scheme"),
                           choices = c(setNames("lineage", t("comp_scheme_lineage")),
                                       setNames("manual",  t("comp_scheme_manual"))),
                           selected = scheme_sel, inline = TRUE)
            )
          )
        )
      ),
      uiOutput("comp_plot_ui")
    )
  })

  output$comp_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    plotOutput("comp_plot",
               height = paste0(input$plot_height, "px"),
               width = paste0(input$plot_width, "px"))
  })

  output$comp_plot <- renderPlot({
    req(seurat_obj(), input$comp_cluster, input$comp_x)
    obj <- seurat_obj()
    pt <- plot_theme()
    meta <- obj@meta.data

    fill_var <- input$comp_cluster
    x_var <- input$comp_x
    facet_vars <- input$comp_facets
    facet_vars <- facet_vars[facet_vars %in% names(meta)]
    group_vars <- unique(c(x_var, facet_vars))
    key_cols <- unique(c(fill_var, group_vars))

    # --- 割合を計算 (group_vars 内で proportion を算出) ---
    df <- meta[key_cols]
    for (cc in key_cols) df[[cc]] <- as.character(df[[cc]])
    counts <- aggregate(list(n = rep(1L, nrow(df))), by = df, FUN = sum)
    grp_key <- do.call(paste, c(counts[group_vars], sep = "\r"))
    totals <- tapply(counts$n, grp_key, sum)
    counts$proportion <- counts$n / as.numeric(totals[grp_key])

    # --- 配色 ---
    levs <- sort(unique(counts[[fill_var]]))
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

    p <- ggplot(counts, aes(x = .data[[x_var]], y = proportion,
                            fill = .data[[fill_var]])) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = fill_colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = t("comp_title"), x = x_var, y = t("comp_yaxis"), fill = fill_var) +
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

    # --- ファセット (ggh4x があれば nested、なければ facet_grid) ---
    if (length(facet_vars) > 0) {
      fct <- stats::as.formula(paste("~", paste(facet_vars, collapse = " + ")))
      if (requireNamespace("ggh4x", quietly = TRUE)) {
        p <- p + ggh4x::facet_nested(fct, scales = "free_x",
                                     space = "free_x", nest_line = TRUE)
      } else {
        p <- p + facet_grid(fct, scales = "free_x", space = "free_x")
      }
    }

    p
  }, bg = "transparent")

  # ==========================================================================
  # マーカー設定UI（Heatmap / Dot plot 共通）
  # ==========================================================================
  # prefix で input ID を分け、tab ごとに状態を保持する
  marker_settings_card <- function(prefix, extra = NULL) {
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

    div(
      class = "card mb-3",
      div(
        class = "card-body",
        h6(t("marker_settings"), class = "card-title text-primary"),
        fluidRow(
          column(6,
            selectInput(paste0(prefix, "_cluster"), t("marker_cluster"),
                        choices = cat_cols, selected = clu_sel)
          ),
          column(6,
            selectInput(paste0(prefix, "_set"), t("marker_set"),
                        choices = names(marker_sets), selected = set_sel)
          )
        ),
        extra
      )
    )
  }

  # --- マーカーセットの平均発現マトリクス (cluster x gene) を計算 ---
  # genes はセット順を維持。z-score 行(=遺伝子)スケーリングは呼び出し側で。
  marker_avg_matrix <- function(obj, cluster_var, genes) {
    mat <- get_expr_matrix(obj, genes)        # genes x cells（存在する遺伝子のみ）
    if (is.null(mat)) return(NULL)
    labels <- as.character(obj@meta.data[[cluster_var]])
    keep <- !is.na(labels)
    mat <- mat[, keep, drop = FALSE]
    labels <- labels[keep]
    # クラスターごとの平均
    cl_levels <- sort(unique(labels))
    avg <- sapply(cl_levels, function(cl) {
      Matrix::rowMeans(mat[, labels == cl, drop = FALSE])
    })
    if (is.null(dim(avg))) avg <- matrix(avg, nrow = nrow(mat),
                                         dimnames = list(rownames(mat), cl_levels))
    avg  # genes x clusters
  }

  # ==========================================================================
  # Heatmap
  # ==========================================================================
  output$heatmap_panel_ui <- renderUI({
    lang <- input$lang
    if (!data_loaded()) return(placeholder_ui())
    col_types <- meta_col_types()
    if (length(col_types$cat) == 0) {
      return(div(class = "text-center text-muted py-4", h5(t("comp_no_cat"))))
    }
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      return(div(class = "text-center text-warning py-4", h5(t("hm_pkg_missing"))))
    }

    extra <- fluidRow(
      column(4,
        checkboxInput("hm_scale", t("hm_scale"),
                      value = isolate(input$hm_scale) %||% TRUE)
      ),
      column(4,
        checkboxInput("hm_cluster_rows", t("hm_cluster_rows"),
                      value = isolate(input$hm_cluster_rows) %||% TRUE)
      ),
      column(4,
        checkboxInput("hm_cluster_cols", t("hm_cluster_cols"),
                      value = isolate(input$hm_cluster_cols) %||% TRUE)
      )
    )

    tagList(
      marker_settings_card("hm", extra),
      uiOutput("heatmap_plot_ui")
    )
  })

  output$heatmap_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    plotOutput("heatmap_plot",
               height = paste0(input$plot_height, "px"),
               width = paste0(input$plot_width, "px"))
  })

  output$heatmap_plot <- renderPlot({
    req(seurat_obj(), input$hm_cluster, input$hm_set,
        requireNamespace("pheatmap", quietly = TRUE))
    obj <- seurat_obj()
    pt <- plot_theme()

    set_df <- marker_set_to_df(marker_sets[[input$hm_set]])
    genes <- unique(set_df$feature)
    avg <- marker_avg_matrix(obj, input$hm_cluster, genes)   # genes x clusters
    validate(need(!is.null(avg) && nrow(avg) > 0, t("hm_no_genes")))

    mat <- avg
    if (isTRUE(input$hm_scale)) {
      mat <- t(scale(t(mat)))          # 遺伝子(行)ごとに Z-score
      mat <- mat[stats::complete.cases(mat), , drop = FALSE]
      mat[mat > 2] <- 2
      mat[mat < -2] <- -2
    }
    validate(need(nrow(mat) > 0, t("hm_no_genes")))

    # 行=クラスター, 列=遺伝子 にして遺伝子名を45度ラベル
    pheatmap::pheatmap(
      mat = t(mat),
      color = grDevices::colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(50),
      border_color = "white",
      show_rownames = TRUE,
      show_colnames = TRUE,
      cluster_rows = isTRUE(input$hm_cluster_rows),
      cluster_cols = isTRUE(input$hm_cluster_cols),
      treeheight_row = 0,
      treeheight_col = 0,
      angle_col = 45,
      fontsize = 11,
      scale = "none",
      silent = FALSE
    )
  }, bg = "white")

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

    extra <- fluidRow(
      column(6,
        sliderInput("dot_scale", t("dot_scale"),
                    min = 2, max = 12,
                    value = isolate(input$dot_scale) %||% 6, step = 1)
      ),
      column(6,
        checkboxInput("dot_facet", t("dot_facet"),
                      value = isolate(input$dot_facet) %||% TRUE)
      )
    )

    tagList(
      marker_settings_card("dot", extra),
      uiOutput("dotplot_plot_ui")
    )
  })

  output$dotplot_plot_ui <- renderUI({
    if (!data_loaded()) return(NULL)
    plotOutput("dotplot_plot",
               height = paste0(input$plot_height, "px"),
               width = paste0(input$plot_width, "px"))
  })

  output$dotplot_plot <- renderPlot({
    req(seurat_obj(), input$dot_cluster, input$dot_set)
    obj <- seurat_obj()
    pt <- plot_theme()

    set_df <- marker_set_to_df(marker_sets[[input$dot_set]])
    genes_all <- unique(set_df$feature)
    mat <- get_expr_matrix(obj, genes_all)     # genes x cells
    validate(need(!is.null(mat) && nrow(mat) > 0, t("hm_no_genes")))

    labels <- as.character(obj@meta.data[[input$dot_cluster]])
    keep <- !is.na(labels)
    mat <- mat[, keep, drop = FALSE]
    labels <- labels[keep]
    cl_levels <- sort(unique(labels))

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

    p <- ggplot(dot, aes(x = feature, y = id)) +
      geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
      scale_radius(range = c(0, input$dot_scale), limits = c(0, 100)) +
      scale_color_gradient2(midpoint = 0, low = "#3C5488FF",
                            mid = "grey90", high = "#DC0000FF", space = "Lab") +
      labs(x = NULL, y = input$dot_cluster,
           size = t("dot_pct"), color = t("dot_avg")) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
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
  }, bg = "transparent")

  # ==========================================================================
  # DEG解析
  # ==========================================================================
  observeEvent(input$run_deg, {
    req(seurat_obj(), input$deg_group_var)

    showNotification(t("notify_deg_running"), type = "message", id = "deg_run")

    tryCatch({
      obj <- seurat_obj()

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
                               logfc.threshold = 0,
                               min.pct = 0.1)
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
                               logfc.threshold = 0,
                               min.pct = 0.1)
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

    tagList(
      # Volcano Plot
      div(class = "mb-3",
        plotOutput("volcano_plot",
                   height = "500px",
                   width = paste0(input$plot_width, "px"))
      ),
      # DEGテーブル
      div(
        DTOutput("deg_table")
      )
    )
  })

  # --- Volcano Plot ---
  output$volcano_plot <- renderPlot({
    req(deg_results())
    res <- deg_results()
    pt <- plot_theme()

    colors <- c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#95a5a6")

    # ラベル用: top 10遺伝子
    top_genes <- res[res$significance != "NS", ]
    top_genes <- top_genes[order(top_genes$p_val_adj), ]
    top_genes <- head(top_genes, 10)

    p <- ggplot(res, aes(x = avg_log2FC, y = neg_log10_pval, color = significance)) +
      geom_point(alpha = 0.6, size = 1.5) +
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

    # トップ遺伝子ラベル
    if (nrow(top_genes) > 0) {
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
  }, bg = "transparent")

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
