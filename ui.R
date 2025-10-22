library(shiny)
library(shinyWidgets)
library(DT)

ui <- fluidPage(
  tags$head(tags$style(HTML("
/* DataTables：禁止横向滚动并让表头表体同宽 */
.dataTables_wrapper .dataTables_scrollBody{ overflow-x:hidden!important; }
.dataTables_scrollHeadInner,
.dataTables_scrollHead table,
.dataTables_scrollBody table{ width:100%!important; }

/* 顶部开关布局 */
.controls-row{display:flex;align-items:center;gap:16px;flex-wrap:nowrap;}

/* 图表容器 */
#charts-row{ margin-bottom:16px; }
.slot-left,.slot-right{
  display:block; margin:0; padding:0 8px;
}
.slot-left .shiny-image-output,
.slot-right .shiny-image-output,
.slot-left .shiny-plot-output,
.slot-right .shiny-plot-output{
  width:100%!important; height:520px!important; margin:0;
}

/* 曼哈顿/QQ 图 */
#manhattan,#qq_plot{margin:0!important;padding:0!important;}
#manhattan img,#qq_plot img{
  width:100%!important; height:100%!important;
  object-fit:contain; margin:0!important; padding:0!important;
}

/* 表格紧凑 */
.match-rows.dataTable{font-size:13px;}
.match-rows.dataTable tbody td{
  padding:6px 8px!important; line-height:16px!important; white-space:nowrap;
  overflow:hidden; text-overflow:ellipsis;
}
.match-rows.dataTable tbody tr{height:28px!important;}

/* 空状态 */
.empty-hint{
  border:1px dashed #d9d9d9; background:#fafafa; border-radius:10px;
  padding:10px 12px; color:#6b7280; font-size:13px; margin-top:8px;
}

/* 标题 */
.chart-title{ margin:8px 0 4px!important; font-weight:500; line-height:1.2; }

/* 间距辅助 */
.after-table{ margin-top:30px!important; }
.plot-gap{ margin-top:16px!important; }

/* Step 1 的附属子面板 */
.step-sub {
  margin: 6px 0 0 16px;         
  padding: 10px 12px;
  border-left: 3px solid #e5e7eb;
  background: #fafafa;
  border-radius: 8px;
}
.step-sub h6 {
  margin: 0 0 6px;
  font-weight: 600;
  font-size: 13px;
  color: #374151;
}
.step-sub .subhelp {
  margin-top: 6px;
  color: #6b7280;
  font-size: 12px;
}

/* 侧栏步骤标题加粗 */
.sidebar-steps h5 {
  font-weight: 600;
  font-size: 14px;
  margin-top: 16px;
  margin-bottom: 6px;
}
"))),
  
  titlePanel("HELLO: Heritability Estimation and Local LD by blOck"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "sidebar-steps",
          
          h5("1. Select/Upload GWAS"),
          fileInput("gwas_file",NULL,accept=c(".gz",".txt",".tsv",".csv")),
          helpText("File must contain columns: SNP, CHR, POS, P (tab- or comma-delimited)"),
          
          div(class = "step-sub",
              h6("Map Columns"),
              div(
                style = "display:grid;grid-template-columns:60px 1fr;gap:6px 8px;align-items:center;",
                div(tags$strong("SNP")), selectInput("map_snp", NULL, choices = character(0)),
                div(tags$strong("CHR")), selectInput("map_chr", NULL, choices = character(0)),
                div(tags$strong("POS")), selectInput("map_pos", NULL, choices = character(0)),
                div(tags$strong("P")),   selectInput("map_p",   NULL, choices = character(0))
              ),
              div(class = "subhelp", "Auto-detected defaults will appear here; adjust if needed.")
          ),
          
          tags$hr(),
          selectInput("gwas_builtin","Example GWAS",choices=character(0)),
          tags$hr(),
          
          h5("2. Select Population"),
          selectInput("ancestry",NULL,c("EUR","AFR","EAS","SAS","Mixed"),"EUR"), tags$hr(),
          
          h5("3. Select GRCh"),
          selectInput("build",NULL,c("37","38"),"37"), tags$hr(),
          
          h5("4. P-value"),
          div(class = "controls-row",
              span("Loci"),
              textInput("p_threshold", NULL, value = "5e-8", width = "180px")
          ),
          div(class = "controls-row",
              span("Block"),
              textInput("block_p_threshold", NULL, value = "5e-8", width = "180px")
          ),
          tags$hr(),
          
          actionButton("apply_block","Apply Selection", class="btn-primary"), tags$hr()
      )
    ),
    
    mainPanel(
      width = 9,
      
      fluidRow(
        id = "charts-row",
        column(
          6, style = "padding-left:8px; padding-right:8px;",
          div(class="slot-left",
              imageOutput("manhattan", width="100%", height="520px")
          )
        ),
        column(
          6, style = "padding-left:8px; padding-right:8px;",
          div(class="slot-right",
              imageOutput("qq_plot", width="100%", height="520px")
          )
        )
      ),
      
      fluidRow(
        column(6,
               h4("Significant Loci"),
               DTOutput("top_table"), uiOutput("top_table_empty"),
               div(class = "after-table"),
               h4(textOutput("locuszoom_title"), class="chart-title"),
               div(class = "plot-gap",
                   plotOutput("locuszoom_plot", height = "340px")
               )
        ),
        column(6,
               h4("Significant Block"),
               DTOutput("signal_block_table"), uiOutput("signal_block_empty"),
               div(class = "after-table"),
               h4(textOutput("blockzoom_title"), class="chart-title"),
               div(class = "plot-gap",
                   plotOutput("blockzoom_plot", height = "520px")
               )
        )
      )
    )
  )
)
