#' Render general outputs for the QC tab
#'
#' Create rendering expressions for the QC tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for QC tab features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny renderUI renderPlot tagList column selectInput isolate fluidRow
#' @importFrom plotly renderPlotly ggplotly layout config plot_ly 
#'   toWebGL add_segments add_annotations %>%
#' @importFrom ggplot2 scale_x_discrete guide_axis scale_y_continuous 
#'   geom_bar theme_minimal labs ggplot scale_fill_manual
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shinyjs js
#' @importFrom dittoSeq dittoColors
#' @importFrom grid grid.newpage grid.text
#' @importFrom DT renderDT datatable formatStyle
#' @importFrom stats cor as.formula
#' @rdname INTERNAL_create_qc_output
.create_qc_outputs <- function(input, output, robjects) {
    # nocov start
    output$pca.comps <- renderUI({
        req(robjects$pc)
        pcs <- robjects$pc

        tagList(
            fluidRow(
                column(4, selectInput("dim1", "Dim1:", choices = pcs$components, selected = "PC1")),
                column(4, selectInput("dim2", "Dim2:", choices = pcs$components, selected = "PC2")),
                if (length(pcs$components) > 2) {
                    column(4, selectInput("dim3", "Dim3:", choices = pcs$components, selected = "PC3"))
                }
            )
        )
    })
    # nocov end

    # nocov start
    output$qc.pca <- renderPlotly({
        req(input$dim1, input$dim2)

        input$pca.update

        pca.res <- robjects$pc

        colorb <- NULL
        shapeb <- NULL

        if (isolate(input$bip.color) != "") {
            colorb <- isolate(input$bip.color)
        }

        if (isolate(input$bip.shape) != "") {
            shapeb <- isolate(input$bip.shape)
        }

        if ((isolate(input$bip.twod) | length(pca.res$components) < 3)) {
            dizzy <- NULL
        } else {
            dizzy <- isolate(input$dim3)
        }

        if (any(unlist(lapply(list(pca.res, isolate(input$dim1), isolate(input$dim2)), is.null)))) {
            fig <- .empty_plot("No counts provided.\nNo PCA performed.", plotly = TRUE)
        } else {
            fig <- plot_pca_biplot(
                pca.res = pca.res,
                dim.x = isolate(input$dim1),
                dim.y = isolate(input$dim2),
                dim.z = dizzy,
                color.by = colorb,
                shape.by = shapeb,
                hover.info = "Label",
                show.loadings = isolate(input$bip.loadings),
                n.loadings = isolate(input$bip.n.loadings),
                pt.size = isolate(input$pca.pt.size)
            )
        }

        robjects$plot.qc.pca <- fig

        fig
    })
    # nocov end

    # nocov start
    output$qc.gini <- renderPlotly({
        if (!is.null(robjects$count.summary)) {
            fig <- plot_bar(robjects$count.summary)
        } else {
            fig <- .empty_plot("No summary provided.\nNo gini values available.", plotly = TRUE)
        }
        robjects$plot.qc.gini <- fig
        fig
    })
    # nocov end

    # nocov start
    output$qc.missed <- renderPlotly({
        if (!is.null(robjects$count.summary)) {
            fig <- plot_bar(robjects$count.summary,
                x = "Label", y = "Zerocounts", fill = "#394E80",
                ylab = "Zero Count sgRNAs", title = "Fully Depleted sgRNAs", yaxis.addition = 10
            )
        } else {
            fig <- .empty_plot("No summary provided.\nNo sgRNA depletion metrics available.", plotly = TRUE)
        }

        robjects$plot.qc.missed <- fig
        fig
    })
    # nocov end

    # nocov start
    output$qc.map <- renderPlotly({
        if (!is.null(robjects$count.summary)) {
            df <- robjects$count.summary

            #TODO: Move this to function, add controls for colors.
            df$Type <- "Mapped"
            df <- df[, c("Label", "Mapped", "Reads", "Percentage", "Type")]

            df <- rbind(
                df,
                data.frame(
                    Label = df$Label,
                    Mapped = df$Mapped,
                    Reads = df$Reads,
                    Percentage = 1 - df$Mapped / df$Reads,
                    Type = rep("Unmapped", nrow(df))
                )
            )

            fig <- ggplot(df, aes(.data[["Label"]], .data[["Percentage"]], fill = .data[["Type"]])) +
                geom_bar(stat = "identity") +
                theme_minimal() +
                labs(title = "Mapping Rates", x = NULL, y = "Mapped Ratio") +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                    axis.text.y = element_text(size = 12),
                    legend.title = element_blank()
                ) +
                scale_y_continuous(expand = c(0, 0)) +
                scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))

            fig <- ggplotly(fig) %>%
                layout(
                    xaxis = list(tickangle = 315),
                    legend=list(title=list(text=NULL))
                ) %>%
                config(
                    toImageButtonOptions = list(format = "svg"),
                    displaylogo = FALSE,
                    plotGlPixelRatio = 7,
                    edits = list(
                        axisTitleText = TRUE,
                        titleText = TRUE
                    )
                )
        } else {
            fig <- .empty_plot("No summary provided.\nNo mapping metrics available.",
                plotly = TRUE
            )
        }
        robjects$plot.qc.map <- fig
        fig
    })
    # nocov end

    # nocov start
    output$qc.histplot <- renderPlotly({
        n.counts <- robjects$norm.counts
        n.counts.log <- as.matrix(log2(n.counts[, c(-1, -2)] + 1))
        colnames(n.counts.log) <- colnames(n.counts)[c(-1, -2)]

        # TODO: Add input to control gridlines.
        if (!is.null(n.counts.log) & nrow(n.counts.log) > 0) {
            fig <- plot_hist(n.counts.log,
                title = "Distribution of read counts",
                xlab = "log2(counts + 1)", ylab = "Frequency", show.grid = FALSE
            )
        } else {
            fig <- .empty_plot("No counts provided.\nCannot make distribution.", plotly = TRUE)
        }

        robjects$plot.qc.hist <- fig
        fig
    })
    # nocov end

    # nocov start
    output$qc.corr <- renderPlot({
        input$corr.update

        n.counts <- robjects$norm.counts
        n.counts.log <- as.matrix(log2(n.counts[, c(-1, -2)] + 1))

        if (!is.null(isolate(input$count.summary_rows_all)) & isolate(input$meta.filt)) {
            n.counts.log <- as.matrix(n.counts.log[, isolate(input$count.summary_rows_all)])
        }

        if (ncol(n.counts.log) > 1) {
            cor.mat <- cor(n.counts.log)
            fig <- plot_correlation_heatmap(cor.mat,
                min.color = isolate(input$corr.min.col),
                max.color = isolate(input$corr.max.col)
            )
        } else {
            fig <- .empty_plot("No correlation possible.\nNo counts or only 1 sample provided.")
        }

        robjects$plot.qc.corr <- fig

        fig
    })
    # nocov end

    # nocov start
    output$count.summary <- renderDT(server = FALSE, {
        datatable(robjects$count.summary,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons", "Scroller"),
            options = list(
                search = list(regex = TRUE),
                lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "all")),
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print"),
                scrollX = TRUE,
                deferRender = TRUE,
                scrollY = 600,
                scroller = TRUE
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "80%")
    })
    # nocov end
    invisible(NULL)
}
