library(tidyverse)
library(data.table)
library(magrittr)
library(ggridges)
library(ggthemes)
library(harmony)
library(singlecellmethods)
library(ggrepel)
library(Matrix)
library(patchwork)
library(sctransform)
library(scales)
library(uwot)
library(pheatmap)


fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}


do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE, do_labels = TRUE, nice_names, palette_use,
                       pt_size = 4, point_size = .5, base_size = 12, do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
    plt_df <- umap_use %>% data.frame() %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L) 
    plt_df$given_name <- plt_df[[label_name]]
    
    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))

        plt_df[[label_name]] <- plt_df$nice_name        
    }
        
    plt <- plt_df %>% 
        ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
#         theme_tufte(base_size = base_size) + 
        theme_test(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4)), alpha = FALSE) +
        scale_color_manual(values = palette_use) + 
        scale_fill_manual(values = palette_use) +    
        theme(plot.title = element_text(hjust = .5)) + 
        labs(x = "UMAP 1", y = "UMAP 2") 
    
    if (do_points) 
#         plt <- plt + geom_point_rast(dpi = 300, width = w, height = h, size = point_size) 
        plt <- plt + geom_point(shape = '.')
    if (do_density) 
        plt <- plt + geom_density_2d()    
        

    if (no_guides)
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    
    if (do_labels) 
        plt <- plt + geom_label_repel(data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], label.size = NA,
                                      aes_string(label = label_name), color = "white", size = pt_size, alpha = 1, segment.size = 0) + 
        guides(col = FALSE, fill = FALSE)
    return(plt)
}



setupVals <- function(data_mat, feature, qlo, qhi) {
    .x <- data_mat[feature, , drop = FALSE] %>% as("dgTMatrix")
    cutoffs <- quantileSparse(.x, c(qlo, qhi))
    cutoffs[2] <- max(cutoffs[2], min(.x@x))
    if (qlo == 0 & qhi == 1) {
        return(.x)
    } 
    
    if (qlo > 0) {
        .x@x[.x@x < cutoffs[1]] <- cutoffs[1]
#         message(sprintf("For %s, lo = %.3f", feature, ifelse(length(.x@x) == ncol(.x), cutoffs[1], NA)))
    }
    if (qhi < 1) {
        .x@x[.x@x > cutoffs[2]] <- cutoffs[2]
#         message(sprintf("For %s, hi = %.3f", feature, cutoffs[2]))
        
    }
    return(.x)
}


quantileSparse <- function(.x, qlist) {
    ratio_zero <- 1 - (length(.x@x) / ncol(.x))
    q_nz <- which(qlist > ratio_zero)
    q_adj <- (qlist[q_nz] - ratio_zero) / (1 - ratio_zero)
    res <- rep(0, length(qlist))
    res[q_nz] <- quantile(.x@x, q_adj)
    res
}

## TODO: test is feature is present
## TODO: allow for different cutoffs, for each marker
## TODO: somehow draw canvas first, then do plotting? 
library(patchwork)
library(ggthemes)

plotFeatures <- function(data_mat, dim_df, features, nrow = 1, 
                         qlo = 0.05, qhi = 1, order_by_expression = FALSE, 
                         pt_shape = 16, pt_size = .5, no_guide = FALSE,
                         .xlim = c(NA, NA), .ylim = c(NA, NA)) {
    plt_df <- data.frame(dim_df[, 1:2])
    colnames(plt_df) <- c("X1", "X2")


    plt_list <- lapply(features, function(feature) {
        .x <- setupVals(data_mat, feature, qlo, qhi)
        plt_df$value <- 0
        plt_df[.x@j + 1, "value"] <- .x@x
        if (order_by_expression) {
            plt_df %<>% dplyr::arrange(value)             
        } else {
            plt_df %<>% dplyr::sample_frac(1L)
        }

        plt <- plt_df %>% 
            ggplot(aes(X1, X2, color = value)) + 
            geom_point(size = .5, shape = pt_shape) + 
#             geom_point_rast(dpi = 300, width = 6, height = 4, size = .5, shape = pt_shape) + 
#             geom_point(shape = ".") + 
            scale_color_gradient2(na.value = "lightgrey", mid = "lightgrey", midpoint = 0, high = muted("blue")) + 
            theme_tufte(base_size = 14, base_family = "Helvetica") + 
            theme(panel.background = element_rect(), plot.title = element_text(hjust = .5)) +
            labs(x = "UMAP 1", y = "UMAP 2", title = feature) + 
            NULL
        if (no_guide) {
            plt <- plt + 
            guides(color = FALSE) 
        }
        
        if (sum(is.na(.xlim)) < 2) 
            plt <- plt + xlim(.xlim)
        if (sum(is.na(.ylim)) < 2) 
            plt <- plt + ylim(.ylim)
        plt

    })
    if (length(plt_list) > 1) {
        Reduce(`+`, plt_list) + patchwork::plot_layout(nrow = nrow)
    } else {
        plt_list[[1]]
    }
}


read10x <- function(run, suffix) {
#     barcode.loc <- file.path(run, "barcodes.tsv")
#     gene.loc <- file.path(run, "genes.tsv")
#     matrix.loc <- file.path(run, "matrix.mtx")
    barcode.loc <- list.files(run, pattern = 'barcodes.tsv', full.names = TRUE)
    gene.loc <- list.files(run, pattern = 'features.tsv|genes.tsv', full.names = TRUE)
    matrix.loc <- list.files(run, pattern = 'matrix.mtx', full.names = TRUE)

    data <- readMM(file = matrix.loc) %>% as("dgCMatrix")
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)
    if (!missing(suffix)) {
        cell.names %<>% paste(suffix, sep = "_")
    }
    
    gene.names <- fread(gene.loc, header = FALSE)$V2
    row.names(data) <- gene.names
    colnames(data) <- cell.names

    
    return(as(data, "dgCMatrix"))
#     return(as(sumOverRowNames(data), "dgCMatrix"))
}
