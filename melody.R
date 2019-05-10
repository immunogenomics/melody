melody <- function(exprs, meta_data, R, R_thresh=0.05, family='poisson') {
    models <- get_model_pars_mix(exprs, 'y ~ depth', meta_data, R, R_thresh, family)
    residuals <- compute_residuals(exprs, models, meta_data)
    exprs_corr <- correct_mix(residuals, models, meta_data, R)
    return(exprs_corr)
}


get_model_pars_mix <- function (umi, model_str, model_data, R, R_thresh=.05, family='poisson') { 
    apply(R, 1, function(cluster_probs) {        
        idx <- which(cluster_probs >= R_thresh)
        if (family == 'lm') {
            message('TODO: design matrix is hard coded right now. Fix this in next iteration.')
            y <- umi[, idx] %>% as.matrix %>% t %>% log1p
            fit <- lm(formula = y ~ model_data$depth[idx], weights = cluster_probs[idx]) 
            model_pars <- data.frame(t(fit$coefficients))
        } else {
            model_pars <- Reduce(rbind, apply(umi, 1, function(y) {
                y <- y[idx]

                ## TODO: replace y=0 response with appropriately sized vector
                if (family == 'poisson') {
                    fit <- glm(
                        as.formula(model_str), 
                        data = model_data[idx, ], 
                        family = poisson, 
                        weights = cluster_probs[idx]
                    )   
                } else if (family == 'gaussian') {
                    y <- log1p(y)
                    fit <- glm(
                        as.formula(model_str), 
                        data = model_data[idx, ], 
                        family = gaussian, 
                        weights = cluster_probs[idx]
                    )                
                } 
                return(data.frame(t(fit$coefficients)))
            }))
        }
        
        row.names(model_pars) <- row.names(umi)
        return(model_pars)
    })
}

compute_residuals <- function (exprs, models, cell_attr) {  
    design <- cbind(1, data.frame(cell_attr)[, c('depth')])
    lapply(models, function(model) {
        mu <- exp(tcrossprod(as.matrix(model), design))
        exprs - mu
    })
}

## add expected values to residuals
correct_mix <- function (residuals, models, cell_attr, R, force_pos = TRUE) {  
    ## TODO: (quantile) cut cells with outlier depth. Mean is sensitive to these
#     quantile(meta_data$depth, c(.01, .99))
    
    message('TODO: when correcting counts, do expm1 since you already used logp1')
    
    design <- cbind(1, data.frame(cell_attr)[, c('depth')])
    if (nrow(R) == 1) {
        expected_depth <- tcrossprod(cell_attr$depth, R) / sum(R)
    } else {
        expected_depth <- tcrossprod(cell_attr$depth, R) %*% diag(1 / rowSums(R))        
    }
    
    exprs_corrected <- Reduce(`+`, lapply(seq_len(nrow(R)), function(k) {
        ## Set depth to the cluster specific depth 
        ## CAUTION: assuming that second column is depth
        design[, 2] <- rep(expected_depth[k], nrow(design))
        mu <- exp(tcrossprod(as.matrix(models[[k]]), design))
        exprs_cluster <- mu + residuals[[k]]
        sweep(exprs_cluster, 2, R[k, ], '*')
    }))
    
    if (force_pos) {
        exprs_corrected[exprs_corrected < 0] <- 0
    }
    return(exprs_corrected)
}




model.matrix.full <- function(vars_use, data_df, intercept = 1, 
                              numerical_mode = c('none', 'scale', 'bin')[1],
                              binning.method = "equal_frequency", num.bin = 10) {    
    discrete_vars <- c()
    for (varname in vars_use) {
        if ("character" %in% class(data_df[[varname]])) {
            data_df[[varname]] <- factor(data_df[[varname]])
            discrete_vars <- c(discrete_vars, varname)
        } else if ("factor" %in% class(data_df[[varname]])) {
            discrete_vars <- c(discrete_vars, varname)
        } else {
            ## handle continuous covariates
            if (numerical_mode == 'bin') {
                ## METHOD 1: bin it into categorical variable
                data_df[[varname]] <- factor(do_binning(data_df[[varname]], binning.method, num.bin))
                discrete_vars <- c(discrete_vars, varname)
            } else if (numerical_mode == 'scale') {
                ## METHOD 2: transform it into z-score 
                data_df[[varname]] <- scale(data_df[[varname]], center = TRUE, scale = TRUE)
            } else {
                ## METHOD 3: keep it as is
                data_df[[varname]] <- data_df[[varname]]                
            }
        }
    }
    
    
    if (length(discrete_vars) > 0) {
        contrasts_df <- lapply(discrete_vars, function(x) {
            diag(nlevels(data.frame(data_df)[, x, drop = TRUE]))
        })

        names(contrasts_df) <- discrete_vars
        res <- model.matrix(as.formula(sprintf("~ %d + %s", intercept, paste(vars_use, collapse = "+"))), 
                     data=data_df, contrasts.arg=contrasts_df)    
        
    } else {
        res <- model.matrix(as.formula(sprintf("~ %d + %s", intercept, paste(vars_use, collapse = "+"))), 
                     data=data_df)        
    }
    if (intercept) {
        ## don't keep parentheses around intercept name
        colnames(res) <- c('Intercept', tail(colnames(res), -1))
    }
    return(res)
}

do_binning <- function(values, binning.method = "equal_frequency", num.bin = 10) {
    if (binning.method == "equal_width") {
        .breaks <- num.bin
    }
    else if (binning.method == "equal_frequency") {
        .breaks <- quantile(values, probs = seq(0, 1, length.out = num.bin + 1))
    }
    else {
        stop(paste0("Invalid selection: '", binning.method, "' for 'binning.method'."))
    }
    cut(values, .breaks, include.lowest = TRUE) %>% as.integer
}
