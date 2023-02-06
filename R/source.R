.add_predictor <-
function (mri_data, data_table, pred_transf) 
{
    out = mri_data
    for (i in 1:length(pred_transf)) {
        pred_transf_i = pred_transf[[i]]
        data_col = data_table[, match(pred_transf_i[1], colnames(data_table))]
        out[[length(out) + 1]] = as.numeric(switch(pred_transf_i[2], 
            factor = data_col == pred_transf_i[3], numeric = data_col))
    }
    out
}
.aprior <-
function (gamma.hat) 
{
    m <- mean(gamma.hat)
    s2 <- var(gamma.hat)
    (2 * s2 + m^2)/s2
}
.assign.folds_one_site <-
function (y, family, nfolds) 
{
    folds = rep(NA, length(y))
    folds = switch(family, binomial = {
        for (i in 0:1) {
            indexs_i = which(y == i)
            folds[indexs_i] = .assign.folds_simple(indexs_i, 
                nfolds)
        }
        folds
    }, cox = {
        sorted_idx = sort.int(y[, 1], index.return = TRUE)$ix
        y_w_indx = cbind(y, sorted_idx)
        idx_g0 = y_w_indx[y[, 2] == 0, 3]
        idx_g1 = y_w_indx[y[, 2] == 1, 3]
        n_g0 = length(idx_g0)
        if (n_g0 > 0) {
            for (i in seq(0, n_g0 - 1, by = nfolds)) {
                n_i = min(c(nfolds, n_g0 - i))
                indexs_i = idx_g0[i + (1:n_i)]
                folds[indexs_i] = sample(n_i)
            }
        }
        n_g1 = length(idx_g1)
        if (n_g1 > 0) {
            for (i in seq(0, n_g1 - 1, by = nfolds)) {
                n_i = min(c(nfolds, n_g1 - i))
                indexs_i = idx_g1[i + (1:n_i)]
                folds[indexs_i] = sample(n_i)
            }
        }
        y_w_indx = data.frame(time = y_w_indx[, 1], status = y_w_indx[, 
            2], sorted_idx = y_w_indx[, 3], folds = folds)
        y_w_indx_original_order = y_w_indx[order(sorted_idx), 
            ]
        folds = y_w_indx_original_order$folds
        folds
    }, gaussian = {
        sorted_idx = order(y)
        n = length(y)
        for (i in seq(0, n - 1, by = nfolds)) {
            n_i = min(c(nfolds, n - i))
            indexs_i = sorted_idx[i + (1:n_i)]
            folds[indexs_i] = .assign.folds_simple(indexs_i, 
                nfolds)
        }
        table_y = table(y)
        if (length(table_y) > 2 || (length(table_y) == 2 && all(table_y) >= 
            2)) {
            fold_with_constant_training_sample = NA
            for (i in 1:nfolds) {
                if (var(y[which(folds != i)]) == 0) {
                  fold_with_constant_training_sample = i
                }
            }
            if (!is.na(fold_with_constant_training_sample)) {
                index_of_one_repeated_y_from_other_folds = sample(which(folds != 
                  fold_with_constant_training_sample), 1)
                index_of_one_different_y_from_the_fold = sample(which(y != 
                  y[index_of_one_repeated_y_from_other_folds]), 
                  1)
                tmp = folds[index_of_one_different_y_from_the_fold]
                folds[index_of_one_different_y_from_the_fold] = folds[index_of_one_repeated_y_from_other_folds]
                folds[index_of_one_repeated_y_from_other_folds] = tmp
            }
        }
        folds
    })
    1 + (folds + sample(1:nfolds, 1))%%nfolds
}
.assign.folds_simple <-
function (y, nfolds) 
{
    n = length(y)
    if (n < nfolds) {
        folds = sample(1:nfolds, n)
    }
    else {
        folds = sample(1 + (1:n%%nfolds))
    }
    1 + (folds + sample(1:nfolds, 1))%%nfolds
}
.bprior <-
function (gamma.hat) 
{
    m <- mean(gamma.hat)
    s2 <- var(gamma.hat)
    (m * s2 + m^3)/s2
}
.calculate_ensemble_accuracy <-
function (cv_results, family) 
{
    acc_per_fold = c()
    for (i in 1:18) {
        fold_results = mp$cv_results[mp$cv_results$iteration == 
            i, ]
        if (family == "binomial") {
            print(.metrics_binary(real = fold_results$real, predictions = fold_results$pred > 
                0.5))
        }
        else if (family == "gaussian") {
            print(sqrt(mean((fold_results$pred - fold_results$real)^2)))
        }
        else {
        }
    }
}
.check_response_var <-
function (data_table, response_var, response_settings) 
{
    response_col = match(response_var, colnames(data_table))
    if (anyNA(response_col)) {
        stop(paste("Response variable", response_var, "not found"))
    }
    data_table[, response_col]
}
.check_stop <-
function (message) 
{
    stop(paste("\n", message, sep = ""))
}
.cite_mripredict <-
function () 
{
    cat("Please cite this software as:\n\n")
    cat("Solanes A, Mezquida G, Janssen J, Amoretti S, Lobo A, González-Pinto A, Arango C, Vieta E, Castro-Fornieles J, Bergé D, Albacete A, Giné E, Parellada M, Bernardo M; PEPs group (collaborators), Pomarol-Clotet E, Radua J. Combining MRI and clinical data to detect high relapse risk after the first episode of psychosis. Schizophrenia (Heidelb). 2022 Nov 17;8(1):100. doi: 10.1038/s41537-022-00309-w. PMID: 36396933; PMCID: PMC9672064.\n")
}
.combat_tmp1 <-
function (dat, batch, levels_batch, mod) 
{
    batchmod <- model.matrix(~-1 + batch)
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels_batch[i])
    }
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    batch.design <- design[, 1:n.batch]
    return(list(dat = dat, batchmod = batchmod, n.batch = n.batch, 
        batches = batches, n.batches = n.batches, n.array = n.array, 
        design = design, batch.design = batch.design))
}
.combat_tmp2 <-
function (tmp1, verbose = TRUE) 
{
    if (verbose) {
        cat("[combat] Adjusting for", ncol(tmp1$design) - ncol(tmp1$batchmod), 
            "covariate(s) or covariate level(s)\n")
    }
    if (qr(tmp1$design)$rank < ncol(tmp1$design)) {
        if (ncol(tmp1$design) == (tmp1$n.batch + 1)) {
            stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
        }
        if (ncol(tmp1$design) > (tmp1$n.batch + 1)) {
            if ((qr(tmp1$design[, -c(1:tmp1$n.batch)])$rank < 
                ncol(tmp1$design[, -c(1:tmp1$n.batch)]))) {
                stop("The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.")
            }
            else {
                stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
            }
        }
    }
    B.hat <- solve(t(tmp1$design) %*% tmp1$design) %*% t(tmp1$design) %*% 
        t(as.matrix(tmp1$dat))
    grand.mean <- t(tmp1$n.batches/tmp1$n.array) %*% B.hat[1:tmp1$n.batch, 
        ]
    var.pooled <- ((tmp1$dat - t(tmp1$design %*% B.hat))^2) %*% 
        rep(1/tmp1$n.array, tmp1$n.array)
    return(list(B.hat = B.hat, grand.mean = grand.mean, var.pooled = var.pooled))
}
.combat_tmp3 <-
function (dat, tmp1, tmp2, verbose = TRUE) 
{
    if (verbose) {
        cat("[combat] Standardizing data across features\n")
    }
    stand.mean <- t(tmp2$grand.mean) %*% t(rep(1, tmp1$n.array))
    if (!is.null(tmp1$design)) {
        tmp <- tmp1$design
        tmp[, c(1:tmp1$n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% tmp2$B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(tmp2$var.pooled) %*% t(rep(1, 
        tmp1$n.array)))
    return(list(stand.mean = stand.mean, s.data = s.data))
}
.combat_tmp4 <-
function (tmp1, tmp2, tmp3, eb = TRUE, verbose = TRUE) 
{
    if (eb) {
        if (verbose) {
            cat("[combat] Fitting L/S model and finding priors\n")
        }
    }
    else {
        if (verbose) {
            cat("[combat] Fitting L/S model\n")
        }
    }
    gamma.hat <- solve(t(tmp1$batch.design) %*% tmp1$batch.design) %*% 
        t(tmp1$batch.design) %*% t(as.matrix(tmp3$s.data))
    delta.hat <- NULL
    for (i in tmp1$batches) {
        delta.hat <- rbind(delta.hat, apply(tmp3$s.data[, i], 
            1, var, na.rm = T))
    }
    gamma.star <- delta.star <- NULL
    gamma.bar <- t2 <- a.prior <- b.prior <- NULL
    if (eb) {
        gamma.bar <- apply(gamma.hat, 1, mean)
        t2 <- apply(gamma.hat, 1, var)
        a.prior <- apply(delta.hat, 1, .aprior)
        b.prior <- apply(delta.hat, 1, .bprior)
        if (verbose) {
            cat("[combat] Finding parametric adjustments\n")
        }
        for (i in 1:tmp1$n.batch) {
            temp <- .it.sol(tmp3$s.data[, tmp1$batches[[i]]], 
                gamma.hat[i, ], delta.hat[i, ], gamma.bar[i], 
                t2[i], a.prior[i], b.prior[i])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    return(list(gamma.hat = gamma.hat, delta.hat = delta.hat, 
        gamma.star = gamma.star, delta.star = delta.star, gamma.bar = gamma.bar, 
        t2 = t2, a.prior = a.prior, b.prior = b.prior))
}
.combat_tmp5 <-
function (tmp1, tmp2, tmp3, tmp4, eb = TRUE, verbose = TRUE) 
{
    if (verbose) {
        cat("[combat] Adjusting the data\n")
    }
    bayesdata <- tmp3$s.data
    j <- 1
    for (i in tmp1$batches) {
        if (eb) {
            bayesdata[, i] <- (bayesdata[, i] - t(tmp1$batch.design[i, 
                ] %*% tmp4$gamma.star))/(sqrt(tmp4$delta.star[j, 
                ]) %*% t(rep(1, tmp1$n.batches[j])))
        }
        else {
            bayesdata[, i] <- (bayesdata[, i] - t(tmp1$batch.design[i, 
                ] %*% tmp4$gamma.hat))/(sqrt(tmp4$delta.hat[j, 
                ]) %*% t(rep(1, tmp1$n.batches[j])))
        }
        j <- j + 1
    }
    return((bayesdata * (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))) + 
        tmp3$stand.mean)
}
.correc <-
function (i, n) 
{
    c1 = c(9.5, 28.699999999999999, 1.8999999999999999, 0, -7, 
        -6.2000000000000002, -1.6000000000000001)
    c2 = c(-6195, -9569, -6728, -17614, -8278, -3570, 1075)
    c3 = c(93380, 175160, 410400, 2157600, 2376000, 2065000, 
        2065000)
    mic = 9.9999999999999995e-07
    c14 = 1.9000000000000001e-05
    if (i * n == 4) {
        return(c14)
    }
    if (i < 1 || i > 7) {
        return(0)
    }
    if (i != 4 && n > 20) {
        return(0)
    }
    if (i == 4 && n > 40) {
        return(0)
    }
    n = 1/n^2
    i = i - 1
    out = (c1[i] + n * (c2[i] + n * c3[i])) * mic
}
.coxph_RD2 <-
function (predictor, stime, sevent) 
{
    surv_obj <- Surv(stime, sevent, type = "right")
    coxphfit <- summary(coxph(surv_obj ~ predictor))
    b <- coxphfit$coefficients[1]
    PI = sort.int((predictor - mean(predictor)) * b, index.return = TRUE)
    rankits = .normOrder(length(PI$x))
    kappa = sqrt(8/pi)
    rankits = rankits/kappa
    surv_obj <- Surv(stime[PI$ix], sevent[PI$ix], type = "right")
    coxphfit <- summary(coxph(surv_obj ~ rankits))
    b <- coxphfit$coefficients[1]
    D = exp(b)
    sigma2 = pi^2/6
    RD2 = (D^2/kappa^2)/(sigma2 + D^2/kappa^2)
    list(RD2 = RD2, b = b)
}
.create_covX <-
function (data_table, covX_transf, nrows = nrow(data_table)) 
{
    covX = matrix(1, nrows)
    colnames_covX = c()
    if (length(covX_transf) > 0) {
        for (i in 1:length(covX_transf)) {
            covX_transf_i = covX_transf[[i]]
            data_col = data_table[, match(covX_transf_i[1], colnames(data_table))]
            covX = cbind(covX, as.numeric(switch(covX_transf_i[2], 
                factor = data_col == covX_transf_i[3], numeric = data_col)))
            colnames_covX = c(colnames_covX, switch(covX_transf_i[2], 
                factor = sprintf("%s_%s", covX_transf_i[1], covX_transf_i[3]), 
                numeric = covX_transf_i[1]))
        }
    }
    colnames(covX) <- c("intercept", colnames_covX)
    covX
}
.create_covX_old <-
function (data_table, covX_transf, nrows = nrow(data_table)) 
{
    covX = matrix(1, nrows)
    colnames_covX = c()
    if (length(covX_transf) > 0) {
        for (i in 1:length(covX_transf)) {
            covX_transf_i = covX_transf[[i]]
            data_col = data_table[, match(covX_transf_i[1], colnames(data_table))]
            covX = cbind(covX, as.numeric(switch(covX_transf_i[2], 
                factor = data_col == covX_transf_i[3], numeric = data_col)))
            colnames_covX = c(colnames_covX, switch(covX_transf_i[2], 
                factor = sprintf("%s_%s", covX_transf_i[1], covX_transf_i[3]), 
                numeric = covX_transf_i[1]))
        }
    }
    colnames(covX) <- c("intercept", colnames_covX)
    covX
}
.create_undummy_table <-
function (x, transf) 
{
    dummy_colnames <- sapply(strsplit(colnames(x), ":"), function(tmp) {
        tmp[1]
    })
    unique_dummy_colnames <- unique(dummy_colnames)
    xp <- setNames(data.frame(matrix(data = NA, ncol = length(unique_dummy_colnames), 
        nrow = nrow(x))), unique_dummy_colnames)
    dummy_values <- sapply(strsplit(colnames(x), ":"), function(tmp) {
        tmp[2]
    })
    for (i in 1:length(dummy_values)) {
        if (is.na(dummy_values[i])) {
            xp[, dummy_colnames[i]] <- x[, i]
        }
        else {
            xp[which(x[, i] == 1), dummy_colnames[i]] <- dummy_values[i]
        }
    }
    for (i in 1:nrow(transf$factor_ref)) {
        xp[which(is.na(xp[, transf$factor_ref[i, 1]])), transf$factor_ref[i, 
            1]] <- transf$factor_ref[i, 2]
    }
    xp
}
.create_Y <-
function (data_table, response_var, response_event) 
{
    out = c()
    if (length(unique(data_table[, response_var])) > 2) {
        out = as.numeric(data_table[, response_var])
    }
    else {
        out = as.numeric(data_table[, match(response_var, colnames(data_table))] == 
            response_event)
    }
    out
}
.find_best_time <-
function (time, status) 
{
    best.x = NA
    best.err2 = Inf
    for (x in unique(sort(time))) {
        err2 = (sum(time < x & status == 1) - sum(time >= x))^2
        if (err2 < best.err2) {
            best.x = x
            best.err2 = err2
        }
    }
    best.x
}
.get_predictor_matrix <-
function (data_table, pred_transf) 
{
    out = c()
    for (i in 1:length(pred_transf)) {
        pred_transf_i = pred_transf[[i]]
        data_col = data_table[, match(pred_transf_i[1], colnames(data_table))]
        out = c(out, as.numeric(switch(pred_transf_i[2], factor = data_col == 
            pred_transf_i[3], numeric = data_col)))
    }
    matrix(out, ncol = length(pred_transf))
}
.indx_equivalent_columns_from_transf <-
function (data_table_transf, columns, has_intercept = TRUE) 
{
    out = c()
    for (i in 1:length(data_table_transf)) {
        if (data_table_transf[[i]][2] == "factor") {
            if (!is.null(data_table_transf) && sprintf("%s_%s", 
                data_table_transf[[i]][1], data_table_transf[[i]][3]) %in% 
                columns) {
                out = c(out, i)
            }
        }
        else {
            if (!is.null(data_table_transf) && data_table_transf[[i]][1] %in% 
                columns) {
                out = c(out, i)
            }
        }
    }
    if (has_intercept) {
        out = out + 1
    }
    out
}
.it.sol <-
function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) 
{
    n <- apply(!is.na(sdat), 1, sum)
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while (change > conv) {
        g.new <- .postmean(g.hat, g.bar, n, d.old, t2)
        sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
            1, sum, na.rm = T)
        d.new <- .postvar(sum2, n, a, b)
        change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
        g.old <- g.new
        d.old <- d.new
        count <- count + 1
    }
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g.star", "d.star")
    adjust
}
.load_mri <-
function (mri_paths, mask = NULL, space = "MNI") 
{
    mri = c()
    if (is.numeric(mri_paths) && length(dim(mri_paths)) == 4) {
        mri$data <- mri_paths
        mri$n <- dim(mri_paths)[4]
        if (!is.null(mask)) {
            for (i in 1:mri$n) {
                if (length(dim(mask)) == 4) {
                  mask = mask[, , , 1]
                }
                mri$data[, , , i] = mri$data[, , , i] * mask
            }
        }
    }
    else {
        n = length(mri_paths)
        if (!n) {
            stop("No MRI data have been specified")
        }
        mri <- .read_mri(as.character(mri_paths[1]), read_data = 0)
        dim <- mri$dim
        mri$data <- array(numeric(), c(dim, n))
        mat <- mri$sto.xyz
        mri$n = n
        for (i in 1:n) {
            mri_i = .read_mri(as.character(mri_paths[i]))
            if (length(mri_i$dim) != 3) {
                stop(paste(mri_paths[i], "should be 3D"))
            }
            if (any(mri_i$dim != dim)) {
                stop(paste(mri_paths[i], "should have the same dimensions than other MRI data"))
            }
            if ((mri_i$sform.code != "NIFTI.XFORM.MNI.152") && 
                (space != "NO_CHECK")) {
                stop(paste(mri_paths[i], "should be in MNI space"))
            }
            if (any(mri_i$sto.xyz != mat)) {
                stop(paste(mri_paths[i], "should have the same MNI transformation matrix than other MRI data"))
            }
            if (mri_i$scl.slope != 0) {
                if (length(mask) == 0) {
                  mri$data[, , , i] = mri_i$scl.inter + mri_i$scl.slope * 
                    mri_i$data
                }
                else {
                  mri_tmp = mri_i$scl.inter + mri_i$scl.slope * 
                    mri_i$data
                  mri$data[, , , i] = mri_tmp * mask[, , , 1]
                }
            }
            else {
                if (length(mask) == 0) {
                  mri$data[, , , i] = mri_i$data
                }
                else {
                  mri_tmp = mri_i$scl.inter + mri_i$scl.slope * 
                    mri_i$data
                  mri$data[, , , i] = mri_tmp * mask[, , , 1]
                }
            }
        }
    }
    class(mri) = "mripredict_data"
    mri
}
.load_mri_old <-
function (mri_paths, space = "MNI") 
{
    n = length(mri_paths)
    if (!n) {
        stop("No MRI data have been specified")
    }
    mri <- .read_mri(mri_paths[1], read_data = 0)
    dim <- mri$dim
    mat <- mri$sto.xyz
    mri$n = n
    mri$data <- array(numeric(), c(dim, n))
    for (i in 1:n) {
        mri_i = .read_mri(mri_paths[i])
        if (length(mri_i$dim) != 3) {
            stop(paste(mri_paths[i], "should be 3D"))
        }
        if (any(mri_i$dim != dim)) {
            stop(paste(mri_paths[i], "should have the same dimensions than other MRI data"))
        }
        if ((mri_i$sform.code != "NIFTI.XFORM.MNI.152") && (space != 
            "NO_CHECK")) {
            stop(paste(mri_paths[i], "should be in MNI space"))
        }
        if (any(mri_i$sto.xyz != mat)) {
            stop(paste(mri_paths[i], "should have the same MNI transformation matrix than other MRI data"))
        }
        if (mri_i$scl.slope || (mri_i$scl.inter || mri_i$scl.slope != 
            1)) {
            mri$data[, , , i] = mri_i$scl.inter + mri_i$scl.slope * 
                mri_i$data
        }
        else {
            mri$data[, , , i] = mri_i$data
        }
    }
    mri
}
.mask_from_which <-
function (which_array, mask_dim) 
{
    mask = array(F, mask_dim)
    mask[which_array] = T
}
.mean_RD2 <-
function (betes) 
{
    betes_mean = mean(betes)
    D = exp(betes_mean)
    kappa = sqrt(8/pi)
    sigma2 = pi^2/6
    RD2 = (D^2/kappa^2)/(sigma2 + D^2/kappa^2)
    RD2
}
.metrics_binary <-
function (real, predictions, folder = "", save = FALSE) 
{
    classlabels = c(1, 0)
    if (length(classlabels) > 2) {
        stop("Only binary classification supported")
    }
    class1 = classlabels[1]
    class2 = classlabels[2]
    tp = sum(predictions == class1 & real == class1)
    fp = sum(predictions == class1 & real == class2)
    fn = sum(predictions == class2 & real == class1)
    tn = sum(predictions == class2 & real == class2)
    acc_class1 = tp/(tp + fn)
    acc_class2 = tn/(tn + fp)
    bac = 0.5 * (acc_class1 + acc_class2)
    metric <- data.frame(class1_label = class1, class2_label = class2, 
        n_class1 = sum(real == class1), n_class2 = sum(real == 
            class2), tp = tp, fp = fp, fn = fn, tn = tn, ppv = round(tp/(tp + 
            fp), 3), npv = round(tn/(tn + fn), 3), sensitivity = round(acc_class1, 
            3), specificity = round(acc_class2, 3), bac = round(bac, 
            3))
    if (save) 
        write.csv(metric, sprintf("%s/binary_results_fold%d.csv", 
            folder))
    metric
}
.metrics_cox <-
function (results, frontier_time, iteration = 1, folder = "", 
    save = TRUE) 
{
    sorted_results = results[order(results$linear_predictor), 
        ]
    sorted_results_desc = results[order(results$linear_predictor, 
        decreasing = TRUE), ]
    idx_g1 = which(results$time < frontier_time & !is.na(results$linear_predictor) & 
        results$status == 1)
    idx_g2 = which(results$time >= frontier_time & !is.na(results$linear_predictor))
    if (length(idx_g1) > 0 && length(idx_g2) > 0) {
        g1 = cbind(results[idx_g1, ], risk = 1)
        g2 = cbind(results[idx_g2, ], risk = 0)
        g1_th_linPred = sorted_results_desc$linear_predictor[dim(g1)[1]]
        g2_th_linPred = sorted_results$linear_predictor[dim(g2)[1]]
        th_linPred = (g1_th_linPred + g2_th_linPred)/2
        g1_predicted = g1$linear_predictor >= th_linPred
        g2_predicted = !(g2$linear_predictor < th_linPred)
        g12 = rbind(g1, g2)
        g12_predicted = c(g1_predicted, g2_predicted)
        perf = .metrics_binary(g12[, ncol(g12)], g12_predicted)
    }
    else {
        message("Warning: not enough samples with status equal to 1 to calculate the performance of the model.")
        perf = NULL
    }
    results$iteration = iteration
    if (save) {
        if (!is.null(results$fold)) 
            write.csv(results, sprintf("%s_cox_results_fold%d.csv", 
                folder, results$fold[1]), row.names = FALSE)
        else write.csv(results, sprintf("%s_cox_results.csv", 
            folder), row.names = FALSE)
    }
    perf
}
.most_frequent_variables <-
function (model_list, mp = NULL, file = NULL) 
{
    mnis = list(gm = NULL, gm_mod = NULL, wm = NULL, wm_mod = NULL)
    out = data.frame(mni_or_variable = NULL, ijk = NULL, modality = NULL, 
        beta = NULL)
    for (i in seq_len(length(model_list))) {
        tmp <- data.frame(mni_or_variable = NULL, modality = NULL, 
            betas = NULL)
        if (length(model_list[[i]]$lasso_mni$gm) > 0) {
            model_list[[i]]$lasso_mni$gm <- matrix(model_list[[i]]$lasso_mni$gm, 
                nrow = 3)
            gm <- paste(model_list[[i]]$lasso_mni$gm[1, ], model_list[[i]]$lasso_mni$gm[2, 
                ], model_list[[i]]$lasso_mni$gm[3, ], sep = "_")
            mnis$gm <- c(mnis$gm, gm)
            tmp <- rbind(data.frame(mni_or_variable = gm, modality = "gm", 
                betas = model_list[[i]]$lasso_mni$gm_betas))
        }
        if (length(model_list[[i]]$lasso_mni$gm_mod) > 0) {
            model_list[[i]]$lasso_mni$gm_mod <- matrix(model_list[[i]]$lasso_mni$gm_mod, 
                nrow = 3)
            gm_mod <- paste(model_list[[i]]$lasso_mni$gm_mod[1, 
                ], model_list[[i]]$lasso_mni$gm_mod[2, ], model_list[[i]]$lasso_mni$gm_mod[3, 
                ], sep = "_")
            mnis$gm_mod <- c(mnis$gm_mod, gm_mod)
            tmp <- rbind(tmp, data.frame(mni_or_variable = gm_mod, 
                modality = "gm_mod", betas = model_list[[i]]$lasso_mni$gm_mod_betas))
        }
        if (length(model_list[[i]]$lasso_mni$wm) > 0) {
            model_list[[i]]$lasso_mni$wm <- matrix(model_list[[i]]$lasso_mni$wm, 
                nrow = 3)
            wm <- paste(model_list[[i]]$lasso_mni$wm[1, ], model_list[[i]]$lasso_mni$wm[2, 
                ], model_list[[i]]$lasso_mni$wm[3, ], sep = "_")
            mnis$wm <- c(mnis$wm, wm)
            tmp <- rbind(tmp, data.frame(mni_or_variable = wm, 
                modality = "wm", betas = model_list[[i]]$lasso_mni$wm_betas))
        }
        if (length(model_list[[i]]$lasso_mni$wm_mod) > 0) {
            model_list[[i]]$lasso_mni$wm_mod <- matrix(model_list[[i]]$lasso_mni$wm_mod, 
                nrow = 3)
            wm_mod <- paste(model_list[[i]]$lasso_mni$wm_mod[1, 
                ], model_list[[i]]$lasso_mni$wm_mod[2, ], model_list[[i]]$lasso_mni$wm_mod[3, 
                ], sep = "_")
            mnis$wm_mod <- c(mnis$wm_mod, wm_mod)
            tmp <- rbind(tmp, data.frame(mni_or_variable = wm_mod, 
                modality = "wm_mod", betas = model_list[[i]]$lasso_mni$wm_mod_betas))
        }
        if (nrow(tmp) > 0) 
            tmp <- cbind(tmp, data.frame(fold = i))
        if (length(model_list[[i]]$lasso_predX_indx) > 0) {
            predictor_vars = mp$pred_var
            betas_i = model_list[[i]]$lasso$i > model_list[[i]]$n_voxels_mask
            tmp <- rbind(tmp, data.frame(mni_or_variable = predictor_vars[model_list[[i]]$lasso_predX_indx], 
                modality = "predictor_variable", betas = model_list[[i]]$lasso$beta[betas_i], 
                fold = i))
        }
        if (mp$response_family != "cox") 
            tmp <- rbind(data.frame(mni_or_variable = "intercept", 
                modality = "intercept", betas = model_list[[i]]$lasso$a0, 
                fold = i), tmp)
        out <- rbind(out, tmp)
    }
    out$name = paste0(out$mni_or_variable, ":", out$modality)
    betes = out
    resum_betes = do.call(data.frame, aggregate(betes$betas ~ 
        betes$name, data = betes, FUN = function(x) c(sum = sum(x), 
        n = length(x))))
    colnames(resum_betes) = c("variable", "beta_sum", "beta_count")
    if (!is.null(file)) {
        write.csv(resum_betes, file = gsub(".csv", "_summary.csv", 
            file))
        write.csv(out, file = file, row.names = F)
    }
    rownames(out) <- NULL
    resum_betes
}
.most_frequent_variables_ijk <-
function (model_list, mp = NULL, file = NULL) 
{
    ijks = c()
    out = data.frame(mni_or_variable = NULL, ijk = NULL, modality = NULL, 
        beta = NULL, covB_intercept = NULL, covB_age = NULL, 
        covB_sex = NULL)
    for (i in seq_len(length(model_list))) {
        tmp <- data.frame(mni_or_variable = NULL, modality = NULL)
        if (length(model_list[[i]]$lasso_ijk) > 0) {
            model_list[[i]]$lasso_ijk <- matrix(model_list[[i]]$lasso_ijk, 
                nrow = 3)
            coords <- paste(model_list[[i]]$lasso_ijk[1, ], model_list[[i]]$lasso_ijk[2, 
                ], model_list[[i]]$lasso_ijk[3, ], sep = "_")
            ijks <- c(ijks, coords)
            tmp <- data.frame(mni_or_variable = coords, modality = "voxel")
        }
        if (nrow(tmp) > 0) 
            model_list[[i]]$lasso_covB = as.matrix(model_list[[i]]$lasso_covB)
        if (length(model_list[[i]]$lasso_covB[1, ]) > 0) {
            tmp <- cbind(tmp, data.frame(betas = model_list[[i]]$lasso$beta[c(1:nrow(tmp))], 
                fold = i, covB_intercept = model_list[[i]]$lasso_covB[1, 
                  ], covB_age = model_list[[i]]$lasso_covB[2, 
                  ], covB_sex = model_list[[i]]$lasso_covB[3, 
                  ]))
        }
        if (length(model_list[[i]]$lasso_predX_indx) > 0) {
            predictor_vars = mp$pred_var
            betas_i = model_list[[i]]$lasso$i > model_list[[i]]$n_voxels_mask
            tmp <- rbind(tmp, data.frame(mni_or_variable = predictor_vars[model_list[[i]]$lasso_predX_indx], 
                modality = "predictor_variable", betas = model_list[[i]]$lasso$beta[betas_i], 
                fold = i, covB_intercept = NA, covB_age = NA, 
                covB_sex = NA))
        }
        if (mp$response_family != "cox") 
            tmp <- rbind(data.frame(mni_or_variable = "intercept", 
                modality = "intercept", betas = model_list[[i]]$lasso$a0, 
                fold = i), tmp)
        out <- rbind(out, tmp)
    }
    if (!is.null(file)) 
        write.csv(out, file = file)
    out
}
.normOrder <-
function (N) 
{
    if (length(N) > 1) {
        n = length(N)
    }
    else {
        n = N
    }
    s = N/2
    n2 = floor(s)
    eps = c(0.41988500000000001, 0.45053599999999999, 0.45693600000000001, 
        0.46848800000000002)
    dl1 = c(0.112063, 0.12177, 0.23929900000000001, 0.21515899999999999)
    dl2 = c(0.080121999999999999, 0.111348, -0.211867, -0.115049)
    gam = c(0.474798, 0.469051, 0.208597, 0.25978400000000001)
    lam = c(0.28276499999999999, 0.30485600000000002, 0.40770800000000001, 
        0.41409299999999999)
    bb = -0.283833
    d = -0.10613599999999999
    b1 = 0.56418959999999996
    s[1] = b1
    k = 3
    len_s = k
    if (n2 < k) {
        k = n2
    }
    s = vector(mode = "numeric", length = n2)
    for (i in 1:k) {
        e1 = (i - eps[i])/(n + gam[i])
        e2 = e1^(lam[i])
        s[i] = e1 + e2 * (dl1[i] + e2 * dl2[i])/n - .correc(i + 
            1, n)
    }
    if (n2 > k) {
        for (i in 4:n2) {
            e1 = (i - eps[4])/(n + gam[4])
            e2 = e1^(lam[4] + bb/(i + d))
            s[i] = e1 + e2 * (dl1[4] + e2 * dl2[4])/n - .correc(i + 
                1, n)
        }
    }
    for (i in 1:n2) {
        s[i] = -qnorm(s[i], 0, 1)
    }
    if (n%%2 == 0) {
        out = c(-s, s[length(s):1])
    }
    else {
        out = c(-s, 0, s[length(s):1])
    }
    out
}
.postmean <-
function (g.hat, g.bar, n, d.star, t2) 
{
    (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}
.postvar <-
function (sum2, n, a, b) 
{
    (0.5 * sum2 + b)/(n/2 + a - 1)
}
.print_action <-
function (message) 
{
    cat(message, "... ", sep = "")
    flush.console()
}
.print_metrics <-
function (linPred, family, frontier_time = 0, testY, id) 
{
    cat("\n")
    pred = c()
    if (family == "binomial") {
        prob = 1/(1 + exp(-linPred))
        pred = prob
        bac = .metrics_binary(testY, pred > mean(testY))$bac
        results = data.frame(id = id, linear_predictor = pred[, 
            1], real = testY[, 1])
        print(results)
    }
    else if (family == "gaussian") {
        pred = as.matrix(linPred)
        results = data.frame(id = id, linear_predictor = pred[, 
            1], real = testY[, 1])
        print(results)
    }
    else if (family == "cox") {
        pred = linPred
        rd2 = .coxph_RD2(pred, testY[, 1], testY[, 2])
        results = data.frame(id = id, linear_predictor = pred[, 
            1], time = testY[, 1], status = testY[, 2])
        print(results)
        print(.metrics_cox(results, frontier_time, rd2$b, save = FALSE))
    }
}
.print_ok <-
function (prefix = "") 
{
    cat(prefix, "Ok\n", sep = "")
    flush.console()
}
.read_1col_file <-
function (mri_paths_file) 
{
    read.table(mri_paths_file, as.is = TRUE)$V1
}
.read_2col_file <-
function (mri_paths_file) 
{
    read.table(mri_paths_file, as.is = TRUE)$V2
}
.read_3col_file <-
function (mri_paths_file) 
{
    read.table(mri_paths_file, as.is = TRUE)$V3
}
.read_4col_file <-
function (mri_paths_file) 
{
    read.table(mri_paths_file, as.is = TRUE)$V4
}
.read_covars_file <-
function (path_file, split_txt = ",") 
{
    con = file(path_file, open = "r")
    linn = list()
    lines = readLines(con)
    for (i in 1:4) {
        linn[i] = list(unlist(strsplit(lines[i], split = split_txt)))
        if (is.na(linn[i])) 
            linn[i] = ""
    }
    close(con)
    linn
}
.read_data_table <-
function (data_table_file) 
{
    text = gsub(",", "\t", readLines(data_table_file))
    text = read.csv(textConnection(text), sep = "\t", check.names = FALSE)
    as.matrix(text)
}
.read_folds_file <-
function (folds_file_path) 
{
    con = file(folds_file_path, open = "r")
    lines = readLines(con)
    linn = list()
    long = length(lines)
    for (i in 1:long) {
        linn[i] = list(as.numeric(unlist(strsplit(lines[i], split = ","))))
    }
    close(con)
    linn
}
.read_mri <-
function (mri_paths, read_data = TRUE) 
{
    mri_tmp <- readNIfTI(mri_paths, read_data = read_data)
    dim <- mri_tmp@dim_[2:(2 + mri_tmp@dim_[1] - 1)]
    sform_code = ""
    if (mri_tmp@sform_code == 4) {
        sform_code = "NIFTI.XFORM.MNI.152"
    }
    mri <- list(dim = dim, sto.xyz = matrix(c(mri_tmp@srow_x, 
        mri_tmp@srow_y, mri_tmp@srow_z, c(0, 0, 0, 1)), nrow = 4, 
        ncol = 4, byrow = TRUE), scl.slope = mri_tmp@scl_slope, 
        scl.inter = mri_tmp@scl_inter, sform.code = sform_code, 
        data = mri_tmp@.Data)
    if (sum(mri$sto.xyz) == 1) {
        mri$sto.ijk <- solve(matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 
            0, 0, 1, 0, 1, 1, 1, 1), nrow = 4))
    }
    else {
        mri$sto.ijk <- (matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
            1, 0, 1, 1, 1, 1), nrow = 4) %*% solve(mri$sto.xyz))
    }
    mri
}
.require <-
function (package) 
{
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
        install.packages(package, repos = "http://cran.uk.r-project.org")
        library(package, character.only = TRUE, quietly = TRUE)
    }
}
apply_model <-
function (mp, mri_data, mri_fu_data, mri_wm_data = NULL, mri_wm_fu_data = NULL, 
    covX_test = NULL, signif_indx, lasso_covB, lasso_covB_fu = NULL, 
    mask, predX_test, scale_clinical = NULL, scale_mri = NULL, 
    lasso, lasso_predX_indx, tipett_take_un = NULL, img_kappa = NULL, 
    use_significant_voxels = FALSE, covX_site = c(), n_voxels_mask, 
    combat = NULL, standardize_images = FALSE, masks_3d) 
{
    X = NULL
    X_signif = c()
    if (mp$modulation == "op") 
        TIPETT = TRUE
    else TIPETT = FALSE
    if (standardize_images) {
        mri_data = mri_data - scale_mri$mean_mri
        if (!is.null(mri_fu_data)) {
            mri_fu_data = (mri_fu_data - scale_mri$mean_mri_fu)/scale_mri$sd_mri_fu * 
                scale_mri$sd_mri
        }
        if (!is.null(mri_wm_data)) {
            mri_wm_data = (mri_wm_data - scale_mri$mean_mri_wm)/scale_mri$sd_mri_wm * 
                scale_mri$sd_mri
        }
        if (!is.null(mri_wm_fu_data)) {
            mri_wm_fu_data = (mri_wm_fu_data - scale_mri$mean_mri_wm_fu)/scale_mri$sd_mri_wm_fu * 
                scale_mri$sd_mri
        }
    }
    n_image_modalities = 0
    if (!is.null(n_voxels_mask) && any(lasso$i <= n_voxels_mask)) {
        if (!is.null(mri_data)) {
            n_image_modalities = n_image_modalities + 1
            n_subj = nrow(covX_test)
            mri_test = lapply(seq_len(n_subj), function(i) {
                mri_data[, , , i][masks_3d$un_gm]
            })
            mri_test = matrix(unlist(mri_test), nrow = n_subj, 
                byrow = TRUE)
            if (!is.null(covX_site)) {
                mri_test = combat_apply(combat$gm, mri_test, 
                  covX_site, mod = covX_test, verbose = FALSE)$dat.combat
            }
            img_3d = mri_data[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_data)[4]) {
                img_3d[masks_3d$un_gm] = mri_test[img_i, ]
                mri_data[, , , img_i] = img_3d
            }
        }
        if (!is.null(mri_fu_data)) {
            n_image_modalities = n_image_modalities + 1
            mri_test = lapply(seq_len(n_subj), function(i) {
                mri_fu_data[, , , i][masks_3d$fu_gm]
            })
            mri_test = matrix(unlist(mri_test), nrow = n_subj, 
                byrow = TRUE)
            if (!is.null(covX_site)) {
                mri_test = combat_apply(combat$gm_mod, mri_test, 
                  covX_site, mod = covX_test, verbose = FALSE)$dat.combat
            }
            img_3d = mri_fu_data[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_fu_data)[4]) {
                img_3d[masks_3d$fu_gm] = mri_test[img_i, ]
                mri_fu_data[, , , img_i] = img_3d
            }
        }
        if (!is.null(mri_wm_data)) {
            n_image_modalities = n_image_modalities + 1
            mri_test = lapply(seq_len(n_subj), function(i) {
                mri_wm_data[, , , i][masks_3d$un_wm]
            })
            mri_test = matrix(unlist(mri_test), nrow = n_subj, 
                byrow = TRUE)
            if (!is.null(covX_site)) {
                mri_test = combat_apply(combat$wm, mri_test, 
                  covX_site, mod = covX_test, verbose = FALSE)$dat.combat
            }
            img_3d = mri_wm_data[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_wm_data)[4]) {
                img_3d[masks_3d$un_wm] = mri_test[img_i, ]
                mri_wm_data[, , , img_i] = img_3d
            }
        }
        if (!is.null(mri_wm_fu_data)) {
            n_image_modalities = n_image_modalities + 1
            mri_test = lapply(seq_len(n_subj), function(i) {
                mri_wm_fu_data[, , , i][masks_3d$fu_wm]
            })
            mri_test = matrix(unlist(mri_test), nrow = n_subj, 
                byrow = TRUE)
            if (!is.null(covX_site)) {
                mri_test = combat_apply(combat$wm_mod, mri_test, 
                  covX_site, mod = covX_test, verbose = FALSE)$dat.combat
            }
            img_3d = mri_wm_fu_data[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_wm_fu_data)[4]) {
                img_3d[masks_3d$fu_wm] = mri_test[img_i, ]
                mri_wm_fu_data[, , , img_i] = img_3d
            }
        }
        if (!is.null(mri_data)) {
            n_subj = dim(mri_data)[4]
            if (mp$modulation != "un" && mp$modulation != "fu") {
                orig_dims = dim(mri_data)
                new_mri = array(NA, dim = dim(mri_data) * c(n_image_modalities, 
                  1, 1, 1))
                for (img_i in 1:dim(mri_data)[4]) {
                  if (n_image_modalities > 0) {
                    img_j = 1
                    first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                      1)
                    new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                      1), , , img_i] = mri_data[, , , img_i]
                  }
                  if (n_image_modalities > 1) {
                    img_j = 2
                    first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                      1)
                    new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                      1), , , img_i] = mri_fu_data[, , , img_i]
                  }
                  if (n_image_modalities > 2) {
                    img_j = 3
                    first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                      1)
                    new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                      1), , , img_i] = mri_wm_data[, , , img_i]
                  }
                  if (n_image_modalities > 3) {
                    img_j = 4
                    first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                      1)
                    new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                      1), , , img_i] = mri_wm_fu_data[, , , img_i]
                  }
                }
                mri_data = new_mri
            }
            mri_test = lapply(seq_len(n_subj), function(i) {
                mri_data[, , , i][which(mask)]
            })
            mri_test = matrix(unlist(mri_test), nrow = n_subj, 
                byrow = TRUE)
            n_voxels = sum(mask)
            if (use_significant_voxels) {
                mri_test = mri_test[, signif_indx]
                n_voxels = sum(signif_indx)
            }
            mri_test = mri_test[, lasso$i[which(lasso$i <= n_voxels)]]
            if (mp$modulation %in% c("fu", "op")) {
                mri_fu_test = lapply(seq_len(n_subj), function(i) {
                  mri_fu_data[, , , i][which(mask)]
                })
                mri_fu_test = matrix(unlist(mri_fu_test), nrow = n_subj, 
                  byrow = TRUE)
                if (use_significant_voxels) {
                  mri_fu_test = mri_fu_test[, signif_indx]
                }
                mri_fu_test = mri_fu_test[, lasso$i[which(lasso$i <= 
                  n_voxels)]]
            }
            if (!is.null(covX_test)) {
                X_un = mri_test - covX_test %*% lasso_covB
            }
            else {
                X_un = mri_test
            }
            if (!is.null(lasso_covB_fu)) {
                if (!is.null(covX_test)) {
                  X_fu = mri_fu_test - covX_test %*% lasso_covB_fu
                }
                else {
                  X_fu = mri_fu_test
                }
            }
            else {
                X_fu = NULL
            }
            if (TIPETT) {
                if (use_significant_voxels) {
                  tipett_take_un = tipett_take_un[which(signif_indx)]
                }
                take_X_un_val = tipett_take_un[lasso$i]
                X_op = matrix(0, nrow = nrow(X_un), ncol = ncol(X_un))
                X_op[, take_X_un_val] = X_un[, take_X_un_val]
                X_op[, !take_X_un_val] = X_fu[, !take_X_un_val]
            }
            else if (mp$modulation == "op") {
                for (i in 1:ncol(lasso_covB)) {
                  X_op[, i] = (1 - img_kappa[i]) * X_un[, i] + 
                    img_kappa[i] * X_fu[, i]
                }
            }
            X_signif = switch(mp$modulation, un = X_un, fu = X_fu, 
                op = X_op, all = X_un)
        }
    }
    if (!is.null(mp$pred_transf)) {
        predX_test = as.matrix(predX_test)
        if (!is.null(scale_clinical) & !any(is.na(scale_clinical$MEANS))) {
            predX_test = t(apply(predX_test, 1, function(x) {
                x - scale_clinical$MEANS
            }))
            predX_test = t(apply(predX_test, 1, function(x) {
                x/scale_clinical$SDS * scale_clinical$MODE_SD
            }))
        }
        predX_test = predX_test[, lasso_predX_indx]
        X = cbind(X_signif, predX_test)
    }
    else {
        X = X_signif
    }
    linPred = as.matrix(X) %*% lasso$beta
    if (mp$response_family != "cox") {
        linPred = lasso$a0 + linPred
    }
    linPred
}
aprior <-
function (gamma.hat) 
{
    m = mean(gamma.hat)
    s2 = var(gamma.hat)
    (2 * s2 + m^2)/s2
}
assign.folds <-
function (y, family = c("binomial", "cox", "gaussian"), site = NULL, 
    nfolds = 10) 
{
    if (!((is.vector(family) || is.factor(family)) && length(family) == 
        1 && family %in% c("binomial", "cox", "gaussian"))) {
        stop("family must be \"binomial\", \"cox\", or \"gaussian\"")
    }
    n = switch(family, binomial = {
        if (!(is.vector(y) && all(y %in% 0:1))) {
            stop("for \"binomial\", y must be a binary vector")
        }
        length(y)
    }, cox = {
        if (class(y) != "Surv") {
            stop("for \"cox\", y must be a \"Surv\" object")
        }
        nrow(y)
    }, gaussian = {
        if (!(is.vector(y) && is.numeric(y))) {
            stop("for \"gaussian\", y must be a numeric vector")
        }
        length(y)
    })
    if (!(is.null(site) || ((is.vector(site) || is.factor(site)) && 
        length(site) == n))) {
        stop("site must be a vector with the same length as y, or NULL")
    }
    if (!(is.vector(nfolds) && is.numeric(nfolds) && length(nfolds) == 
        1 && nfolds > 0)) {
        stop("nfolds must be a positive number")
    }
    if (is.null(site)) {
        return(.assign.folds_one_site(y, family, nfolds))
    }
    folds = rep(NA, length(y))
    for (site_i in site) {
        i = which(site == site_i)
        folds[i] = .assign.folds_one_site(y[i], family, nfolds)
    }
    folds
}
bprior <-
function (gamma.hat) 
{
    m = mean(gamma.hat)
    s2 = var(gamma.hat)
    (m * s2 + m^3)/s2
}
calculate_effects_controls <-
function (mri, covX, path) 
{
    B = solve(t(covX) %*% covX) %*% t(covX)
    beta = apply(mri, 1:3, function(x) {
        B %*% x
    })
    paths = c()
    for (i in 1:ncol(covX)) {
        mri_example[] = nifti(beta[i, , , ])
        writeNIfTI(mri_example, paste(path, "_beta", i, sep = ""))
        paths = c(paths, paste(path, "_beta", i, sep = ""))
    }
    paths
}
combat_apply <-
function (tmp, dat, batch, mod = NULL, verbose = TRUE) 
{
    if (!is.factor(batch) || nlevels(batch) != length(tmp$levels_batch) || 
        any(levels(batch) != tmp$levels_batch)) {
        stop("batch must be a factor with the same levels than when fitting combat")
    }
    if (is.data.frame(dat)) {
        dat <- as.matrix(dat)
    }
    else if (!is.matrix(dat)) {
        stop("dat must be a matrix")
    }
    if (tmp$transpose) {
        if (nrow(dat) != length(batch)) {
            stop("dat must have the same number of rows than the length of batch")
        }
        dat <- t(dat)
    }
    else {
        if (ncol(dat) != length(batch)) {
            stop("dat must have the same number of columns than the length of batch")
        }
    }
    dat.combat <- dat
    dat <- dat[tmp$not_constant, ]
    if (is.data.frame((mod))) {
        mod <- as.matrix(mod)
    }
    else if (!(is.matrix(mod) || is.null(mod))) {
        stop("mod must be a matrix or NULL")
    }
    tmp1 <- .combat_tmp1(dat, batch, tmp$levels_batch, mod)
    tmp3 <- .combat_tmp3(dat, tmp1, tmp$tmp2, verbose)
    tmp5 <- .combat_tmp5(tmp1, tmp$tmp2, tmp3, tmp$tmp4, tmp$eb, 
        verbose)
    dat.combat[tmp$not_constant, ] <- tmp5
    if (tmp$transpose) {
        dat.combat <- t(dat.combat)
    }
    return(list(dat.combat = dat.combat, gamma.hat = tmp$tmp4$gamma.hat, 
        delta.hat = tmp$tmp4$delta.hat, gamma.star = tmp$tmp4$gamma.star, 
        delta.star = tmp$tmp4$delta.star, gamma.bar = tmp$tmp4$gamma.bar, 
        t2 = tmp$tmp4$t2, a.prior = tmp$tmp4$a.prior, b.prior = tmp$tmp4$b.prior, 
        batch = tmp1$batch, mod = tmp1$mod, stand.mean = tmp3$stand.mean, 
        stand.sd = sqrt(tmp$tmp2$var.pooled)[, 1]))
}
combat_fit <-
function (dat, batch, mod = NULL, eb = TRUE, verbose = TRUE) 
{
    if (!is.factor(batch)) {
        stop("batch must be a factor")
    }
    if (is.data.frame(dat)) {
        dat <- as.matrix(dat)
    }
    else if (!is.matrix(dat)) {
        stop("dat must be a matrix")
    }
    if (ncol(dat) == length(batch)) {
        transpose <- FALSE
        if (verbose) {
            cat("[combat] Subjects are COLUMNS\n")
        }
    }
    else if (nrow(dat) == length(batch)) {
        transpose <- TRUE
        dat <- t(dat)
        if (verbose) {
            cat("[combat] Subjects are ROWS\n")
        }
    }
    else {
        stop("dat must have the same number of columns or rows than the length of batch")
    }
    if (is.data.frame((mod))) {
        mod <- as.matrix(mod)
    }
    else if (!(is.matrix(mod) || is.null(mod))) {
        stop("mod must be a matrix or NULL")
    }
    if (any(is.na(dat))) {
        if (verbose) {
            cat("[combat] Imputing missing data (only for fit)\n")
        }
        for (batch_i in sort(unique(batch))) {
            i <- which(batch == batch_i)
            dat_i <- dat[, i]
            if (any(is.na(dat_i))) {
                for (j in 1:nrow(dat)) {
                  dat_ji <- dat_i[j, ]
                  is_na <- which(is.na(dat_ji))
                  if (length(is_na) > 0 && length(is_na) < length(i)) {
                    if (length(is_na) == 1) {
                      mod_i_is_na <- matrix(mod[i[is_na], ], 
                        nrow = 1)
                    }
                    else {
                      mod_i_is_na <- mod[i[is_na], ]
                    }
                    beta <- matrix(coef(lm(dat_ji ~ mod[i, ])))
                    beta[which(is.na(beta))] <- 0
                    dat[j, i[is_na]] <- cbind(1, mod_i_is_na) %*% 
                      beta
                  }
                  else {
                    dat[j, i[is_na]] <- mean(dat_ji, na.rm = TRUE)
                  }
                }
            }
        }
        for (batch_i in sort(unique(batch))) {
            i <- which(batch == batch_i)
            if (any(is.na(dat[, i]))) {
                for (j in 1:nrow(dat)) {
                  dat_j <- dat[j, ]
                  if (is.na(dat_j[i[1]])) {
                    if (!is.null(mod)) {
                      beta <- matrix(coef(lm(dat_j ~ mod)))
                      beta[which(is.na(beta))] <- 0
                      dat[j, i] <- cbind(1, mod[i, ]) %*% beta
                    }
                    else {
                      dat[j, i] <- mean(dat_j, na.rm = TRUE)
                    }
                  }
                }
            }
        }
    }
    not_constant <- which(apply(dat, 1, function(x) {
        var(x) > 0
    }))
    dat <- dat[not_constant, ]
    if (eb) {
        if (verbose) {
            cat("[combat] Performing ComBat with empirical Bayes\n")
        }
    }
    else {
        if (verbose) {
            cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
        }
    }
    if (verbose) {
        cat("[combat] Found", nlevels(batch), "batches (e.g., sites)\n")
    }
    levels_batch <- levels(batch)
    tmp1 <- .combat_tmp1(dat, batch, levels_batch, mod)
    tmp2 <- .combat_tmp2(tmp1, verbose)
    tmp3 <- .combat_tmp3(dat, tmp1, tmp2, verbose)
    tmp4 <- .combat_tmp4(tmp1, tmp2, tmp3, eb, verbose)
    return(list(levels_batch = levels_batch, transpose = transpose, 
        not_constant = not_constant, eb = eb, tmp2 = tmp2, tmp4 = tmp4))
}
createMatchingIndices <-
function (x, batch, xmin = NULL, xmax = NULL, step = 1) 
{
    stopifnot(length(x) == length(batch))
    batches <- unique(batch)
    n.batches <- length(batches)
    x_per_batch <- split(x, f = batch)[batches]
    if (is.null(xmin)) 
        xmin <- min(x)
    if (is.null(xmax)) 
        xmax <- max(x)
    grid <- seq(xmin, xmax, step)
    n.bins <- length(grid) - 1
    counts <- matrix(0, n.bins, n.batches)
    for (i in 1:n.bins) {
        counts[i, ] <- unlist(lapply(x_per_batch, function(temp) {
            sum(temp >= grid[i] & temp < grid[i + 1])
        }))
    }
    mins <- unlist(apply(counts, 1, min))
    indices <- c()
    for (i in 1:n.bins) {
        for (j in 1:n.batches) {
            min <- mins[i]
            if (min != 0) {
                cand <- which(x >= grid[i] & x < grid[i + 1] & 
                  batch == batches[j])
                if (length(cand) != 1) {
                  cand <- sample(cand, min)
                }
                indices <- c(indices, cand)
            }
        }
    }
    return(indices)
}
cv <-
function (x, y, family = c("binomial", "cox", "gaussian"), fit_fun, 
    predict_fun, site = NULL, covar = NULL, nfolds = 10, ...) 
{
    if (!is.null(site)) {
        if (!is.null(covar)) {
            cat("[cv] Cross-validation with sites and covariates\n")
            type = "site+covar"
        }
        else {
            cat("[cv] Cross-validation with sites\n")
            type = "site"
        }
        site = as.factor(site)
    }
    else {
        if (!is.null(covar)) {
            cat("[cv] Cross-validation with covariates\n")
            type = "covar"
        }
        else {
            cat("[cv] Simple cross-validation\n")
            type = "simple"
        }
    }
    folds = assign.folds(y, family, site = site, nfolds = nfolds)
    models = list()
    output = data.frame(fold = folds, y, y.pred = NA)
    for (fold in 1:nfolds) {
        cat("[cv] Fold", fold, " - Training\n")
        training = which(folds != fold)
        x_training = x[training, ]
        y_training = switch(family, cox = y[training, ], y[training])
        model = switch(type, simple = fit_fun(x_training, y_training, 
            NULL, NULL, ...), site = fit_fun(x_training, y_training, 
            site[training], NULL, ...), covar = fit_fun(x_training, 
            y_training, NULL, covar[training, ], ...), `site+covar` = fit_fun(x_training, 
            y_training, site[training], covar[training, ], ...))
        models[[fold]] = model
        cat("[cv] Fold", fold, " - Test\n")
        test = which(folds == fold)
        x_test = matrix(x[test, ], ncol = ncol(x))
        output$y.pred[test] = switch(type, simple = predict_fun(model, 
            x_test, NULL, NULL, ...), site = predict_fun(model, 
            x_test, site[test], NULL, ...), covar = predict_fun(model, 
            x_test, NULL, covar[test, ], ...), `site+covar` = predict_fun(model, 
            x_test, site[test], covar[test, ], ...), )
    }
    list(predictions = output, models = models)
}
data.frame2glmnet.matrix <-
function (m, x) 
{
    if (!is.data.frame(x)) {
        stop("x must be a data.frame")
    }
    if (class(m) != "data.frame2glmnet.matrix_fit") {
        stop("m must be a \"data.frame2glmnet.matrix\" object")
    }
    xp = NULL
    if (length(m) > 0) {
        for (i in 1:length(m)) {
            transf_i = m[[i]]
            xj = x[, match(transf_i[1], colnames(x))]
            xp = cbind(xp, switch(transf_i[2], factor = {
                if (length(transf_i) == 4) {
                  xpj = matrix(as.numeric(xj == transf_i[4]))
                  colnames(xpj) = paste0(transf_i[1], ":", transf_i[4])
                } else {
                  xpj = NULL
                  for (k in 3:length(transf_i)) {
                    xpj = cbind(xpj, as.numeric(xj == transf_i[k]))
                  }
                  colnames(xpj) = paste0(transf_i[1], ":", transf_i[3:length(transf_i)])
                }
                xpj
            }, numeric = {
                xpj = matrix(xj)
                colnames(xpj) = transf_i[1]
                xpj
            }))
        }
    }
    xp
}
data.frame2glmnet.matrix_fit <-
function (x) 
{
    if (!is.data.frame(x)) {
        stop("x must be a data.frame")
    }
    m = list()
    if (ncol(x) > 0) {
        for (j in 1:ncol(x)) {
            xj_name = colnames(x)[j]
            xj_char = as.character(x[, j])
            xj_not_na_num = suppressWarnings(as.numeric(xj_char[which(!is.na(xj_char))]))
            xj_levels = sort(unique(xj_char))
            if (length(xj_levels) < 2) {
                stop(paste("variable", xj_name, "has no different values"))
            }
            if (any(is.na(xj_not_na_num))) {
                m[[length(m) + 1]] = c(xj_name, "factor", xj_levels)
            }
            else {
                m[[length(m) + 1]] = c(xj_name, "numeric")
            }
        }
    }
    class(m) = "data.frame2glmnet.matrix_fit"
    m
}
fit_model <-
function (mp, data_informative_table, Y, mri, mri_fu, preloaded_covB = NULL, 
    preloaded_covB_fu = NULL, SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998), 
    iter = -1, internal_folds = c(), n_cores = 1, use_significant_voxels = FALSE, 
    covX_site = c(), mri_wm = NULL, mri_fu_wm = NULL, name_combat = "_combat", 
    mask = NULL, standardize_images = FALSE) 
{
    tipett_take_un = c()
    if (!is.null(data_informative_table)) {
        if (ncol(data_informative_table) > 1) {
            covX_training = as.matrix(cbind(1, data_informative_table[, 
                mp$covX_var]))
            predX_train = as.matrix(data_informative_table[, 
                mp$pred_var])
        }
        else {
            covX_training = as.matrix(data_informative_table)
        }
        data_informative_table = as.matrix(data_informative_table)
    }
    else {
        covX_training = matrix(1, nrow = nrow(as.matrix(Y)))
    }
    .print_ok()
    masked_X = c()
    masked_X_sign = c()
    scale_mri = list()
    scale_mri$sd_mri = NA
    scale_mri$sd_mri_fu = NA
    scale_mri$sd_mri_wm = NA
    scale_mri$sd_mri_wm_fu = NA
    scale_mri$mean_mri = NA
    scale_mri$mean_mri_fu = NA
    scale_mri$mean_mri_wm = NA
    scale_mri$mean_mri_wm_fu = NA
    if (standardize_images) {
        scale_mri$sd_mri = sd(mri)
        scale_mri$mean_mri = mean(mri)
        mri = mri - scale_mri$mean_mri
        if (!is.null(mri_fu)) {
            scale_mri$sd_mri_fu = sd(mri_fu)
            scale_mri$mean_mri_fu = mean(mri_fu)
            mri_fu = (mri_fu - scale_mri$mean_mri_fu)/scale_mri$sd_mri_fu * 
                scale_mri$sd_mri
        }
        if (!is.null(mri_wm)) {
            scale_mri$sd_mri_wm = sd(mri_wm)
            scale_mri$mean_mri_wm = mean(mri_wm)
            mri_wm = (mri_wm - scale_mri$mean_mri_wm)/scale_mri$sd_mri_wm * 
                scale_mri$sd_mri
        }
        if (!is.null(mri_fu_wm)) {
            scale_mri$sd_mri_wm_fu = sd(mri_wm_fu)
            scale_mri$mean_mri_wm_fu = mean(mri_wm_fu)
            mri_fu_wm = (mri_fu_wm - scale_mri$mean_mri_wm_fu)/scale_mri$sd_mri_wm_fu * 
                scale_mri$sd_mri
        }
    }
    signif_indx = c()
    masks_3d = list()
    combat = list()
    if (!is.null(mp$mri_paths) || !is.null(mp$mri_un)) {
        if (!is.null(mask)) {
            mask = mask
        }
        else {
            mask = c()
            mask = array(TRUE, dim = dim(mri[, , , 1]))
            mask_o = mask
        }
        mask = mask & (apply(mri, 1:3, var) > 0)
        mask_ijk = which(mask, arr.ind = TRUE)
        if (iter != -1) {
            .print_action("Training sample: cutting brain")
            mask = rotate_coordinates(list_coordinates = mask_ijk, 
                iter = iter, return_3d = TRUE, img_3d = mask_o) > 
                0
        }
        n_image_modalities = 0
        if (!is.null(mri)) {
            n_image_modalities = n_image_modalities + 1
            masks_3d = list()
            if (!is.null(covX_site)) {
                .print_action("Applying combat")
            }
            list_masked = remove_effects(mask_3d = mask, mri = mri, 
                covX_training = covX_training, trainY = Y, response_family = mp$response_family, 
                REMOVE_EFFECT_TO_X = FALSE, n_cores = n_cores, 
                SIGNIFICANT = (use_significant_voxels || mp$modulation == 
                  "op"), covX_site = covX_site, modality = paste0("_un_gm", 
                  name_combat))
            combat$gm = list_masked$combat
            img_3d = mri[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri)[4]) {
                img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i, 
                  ]
                mri[, , , img_i] = img_3d
            }
            masks_3d$un_gm = list_masked$mask_3d
        }
        if (!is.null(mri_fu)) {
            n_image_modalities = n_image_modalities + 1
            list_masked = remove_effects(mask_3d = mask, mri = mri_fu, 
                covX_training = covX_training, trainY = Y, response_family = mp$response_family, 
                REMOVE_EFFECT_TO_X = FALSE, n_cores = n_cores, 
                SIGNIFICANT = (use_significant_voxels || mp$modulation == 
                  "op"), covX_site = covX_site, modality = paste0("_fu_gm", 
                  name_combat))
            combat$gm_mod = list_masked$combat
            img_3d = mri_fu[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_fu)[4]) {
                img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i, 
                  ]
                mri_fu[, , , img_i] = img_3d
            }
            masks_3d$fu_gm = list_masked$mask_3d
        }
        if (!is.null(mri_wm)) {
            n_image_modalities = n_image_modalities + 1
            list_masked = remove_effects(mask_3d = mask, mri = mri_wm, 
                covX_training = covX_training, trainY = Y, response_family = mp$response_family, 
                REMOVE_EFFECT_TO_X = FALSE, n_cores = n_cores, 
                SIGNIFICANT = (use_significant_voxels || mp$modulation == 
                  "op"), covX_site = covX_site, modality = paste0("_un_wm", 
                  name_combat))
            combat$wm = list_masked$combat
            img_3d = mri_wm[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_wm)[4]) {
                img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i, 
                  ]
                mri_wm[, , , img_i] = img_3d
            }
            masks_3d$un_wm = list_masked$mask_3d
        }
        if (!is.null(mri_fu_wm)) {
            n_image_modalities = n_image_modalities + 1
            list_masked = remove_effects(mask_3d = mask, mri = mri_fu_wm, 
                covX_training = covX_training, trainY = Y, response_family = mp$response_family, 
                REMOVE_EFFECT_TO_X = FALSE, n_cores = n_cores, 
                SIGNIFICANT = (use_significant_voxels || mp$modulation == 
                  "op"), covX_site = covX_site, modality = paste0("_fu_wm", 
                  name_combat))
            combat$wm_mod = list_masked$combat
            img_3d = mri_fu_wm[, , , 1]
            img_3d = img_3d * 0
            for (img_i in 1:dim(mri_fu_wm)[4]) {
                img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i, 
                  ]
                mri_fu_wm[, , , img_i] = img_3d
            }
            masks_3d$fu_wm = list_masked$mask_3d
        }
        mask_single_image <- mask
        if (mp$modulation != "un" && mp$modulation != "fu") {
            orig_dims = dim(mask)
            new_dims = orig_dims * c(n_image_modalities, 1, 1)
            new_mask = array(NA, dim = new_dims)
            for (img_i in 1:n_image_modalities) {
                mask_i <- switch(as.character(img_i), `1` = masks_3d$un_gm, 
                  `2` = masks_3d$fu_gm, `3` = masks_3d$un_wm, 
                  `4` = masks_3d$fu_wm)
                first_index_x = (((img_i - 1) * (orig_dims[1])) + 
                  1)
                new_mask[first_index_x:(first_index_x + orig_dims[1] - 
                  1), , ] = mask_i
            }
            mask = new_mask
            new_mask <- NULL
            new_mri = array(NA, dim = dim(mri) * c(n_image_modalities, 
                1, 1, 1))
            if (n_image_modalities > 0) {
                img_j = 1
                first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                  1)
                new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                  1), , , ] = mri
            }
            if (n_image_modalities > 1) {
                img_j = 2
                first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                  1)
                new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                  1), , , ] = mri_fu
            }
            if (n_image_modalities > 2) {
                img_j = 3
                first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                  1)
                new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                  1), , , ] = mri_wm
            }
            if (n_image_modalities > 3) {
                img_j = 4
                first_index_x = (((img_j - 1) * (orig_dims[1])) + 
                  1)
                new_mri[first_index_x:(first_index_x + orig_dims[1] - 
                  1), , , ] = mri_fu_wm
            }
            mri <- new_mri
            mri_wm <- NULL
            mri_fu_wm <- NULL
            new_mri <- NULL
        }
        mask_ijk = cbind(which(mask, arr.ind = TRUE), 1)
        mask_ijk = t(mask_ijk[which(mask_ijk[, 4] == 1), ])
        .print_ok()
        .print_action("Training sample: removing effects of covariates\n")
        if (is.null(preloaded_covB)) {
            list_masked = remove_effects(mask_3d = mask, mri = mri, 
                covX_training = covX_training, trainY = Y, response_family = mp$response_family, 
                n_cores = n_cores, SIGNIFICANT = (use_significant_voxels || 
                  mp$modulation == "op"))
        }
        else {
            list_masked = remove_effects_precalculated(mask = mask, 
                preloaded_covB = preloaded_covB, mri = mri, covX = covX_training, 
                trainY = Y, response_family = mp$response_family, 
                n_cores = n_cores, SIGNIFICANCE = (use_significant_voxels || 
                  mp$modulation == "op"), covX_site = covX_site)
        }
        mask_ijk = list_masked$mask_ijk
        mask = list_masked$mask_3d
        .print_ok()
        if (mp$modulation == "fu" || mp$modulation == "op") {
            .print_action("Training sample: Removing effect of covariates on fully modulated data")
            if (is.null(preloaded_covB)) 
                list_fu_masked = remove_effects(mask_3d = mask, 
                  mri = mri_fu, covX_training = covX_training, 
                  trainY = Y, n_cores = n_cores, response_family = mp$response_family)
            else list_fu_masked = remove_effects_precalculated(mask = mask, 
                preloaded_covB = preloaded_covB_fu, mri = mri_fu, 
                covX = covX_training, trainY = Y, response_family = mp$response_family, 
                n_cores = n_cores, SIGNIFICANCE = (use_significant_voxels || 
                  mp$modulation == "op"))
            .print_ok()
        }
        mri_fu <- NULL
        if (mp$modulation == "op") {
            .print_action("Training sample: optimal modulation")
            list_op_masked = tipett_modulation(mask_X = list_masked$X, 
                t_or_z_vals = list_fu_masked$t_or_z_vals, mask_fu_X = list_fu_masked$X, 
                t_or_z_fu_vals = list_fu_masked$t_or_z_vals, 
                covX_training = covX_training, trainY = Y, response_family = mp$response_family, 
                SIGNIFICANCE = use_significant_voxels)
            tipett_take_un = list_op_masked$tipett_take_un
            .print_ok()
        }
        if (mp$modulation == "un") {
            masked_X = list_masked$X
            mask_X = ""
            tvals = list_masked$t_or_z_vals
            masked_covB = list_masked$covB
            masked_covB_fu = c()
        }
        else if (mp$modulation == "fu") {
            masked_X = list_fu_masked$X
            mask_fu_X = ""
            tvals = list_fu_masked$t_or_z_vals
            masked_covB = list_fu_masked$covB
            masked_covB_fu = list_fu_masked$covB
        }
        else if (mp$modulation == "op") {
            masked_X = list_op_masked$mask_op_X
            mask_op_X = ""
            tvals = list_op_masked$t_or_z_op_vals
            masked_covB = list_masked$covB
            masked_covB_fu = list_fu_masked$covB
        }
        else {
            masked_X = list_masked$X
            mask_X = ""
            tvals = list_masked$t_or_z_vals
            masked_covB = list_masked$covB
            masked_covB_fu = c()
        }
        signif_indx = c()
        if (use_significant_voxels) {
            signif_indx = abs(tvals) > SIGNIFICANCE_THRESHOLD
            if (sum(signif_indx) < 2) 
                signif_indx = which(abs(tvals) == sort(abs(tvals), 
                  decreasing = TRUE)[1:2])
            mask_ijk = mask_ijk[, signif_indx]
            if (ncol(as.matrix(masked_covB)) == 1) {
                masked_covB = masked_covB[signif_indx]
            }
            else {
                masked_covB = masked_covB[, signif_indx]
            }
            if (mp$modulation == "op" || mp$modulation == "fu") {
                if (ncol(as.matrix(masked_covB_fu)) == 1) {
                  masked_covB_fu = masked_covB_fu[signif_indx]
                }
                else {
                  masked_covB_fu = masked_covB_fu[, signif_indx]
                }
            }
            masked_X = masked_X[, signif_indx]
        }
    }
    n_voxels_mask = ifelse(is.null(masked_X), 0, ncol(masked_X))
    scale_clinical = list()
    scale_clinical$MODE_SD = NA
    scale_clinical$MEANS = NA
    scale_clinical$SDS = NA
    if (!is.null(mp$pred_transf)) {
        if (!is.null(masked_X)) {
            scale_clinical = list(MODE_SD = c(), MEANS = c(), 
                SDS = c())
            DF = nrow(masked_X) - 1
            SD = apply(masked_X, 2, function(x) {
                sqrt(sum(x^2)/DF)
            })
            HIST = hist(SD, plot = FALSE, breaks = 100)
            LOWESS = lowess(HIST$mids, HIST$density)
            scale_clinical$MODE_SD = LOWESS$x[which(LOWESS$y == 
                max(LOWESS$y[-1:-10]))]
            scale_clinical$MEANS = apply(predX_train, 2, mean)
            predX_train = as.matrix(apply(predX_train, 1, function(x) {
                x - scale_clinical$MEANS
            }))
            if (length(scale_clinical$MEANS) > 1) 
                predX_train = t(predX_train)
            scale_clinical$SDS = apply(predX_train, 2, sd)
            predX_train = as.matrix(apply(predX_train, 1, function(x) {
                x/scale_clinical$SDS * scale_clinical$MODE_SD
            }))
            if (length(scale_clinical$MEANS) > 1) 
                predX_train = t(predX_train)
            columnsNaN <- which(apply(predX_train, 2, function(x) {
                all(is.na(x))
            }))
            predX_train[, columnsNaN] <- 0
            masked_X = cbind(masked_X, predX_train)
        }
        else {
            columnsNaN <- which(apply(predX_train, 2, function(x) {
                all(is.na(x))
            }))
            predX_train[, columnsNaN] <- 0
            masked_X = predX_train
        }
    }
    .print_ok()
    .print_action("Training sample: fitting lasso regression")
    lasso = glmnet_fit(x = masked_X, y = Y, family = mp$response_family, 
        foldid = internal_folds, nfolds = length(unique(internal_folds)), 
        standardize = FALSE)
    lassoB0 = lasso$a0
    cat("\nVariables used in the model: \n")
    voxels_used = lasso$i <= n_voxels_mask
    predictors_used = lasso$i > n_voxels_mask
    if (!all(voxels_used == FALSE)) 
        cat(sprintf("\nNumber of voxels: %s; Lasso index (voxel variable): %s", 
            n_voxels_mask, lasso$i[voxels_used]))
    if (!all(predictors_used == FALSE)) 
        cat(sprintf("\nNumber of voxels: %s; Lasso index: %s; Predictor variable: %s", 
            n_voxels_mask, lasso$i[predictors_used], names(predictors_used[predictors_used])))
    cat("\n")
    if (length(lasso$beta) == 1) {
        lassoB = t(lasso$beta)
    }
    else {
        lassoB = lasso$beta
    }
    lasso_covB = c()
    lasso_covB_fu = c()
    lasso_ijk = c()
    if (!is.null(mp$mri_paths) || !is.null(mp$mri_un)) {
        if (ncol(as.matrix(masked_covB)) == 1) {
            lasso_covB = as.matrix(masked_covB)[lasso$i[which(lasso$i <= 
                n_voxels_mask)]]
            if (mp$modulation %in% c("fu", "op")) 
                lasso_covB_fu = as.matrix(masked_covB_fu)[lasso$i[which(lasso$i <= 
                  n_voxels_mask)]]
        }
        else {
            lasso_covB = as.matrix(masked_covB)[, lasso$i[which(lasso$i <= 
                n_voxels_mask)]]
            if (mp$modulation %in% c("fu", "op")) 
                lasso_covB_fu = as.matrix(masked_covB_fu[, lasso$i[which(lasso$i <= 
                  n_voxels_mask)]])
        }
        lasso_ijk = matrix(mask_ijk[, lasso$i[which(lasso$i <= 
            n_voxels_mask)]], nrow = 4)[1:3, ]
        lasso_predX_indx = lasso$i[which(lasso$i > n_voxels_mask)] - 
            n_voxels_mask
    }
    else {
        lasso_predX_indx = lasso$i
    }
    lasso_mni = NULL
    if (!mp$modulation %in% c("cl", "clinical")) {
        lasso_ijk = matrix(lasso_ijk, nrow = 3)
        dim_x = mp$mri_params$dim[1]
        if (mp$modulation == "fu") {
            idx_gm_mod = which(lasso_ijk[1, ] > dim_x * 0 & lasso_ijk[1, 
                ] <= dim_x * 1)
            idx_gm = which(lasso_ijk[1, ] > dim_x * 1 & lasso_ijk[1, 
                ] <= dim_x * 2)
            if (length(idx_gm) > 0) {
                ijk_gm = lasso_ijk[, idx_gm] - c(dim_x * 1, 0, 
                  0)
                lasso_mni$gm = ijk2mni(rbind(matrix(ijk_gm, nrow = 3), 
                  1), mp$mri_params$sto.ijk)[, 1:length(idx_gm)]
                lasso_mni$gm_betas = lasso$beta[idx_gm]
            }
            else {
                lasso_mni$gm = NULL
                lasso_mni$gm_betas = NULL
            }
            if (length(idx_gm_mod) > 0) {
                ijk_gm_mod = lasso_ijk[, idx_gm_mod] - c(dim_x * 
                  0, 0, 0)
                lasso_mni$gm_mod = ijk2mni(rbind(matrix(ijk_gm_mod, 
                  nrow = 3), 1), mp$mri_params$sto.ijk)[, 1:length(idx_gm_mod)]
                lasso_mni$gm_mod_betas = lasso$beta[idx_gm_mod]
            }
            else {
                lasso_mni$gm_mod = NULL
                lasso_mni$gm_mod_betas = NULL
            }
        }
        else {
            idx_gm = which(lasso_ijk[1, ] > dim_x * 0 & lasso_ijk[1, 
                ] <= dim_x * 1)
            idx_gm_mod = which(lasso_ijk[1, ] > dim_x * 1 & lasso_ijk[1, 
                ] <= dim_x * 2)
            ijk_gm = lasso_ijk[, idx_gm] - c(dim_x * 0, 0, 0)
            if (length(idx_gm) > 0) {
                lasso_mni$gm = ijk2mni(rbind(matrix(ijk_gm, nrow = 3), 
                  1), mp$mri_params$sto.ijk)[, 1:length(idx_gm)]
                lasso_mni$gm_betas = lasso$beta[idx_gm]
            }
            else {
                lasso_mni$gm = NULL
                lasso_mni$gm_betas = NULL
            }
            if (length(idx_gm_mod) > 0) {
                ijk_gm_mod = lasso_ijk[, idx_gm_mod] - c(dim_x * 
                  1, 0, 0)
                lasso_mni$gm_mod = ijk2mni(rbind(matrix(ijk_gm_mod, 
                  nrow = 3), 1), mp$mri_params$sto.ijk)[, 1:length(idx_gm_mod)]
                lasso_mni$gm_mod_betas = lasso$beta[idx_gm_mod]
            }
            else {
                lasso_mni$gm_mod = NULL
                lasso_mni$gm_mod_betas = NULL
            }
        }
        idx_wm = which(lasso_ijk[1, ] > dim_x * 2 & lasso_ijk[1, 
            ] <= dim_x * 3)
        if (length(idx_wm) > 0) {
            ijk_wm = lasso_ijk[, idx_wm] - c(dim_x * 2, 0, 0)
            lasso_mni$wm = ijk2mni(rbind(matrix(ijk_wm, nrow = 3), 
                1), mp$mri_params$sto.ijk)[, 1:length(idx_wm)]
            lasso_mni$wm_betas = lasso$beta[idx_wm]
        }
        else {
            lasso_mni$wm = NULL
            lasso_mni$wm_betas = NULL
        }
        idx_wm_mod = which(lasso_ijk[1, ] > dim_x * 3 & lasso_ijk[1, 
            ] <= dim_x * 4)
        if (length(idx_wm_mod) > 0) {
            ijk_wm_mod = lasso_ijk[, idx_wm_mod] - c(dim_x * 
                3, 0, 0)
            lasso_mni$wm_mod = ijk2mni(rbind(matrix(ijk_wm_mod, 
                nrow = 3), 1), mp$mri_params$sto.ijk)[, 1:length(idx_wm_mod)]
            lasso_mni$wm_mod_betas = lasso$beta[idx_wm_mod]
        }
        else {
            lasso_mni$wm_mod = NULL
            lasso_mni$wm_mod_betas = NULL
        }
        signif_indx = signif_indx
        lasso_covB = lasso_covB
    }
    else {
        n_voxels_mask = 0
    }
    list(data_table_imputed_train = data_informative_table, signif_indx = signif_indx, 
        lasso_ijk = lasso_ijk, lasso_mni = lasso_mni, lasso_covB = lasso_covB, 
        lasso_fu_covB = lasso_covB_fu, mask = mask, lasso = lasso, 
        lasso_predX_indx = lasso_predX_indx, scale_predictors = scale_clinical, 
        scale_mri = scale_mri, tipett_take_un = tipett_take_un, 
        take_significant = use_significant_voxels, masks_3d = masks_3d, 
        n_voxels_mask = n_voxels_mask, combat = combat)
}
glmnet_fit <-
function (x, y, family = "binomial", foldid = NULL, nfolds = 10, 
    standardize = TRUE, min.beta = 9.9999999999999998e-13) 
{
    if (!((is.vector(family) || is.factor(family)) && length(family) == 
        1 && family %in% c("binomial", "cox", "gaussian"))) {
        stop("family must be \"binomial\", \"cox\", or \"gaussian\"")
    }
    if (family == "binomial") {
        if (!(is.vector(y) && all(y %in% 0:1))) {
            stop("for \"binomial\", y must be a binary vector")
        }
        n = length(y)
    }
    if (family == "cox") {
        if (class(y) != "Surv") {
            stop("for \"cox\", y must be a Surv object")
        }
        n = nrow(y)
    }
    if (family == "gaussian") {
        if (!(is.vector(y) && is.numeric(y))) {
            stop("for \"gaussian\", y must be a numeric vector")
        }
        n = length(y)
    }
    if (!(is.matrix(x) && nrow(x) == n)) {
        stop("x must be a matrix with the same height as y")
    }
    if (!(is.null(foldid) || (is.vector(foldid) && is.numeric(foldid) && 
        length(foldid) == n))) {
        stop("foldid must be a numeric vector with the same length as y, or NULL")
    }
    if (!(is.vector(nfolds) && is.numeric(nfolds) && length(nfolds) == 
        1 && nfolds > 0)) {
        stop("nfolds must be a positive number")
    }
    if (!(is.vector(standardize) && is.logical(standardize) && 
        length(standardize) == 1)) {
        stop("standardize must be TRUE or FALSE")
    }
    if (!(is.vector(min.beta) && is.numeric(min.beta) && length(min.beta) == 
        1 && min.beta > 0)) {
        stop("min.beta must be a positive number")
    }
    if (ncol(x) == 1) {
        coef = switch(family, binomial = coef(glm(y ~ x, family = binomial)), 
            cox = coef(coxph(Surv(time = y[, 1], event = y[, 
                2]) ~ x)), gaussian = coef(lm(y ~ x)))
        if (family == "cox") {
            a0 = NULL
            betas = coef
        }
        else {
            a0 = coef[1]
            betas = coef[2]
        }
        i = 1
    }
    else {
        type_measure = switch(family, binomial = "class", cox = "deviance", 
            gaussian = "mse")
        if (family == "cox") {
            colnames(y) <- c("time", "status")
        }
        if (is.null(foldid)) {
            foldid = assign.folds(y, family, nfolds = nfolds)
        }
        if (family == "binomial" && min(table(y)) < 3) {
            stop("too few subjects of one group")
        }
        cv = cv.glmnet(x, y, type.measure = type_measure, family = family, 
            foldid = foldid, nfolds = nfolds, standardize = standardize)
        idx_lambda = match(cv$lambda.min, cv$lambda)
        if (cv$glmnet.fit$df[idx_lambda] == 0) {
            idx_lambda = which(cv$glmnet.fit$df > 0)[1]
        }
        glmnet.control(fdev = 0)
        lasso = glmnet(x, y, family, lambda = cv$lambda, standardize = standardize)
        a0 = lasso$a0[idx_lambda]
        betas = lasso$beta[, idx_lambda]
        betas[which(abs(betas) < min.beta)] = 0
        i = which(betas != 0)
    }
    m = list(family = family, a0 = a0, i = i, beta = betas[i])
    class(m) = "glmnet_fit"
    m
}
glmnet_predict <-
function (m, x) 
{
    if (class(m) != "glmnet_fit") {
        stop("m must be a \"glmnet_fit\" object")
    }
    if (!is.matrix(x)) {
        stop("x must be a matrix")
    }
    y = matrix(x[, m$i], ncol = length(m$i)) %*% m$beta
    if (m$family != "cox") {
        y = m$a0 + y
    }
    if (m$family == "binomial") {
        y = 1/(1 + exp(-y))
    }
    y
}
glmnet_select <-
function (x) 
{
    if (!(is.list(x) && length(x) > 0 && class(x[[1]]) == "glmnet_fit")) {
        stop("x must be a list of objects of class \"glmnet_fit\" objects")
    }
    if (length(x) == 1) {
        warning("The list contains only one model")
        return(x[[1]])
    }
    vars = c()
    for (i in 1:length(x)) {
        vars = c(vars, x[[i]]$i)
    }
    vars = unique(sort(vars))
    models = matrix(0, ncol = length(vars), nrow = length(x))
    for (i in 1:length(x)) {
        models[i, match(x[[i]]$i, vars)] = 1
    }
    cat("[glmnet_select] - Calculating Dice coefficients...\n")
    dice = matrix(NA, ncol = nrow(models), nrow = nrow(models))
    for (i1 in 1:(nrow(models) - 1)) {
        for (i2 in (i1 + 1):nrow(models)) {
            dice_ij = 2 * sum(models[i1, ] * models[i2, ])/(sum(models[i1, 
                ]) + sum(models[i2, ]))
            dice[i1, i2] = dice_ij
            dice[i2, i1] = dice_ij
        }
    }
    mean_dice = apply(dice, 1, mean, na.rm = TRUE)
    cat("[glmnet_select] - Selecting the model with the highest Dice coefficient...\n")
    selected = which(mean_dice == max(mean_dice))
    y = x[[selected[1]]]
    if (length(selected) > 1) {
        for (i in 2:length(selected)) {
            y$a0 = c(y$a0, x[[selected[i]]]$a0)
            y$beta = rbind(y$beta, x[[selected[i]]]$beta)
        }
        if (!is.null(y$a0)) {
            y$a0 = mean(y$a0)
        }
        y$beta = apply(y$beta, 2, mean)
    }
    y
}
ijk2mni <-
function (ijk, sto_ijk) 
{
    ijk <- round(ijk)
    mni <- solve(sto_ijk) %*% rbind(matrix(ijk[1:3, ], nrow = 3), 
        1)
    as.matrix(mni[1:3, ])
}
impute.glmnet.matrix <-
function (m, x, nimp = 20) 
{
    if (class(m) != "impute.glmnet.matrix_fit") {
        stop("m must be a \"impute.glmnet.matrix_fit\" object")
    }
    if (!is.matrix(x)) {
        stop("x must be a matrix")
    }
    if (!(is.vector(nimp) && is.numeric(nimp) && length(nimp) == 
        1 && nimp > 0)) {
        stop("nimp must be a positive number")
    }
    start.time <- Sys.time()
    X_na = is.na(x)
    x.imp = list()
    for (imp in seq_len(nimp)) {
        x.imp[[imp]] = x
    }
    for (j in seq_len(ncol(x))) {
        X_na_j = X_na[, j]
        mj = m[[j]]
        family = mj$family
        for (i in 1:length(mj$imp.models)) {
            if (any(X_na_j)) {
                imp.model = mj$imp.models[[i]]
                if (!is.null(imp.model)) {
                  x.x = matrix(x[, -j], ncol = ncol(x) - 1)
                  cols.used = imp.model$i
                  rows.to_predict.complete = which(X_na_j & apply(x.x, 
                    1, function(tmp) {
                      all(cols.used %in% which(!is.na(tmp)))
                    }))
                  x.x.complete = matrix(x.x[rows.to_predict.complete, 
                    ], ncol = ncol(x) - 1)
                  x.y.to_predict = glmnet_predict(imp.model, 
                    x.x.complete)
                  for (imp in 1:nimp) {
                    if (family == "binomial") {
                      Ximp = as.numeric(x.y.to_predict > runif(length(rows.to_predict.complete)))
                    }
                    else {
                      Ximp = x.y.to_predict + rnorm(length(rows.to_predict.complete), 
                        0, mj$errors[i])
                    }
                    x.imp[[imp]][rows.to_predict.complete, j] = Ximp
                  }
                  X_na_j[rows.to_predict.complete][which(!is.na(x.y.to_predict))] = FALSE
                }
            }
        }
        if (any(X_na_j)) {
            for (imp in 1:nimp) {
                Ximp = sample(mj$data, sum(X_na_j), replace = TRUE)
                x.imp[[imp]][which(X_na_j), j] = Ximp
            }
        }
    }
    cat("[impute.glmnet.matrix] Running time:", difftime(Sys.time(), 
        start.time, units = "mins"), "minutes\n")
    x.imp
}
impute.glmnet.matrix_fit <-
function (x, n_cores = 1) 
{
    if (!is.matrix(x)) {
        stop("x must be a matrix")
    }
    cat("[impute.glmnet.matrix_fit] Estimating imputation models...\n")
    library(doSNOW)
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    iterations = ncol(x)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    start.time <- Sys.time()
    X_na = is.na(x)
    m = list()
    m = foreach(j = 1:ncol(x), .export = c("assign.folds", ".assign.folds_one_site", 
        "cv.glmnet", "glmnet.control", "glmnet", "glmnet_predict", 
        "glmnet_fit", ".assign.folds_simple"), .options.snow = opts) %dopar% 
        {
            X_na_j = X_na[, j]
            x.x = matrix(x[which(!X_na_j), -j], ncol = ncol(x) - 
                1)
            x.y = x[which(!X_na_j), j]
            family = ifelse(setequal(x.y, c(0, 1)), "binomial", 
                "gaussian")
            completeCols = apply(x.x, 1, function(tmp) {
                complete = which(!is.na(tmp))
                ifelse(length(complete) > 0, paste(complete, 
                  collapse = ","), NA)
            })
            completeCols = setdiff(unique(completeCols), NA)
            imp.models = list()
            errors = c()
            for (k in seq_len(length(completeCols))) {
                name = paste("col", j, "-set", k, sep = "")
                cols.complete = as.numeric(strsplit(completeCols[k], 
                  ",", fixed = TRUE)[[1]])
                rows.complete = which(apply(x.x, 1, function(tmp) {
                  all(cols.complete %in% which(!is.na(tmp)))
                }))
                x.x.complete = matrix(x.x[rows.complete, cols.complete], 
                  ncol = length(cols.complete))
                x.y.complete = x.y[rows.complete]
                if (length(unique(x.y.complete)) > 1 && ((family == 
                  "binomial" && length(table(x.y.complete, exclude = NULL)) == 
                  2 && min(table(x.y.complete, exclude = NULL)) > 
                  2) || (family == "gaussian" && (length(table(x.y.complete, 
                  exclude = NULL)) > 2 || (length(table(x.y.complete, 
                  exclude = NULL)) == 2 && min(table(x.y.complete, 
                  exclude = NULL)) > 2))))) {
                  imp.model = glmnet_fit(x.x.complete, x.y.complete, 
                    family)
                  imp.model$name = name
                  x.y.complete.pred = glmnet_predict(imp.model, 
                    x.x.complete)
                  imp.model$i = cols.complete[imp.model$i]
                  error = ifelse(family == "binomial", mean(((x.y.complete.pred > 
                    0.5) != x.y.complete)), sqrt(mean((x.y.complete.pred - 
                    x.y.complete)^2)))
                  names(error) = name
                }
                else {
                  imp.model = NULL
                  error = Inf
                }
                imp.models[[k]] = imp.model
                errors = c(errors, error)
            }
            setTxtProgressBar(pb, j)
            list(family = family, data = x.y, imp.models = imp.models[order(errors)], 
                errors = errors[order(errors)])
        }
    close(pb)
    stopCluster(cl)
    cat("[impute.glmnet.matrix_fit] Running time:", difftime(Sys.time(), 
        start.time, units = "mins"), "minutes\n")
    class(m) = "impute.glmnet.matrix_fit"
    m
}
it.sol <-
function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) 
{
    n <- apply(!is.na(sdat), 1, sum)
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while (change > conv) {
        g.new <- postmean(g.hat, g.bar, n, d.old, t2)
        sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
            1, sum, na.rm = T)
        d.new <- postvar(sum2, n, a, b)
        change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
        g.old <- g.new
        d.old <- d.new
        count <- count + 1
    }
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g.star", "d.star")
    adjust
}
launchApp <-
function () 
{
    shiny::runApp(system.file("shinyapp", package = "mripredict"))
}
load_preloaded_covB <-
function (path) 
{
    preloaded_covB = read.csv(file = path)
    ijk = preloaded_covB[, 1:3]
    covB = preloaded_covB[, 4:ncol(preloaded_covB)]
    list(covB = covB, ijk = ijk)
}
loadSampleMRI <-
function () 
{
}
mni2ijk <-
function (mni, sto_ijk) 
{
    coordinate <- round(sto_ijk %*% rbind(matrix(mni[1:3, ], 
        nrow = 3), 1))
    as.matrix(coordinate[1:3, ])
}
mripredict <-
function (mri_data = NULL, clinical_data, response_var, covariates = NULL, 
    predictor = NULL, response_family, modulation, information_variables = NULL, 
    mask_files_path = NULL) 
{
    if ((is.null(mri_data)) && (is.null(predictor))) {
        .check_stop("No MRI files or predictor variables found.")
    }
    mri_data_format = "character"
    if (!is.character(mri_data)) {
        if (is.data.frame(mri_data) && is.character(mri_data[1, 
            1])) {
            mri_data_format <- "dataframe_character"
        }
        else if (is.list(mri_data) && inherits(mri_data[[1]], 
            "mripredict_data")) {
            mri_data_format <- "mripredict_data_list"
        }
        else if (inherits(mri_data, "mripredict_data")) {
            mri_data_format <- "mripredict_data"
        }
    }
    .print_action("Setting new MRIPredict model")
    mri_paths = NULL
    mri_fu_paths = NULL
    mri_wm_paths = NULL
    mri_wmmod_paths = NULL
    mask_path = NULL
    mask_fu_path = NULL
    mask_wm_path = NULL
    mask_wmmod_path = NULL
    mri_un = NULL
    mri_fu = NULL
    mri_wm_un = NULL
    mri_wm_fu = NULL
    if (is.null(mri_data) & modulation != "clinical") {
        stop("No MRI data is selected, and the modulation is not set to 'clinical'. Please introduce MRI data or change modulation to 'clinical'.")
    }
    if (modulation != "clinical" && mri_data_format != "mripredict_data" && 
        mri_data_format != "mripredict_data_list") {
        if (!is.null(mri_data)) 
            mri_data = as.data.frame(mri_data)
        if (length(mri_data) == 1 && !all(grepl(".nii.gz|.nii", 
            mri_data))) {
            if (all(mri_data != "")) {
                if (file.exists(mri_data)) {
                  info = file.info(mri_data)
                  if (info$size > 0) {
                    mri_paths = .read_1col_file(mri_data)
                    if (!is.null(mask_files_path)) 
                      mask_path = mask_files_path[1]
                    if (modulation != "un") {
                      mri_fu_paths = .read_2col_file(mri_data)
                      if (!is.null(mask_files_path)) 
                        mask_fu_path = mask_files_path[2]
                    }
                    if (modulation == "all") {
                      mri_wm_paths = .read_3col_file(mri_data)
                      mri_wmmod_paths <- .read_4col_file(mri_data)
                      if (!is.null(mask_files_path)) {
                        mask_wm_path = mask_files_path[3]
                        mask_wmmod_path = mask_files_path[4]
                      }
                    }
                  }
                  else {
                    mri_data = NULL
                  }
                }
                else {
                  stop("MRI paths file does not exist.")
                }
            }
        }
        else if (dim(mri_data)[2] <= 4) {
            mri_data <- as.data.frame(mri_data)
            if (all(mri_data != "")) {
                mri_paths = mri_data[, 1]
                if (ncol(mri_data) > 1 && modulation != "un") 
                  mri_fu_paths = mri_data[, 2]
                if (ncol(mri_data) > 2 && modulation == "all") {
                  mri_wm_paths = mri_data[, 3]
                  mri_wmmod_paths = mri_data[, 4]
                }
            }
            if (!is.null(mask_files_path)) 
                mask_path = mask_files_path[1]
            if (!is.null(mask_files_path) && length(mask_files_path) > 
                1) 
                mask_fu_path = mask_files_path[2]
            if (!is.null(mask_files_path) && length(mask_files_path) > 
                2) 
                mask_wm_path = mask_files_path[3]
            if (!is.null(mask_files_path) && length(mask_files_path) > 
                3) 
                mask_wm_fu_path = mask_files_path[4]
        }
        else if (is.numeric(mri_data) && dim(mri_data)[2] > 4) {
            mri_paths <- mri_data
        }
        else if (is.list(mri_data)) {
            mri_paths <- mri_data$gm_un
            mri_fu_paths <- mri_data$gm_fu
            mri_wm_paths <- mri_data$wm_un
            mri_wmmod_paths <- mri_data$wm_fu
        }
        else {
            stop("Wrong MRI data format.")
        }
    }
    if (mri_data_format == "mripredict_data_list") {
        n_modalities = length(mri_data)
        mri_un = mri_data[[1]]
        if (n_modalities > 1) 
            mri_fu = mri_data[[2]]
        if (n_modalities > 2) 
            mri_wm_un = mri_data[[3]]
        if (n_modalities > 3) 
            mri_wm_fu = mri_data[[4]]
    }
    else if (mri_data_format == "mripredict_data") {
        mri_un = mri_data
    }
    if (!is.null(clinical_data) && any(clinical_data != "")) {
        if (length(clinical_data) == 1) {
            data_table = .read_data_table(clinical_data)
        }
        else {
            data_table = clinical_data
        }
    }
    else {
        data_table = NULL
    }
    if (!is.null(clinical_data)) 
        response = .check_response_var(data_table, response_var)
    else response = NULL
    response_levels = c()
    if (response_family == "cox") {
        response_levels = unique(response[, 2])
        if (length(response_levels) != 2) {
        }
    }
    else if ((response_family != "gaussian")) {
        response_levels = unique(response)
        if (length(response_levels) != 2) {
            mp = list(mri_paths = mri_paths, mri_fu_paths = mri_fu_paths, 
                mri_wm_paths = mri_wm_paths, mri_wmmod_paths = mri_wmmod_paths, 
                data_table = data_table, response_var = response_var, 
                response_ref = response_levels[1], response_event = response_levels[2], 
                modulation = modulation)
            attr(mp, "status") <- paste("Response variable", 
                response_var, "should have two levels")
            return(mp)
        }
    }
    mp = list(mri_paths = mri_paths, mri_fu_paths = mri_fu_paths, 
        mri_wm_paths = mri_wm_paths, mri_wmmod_paths = mri_wmmod_paths, 
        mri_loaded = F, mask_path = mask_path, mask_fu_path = mask_fu_path, 
        mask_wm_path = mask_wm_path, mask_wmmod_path = mask_wmmod_path, 
        data_table = data_table, response_var = response_var, 
        response_ref = response_levels[1], response_event = response_levels[2], 
        modulation = modulation)
    if (mri_data_format != "character" && mri_data_format != 
        "dataframe_character" && modulation != "clinical") {
        mp$mri_un = mri_un
        mp$mri_fu = mri_fu
        mp$mri_wm_un = mri_wm_un
        mp$mri_wm_fu = mri_wm_fu
        mp$mri_loaded = T
    }
    if (!is.null(covariates) & length(covariates) != 0) {
        data_cov = seleccio_covariables_data.frame2glmnet.matrix_fit(data_table, 
            covariates)
        covX_info_fit = data.frame2glmnet.matrix_fit(data_cov)
        covX_info = data.frame2glmnet.matrix(covX_info_fit, data_cov)
        mp$covX_transf = covX_info_fit
        mp$covX_var = colnames(covX_info)
    }
    else {
        mp$covX_transf = NULL
        mp$covX_var = NULL
    }
    if (is.null(information_variables) & (!is.null(covariates) | 
        !is.null(predictor))) {
        information_variables = c(covariates, predictor)
    }
    if (!is.null(information_variables) & (length(information_variables) != 
        0)) {
        data_info = seleccio_covariables_data.frame2glmnet.matrix_fit(data_table, 
            information_variables)
        data_info_fit = data.frame2glmnet.matrix_fit(data_info)
        data_info = data.frame2glmnet.matrix(data_info_fit, data_info)
        mp$data_table_transf = data_info_fit
        mp$data_table_var = information_variables
    }
    else {
        mp$data_table_transf = NULL
        mp$data_table_var = NULL
    }
    if (!is.null(predictor) & length(predictor) != 0) {
        data_pred = seleccio_covariables_data.frame2glmnet.matrix_fit(data_table, 
            predictor)
        pred_info_fit = data.frame2glmnet.matrix_fit(data_pred)
        pred_info = data.frame2glmnet.matrix(pred_info_fit, data_pred)
        mp$pred_transf = pred_info_fit
        mp$pred_var = colnames(pred_info)
    }
    else {
        mp$pred_transf = NULL
        mp$pred_var = NULL
    }
    mp$response_family = response_family
    if (!is.null(mri_paths)) {
        mp$mri_params = .load_mri(mri_paths = mri_paths[1])
        mp$mri_params$data <- NULL
    }
    class(mp) = "mripredict"
    .print_ok()
    attr(mp, "status") <- "OK"
    mp
}
mripredict_cv <-
function (mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, 
    preloaded_covB_fu_path = NULL, folds_file = NULL, n_cores = 1, 
    use_significant_voxels = FALSE, use_ensemble_voxels = FALSE, 
    use_ensemble_subjects = FALSE, n_folds = 10, ide_shiny = FALSE, 
    standardize_images = FALSE) 
{
    EXPERIMENTAL_THRESHOLD = F
    .require("glmnet")
    .require("oro.nifti")
    .require("survival")
    .require("doParallel")
    .require("logistf")
    .require("doSNOW")
    SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998)
    if (use_ensemble_voxels || use_ensemble_subjects) {
        N_ITERATIONS = 18
    }
    else {
        N_ITERATIONS = 1
    }
    N_M_IMPUTATIONS = 20
    new_folds = F
    if (is.null(folds_file)) {
        new_folds = T
    }
    if (!grepl("/", save_name, fixed = TRUE)) {
        if (!dir.exists(paste0("output/", save_name))) 
            dir.create(paste0("output/", save_name), recursive = T)
        save_name = paste0("output/", save_name, "/", save_name)
    }
    if (n_cores == "auto") {
        n_cores = max(round(detectCores()/2), 1)
    }
    mp$mask$data <- mp$mask_fu$data <- mp$mask_wm$data <- mp$mask_wmmmod$data <- NULL
    if (!is.null(mp$mask_path) && mp$mask_path != "") 
        mp$mask = .load_mri(mp$mask_path, space = "NO_CHECK")
    if (!is.null(mp$mask_fu_path) && mp$mask_fu_path != "") 
        mp$mask_fu = .load_mri(mp$mask_fu_path, space = "NO_CHECK")
    if (!is.null(mp$mask_wm_path) && mp$mask_wm_path != "") 
        mp$mask_wm = .load_mri(mp$mask_wm_path, space = "NO_CHECK")
    if (!is.null(mp$mask_wmmod_path) && mp$mask_wmmod_path != 
        "") 
        mp$mask_wmmod = .load_mri(mp$mask_wmmod_path, space = "NO_CHECK")
    .print_action("Loading MRI data")
    a = Sys.time()
    if (!mp$mri_loaded) {
        mri = NULL
        if (!is.null(mp$mri_paths)) {
            mri = .load_mri(mp$mri_paths, mask = mp$mask$data, 
                space = space)
        }
        mri_fu = NULL
        if (!mp$modulation == "un" && !is.null(mp$mri_fu_paths)) {
            mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, 
                space = space)
        }
        mri_wm = NULL
        mri_wm_fu = NULL
        if (!is.null(mp$mri_wmmod_paths) && mp$mri_wmmod_paths != 
            "") {
            mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, 
                space = space)
            mri_wm = .load_mri(mp$mri_wm_paths, mp$mask_wm$data, 
                space = space)
            mri_wm_fu = .load_mri(mp$mri_wmmod_paths, mp$mask_wmmod$data, 
                space = space)
        }
    }
    else {
        mri = mp$mri_un
        mri_fu = mp$mri_fu
        mri_wm = mp$mri_wm
        mri_wm_fu = mp$mri_wm_fu
    }
    cat("Time loading data:", round(difftime(Sys.time(), a, units = "mins"), 
        digits = 2), "mins.\n")
    .print_ok()
    .print_action("Creating response vector and covariate matrix")
    Y = c()
    if (mp$response_family == "cox") {
        Y = mp$data_table[, match(mp$response_var, colnames(mp$data_table))]
        Y = Surv(time = as.numeric(Y[, 1]), event = as.numeric(Y[, 
            2]))
    }
    else {
        Y = .create_Y(mp$data_table, mp$response_var, mp$response_event)
    }
    if (mp$response_family == "binomial" && min(table(Y)) < n_folds) {
        message("Warning: Too few samples in one outcome group. The program may be unable to perform the cross-validation.")
    }
    n_subjects = nrow(mp$data_table)
    covX = .create_covX(mp$data_table, mp$covX_transf)
    if (!is.null(mp$data_table_transf)) 
        data_informative_table = data.frame2glmnet.matrix(mp$data_table_transf, 
            mp$data_table)
    else data_informative_table = NULL
    .print_ok()
    if ("site" %in% colnames(mp$data_table)) {
        sites = factor(x = mp$data_table[, "site"])
    }
    else {
        sites = NULL
    }
    if (is.null(folds_file)) 
        folds_file = sprintf("%s_list_folds_n%s.txt", save_name, 
            n_folds)
    if (file.exists(folds_file) & !new_folds) {
        warning(sprintf("Previous folds file found: %s. The same previous folds distribution will be used.\nNote: To not use this folds file, provide a different name using the parameter save_name or remove the previous file.", 
            folds_file))
        list_folds = .read_folds_file(folds_file)
    }
    else {
        list_folds = NULL
    }
    if (is.null(list_folds)) {
        message("No folds file preloaded. New folds distribution will be created.")
        assigned_fold = assign.folds(y = Y, family = mp$response_family, 
            nfolds = n_folds, site = sites)
        write(assigned_fold, file = sprintf("%s_list_folds_n%s.txt", 
            save_name, n_folds), sep = ",", ncolumns = length(assigned_fold))
        is_folds_loaded = FALSE
    }
    else if (n_folds != length(unique(list_folds[[1]]))) {
        message("The folds file specified has different number of folds. New folds distribution will be created.")
        assigned_fold = assign.folds(y = Y, family = mp$response_family, 
            nfolds = n_folds, site = sites)
        write(assigned_fold, file = sprintf("%s_list_folds_n%s.txt", 
            save_name, n_folds), sep = ",", ncolumns = length(assigned_fold))
        is_folds_loaded = FALSE
    }
    else {
        assigned_fold = list_folds[[1]]
        if (length(list_folds) < (n_folds + 1)) {
            is_folds_loaded = FALSE
        }
        else {
            is_folds_loaded = TRUE
        }
        if (length(assigned_fold) != n_subjects) {
            stop("The number of subjects in the fold file, and the subjects selected are not the same. Please check that you selected the correct fold file, or the correct subjects.")
        }
    }
    if (!is.null(preloaded_covB_path) && preloaded_covB_path != 
        "") {
        preloaded_covB = .load_mri(preloaded_covB_path, space = "NO_CHECK")
        if (!mp$modulation == "un" && mp$mri_fu_paths != "") {
            preloaded_covB_fu = .load_mri(preloaded_covB_fu_path, 
                space = "NO_CHECK")
        }
        if (!(all(dim(preloaded_covB$data[, , , 1]) == dim(mri$data[, 
            , , 1])))) {
            stop("Preloaded effects do not have the same dimensions that the images.")
        }
    }
    else {
        preloaded_covB = NULL
    }
    cv_table = matrix(nrow = nrow(mp$data_table), ncol = 2)
    cv_table[, 1] = assigned_fold
    cv_accuracy = c()
    cv_betas = c()
    cv_table_predictions = data.frame()
    results = c()
    n_subjects = min(length(Y), nrow(Y))
    subjects_used = matrix(nrow = n_subjects, ncol = N_ITERATIONS)
    subjects_used = list()
    z_cuts = c()
    train_thresholds_all = c()
    mp$models = list()
    model_counter = 1
    time_points = c(30, 60, 180, 360, 720)
    time_folds = c()
    for (fold in 1:n_folds) {
        rdata_file = Sys.glob(sprintf("%s_fold*.Rdata", save_name))
        start_fold.time <- Sys.time()
        end_time = 0
        if (length(rdata_file) == 0 | new_folds) {
            cat(paste("Starting fold", fold, "of", n_folds), 
                "...\n")
            training = which(assigned_fold != fold)
            test = which(assigned_fold == fold)
            if (N_ITERATIONS > 1 && use_ensemble_subjects) {
                subjects_used[[iter]] = training[sample(1:length(training), 
                  replace = TRUE)]
                training = subjects_used[[iter]]
            }
            if (mp$response_family == "cox") {
                trainY = Y[training, ]
                testY = Y[test, ]
            }
            else {
                trainY = Y[training]
                testY = Y[test]
            }
            sites_training = sites[training]
            sites_test = sites[test]
            internal_folds = c()
            if (is_folds_loaded) {
                internal_folds = list_folds[[fold + 1]]
                if (length(internal_folds) != length(training)) {
                  stop("The number of folds in the provided folds file seems to not correspond to the current data. Please specify a different save_name or delete the previous folds file.")
                }
            }
            else {
                internal_folds = assign.folds(y = trainY, family = mp$response_family, 
                  nfolds = n_folds, site = sites_training)
                write(assigned_fold, file = sprintf("%s_list_folds_n%s.txt", 
                  save_name, n_folds), sep = ",", ncolumns = length(assigned_fold))
                write(internal_folds, file = sprintf("%s_list_folds_n%s.txt", 
                  save_name, n_folds), sep = ",", ncolumns = length(internal_folds), 
                  append = TRUE)
            }
            n_multiple_imputations = 1
            if (!is.null(data_informative_table)) {
                data_table_imputed_train = data_informative_table[training, 
                  , drop = F]
                data_table_imputed_test = data_informative_table[test, 
                  , drop = F]
                if (any(is.na(data_informative_table))) {
                  n_multiple_imputations = N_M_IMPUTATIONS
                  if (length(Sys.glob(sprintf("%s_fold%d_iteration*.txt", 
                    save_name, fold))) != N_ITERATIONS) {
                    impute_obj = impute.glmnet.matrix_fit(x = data_table_imputed_train, 
                      n_cores = n_cores)
                    imp.data_informative_table_train = impute.glmnet.matrix(m = impute_obj, 
                      x = data_table_imputed_train, nimp = N_M_IMPUTATIONS)
                    imp.data_informative_table_test = impute.glmnet.matrix(m = impute_obj, 
                      x = data_table_imputed_test, nimp = N_M_IMPUTATIONS)
                  }
                }
            }
            else {
                data_table_imputed_train = NULL
                imp.data_informative_table_train = NULL
                data_table_imputed_test = NULL
                imp.data_informative_table_test = NULL
            }
            train_linear_predictor_to_find_the_threshold_all = vector(mode = "list", 
                length = length(time_points))
            for (iter in 1:N_ITERATIONS) {
                name_base <- sprintf("%s_model_FOLD_%d", save_name, 
                  fold)
                fold_rds_name = sprintf("%s_mp.rds", name_base)
                if (!file.exists(fold_rds_name) | new_folds) {
                  linPreds = c()
                  if (mp$response_family == "cox") {
                    train_thresholds = NULL
                  }
                  start.time <- Sys.time()
                  if (N_ITERATIONS > 1) 
                    cat("Fold:", fold, ". Iteration ", iter, 
                      "of", N_ITERATIONS, "\n", save_name)
                  sprintf("/nIteration: %d", iter)
                  preds = c()
                  time_imp = c()
                  start_imputation.time = Sys.time()
                  last_imputation.time = start_imputation.time
                  for (iter_imputation in 1:n_multiple_imputations) {
                    gc()
                    if (n_multiple_imputations > 1) {
                      data_table_imputed_train = imp.data_informative_table_train[[iter_imputation]]
                    }
                    if (n_multiple_imputations > 1) {
                      .print_action(paste("Fold:", fold, "Imputation", 
                        iter_imputation, "of", n_multiple_imputations, 
                        ". Iteration:", iter, "/", N_ITERATIONS, 
                        "\n", save_name))
                      if (iter_imputation > 1) {
                        time_per_imp = round(difftime(Sys.time(), 
                          start_imputation.time), digits = 2)
                        time_imp = c(time_imp, difftime(Sys.time(), 
                          start_imputation.time, units = "secs"))
                        seconds_remaining_imput = (mean(time_imp)/iter_imputation) * 
                          (n_multiple_imputations - iter_imputation)
                        seconds_all = (mean(time_imp)/iter_imputation) * 
                          n_multiple_imputations * N_ITERATIONS * 
                          n_folds
                        units = attr(time_per_imp, "units")
                        message("\n[Imputation progress] - Time per imputation: ", 
                          time_per_imp, " ", units, ".")
                        remaining_imp_time = (n_multiple_imputations - 
                          iter_imputation) * mean(time_imp)
                        message("[Imputation progress] - Estimated remaining imputations time for this fold: ", 
                          round(remaining_imp_time/60, digits = 2), 
                          " mins. End at: ", Sys.time() + seconds_remaining_imput)
                        if (fold == 1) {
                          message("[CV] - Estimated end time:", 
                            Sys.time() + seconds_all)
                        }
                        else {
                          message("[CV] - Estimated end time:", 
                            end_time)
                        }
                      }
                    }
                    model_list = fit_model(mp = mp, data_informative_table = data_table_imputed_train, 
                      Y = trainY, mri = mri$data[, , , training], 
                      mri_fu = mri_fu$data[, , , training], mri_wm = mri_wm$data[, 
                        , , training], mri_fu_wm = mri_wm_fu$data[, 
                        , , training], preloaded_covB = preloaded_covB, 
                      preloaded_covB_fu = preloaded_covB_fu, 
                      iter = ifelse(use_ensemble_voxels == FALSE, 
                        -1, iter), internal_folds = internal_folds, 
                      n_cores = n_cores, use_significant_voxels = use_significant_voxels, 
                      covX_site = sites_training, standardize_images = standardize_images)
                    mp$combat = model_list$combat
                    scale_clinical = model_list$scale_predictors
                    scale_mri = model_list$scale_mri
                    if (mp$response_family == "cox" & EXPERIMENTAL_THRESHOLD) {
                      predX_training = matrix(data_table_imputed_train[, 
                        mp$pred_var], nrow = nrow(data_table_imputed_train))
                      covX_training = cbind(1, as.matrix(data_table_imputed_train[, 
                        mp$covX_var]))
                      for (i_time_point in 1:length(time_points)) {
                        status0_at_time_point = which((trainY[, 
                          1] == time_points[i_time_point] & trainY[, 
                          2] == 0) | (trainY[, 1] > time_points[i_time_point]))
                        status1_at_time_point = which(trainY[, 
                          1] <= time_points[i_time_point] & trainY[, 
                          2] == 1)
                        if (length(status0_at_time_point) > 0 && 
                          length(status1_at_time_point) > 0) {
                          any_status_at_time_point = c(status0_at_time_point, 
                            status1_at_time_point)
                          training_for_threshold = training[any_status_at_time_point]
                          train_linear_predictor_to_find_the_threshold = apply_model(mp = mp, 
                            mri_data = mri$data[, , , training_for_threshold], 
                            mri_fu_data = mri_fu$data[, , , training_for_threshold], 
                            mri_wm_data = mri_wm$data[, , , training_for_threshold], 
                            mri_wm_fu_data = mri_wm_fu$data[, 
                              , , training_for_threshold], covX_test = covX_training[any_status_at_time_point, 
                              ], signif_indx = signif_indx, lasso_covB = model_list$lasso_covB, 
                            lasso_covB_fu = model_list$lasso_covB_fu, 
                            mask = model_list$mask, predX_test = predX_training[any_status_at_time_point, 
                              ], scale_clinical = scale_clinical, 
                            scale_mri = scale_mri, lasso = model_list$lasso, 
                            lasso_predX_indx = model_list$lasso_predX_indx, 
                            tipett_take_un = tipett_take_un, 
                            img_kappa = NULL, use_significant_voxels = use_significant_voxels, 
                            covX_site = sites_training[any_status_at_time_point], 
                            masks_3d = model_list$masks_3d, n_voxels_mask = model_list$n_voxels_mask, 
                            combat = model_list$combat, standardize_images = standardize_images)
                          train_linear_predictor_to_find_the_threshold_all[[i_time_point]] = cbind(train_linear_predictor_to_find_the_threshold_all[[i_time_point]], 
                            train_linear_predictor_to_find_the_threshold)
                        }
                      }
                    }
                    mp$models[[model_counter]] <- model_list
                    model_list <- NULL
                    model_counter <- model_counter + 1
                    .print_ok()
                    .print_action("Test sample: applying the model")
                    test_preds = c()
                    for (iter_imputation_test in 1:n_multiple_imputations) {
                      if (n_multiple_imputations > 1) {
                        data_table_imputed_test = as.matrix(imp.data_informative_table_test[[iter_imputation_test]])
                      }
                      if (!is.null(mp$covX_var)) {
                        covX_test = cbind(1, matrix(data_table_imputed_test[, 
                          mp$covX_var], nrow = nrow(data_table_imputed_test)))
                      }
                      else {
                        covX_test = matrix(1, nrow = length(test))
                      }
                      if (!is.null(mp$pred_var)) {
                        predX_test = matrix(data_table_imputed_test[, 
                          mp$pred_var], nrow = nrow(data_table_imputed_test))
                      }
                      else {
                        predX_test = NULL
                      }
                      preds = apply_model(mp = mp, mri_data = mri$data[, 
                        , , test, drop = FALSE], mri_fu_data = mri_fu$data[, 
                        , , test, drop = FALSE], mri_wm_data = mri_wm$data[, 
                        , , test, drop = FALSE], mri_wm_fu_data = mri_wm_fu$data[, 
                        , , test, drop = FALSE], covX_test = covX_test, 
                        signif_indx = signif_indx, lasso_covB = mp$models[[model_counter - 
                          1]]$lasso_covB, lasso_covB_fu = mp$models[[model_counter - 
                          1]]$lasso_fu_covB, mask = mp$models[[model_counter - 
                          1]]$mask, predX_test = predX_test, 
                        scale_clinical = mp$models[[model_counter - 
                          1]]$scale_clinical, scale_mri = mp$models[[model_counter - 
                          1]]$scale_mri, lasso = mp$models[[model_counter - 
                          1]]$lasso, lasso_predX_indx = mp$models[[model_counter - 
                          1]]$lasso_predX_indx, tipett_take_un = mp$models[[model_counter - 
                          1]]$tipett_take_un, img_kappa = NULL, 
                        use_significant_voxels = use_significant_voxels, 
                        covX_site = sites_test, masks_3d = mp$models[[model_counter - 
                          1]]$masks_3d, n_voxels_mask = mp$models[[model_counter - 
                          1]]$n_voxels_mask, standardize_images = standardize_images, 
                        combat = mp$models[[model_counter - 1]]$combat)
                      test_preds = cbind(test_preds, preds)
                    }
                    preds = rowMeans(test_preds)
                    linPreds = cbind(linPreds, preds)
                    if (ide_shiny) 
                      incProgress(1/(n_folds * N_ITERATIONS * 
                        n_multiple_imputations), detail = paste("Doing crossvalidation. Fold", 
                        fold, ". Iteration:", iter, ". Imputation:", 
                        iter_imputation))
                  }
                  linPred = matrix(rowMeans(linPreds))
                  .print_ok()
                  .print_action("Saving the predictions")
                  cat("\n")
                  pred = c()
                  if (mp$response_family == "binomial") {
                    prob = 1/(1 + exp(-linPred))
                    pred = prob
                    bac = .metrics_binary(testY, pred > (sum(trainY == 
                      1)/length(trainY)))$bac
                    cv_accuracy[fold] = bac
                    cv_table[test, 2] = prob
                  }
                  else if (mp$response_family == "gaussian") {
                    pred = as.matrix(linPred)
                    cv_accuracy[fold] = sqrt(mean((pred - testY)^2))
                    cv_table[test, 2] = pred
                  }
                  else if (mp$response_family == "cox") {
                    pred = linPred
                    cv_accuracy[fold] = NA
                    cv_betas[fold] = NA
                    cv_table[test, 2] = pred
                    results = data.frame(id = test, linear_predictor = pred[, 
                      1], time = testY[, 1], status = testY[, 
                      2])
                    if (n_multiple_imputations > 1) 
                      print("Mean results for multiple imputations:")
                  }
                  if (mp$response_family == "cox") {
                    cv_table_predictions = rbind(cv_table_predictions, 
                      cbind(data.frame(id = test, linear_predictor = pred[, 
                        1], rd2_beta = NA, time = testY[, 1], 
                        status = testY[, 2], fold = fold, iteration = iter)))
                  }
                  else if (mp$response_family == "gaussian") {
                    cv_table_predictions = rbind(cv_table_predictions, 
                      data.frame(id = test, linear_predictor = pred[, 
                        1], response = testY, fold = fold, iteration = iter))
                  }
                  else {
                    cv_table_predictions = rbind(cv_table_predictions, 
                      data.frame(id = test, linear_predictor = linPred[, 
                        1], prob = pred[, 1], response = testY, 
                        fold = fold, iteration = iter))
                  }
                  .print_ok()
                  switch(mp$response_family, binomial = message("BAC: ", 
                    round(cv_accuracy[fold], digits = 2)), gaussian = message("RMSE: ", 
                    cv_accuracy[fold]), cox = {
                  })
                  end.time <- Sys.time()
                  if (iter > 1) 
                    message("Time per iteration: ", round(difftime(end.time, 
                      start.time, units = "mins"), digits = 2), 
                      " mins.")
                }
                else {
                  cat("Skipping fold: ", fold, " Iteration: ", 
                    iter, "\n")
                  res_csv = read.csv(sprintf("%s_iteration%d_%s_results_fold%d.csv", 
                    save_name, iter, mp$response_family, fold))
                  cv_table_predictions = rbind(cv_table_predictions, 
                    res_csv)
                  name_base <- sprintf("%s_model_FOLD_%d", save_name, 
                    fold)
                  mp <- readRDS(file = sprintf("%s_mp.rds", name_base))
                }
            }
            cat("\n[Saving] - Saving fold model to", sprintf("%s_model_FOLD_%d", 
                save_name, fold))
            name_base <- sprintf("%s_model_FOLD_%d", save_name, 
                fold)
            if (fold > 1) {
                name_base_previous <- sprintf("%s_model_FOLD_%d", 
                  save_name, fold - 1)
                file.remove(sprintf("%s_mp.rds", name_base_previous))
            }
            saveRDS(object = mp, file = sprintf("%s_mp.rds", 
                name_base))
            if (mp$response_family == "cox") {
                trainY = Y[training, ]
                train_thresholds_row = c()
                for (i_time_point in 1:length(time_points)) {
                  status0_at_time_point = which((trainY[, 1] == 
                    time_points[i_time_point] & trainY[, 2] == 
                    0) | (trainY[, 1] > time_points[i_time_point]))
                  status1_at_time_point = which(trainY[, 1] <= 
                    time_points[i_time_point] & trainY[, 2] == 
                    1)
                  if (!is.null(train_linear_predictor_to_find_the_threshold_all[[i_time_point]])) {
                    train_linear_predictor_to_find_the_threshold = rowMeans(train_linear_predictor_to_find_the_threshold_all[[i_time_point]])
                  }
                  else {
                    train_linear_predictor_to_find_the_threshold = NA
                  }
                  sorted_train_linear_predictor_to_find_the_threshold = sort(unique(train_linear_predictor_to_find_the_threshold))
                  thresholds = (sorted_train_linear_predictor_to_find_the_threshold[-length(sorted_train_linear_predictor_to_find_the_threshold)] + 
                    sorted_train_linear_predictor_to_find_the_threshold[-1])/2
                  bacs = c()
                  for (threshold in thresholds) {
                    true_status = c(rep(0, length(status0_at_time_point)), 
                      rep(1, length(status1_at_time_point)))
                    predicted_status = 1 * (train_linear_predictor_to_find_the_threshold > 
                      threshold)
                    sensitivity = sum(predicted_status == 1 & 
                      true_status == 1)/sum(true_status == 1)
                    specificity = sum(predicted_status == 0 & 
                      true_status == 0)/sum(true_status == 0)
                    bac = (sensitivity + specificity)/2
                    bacs = c(bacs, bac)
                  }
                  train_thresholds_row = c(train_thresholds_row, 
                    thresholds[which.max(bacs)])
                }
                train_thresholds_all = rbind(train_thresholds_all, 
                  train_thresholds_row)
            }
            time_per_fold = round(difftime(Sys.time(), start_fold.time), 
                digits = 2)
            time_folds = c(time_folds, difftime(Sys.time(), start_fold.time, 
                units = "secs"))
            units = attr(time_per_fold, "units")
            message("\n[Fold ended] - Time per fold: ", time_per_fold, 
                " ", units, ".")
            remaining_time = (n_folds - fold) * mean(time_folds)
            end_time = Sys.time() + remaining_time
            message("[Fold ended] - Estimated remaining time: ", 
                round(remaining_time/60, digits = 2), " mins. End at: ", 
                end_time)
        }
        else {
            fold = gsub(sprintf("%s_fold", save_name), "", rdata_file)
            fold = max(as.numeric(gsub(".Rdata", "", fold)))
            load(sprintf("%s_fold%s.Rdata", save_name, fold))
        }
    }
    if (!is.null(train_thresholds_all)) 
        write.csv(train_thresholds_all, sprintf("%s_train_thresholds.csv", 
            save_name), row.names = FALSE)
    if (mp$response_family == "cox") {
        mean_pred_subjs = data.frame(id = sort(unique(cv_table_predictions$id)), 
            lin_pred = tapply(cv_table_predictions$linear_predictor, 
                cv_table_predictions$id, mean), times = tapply(cv_table_predictions$time, 
                cv_table_predictions$id, unique), status = tapply(cv_table_predictions$status, 
                cv_table_predictions$id, unique), fold = tapply(cv_table_predictions$fold, 
                cv_table_predictions$id, unique))
        n_pred_subjs_cols = ncol(mean_pred_subjs)
        th_i = 1
        for (time_point in time_points) {
            for (fold in 1:n_folds) {
                mean_pred_subjs[mean_pred_subjs$fold == fold, 
                  n_pred_subjs_cols + th_i] <- NA
                if (ncol(train_thresholds_all) >= th_i) {
                  mean_pred_subjs[mean_pred_subjs$fold == fold, 
                    n_pred_subjs_cols + th_i] <- mean_pred_subjs$lin_pred[mean_pred_subjs$fold == 
                    fold] > train_thresholds_all[fold, th_i]
                }
                colnames(mean_pred_subjs)[n_pred_subjs_cols + 
                  th_i] <- sprintf("threshold_%s", time_point)
            }
            th_i = th_i + 1
        }
        mp$cv_results = mean_pred_subjs
    }
    else {
        mp$cv_results = cv_table_predictions
    }
    mp$cv_table = cv_table
    mp$cv_accuracy = cv_accuracy
    mp$cv_betas = cv_betas
    mp$subjects_used = subjects_used
    cat("\n[End] - FOLDS performance:", round(cv_accuracy, digits = 2), 
        "\n")
    flush.console()
    switch(mp$response_family, binomial = {
        bin_threshold = sum(cv_table_predictions$response == 
            1)/length(cv_table_predictions$response)
        message("Mean BAC: ", round(.metrics_binary(cv_table_predictions$response, 
            cv_table_predictions$prob > bin_threshold)$bac, digits = 2))
    }, gaussian = message("Mean RMSE: ", sqrt(mean((cv_table_predictions$linear_predictor - 
        cv_table_predictions$response)^2))), cox = {
        final_rd2 = .coxph_RD2(predictor = mp$cv_results$lin_pred, 
            stime = mp$cv_results$time, sevent = mp$cv_results$status)
        message("Final RD2: ", final_rd2$RD2, " Beta:", final_rd2$b)
    })
    if (mp$response_family == "cox") {
        mp$frontier_time = .find_best_time(trainY[, 1], trainY[, 
            2])
        mp$metrics = .metrics_cox(cv_table_predictions, mp$frontier_time, 
            save = FALSE)
        write.csv(mp$metrics, sprintf("%s_custom_metric.txt", 
            save_name))
        write.csv(.coxph_RD2(predictor = mp$cv_results$lin_pred, 
            stime = mp$cv_results$time, sevent = mp$cv_results$status), 
            sprintf("%s_rd2.txt", save_name))
        print(mp$metrics)
    }
    else if (mp$response_family == "gaussian") {
        colnames(mp$cv_results) = c("id", "pred", "real", "fold", 
            "iteration")
        mp$cv_results = aggregate(mp$cv_results[, 2:4], list(mp$cv_results$id), 
            mean)
        colnames(mp$cv_results) = c("id", "prediction", "real", 
            "fold")
        mp$metrics = data.frame(Measure = "RMSE", Value = sqrt(mean((mp$cv_results$prediction - 
            mp$cv_results$real)^2)))
        cat(sprintf("Mean RMSE: %f", mp$metrics), file = sprintf("%s_gaussian.txt", 
            save_name))
    }
    else {
        write.csv(mp$cv_results, sprintf("%s_bin_raw.csv", save_name), 
            row.names = F)
        mp$cv_results = aggregate(mp$cv_results[, 2:5], list(mp$cv_results$id), 
            mean)
        colnames(mp$cv_results) = c("id", "linear_predictor", 
            "probability", "real", "fold")
        mp$metrics <- .metrics_binary(mp$cv_results$real, mp$cv_results$probability > 
            bin_threshold)
        colnames(mp$cv_results)[3] <- sprintf("Prediction (0 = %s; 1 = %s)", 
            mp$response_ref, mp$response_event)
        colnames(mp$cv_results)[4] <- sprintf("Real (0 = %s; 1 = %s)", 
            mp$response_ref, mp$response_event)
        cat(sprintf("Mean BAC: %f", mp$metrics$BAC, file = sprintf("%s_bin.txt", 
            save_name)))
    }
    cat("\n[Saving] - Saving MRIPredict object to:", sprintf("%s_mp.rds", 
        save_name))
    mp$cv_results = mp$cv_results[order(mp$cv_results$id), ]
    saveRDS(mp, file = sprintf("%s_mp.rds", save_name))
    .most_frequent_variables(mp$models, mp = mp, file = sprintf("%s_frequent_variables.csv", 
        save_name))
    write.csv(mp$cv_results, sprintf("%s_subject_linear_predictors.csv", 
        save_name))
    cat("\n[DONE] - CV finished:")
    mp
}
mripredict_fit <-
function (mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, 
    preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1, 
    use_significant_voxels = FALSE, use_ensemble_learning = FALSE, 
    use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, 
    ide_shiny = FALSE, standardize_images = FALSE) 
{
    .require("glmnet")
    .require("oro.nifti")
    .require("survival")
    .require("doParallel")
    .require("logistf")
    data_informative_table = NULL
    imp.data_informative_table_train = NULL
    imp.data_informative_table_test = NULL
    SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998)
    if (use_ensemble_learning) {
        N_ITERATIONS = 18
    }
    else {
        N_ITERATIONS = 1
    }
    n_multiple_imputations = 1
    N_M_IMPUTATIONS = 20
    if (n_cores == "auto") {
        n_cores = max(round(detectCores()/2), 1)
    }
    mp$mask$data <- mp$mask_fu$data <- mp$mask_wm$data <- mp$mask_wmmmod$data <- NULL
    if (!is.null(mp$mask_path)) 
        mp$mask = .load_mri(mp$mask_path, space = "NO_CHECK")
    if (!is.null(mp$mask_fu_path)) 
        mp$mask_fu = .load_mri(mp$mask_fu_path, space = "NO_CHECK")
    if (!is.null(mp$mask_wm_path)) 
        mp$mask_wm = .load_mri(mp$mask_wm_path, space = "NO_CHECK")
    if (!is.null(mp$mask_wmmod_path)) 
        mp$mask_wmmod = .load_mri(mp$mask_wmmod_path, space = "NO_CHECK")
    .print_action("Loading MRI data")
    a = Sys.time()
    mri = NULL
    if (!is.null(mp$mri_paths)) {
        mri = .load_mri(mp$mri_paths, mask = mp$mask$data, space = space)
    }
    mri_fu = NULL
    if (!mp$modulation == "un" && !is.null(mp$mri_fu_paths)) {
        mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, 
            space = space)
    }
    mri_wm = NULL
    mri_wm_fu = NULL
    if (!is.null(mp$mri_wmmod_paths)) {
        mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, 
            space = space)
        mri_wm = .load_mri(mp$mri_wm_paths, mp$mask_wm$data, 
            space = space)
        mri_wm_fu = .load_mri(mp$mri_wmmod_paths, mp$mask_wmmod$data, 
            space = space)
    }
    cat("Time loading data:", difftime(Sys.time(), a, units = "mins"), 
        "mins.\n")
    .print_ok()
    .print_action("Creating response vector and covariate matrix")
    Y = c()
    if (mp$response_family == "cox") {
        Y = mp$data_table[, match(mp$response_var, colnames(mp$data_table))]
        Y = Surv(time = as.numeric(Y[, 1]), event = as.numeric(Y[, 
            2]))
    }
    else {
        Y = .create_Y(mp$data_table, mp$response_var, mp$response_event)
    }
    n_subjects = nrow(mp$data_table)
    if (!is.null(mp$covX_transf)) 
        covX = .create_covX(mp$data_table, mp$covX_transf)
    if (!is.null(mp$data_table_transf)) 
        data_informative_table = data.frame2glmnet.matrix(mp$data_table_transf, 
            mp$data_table)
    .print_ok()
    if ("site" %in% colnames(mp$data_table)) {
        sites = factor(x = mp$data_table[, "site"])
    }
    else {
        sites = NULL
    }
    n_subjects = min(length(Y), nrow(Y))
    subjects_used = matrix(nrow = n_subjects, ncol = N_ITERATIONS)
    subjects_used = list()
    z_cuts = c()
    train_thresholds_all = c()
    mp$models = list()
    model_counter = 1
    time_points = c(30, 60, 180, 360, 720)
    if (!is.null(data_informative_table)) {
        data_table_imputed_train = data_informative_table
        if (any(is.na(data_informative_table))) {
            n_multiple_imputations = N_M_IMPUTATIONS
            impute_obj = impute.glmnet.matrix_fit(x = data_table_imputed_train, 
                n_cores = 4)
            imp.data_informative_table_train = impute.glmnet.matrix(m = impute_obj, 
                x = data_table_imputed_train, nimp = N_M_IMPUTATIONS)
            mp$impute_obj = impute_obj
        }
        else {
            n_multiple_imputations = 1
            mp$impute_obj = NULL
        }
    }
    else {
        data_table_imputed_train = NULL
        imp.data_informative_table_train = NULL
    }
    start_fold.time <- Sys.time()
    mp$n_iterations = N_ITERATIONS
    mp$n_imputations = n_multiple_imputations
    for (iter in 1:N_ITERATIONS) {
        if (N_ITERATIONS > 1 && use_ensemble_subjects) {
            subjects_used[[iter]] = training[sample(1:length(training), 
                replace = TRUE)]
            training = subjects_used[[iter]]
        }
        if (mp$response_family == "cox") {
            trainY = Y
        }
        else {
            trainY = Y
            if (mp$response_family == "binomial") {
                mp$bin_threshold = sum(trainY == 1)/length(trainY)
            }
        }
        linPreds = c()
        if (mp$response_family == "cox") {
            train_thresholds = NULL
        }
        start.time <- Sys.time()
        if (N_ITERATIONS > 1) 
            cat("Iteration ", iter, "of", N_ITERATIONS, "\n", 
                save_name)
        sprintf("/nIteration: %d", iter)
        preds = c()
        for (iter_imputation in 1:n_multiple_imputations) {
            if (n_multiple_imputations > 1 & !is.null(imp.data_informative_table_train)) {
                data_table_imputed_train = imp.data_informative_table_train[[iter_imputation]]
            }
            if (n_multiple_imputations > 1) 
                .print_action(paste("Imputation", iter_imputation, 
                  "of", n_multiple_imputations, "\n", save_name))
            internal_folds = c()
            internal_folds = assign.folds(y = trainY, family = mp$response_family, 
                nfolds = 10, site = sites)
            model_list = fit_model(mp = mp, data_informative_table = data_table_imputed_train, 
                Y = trainY, mri = mri$data, mri_fu = mri_fu$data, 
                mri_wm = mri_wm$data, mri_fu_wm = mri_wm_fu$data, 
                iter = ifelse(use_ensemble_voxels == FALSE, -1, 
                  iter), internal_folds = internal_folds, n_cores = n_cores, 
                use_significant_voxels = use_significant_voxels, 
                covX_site = sites, standardize_images = standardize_images)
            if (mp$modulation %in% c("cl", "clinical")) {
                n_voxels_mask = 0
            }
            tipett_take_un = model_list$tipett_take_un
            scale_clinical = model_list$scale_predictors
            scale_mri = model_list$scale_mri
            mp$models[[model_counter]] <- model_list
            model_counter <- model_counter + 1
            .print_ok()
            if (ide_shiny) 
                incProgress(1/(n_folds * N_ITERATIONS * n_multiple_imputations), 
                  detail = paste("Doing crossvalidation. Fold", 
                    fold, ". Iteration:", iter, ". Imputation:", 
                    iter_imputation))
        }
        .print_ok()
        end.time <- Sys.time()
        if (iter > 1) 
            message("Time per iteration:", difftime(end.time, 
                start.time, units = "mins"), " mins.")
    }
    cat("\n[Saving] - Saving MRIPredict object to:", sprintf("%s_mp.rds", 
        save_name))
    saveRDS(mp, file = sprintf("%s_mp.rds", save_name))
    mp
}
mripredict_predict <-
function (mp, mri_data = NULL, mri_fu_paths_file = NULL, clinical_data = NULL, 
    space, n_cores = 1) 
{
    .require("glmnet")
    .require("oro.nifti")
    .require("survival")
    .require("doParallel")
    .require("parallel")
    .require("logistf")
    SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998)
    N_ITERATIONS = mp$n_iterations
    N_M_IMPUTATIONS = mp$n_imputations
    if (n_cores == "auto") {
        n_cores = 2
    }
    .print_action("Setting new MRIPredict model")
    mp_test = mripredict(mri_data = mri_data, clinical_data = clinical_data, 
        response_var = mp$response_var, covariates = mp$covX_var, 
        predictor = mp$pred_var, response_family = mp$response_family, 
        modulation = mp$modulation)
    .print_action("Loading MRI data")
    a = Sys.time()
    mri = NULL
    if (!is.null(mp_test$mri_paths)) {
        mri = .load_mri(mp_test$mri_paths, mask = mp$mask$data, 
            space = space)
    }
    mri_fu = NULL
    if (!mp_test$modulation == "un" && !is.null(mp_test$mri_fu_paths)) {
        mri_fu = .load_mri(mp_test$mri_fu_paths, mp$mask_fu$data, 
            space = space)
    }
    mri_wm = NULL
    mri_wm_fu = NULL
    if (!is.null(mp_test$mri_wmmod_paths)) {
        mri_fu = .load_mri(mp_test$mri_fu_paths, mp$mask_fu$data, 
            space = space)
        mri_wm = .load_mri(mp_test$mri_wm_paths, mp$mask_wm$data, 
            space = space)
        mri_wm_fu = .load_mri(mp_test$mri_wmmod_paths, mp$mask_wmmod$data, 
            space = space)
    }
    cat("Time loading data:", difftime(Sys.time(), a, units = "mins"), 
        "mins.\n")
    .print_ok()
    .print_action("Creating response vector and covariate matrix")
    n_subjects = nrow(clinical_data)
    if (!is.null(mp$covX_transf)) 
        covX = .create_covX(clinical_data, mp$covX_transf)
    sites = NULL
    if (!is.null(mp$data_table_transf)) {
        data_informative_table = data.frame2glmnet.matrix(m = mp$data_table_transf, 
            x = clinical_data)
        if ("site" %in% colnames(clinical_data)) {
            sites = factor(x = clinical_data[, "site"])
        }
    }
    else {
        data_informative_table = NULL
        imp.data_informative_table_test = NULL
    }
    .print_ok()
    if (!is.null(mp$impute_obj)) {
        imp.data_informative_table_test = impute.glmnet.matrix(m = mp$impute_obj, 
            x = data_informative_table, nimp = N_M_IMPUTATIONS)
    }
    linPreds = c()
    i = 1
    for (iter in 1:N_ITERATIONS) {
        for (iter_imputation in 1:N_M_IMPUTATIONS) {
            if (!is.null(mp$impute_obj)) 
                data_informative_table = imp.data_informative_table_test[[iter_imputation]]
            model_list = mp$models[[i]]
            i = i + 1
            lasso_ijk = model_list$lasso_ijk
            signif_indx = model_list$signif_indx
            lasso_covB = model_list$lasso_covB
            lasso_covB_fu = model_list$lasso_fu_covB
            mask = model_list$mask
            scale_clinical = model_list$predictors_scale
            if (!is.null(data_informative_table)) {
                if (!is.null(mp$covX_var)) {
                  covX_test = cbind(1, matrix(data_informative_table[, 
                    mp$covX_var], nrow = nrow(data_informative_table)))
                }
                else {
                  covX_test = data_informative_table
                }
            }
            else {
                if (is.null(clinical_data)) {
                  covX_test = matrix(1, nrow = mri$n)
                }
                else {
                  covX_test = matrix(1, nrow = nrow(clinical_data))
                }
            }
            if (!is.null(mp$pred_var)) {
                predX_test = matrix(data_informative_table[, 
                  mp$pred_var], nrow = nrow(data_informative_table))
            }
            tipett_take_un = model_list$tipett_take_un
            preds = apply_model(mp = mp, mri_data = mri$data, 
                mri_fu_data = mri_fu$data, mri_wm_data = mri_wm$data, 
                mri_wm_fu_data = mri_wm_fu$data, covX_test = covX_test, 
                signif_indx = model_list$signif_indx, lasso_covB = lasso_covB, 
                lasso_covB_fu = lasso_covB_fu, mask = model_list$mask, 
                masks_3d = model_list$masks_3d, predX_test = predX_test, 
                scale_clinical = model_list$scale_predictors, 
                scale_mri = model_list$scale_mri, lasso = model_list$lasso, 
                lasso_predX_indx = model_list$lasso_predX_indx, 
                covX_site = sites, n_voxels_mask = model_list$n_voxels_mask, 
                combat = model_list$combat)
            linPreds = cbind(linPreds, preds)
        }
    }
    linPred = rowMeans(linPreds)
    pred = c()
    if (mp$response_family == "binomial") {
        prob = 1/(1 + exp(-linPred))
        pred = prob
        label = pred
        label[pred <= mp$bin_threshold] = mp$response_ref
        label[pred > mp$bin_threshold] = mp$response_event
        names(pred) = label
    }
    else if (mp$response_family == "gaussian") {
        pred = as.matrix(linPred)
    }
    else if (mp$response_family == "cox") {
        pred = linPred
    }
    pred
}
optimal_modulation <-
function (mask_X, mask_fu_X, mask_ijk, covX_training, trainY, 
    response_family) 
{
    nx = ceiling(sqrt(ncol(mask_X)))
    tmp_mask = mask_X[1, ]
    length(tmp_mask) <- prod(dim(matrix(mask_X[1, ], ncol = nx)))
    m_tmp_mask = matrix(!is.na(tmp_mask), ncol = nx, byrow = FALSE)
    nim = nifti(array(m_tmp_mask, dim = c(1, nx, nx)), datatype = 16)
    nim@scl_slope = 1
    writeNIfTI(nim, "temp_mask", onefile = TRUE, gzipped = TRUE, 
        verbose = FALSE, warn = -1, compression = 6)
    path_un = c()
    path_fu = c()
    for (subj_i in 1:nrow(covX_training)) {
        subj_x = mask_X[subj_i, ]
        length(subj_x) <- prod(dim(matrix(subj_x, ncol = nx)))
        m_subj_x = matrix(subj_x, ncol = nx, byrow = FALSE)
        m_subj_x[is.na(m_subj_x)] = 0
        nim = nifti(array(m_subj_x, dim = c(1, nx, nx)), datatype = 16)
        nim@scl_slope = 1
        path_un_tmp = paste(getwd(), "/R_tmp/un_", subj_i, sep = "")
        path_un = c(path_un, paste(path_un_tmp, ".nii.gz", sep = ""))
        writeNIfTI(nim, path_un_tmp, onefile = TRUE, gzipped = TRUE, 
            verbose = FALSE, warn = -1, compression = 6)
        subj_x = mask_fu_X[subj_i, ]
        length(subj_x) <- prod(dim(matrix(subj_x, ncol = nx)))
        m_subj_x = matrix(subj_x, ncol = nx, byrow = FALSE)
        m_subj_x[is.na(m_subj_x)] = 0
        nim = nifti(array(m_subj_x, dim = c(1, nx, nx)), datatype = 16)
        nim@scl_slope = 1
        path_fu_tmp = paste(getwd(), "/R_tmp/fu_", subj_i, sep = "")
        path_fu = c(path_fu, paste(path_fu_tmp, ".nii.gz", sep = ""))
        writeNIfTI(nim, path_fu_tmp, onefile = TRUE, gzipped = TRUE, 
            verbose = FALSE, warn = -1, compression = 6)
    }
    write(path_un, file = "temp_un.txt")
    write(path_fu, file = "temp_fu.txt")
    if (ncol(as.matrix(covX_training)) == 1) 
        design_mat = cbind(1, trainY)
    else design_mat = cbind(1, trainY, covX_training[, 2:ncol(covX_training)])
    write.table(design_mat, "design_mat.txt", col.names = FALSE, 
        row.names = FALSE)
    system("../matlab/optimal_modulation_linux64 -d design_mat.txt -i temp_un.txt temp_fu.txt -o temp -m temp_mask.nii.gz")
    img_kappa = readNIfTI("temp_kappa.nii.gz")
    img_kappa = as.vector(img_kappa@.Data[1, , ])
    list_mri_op_training = lapply(seq_len(ncol(mask_ijk)), function(i) {
        (1 - img_kappa[i]) * mri$data[as.numeric(mask_ijk[1, 
            i]), as.numeric(mask_ijk[2, i]), as.numeric(mask_ijk[3, 
            i]), training] + img_kappa[i] * mri_fu$data[as.numeric(mask_ijk[1, 
            i]), as.numeric(mask_ijk[2, i]), as.numeric(mask_ijk[3, 
            i]), test]
    })
    list_mask_all = pbmclapply(list_mri_op_training, function(voxel_training, 
        trainY, covX_training) {
        covm = lm.fit(covX_training, voxel_training)
        covB = coefficients(covm)
        X = residuals(covm)
        if (response_family == "binomial") {
            sig = summary(glm(trainY ~ X, family = binomial()))$coefficients[2, 
                3]
        }
        else if (response_family == "cox") {
            surv = Surv(trainY[, 1], trainY[, 2])
            sig = summary(coxph(surv ~ X))$coefficients[4]
        }
        else {
            sig = summary(lm(trainY ~ X))$coefficients[2, 3]
        }
        c(covB, X, sig)
    }, trainY, covX_training, mc.cores = n_cores)
    list_mri_op_training = ""
    mask_op_all = matrix(unlist(list_mask_all), ncol = n_voxels)
    mask_op_covB = mask_op_all[1:ncol(covX_training), ]
    mask_op_X = mask_op_all[-c(1:ncol(covX_training), nrow(mask_all)), 
        ]
    t_or_z_op_vals = mask_op_all[nrow(mask_op_all), ]
    mask_op_signif_indx = which(abs(t_or_z_op_vals[1:n_voxels]) > 
        SIGNIFICANCE_THRESHOLD)
    list(mask_op_covB = mask_op_covB, mask_op_X = mask_op_X, 
        t_or_z_op_vals = t_or_z_op_vals, mask_op_signif_indx = mask_op_signif_indx)
}
postmean <-
function (g.hat, g.bar, n, d.star, t2) 
{
    (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}
postvar <-
function (sum2, n, a, b) 
{
    (0.5 * sum2 + b)/(n/2 + a - 1)
}
remove_effects <-
function (mask_3d, mri, covX_training, trainY = NULL, n_cores = 1, 
    response_family, SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998), 
    SIGNIFICANT = TRUE, covX_site = c(), modality = "", REMOVE_EFFECT_TO_X = TRUE) 
{
    if (response_family == "gaussian") {
        SIGNIFICANCE_THRESHOLD = qt(0.97499999999999998, nrow(covX_training) - 
            ncol(covX_training))
    }
    if (length(dim(mri)) == 4) {
        mask_constant = mri[, , , 1] > -1
        mask_constant = mask_constant & (apply(mri, 1:3, var) > 
            0)
        mask_3d = mask_3d & mask_constant
        mask_ijk = cbind(which(mask_3d, arr.ind = TRUE), 1)
        mask_ijk = t(mask_ijk[which(mask_ijk[, 4] == 1), ])
        list_mri_training = lapply(seq_len(ncol(mask_ijk)), function(i) {
            mri[as.numeric(mask_ijk[1, i]), as.numeric(mask_ijk[2, 
                i]), as.numeric(mask_ijk[3, i]), ]
        })
    }
    else if (length(dim(mri) == 2)) {
        mask_constant = mri[, 1] > -1 & (apply(mri, 2, var) > 
            0)
        mask_3d = mask_3d & mask_constant
        list_mri_training = lapply(seq_len(ncol(mask_3d)), function(i) {
            mri[, i]
        })
        mask_ijk = NULL
    }
    mri <- NULL
    n_voxels = length(list_mri_training)
    message("n_cores:", n_cores)
    time.start = Sys.time()
    covX_training = as.matrix(covX_training)
    if (!is.null(trainY)) {
        mask_all = lapply(list_mri_training, function(voxel_training, 
            trainY, covX_training, SIGNIFICANT) {
            if (!is.null(covX_site) || !REMOVE_EFFECT_TO_X) {
                covB = rep(0, ncol(covX_training))
                X = voxel_training
            }
            else {
                covm = lm.fit(covX_training, voxel_training)
                covB = coefficients(covm)
                X = residuals(covm)
            }
            if (SIGNIFICANT) {
                if (length(unique(X)) > 1) {
                  if (response_family == "cox") {
                    sig = summary(coxph(trainY ~ X))$coefficients[1, 
                      4]
                    if (is.na(as.numeric(sig))) {
                      sig = 0
                    }
                  }
                  else if (response_family == "binomial") {
                    sig = summary(glm(trainY ~ X, family = binomial))$coefficients[2, 
                      3]
                  }
                  else {
                    sig = summary(lm(trainY ~ X))$coefficients[2, 
                      3]
                  }
                }
                else {
                  sig = 0
                }
            }
            else {
                sig = NA
            }
            c(covB, X, sig)
        }, switch(response_family, cox = Surv(trainY[, 1], trainY[, 
            2]), trainY), covX_training, SIGNIFICANT)
        trainY = switch(response_family, cox = Surv(trainY[, 
            1], trainY[, 2]), trainY)
    }
    list_mri_training = NULL
    mask_all = matrix(unlist(mask_all), ncol = n_voxels)
    mask_covB = mask_all[1:ncol(covX_training), ]
    X = mask_all[-c(1:ncol(covX_training), nrow(mask_all)), ]
    combat = c("")
    if (!is.null(covX_site)) {
        combat = combat_fit(dat = X, batch = covX_site, mod = covX_training, 
            verbose = F)
        X = combat_apply(tmp = combat, dat = X, batch = covX_site, 
            mod = covX_training, verbose = F)
        mask_all = apply(X$dat.combat, 2, function(voxel_training, 
            covX_training) {
            covB = rep(0, ncol(covX_training))
            X = voxel_training
            c(covB, X, NA)
        }, covX_training)
        mask_covB = mask_all[1:ncol(covX_training), ]
        X = mask_all[-c(1:ncol(covX_training), nrow(mask_all)), 
            ]
    }
    t_or_z_vals = mask_all[nrow(mask_all), ]
    mask_all = ""
    mask_signif_indx = which(abs(t_or_z_vals[1:n_voxels]) > SIGNIFICANCE_THRESHOLD)
    list(covB = mask_covB, X = X, t_or_z_vals = t_or_z_vals, 
        signif_indx = mask_signif_indx, mask_ijk = mask_ijk, 
        mask_3d = mask_3d, combat = combat)
}
remove_effects_precalculated <-
function (mask_3d, preloaded_covB, mri, covX_training = "", trainY = "", 
    response_family = "", SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998), 
    SIGNIFICANCE = FALSE, n_cores = 1, covX_site = c()) 
{
    if (response_family == "gaussian") {
        SIGNIFICANCE_THRESHOLD = qt(0.97499999999999998, nrow(covX_training) - 
            ncol(covX_training))
    }
    trainY = switch(response_family, cox = Surv(trainY[, 1], 
        trainY[, 2]), trainY)
    a = Sys.time()
    cat("Removing effects...\n")
    n_subj = length(mri[1, 1, 1, ])
    n_covs = length(preloaded_covB$data[1, 1, 1, ])
    mri_matrix = lapply(seq_len(n_subj), function(i) {
        mri[, , , i][which(mask_3d)]
    })
    mri_matrix = matrix(unlist(mri_matrix), nrow = n_subj, byrow = TRUE)
    covB_matrix = lapply(seq_len(n_covs), function(i) {
        preloaded_covB$data[, , , i][which(mask_3d)]
    })
    covB_matrix = matrix(unlist(covB_matrix), ncol = length(covB_matrix[[1]]), 
        byrow = TRUE)
    X = mri_matrix - covX_training %*% covB_matrix
    combat = c("")
    if (!is.null(covX_site)) {
        combat = combat_fit(dat = X, batch = covX_site, verbose = FALSE)
        X = combat_apply(tmp = combat, dat = X, batch = covX_site, 
            verbose = FALSE)$dat.combat
        save(combat, file = "combat.Rdata")
    }
    n_voxels = dim(mri_matrix)[2]
    message("Running on ", n_cores, " thread(s).")
    if (n_cores > 1) {
        time1 = Sys.time()
        if (SIGNIFICANCE) {
            pb <- txtProgressBar(max = n_voxels, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
            cl <- makeCluster(n_cores, setup_strategy = "sequential", 
                timeout = 0.5)
            registerDoParallel(cl)
            sig <- foreach(i = 1:n_voxels, .combine = "c", .options.snow = opts) %dopar% 
                {
                  if (length(unique(X[, i])) > 1) {
                    if (response_family == "cox") {
                      sig = summary(coxph(trainY ~ X[, i]))$coefficients[1, 
                        4]
                      if (is.na(as.numeric(sig))) 
                        sig = 0
                    }
                    else if (response_family == "binomial") {
                      sig = summary(glm(trainY ~ X[, i], family = binomial))$coefficients[2, 
                        3]
                    }
                    else {
                      sig = summary(lm(trainY ~ X[, i]))$coefficients[2, 
                        3]
                    }
                  }
                  else {
                    sig = 0
                  }
                  sig
                }
            close(pb)
            stopCluster(cl)
            cat("Time removing effects (run with ", n_cores, 
                " cores):", difftime(Sys.time(), a, units = "mins"), 
                "mins.\n")
        }
        else {
            sig <- rep(NA, n_voxels)
        }
    }
    else {
        sig = rep(NA, n_voxels)
        if (SIGNIFICANCE) {
            .require("pbapply")
            pbo = pboptions(type = "timer")
            time1 = Sys.time()
            sig = pbapply::pbsapply(1:n_voxels, function(i) {
                if (length(unique(X[, i])) > 1) {
                  if (response_family == "cox") {
                    sig = summary(coxph(trainY ~ X[, i]))$coefficients[1, 
                      4]
                    if (is.na(as.numeric(sig))) 
                      sig = 0
                  }
                  else if (response_family == "binomial") {
                    sig = summary(glm(trainY ~ X[, i], family = binomial))$coefficients[2, 
                      3]
                  }
                  else {
                    sig = summary(lm(trainY ~ X[, i]))$coefficients[2, 
                      3]
                  }
                }
                else {
                  sig = 0
                }
                sig
            })
        }
        else {
            sig <- rep(NA, n_voxels)
        }
    }
    out_X = X
    out_covB = covB_matrix
    mask_ijk = which(mask_3d, arr.ind = TRUE)
    out_mask_ijk = mask_ijk
    cat("Time removing effects:", difftime(Sys.time(), a, units = "mins"), 
        "mins.\n")
    list(X = out_X, covB = out_covB, t_or_z_vals = sig, mask_ijk = t(out_mask_ijk), 
        mask_3d = mask_3d, combat = combat)
}
rotate_coordinates <-
function (list_coordinates, iter = 1, return_3d = FALSE, img_3d = c()) 
{
    vectors = list(c(0, 1, 0), c(0, 1, 1), c(0, 0, 1), c(0, -1, 
        1), c(0, -1, 0), c(0, -1, -1), c(0, 0, -1), c(0, 1, -1), 
        c(1, 0, 0), c(1, 0, 1), c(-1, 0, 1), c(-1, 0, 0), c(-1, 
            0, -1), c(1, 0, -1), c(1, 1, 0), c(-1, 1, 0), c(-1, 
            -1, 0), c(1, -1, 0))
    vect = vectors[[iter]]
    v1 = vect/sqrt(2)
    U = list_coordinates
    C = apply(U, 1, function(u) {
        sum(u * v1)
    })
    if (runif(1) > 0.5) {
        M = C > median(C)
    }
    else {
        M = C >= median(C)
    }
    if (return_3d) {
        mask_ijk = U
        total = nrow(mask_ijk)
        img_3d = img_3d * 0
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        for (i in seq_len(nrow(mask_ijk))) {
            if (M[i]) 
                img_3d[mask_ijk[i, 1], mask_ijk[i, 2], mask_ijk[i, 
                  3]] = TRUE
            else img_3d[mask_ijk[i, 1], mask_ijk[i, 2], mask_ijk[i, 
                3]] = FALSE
            setTxtProgressBar(pb, i)
        }
        close(pb)
        img_3d
    }
    else {
        cbind(U, M)
    }
}
rotation_matrix_3d <-
function (degrees) 
{
    a = degrees * pi/180
    R.x <- matrix(c(1, 0, 0, 0, 0, cos(a), -sin(a), 0, 0, sin(a), 
        cos(a), 0, 0, 0, 0, 1), 4)
    R.y <- matrix(c(cos(a), 0, -sin(a), 0, 0, 1, 0, 0, sin(a), 
        0, cos(a), 0, 0, 0, 0, 1), 4)
    R.z <- matrix(c(cos(a), -sin(a), 0, 0, sin(a), cos(a), 0, 
        0, 0, 0, 1, 0, 0, 0, 0, 1), 4)
    list(R.x = R.x, R.y = R.y, R.z = R.z)
}
seleccio_covariables_data.frame2glmnet.matrix_fit <-
function (x, covariates) 
{
    j = match(covariates, colnames(x))
    if (any(is.na(j))) {
        stop(paste("Covariate", covariates[which(is.na(j))], 
            "not found"))
    }
    x = as.data.frame(x[, j, drop = FALSE])
    for (j in 1:ncol(x)) {
        x[, j] = factor(x[, j])
    }
    x
}
tipett_modulation <-
function (mask_X, t_or_z_vals, mask_fu_X, t_or_z_fu_vals, mask_ijk, 
    covX_training, trainY, response_family, SIGNIFICANCE = FALSE, 
    SIGNIFICANCE_THRESHOLD = qnorm(0.97499999999999998)) 
{
    indx_max_abs_un = abs(t_or_z_vals) > abs(t_or_z_fu_vals)
    t_or_z_op_vals = t_or_z_fu_vals
    t_or_z_op_vals[indx_max_abs_un] = t_or_z_vals[indx_max_abs_un]
    mask_op_X = mask_fu_X
    mask_op_X[, indx_max_abs_un] = mask_X[, indx_max_abs_un]
    list(mask_op_X = mask_op_X, t_or_z_op_vals = t_or_z_op_vals, 
        tipett_take_un = indx_max_abs_un)
}
translate_matrix_3d <-
function (x, y, z) 
{
    matrix(c(1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1), 
        4)
}
