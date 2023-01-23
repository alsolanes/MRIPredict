apply_model = function(mp, mri_data, mri_fu_data, mri_wm_data=NULL, mri_wm_fu_data=NULL, 
                       covX_test = NULL, signif_indx, lasso_covB, lasso_covB_fu = NULL, mask, 
                       predX_test, scale_clinical=NULL, scale_mri=NULL,lasso, lasso_predX_indx, 
                       tipett_take_un=NULL,img_kappa=NULL,use_significant_voxels=FALSE, covX_site=c(),
                       n_voxels_mask, combat=NULL, standardize_images = FALSE, masks_3d
) {
  
  X = NULL
  X_signif = c()
  if (mp$modulation=='op')
    TIPETT=TRUE
  else
    TIPETT=FALSE
  if (standardize_images){
    mri_data = mri_data - scale_mri$mean_mri
    if(!is.null(mri_fu_data)){
      mri_fu_data = (mri_fu_data - scale_mri$mean_mri_fu)/scale_mri$sd_mri_fu * scale_mri$sd_mri
    }
    if(!is.null(mri_wm_data)){
      mri_wm_data = (mri_wm_data - scale_mri$mean_mri_wm)/scale_mri$sd_mri_wm * scale_mri$sd_mri
    }
    if(!is.null(mri_wm_fu_data)){
      mri_wm_fu_data = (mri_wm_fu_data - scale_mri$mean_mri_wm_fu)/scale_mri$sd_mri_wm_fu * scale_mri$sd_mri
    }
  }

  # remove effect of covariates 
  # passar mri_test i mri_fu_test a 2D
  # 
  n_image_modalities = 0
  if (!is.null(n_voxels_mask) && any(lasso$i <= n_voxels_mask)) { # if the model contains at least one voxel
    if (!is.null(mri_data) ){
      n_image_modalities = n_image_modalities + 1
      n_subj = nrow(covX_test)
      mri_test = lapply(seq_len(n_subj),function(i){mri_data[,,,i][masks_3d$un_gm]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      # Apply combat
      if(!is.null(covX_site)) {
        mri_test = combat_apply(combat$gm, mri_test, covX_site, mod=covX_test, verbose = FALSE)$dat.combat 
      }
      # FI combat
      img_3d = mri_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_data)[4]){ #for every image
        img_3d[masks_3d$un_gm] = mri_test[img_i,]
        mri_data[,,,img_i] = img_3d
      }
    }
    if (!is.null(mri_fu_data)) {
      n_image_modalities = n_image_modalities + 1
      mri_test = lapply(seq_len(n_subj),function(i){mri_fu_data[,,,i][masks_3d$fu_gm]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      if(!is.null(covX_site)) {
        mri_test = combat_apply(combat$gm_mod, mri_test, covX_site, mod=covX_test, verbose = FALSE)$dat.combat 
      }
      img_3d = mri_fu_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_fu_data)[4]){ #for every image
        img_3d[masks_3d$fu_gm] = mri_test[img_i,]
        mri_fu_data[,,,img_i] = img_3d
      }
    }
    if (!is.null(mri_wm_data)) {
      n_image_modalities = n_image_modalities + 1
      mri_test = lapply(seq_len(n_subj),function(i){mri_wm_data[,,,i][masks_3d$un_wm]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      if(!is.null(covX_site)) {
        mri_test = combat_apply(combat$wm, mri_test, covX_site, mod=covX_test, verbose = FALSE)$dat.combat 
      }
      img_3d = mri_wm_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_wm_data)[4]){ #for every image
        img_3d[masks_3d$un_wm] = mri_test[img_i,]
        mri_wm_data[,,,img_i] = img_3d
      }
    }
    if (!is.null(mri_wm_fu_data)) {
      n_image_modalities = n_image_modalities + 1
      mri_test = lapply(seq_len(n_subj),function(i){mri_wm_fu_data[,,,i][masks_3d$fu_wm]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      if(!is.null(covX_site)) {
        mri_test = combat_apply(combat$wm_mod, mri_test, covX_site,mod=covX_test, verbose = FALSE)$dat.combat 
      }
      img_3d = mri_wm_fu_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_wm_fu_data)[4]){ #for every image
        img_3d[masks_3d$fu_wm] = mri_test[img_i,]
        mri_wm_fu_data[,,,img_i] = img_3d
      }
    }

    if (!is.null(mri_data)) {
      n_subj = dim(mri_data)[4]
      
      #################################################################################################################
      if ( mp$modulation != 'un' && mp$modulation!='fu' ){
        
        orig_dims = dim(mri_data)
        new_mri = array(NA, dim = dim(mri_data) * c(n_image_modalities,1,1,1))
        for (img_i in 1:dim(mri_data)[4]){ #for every image
          if(n_image_modalities>0){
            img_j = 1
            first_index_x = (((img_j-1)*(orig_dims[1]))+1)
            new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_data[,,,img_i]
          }
          if(n_image_modalities>1){
            img_j = 2
            first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
            new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_fu_data[,,,img_i]
          }
          if(n_image_modalities>2){
            img_j = 3
            first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
            new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_wm_data[,,,img_i]
          }
          if(n_image_modalities>3){
            img_j = 4
            first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
            new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_wm_fu_data[,,,img_i]
          }
        }
        mri_data = new_mri
      }
      #################################################################################################################

      mri_test = lapply(seq_len(n_subj),function(i){mri_data[,,,i][which(mask)]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      
      n_voxels = sum(mask)
      if(use_significant_voxels) {
        mri_test = mri_test[,signif_indx] # signif_indx respect mask
        n_voxels = sum(signif_indx)
      }
      mri_test = mri_test[, lasso$i[which(lasso$i <= n_voxels)]] # lasso$i respect signif_indx
      
      if (mp$modulation %in% c('fu','op')) {
        mri_fu_test = lapply(seq_len(n_subj),function(i){mri_fu_data[,,,i][which(mask)]})
        mri_fu_test = matrix(unlist(mri_fu_test), nrow = n_subj, byrow = TRUE)
        if(use_significant_voxels) {
          mri_fu_test = mri_fu_test[,signif_indx] # afegit 
        }
        mri_fu_test = mri_fu_test[, lasso$i[which(lasso$i <= n_voxels)]]
      }
      # 2. remove effect sex and age
      if(!is.null(covX_test)) {
        X_un = mri_test - covX_test %*% lasso_covB # X_un will be lasso$i indices over significant voxels
      } else {
        X_un = mri_test
      }
      
      if (!is.null(lasso_covB_fu)) {

        if (!is.null(covX_test)) {
          X_fu = mri_fu_test - covX_test %*% lasso_covB_fu
        } else {
          X_fu = mri_fu_test
        }
      } else {
        X_fu = NULL
      }

      ########################################################
      
      if (TIPETT) {
        if(use_significant_voxels) {
          tipett_take_un = tipett_take_un[which(signif_indx)]
        }
        take_X_un_val = tipett_take_un[lasso$i] # select the significant tippet tvals of the mask
        X_op = matrix(0, nrow = nrow(X_un), ncol = ncol(X_un))
        X_op[,  take_X_un_val] = X_un[,take_X_un_val] # which tipett abs(tval_un)>abs(tval_fu)
        X_op[, !take_X_un_val] = X_fu[, !take_X_un_val]
        
      } else if (mp$modulation=='op') {
        for (i in 1:ncol(lasso_covB)) {
          
          X_op[,i] = (1 - img_kappa[i]) * X_un[,i] + img_kappa[i] * X_fu[,i]
        }
      }
      
      
      X_signif = switch(mp$modulation,
                        un=X_un,
                        fu=X_fu,
                        op=X_op,
                        all=X_un)
    }
    ### TEST: PREDICTOR VARIABLES ###
  }
  if(!is.null(mp$pred_transf)) {
    predX_test = as.matrix(predX_test)
    # Scale clinical variables in the test dataset

    if(!is.null(scale_clinical) & !any(is.na(scale_clinical$MEANS))){
      predX_test = t(apply(predX_test, 1, function(x) {x - scale_clinical$MEANS})) # Center the covariates
      predX_test = t(apply(predX_test, 1, function(x) {x / scale_clinical$SDS * scale_clinical$MODE_SD})) # Scale the covariates to MODE_SD
    }
    predX_test = predX_test[,lasso_predX_indx]
    X = cbind(X_signif, predX_test)
  } else {
    X = X_signif
  }
  ### TEST: END PREDICTOR ###
  linPred = as.matrix(X) %*% lasso$beta
  if (mp$response_family != 'cox') {
    linPred = lasso$a0 + linPred
  }

  linPred
}
