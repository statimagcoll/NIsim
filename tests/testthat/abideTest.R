
abideDir = file.path(tempdir(), 'abide')
derivatives = c("alff", "degree_binarize", "degree_weighted", "dual_regression", "eigenvector_binarize", "eigenvector_weighted", "falff", "func_mask", "
                func_mean", "func_preproc", "lfcd", "reho", "rois_aal", "rois_cc200", "rois_cc400", "rois_dosenbach160", "rois_ez", "rois_ho", "rois_tt", "vmhc")
pipelines = c("ccs", "cpac", "dparsf", "niak")
strategies = c("filt_global", "filt_noglobal", "nofilt_global", "nofilt_noglobal")

