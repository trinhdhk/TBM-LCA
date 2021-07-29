 // * PCA for the CSF --------------------------------------------------------
  for (i in 1:5) {
    vector[sum(obs_Xc_all[:,i+2])] X_unscaled = to_vector(Xc_all[which(obs_Xc_all[:,i+2]),i+2]);
    CSF_scaled[which(obs_Xc_all[:,i+2]), i] = to_array_1d((X_unscaled) ./ sd(X_unscaled));
    CSF_scaled[which_not(obs_Xc_all[:,i+2]), i] = rep_array(0.0, N_all - sum(obs_Xc_all[:,i+2]));
  }
