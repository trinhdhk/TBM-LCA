int nXd2 = nXd - 2 + nTd;
int nXc2 = nXc - 1 + nTc; /*GCS*/
int obs_Xd2_all[N_all, nXd2] ;
int obs_Xd2[N, nXd2];
int obs_Xc2_all[N_all, nXc2];
int obs_Xc2[N, nXc2];

{
// Build the new obs_Xd matrix
obs_Xd2_all[:,1] = obs_Xd_all[:,1];
obs_Xd2_all[:,2:4] = obs_Td_all[:,1:3];
obs_Xd2_all[:,5:7] = obs_Td_all[:,4:6];
obs_Xd2_all[:,8:nXd2 - 1] = obs_Xd_all[:,4:nXd];
obs_Xd2_all[:,nXd2] = obs_Td_all[:,nTd];
obs_Xd2 = obs_Xd2_all[which(keptin),:];

//Build the new obs_Xc matrix
obs_Xc2_all[:,1:(nXc - 1)] = obs_Xc_all[:,1:(nXc - 1)];
obs_Xc2_all[:,nXc:nXc2]    = obs_Tc_all[:,1:(nXc2 - nXc)];
obs_Xc2 = obs_Xc2_all[which(keptin),:];
}
