// vector[3] z_cs[N]; // - Clinical symptoms
// vector[3] z_mp[N]; // - Motor palsy
// 
// for (n in 1:N_miss_cs) z_cs[n_miss_cs[n], d_miss_cs[n]] = z_miss_cs[n];
// for (n in 1:N_pos_cs)  z_cs[n_pos_cs[n] , d_pos_cs[n]]  = z_pos_cs[n];
// for (n in 1:N_neg_cs)  z_cs[n_neg_cs[n] , d_neg_cs[n]]  = z_neg_cs[n];
// 
// for (n in 1:N_miss_mp) z_mp[n_miss_mp[n], d_miss_mp[n]] = z_miss_mp[n];
// for (n in 1:N_pos_mp)  z_mp[n_pos_mp[n] , d_pos_mp[n]]  = z_pos_mp[n];
// for (n in 1:N_neg_mp)  z_mp[n_neg_mp[n] , d_neg_mp[n]]  = z_neg_mp[n];

// vector[3] cs_p = append_row(cs_p12, (1 - cs_p12[1])*(1 - cs_p12[2]));
// vector[3] mp_p = append_row(mp_p12, (1 - mp_p12[1])*(1 - mp_p12[2]));
// vector[3] gcs_a0 = [0, 0, 0]';
