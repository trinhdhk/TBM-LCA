// * Clinical symptoms Td[1:3] ----------------------------------------------
int Td_cs[N,3] = Td[:,1:3];
int obs_cs[N,3] = obs_Td[:,1:3];
int<lower=0> N_miss_cs;
int<lower=1,upper=N> n_miss_cs[(N*3) - sum2d(obs_cs)];
int<lower=1,upper=3> d_miss_cs[size(n_miss_cs)];
int<lower=0> N_pos_cs;
int<lower=1,upper=N> n_pos_cs[sum2d_with_missing(Td_cs, obs_cs)];
int<lower=1,upper=3> d_pos_cs[size(n_pos_cs)];
int<lower=0> N_neg_cs;
int<lower=1,upper=N> n_neg_cs[sum2d(obs_cs) - size(n_pos_cs)];
int<lower=1,upper=3> d_neg_cs[size(n_neg_cs)];

// * Motor palsy Td[4:6] ---------------------------------------------------
int Td_mp[N,3] = Td[,4:6];
int obs_mp[N,3] = obs_Td[,4:6];
int<lower=0> N_miss_mp;
int<lower=1,upper=N> n_miss_mp[(N*3) - sum2d(obs_mp)];
int<lower=1,upper=3> d_miss_mp[size(n_miss_mp)];
int<lower=0> N_pos_mp;
int<lower=1,upper=N> n_pos_mp[sum2d_with_missing(Td_mp, obs_mp)];
int<lower=1,upper=3> d_pos_mp[size(n_pos_mp)];
int<lower=0> N_neg_mp;
int<lower=1,upper=N> n_neg_mp[sum2d(obs_mp) - size(n_pos_mp)];
int<lower=1,upper=3> d_neg_mp[size(n_neg_mp)];

// * CSF laboratory test Xc[3:8] --------------------------------------------
int obs_csf = sum2d(obs_Xc[:,2:7]);

// * GCS
int obs_gcs_compartments = sum2d(obs_Tc[:,1:3]);
int obs_gcs = sum(obs_Xc[:,8]);
// vector<lower=0, upper=1>[3] gcs_a0 = [0, 0, 0]';

// * Blood test
int obs_bld = sum2d(obs_Tc[:,4:5]);
vector[2] bld_a0 = [0, 0]';

// * Var assignments --------------------------------------------------------
N_pos_cs  = size(n_pos_cs);
N_neg_cs  = size(n_neg_cs);
N_miss_cs = size(n_miss_cs);
{
  int i = 1;
  int j = 1;
  int k = 1;
  for (n in 1:N) {
    for (d in 1:3) {
      if (!obs_cs[n,d]){
        n_miss_cs[k] = n;
        d_miss_cs[k] = d;
        k += 1;
      } else if (Td_cs[n,d]) {
        n_pos_cs[i] = n;
        d_pos_cs[i] = d;
        i += 1;
      } else {
        n_neg_cs[j] = n;
        d_neg_cs[j] = d;
        j += 1;
      }
    }
  }
}

N_pos_mp = size(n_pos_mp);
N_neg_mp = size(n_neg_mp);
N_miss_mp = size(n_miss_mp);
{
  int i = 1;
  int j = 1;
  int k = 1;
  for (n in 1:N) {
    for (d in 1:3) {
      if (!obs_mp[n,d]){
        n_miss_mp[k] = n;
        d_miss_mp[k] = d;
        k += 1;
      } else if (Td_mp[n,d]) {
        n_pos_mp[i] = n;
        d_pos_mp[i] = d;
        i += 1;
      } else {
        n_neg_mp[j] = n;
        d_neg_mp[j] = d;
        j += 1;
      }
    }
  }
}

for (n in 1:N_all){
  obs_test_all[n] = obs_Smear_all[n] + obs_Mgit_all[n] + obs_Xpert_all[n] > 0;
  test_all[n] = Y_Smear_all[n] + Y_Mgit_all[n] + Y_Xpert_all[n] > 0;
}
obs_test = obs_test_all[which(keptin)];
test = test_all[which(keptin)];
