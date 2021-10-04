{
      int nA_posneg = nA_pos + nA_neg;
      int A_[nA - nA_posneg];
      int A_posneg[nA_posneg] = append_array(A_pos, A_neg);
      int A_posneg_asc[nA_posneg] = sort_asc(A_posneg);
      int J = 1;
      int k = 1;
      
      for (i in 1:nA) {
        int match = 0;
        if (J <= (nA_posneg))
        for (j in J:(nA_posneg)){
          if (A_posneg_asc[j] == i) {
            match = 1;
            J += 1;
            break;
          }
        }
        if (match == 0) {
          A_[k] = i;
          k += 1;
        }
      }
      a_raw[A_posneg] = append_row(a_pos, a_neg);
      a_raw[A_] = a_;
      a = a_raw*SP[1];
    }   
