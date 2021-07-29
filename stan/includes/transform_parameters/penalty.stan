    real SP[2];
    SP[1] = (adapt_penalty[1] == 1) ? sp[1] : penalty_term[1];
    SP[2] = (penalty_term[2] == 0) ? sp[num_elements(sp)] : (penalty_term[2] == -1) ? SP[1] : penalty_term[2];
