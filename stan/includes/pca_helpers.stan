/*
Functions to do Householder reparametrisation
Taken from https://github.com/RSNirwan/HouseholderBPCA/blob/master/py_stan_code/ppca_house_improved.stan
*/
    matrix V_low_tri_plus_diag_unnormed (int D, int Q, vector v) {
        // Put parameters into lower triangular matrix
        matrix[D, Q] V;

        int idx = 1;
        for (d in 1:D) {
            for (q in 1:Q) {
                if (d >= q) {
                    V[d, q] = v[idx];
                    idx += 1;
                } else
                V[d, q] = 0;
            }
        }
        return V;
    }
    matrix V_low_tri_plus_diag (int D, int Q, vector v) {
        matrix[D, Q] V = V_low_tri_plus_diag_unnormed(D, Q, v);
        for (q in 1:Q){
            V[,q] = V[,q]/sqrt( sum(square(V[,q])) );
        }
        return V;
    }
    real sign(real x){
        if (x < 0.0)
            return -1.0;
        else
            return 1.0;
    }
    matrix Householder (int k, matrix V) {
        // Householder transformation corresponding to kth column of V
        int D = rows(V);
        vector[D] v = V[, k];
        matrix[D,D] H;
        real sgn = sign(v[k]);
        
        //v[k] +=  sgn; //v[k]/fabs(v[k]);
        v[k] += v[k]/fabs(v[k]);
        H = diag_matrix(rep_vector(1, D)) - (2.0 / dot_self(v)) * (v * v');
        H[k:, k:] = -sgn*H[k:, k:];
        return H;
    }
    matrix[] H_prod_right (matrix V) {
        // Compute products of Householder transformations from the right, i.e. backwards
        int D = rows(V);
        int Q = cols(V);
        matrix[D, D] H_prod[Q + 1];
        H_prod[1] = diag_matrix(rep_vector(1, D));
        for (q in 1:Q)
            H_prod[q + 1] = Householder(Q - q + 1, V) * H_prod[q];
        return H_prod;    
    }
    matrix orthogonal_matrix (int D, int Q, vector v) {
        matrix[D, Q] V = V_low_tri_plus_diag(D, Q, v);
        // Apply Householder transformations from right
        matrix[D, D] H_prod[Q + 1] = H_prod_right(V);
        return H_prod[Q + 1][, 1:Q];    
    }
