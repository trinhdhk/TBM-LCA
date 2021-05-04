int num_basis = num_knots + spline_degree - 1; // total number of B-splines
matrix[num_basis, num_data] B;  // matrix of B-splines
vector[spline_degree + num_knots] ext_knots_temp;
vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
for (ind in 1:num_basis)
  B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
B[num_knots + spline_degree - 1, num_data] = 1;