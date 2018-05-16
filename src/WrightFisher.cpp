#include "WrightFisher.hpp"

double WrightFisher::psi_diploid(const llong i, const llong N, const double s, const double h, const double u, const double v) {

    llong j = (2 * N) - i;
    double w_11 = 1 + s;
    double w_12 = 1 + (s * h);
    double w_22 = 1;
    double a = w_11 * i * i;
    double b = w_12 * i * j;
    double c = w_22 * j * j;
    double w_bar = a + (2 * b) + c;
    return (((a + b) * (1 - u)) + ((b + c) * v)) / w_bar;
}

WrightFisher::Row WrightFisher::binom_row(const llong i, const llong Nx, const llong Ny, const double s, const double h, const double u, const double v, const double alpha) {

    llong Ny2 = 2 * Ny;

    // binomial sampling probability
    double p = psi_diploid(i, Nx, s, h, u, v);

    // start and end quantiles for covering 1 - alpha weight of the probability mass
    llong start = (llong)binom_tail_cover(alpha / 2, Ny2, p, true);
    llong end = (llong)binom_tail_cover(alpha / 2, Ny2, p, false);

    // make sure we didn't mess up
    assert(start < end);

    // Initialize row
    Row r(start, end);

    // Start iterative binomial calculation (WFES supplementary, eq 18,19)
    double d = ld_binom(start, Ny2, p);
    double lc = log(p) - log(1 - p);
    r.Q(0) = d;

    // Iterative binomial (in log)
    for(llong j = start + 1; j <= end; j++) {
        d += log(Ny2 - j + 1) - log(j) + lc;
        r.Q(j - start) = d;
    }

    // Exponentiate
    r.Q.exp_inplace();
    // Re-weigh to sum to 1
    r.weight = r.Q.normalize_inplace();

    return r;

}

WrightFisher::Matrix WrightFisher::Equilibrium(const llong N, const double s, const double h, const double u, const double v, const double alpha, const llong block_size) {
    llong N2 = 2 * N;
    WrightFisher::Matrix W(N2 + 1, N2 + 1, n_absorbing(WrightFisher::NON_ABSORBING));
    for(llong i = 0; i <= N2; i++) {
        WrightFisher::Row r = binom_row(i, N, N, s, h, u, v, alpha);
        // I - Q
        r.Q.negate_inplace();
        r.Q(i - r.start) += 1;
        if (r.end == N2) {
            r.Q(r.size - 1) += 1;
            W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1);
        } else {
            // allocate one additional cell (slack = 1)
            W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1, false, 1);
            W.Q.values[W.Q.non_zeros - 1] = 1;
            W.Q.columns[W.Q.non_zeros - 1] = N2;
            W.Q.finalize_row();
        }
    }
    return W;
}

WrightFisher::Matrix WrightFisher::Single(const llong Nx, const llong Ny, const absorption_type abs_t, const double s, const double h, const double u, const double v, const double alpha, const llong block_size) {

    llong Nx2 = 2 * Nx; 
    llong Ny2 = 2 * Ny; 

    llong n_abs = n_absorbing(abs_t);

    Matrix W(Nx2 + 1 - n_abs, Ny2 + 1 - n_abs, n_abs);

    if(abs_t == NON_ABSORBING) {
        for(llong i = 0; i <= Nx2; i++) {
            Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);
            W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1);
        }
    }

    else if (abs_t == EXTINCTION_ONLY) {
        for(llong i = 1; i <= Nx2; i++) {
            Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);
            if (r.start == 0) {
                //                   m0       m1         r0 r1
                W.Q.append_data(r.Q, r.start, r.end - 1, 1, r.size - 1);
                W.R(i - 1, 0) = r.Q(0);
            } else {
                W.Q.append_data(r.Q, r.start - 1, r.end - 1, 0, r.size - 1);
            }
            
        }
    }

    else if (abs_t == FIXATION_ONLY) {
        for(llong i = 0; i <= Nx2 - 1; i++) {
            Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);
            if (r.end == Ny2) {
                //                   m0       m1         r0 r1
                W.Q.append_data(r.Q, r.start, r.end - 1, 0, r.size - 2);
                W.R(i, 0) = r.Q(r.size - 1);
            } else {
                W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1);
            }
            
        }
    }

    else if (abs_t == BOTH_ABSORBING) {
        for(llong i = 1; i <= Nx2 - 1; i++) {
            Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);

            if (r.start == 0 && r.end == Ny2) {
                W.Q.append_data(r.Q, r.start, r.end - 2, 1, r.end - 1);
                W.R(i - 1, 0) = r.Q(0);
                W.R(i - 1, 1) = r.Q(r.size - 1);
            } else if (r.start == 0) {
                W.Q.append_data(r.Q, r.start, r.end - 1, 1, r.size - 1);
                W.R(i - 1, 0) = r.Q(0);
            } else if (r.end == Ny2) {
                W.Q.append_data(r.Q, r.start - 1, r.end - 2, 0, r.size - 2);
                W.R(i - 1, 0) = r.Q(r.size - 1);
            } else {
                W.Q.append_data(r.Q, r.start - 1, r.end - 1, 0, r.size - 1);
            }
        }
    }

    else {
        throw std::runtime_error("WrightFisher::Single() unknown absorption_type");
    }

    return W;


    // llong size = (2 * N) + 1 - n_abs;

    // Matrix W(size, n_abs);

    // for(llong r = 0; r < size; r += block_size) {

    //     llong block = (r + block_size) < size ? block_size : size - r;
    //     vector<Row> buffer(block);

    //     #pragma omp parallel for
    //     for(llong b = 0; b < block; b++) {
    //         buffer[b] = make_row(r + b, N, N, n_abs, s, h, u, v, alpha);
    //     }

    //     for(llong b = 0; b < block; b++) {
    //         Row& row = buffer[b];
    //         W.Q.append_data(row.Q, row.col_start);
    //         for(llong j = 0; j < n_abs; j++) {
    //             W.R(r + b, j) = row.R(j);
    //         }
    //     }
    // }

}

WrightFisher::Matrix WrightFisher::SingleAlt(const llong Nx, const llong Ny, const absorption_type abs_t, const double s, const double h, const double u, const double v, const double alpha, const llong block_size) 
{
    llong Nx2 = 2 * Nx; 
    llong Ny2 = 2 * Ny; 

    llong n_abs = n_absorbing(abs_t);

    Matrix W(Nx2 + 1 - n_abs, Ny2 + 1 - n_abs, n_abs);

    for(llong row = 0; row <= Nx2; row++) {
        Row r = binom_row(row, Nx, Ny, s, h, u, v, alpha);
        llong r_last = r.size - 1;

        switch(abs_t) {
            case NON_ABSORBING:
                W.Q.append_data(r.Q, r.start, r.end, 0, r_last);
            break;

            case EXTINCTION_ONLY:
                if (row == 0) continue;
                else {
                    if (r.start == 0) {
                        W.Q.append_data(r.Q, r.start, r.end - 1, 1, r_last);
                        W.R(row - 1, 0) = r.Q(0);
                    } else {
                        W.Q.append_data(r.Q, r.start - 1, r.end - 1, 0, r_last);
                    }
                }
            break;

            case FIXATION_ONLY:
                if (row == Nx2) continue;
                else {
                    if (r.end == Ny2) {
                        W.Q.append_data(r.Q, r.start, r.end - 1, 0, r_last - 1);
                        W.R(row, 0) = r.Q(r_last);
                    } else {
                        W.Q.append_data(r.Q, r.start, r.end, 0, r_last);
                    }
                }
            break;

            case BOTH_ABSORBING:
                if (row == 0 || row == Nx2) continue;
                else {
                    if (r.start == 0 && r.end == Ny2) {
                        W.Q.append_data(r.Q, r.start, r.end - 2, 1, r_last - 1);
                        W.R(row - 1, 0) = r.Q(0);
                        W.R(row - 1, 1) = r.Q(r_last);
                    } else if (r.start == 0) {
                        W.Q.append_data(r.Q, r.start, r.end - 1, 1, r_last);
                        W.R(row - 1, 0) = r.Q(0);
                    } else if (r.end == Ny2) {
                        W.Q.append_data(r.Q, r.start - 1, r.end - 2, 0, r_last - 1);
                        W.R(row - 1, 1) = r.Q(r_last);
                    } else {
                        W.Q.append_data(r.Q, r.start - 1, r.end - 1, 0, r_last);
                    }
                }
            break;
        }
    }
    return W;
}

std::vector<std::pair<llong, llong>> submatrix_indeces(const lvec& sizes) {
    llong i = 0;
    llong j = 0;

    llong size = sizes.sum();

    std::vector<std::pair<llong, llong>>idx (size);

    for(llong r = 0; r < size; r++) {
        if (j == sizes(i)) {
            j = 0;
            i++;
        }
        idx[r].first = i;
        idx[r].second = j;
        j++;
    }
    return idx;
}

WrightFisher::Matrix WrightFisher::Switching(const lvec& N, const absorption_type abs_t, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, double alpha, llong block_size) {
    llong k = N.size;

    if (abs_t == EXTINCTION_ONLY) {
        // backward mutation rate should be above 0
        for(llong i = 0; i < k; i++) {
            assert(u(i) > 0);
        }
    } else if (abs_t == FIXATION_ONLY) {
        // forward mutation rate should be above 0
        for(llong i = 0; i < k; i++) {
            assert(v(i) > 0);
        }
    }

    llong n_abs_total = n_absorbing(abs_t) * k;

    lvec sizes(k);
    for(llong i = 0; i < k; i++) sizes(i) = 2 * N(i) + 1;
    llong size = sizes.sum();

    Matrix W(size - n_abs_total, size - n_abs_total, n_abs_total);
    std::vector<std::pair<llong, llong>> index = submatrix_indeces(sizes);

    for(llong row = 0; row < size; row++) {
        llong i = index[row].first; // model index
        llong im = index[row].second; // current index within model i

        // coordinate of the submodel start
        llong offset = 0;
        // iterate over submodels
        for(llong j = 0; j < k; j++) {
            Row r = binom_row(im, N(i), N(j), s(i), h(i), u(i), v(i), alpha);
            r.Q *= switching(i, j);

            bool row_complete = (j == (k - 1));
            llong m_start = r.start + offset;
            llong m_end = r.end + offset;
            llong r_last = r.size - 1;

            // W.Q.append_data(r.Q, r.start + offset, r.end + offset, 0, r_last, j == (k - 1));
            switch(abs_t) {
                case NON_ABSORBING:
                    W.Q.append_data(r.Q, m_start, m_end, 0, r_last, row_complete);
                break;

                case EXTINCTION_ONLY:
                    if (im == 0) continue;
                    else {
                        if (r.start == 0) {
                            W.Q.append_data(r.Q, m_start, m_end - 1, 1, r_last, row_complete);
                            W.R(row - 1, j) = r.Q(0);
                        } else {
                            W.Q.append_data(r.Q, m_start - 1, m_end - 1, 0, r_last, row_complete);
                        }
                    }
                break;

                case FIXATION_ONLY:
                    if (im == N(i) * 2) continue;
                    else {
                        if (r.end == N(j) * 2) {
                            W.Q.append_data(r.Q, m_start, m_end - 1, 0, r_last - 1, row_complete);
                            W.R(row, j) = r.Q(r_last);
                        } else {
                            W.Q.append_data(r.Q, m_start, m_end, 0, r_last, row_complete);
                        }
                    }
                break;

                case BOTH_ABSORBING:
                    if (im == 0 || im == N(i) * 2) continue;
                    else {
                        if (r.start == 0 && r.end == N(j) * 2) {
                            W.Q.append_data(r.Q, m_start, m_end - 2, 1, r_last - 1, row_complete);
                            W.R(row - 1, 2 * j) = r.Q(0);
                            W.R(row - 1, (2 * j) + 1) = r.Q(r_last);
                        } else if (r.start == 0) {
                            W.Q.append_data(r.Q, m_start, m_end - 1, 1, r_last, row_complete);
                            W.R(row - 1, 2 * j) = r.Q(0);
                        } else if (r.end == N(j) * 2) {
                            W.Q.append_data(r.Q, m_start - 1, m_end - 2, 0, r_last - 1, row_complete);
                            W.R(row - 1, (2 * j) + 1) = r.Q(r_last);
                        } else {
                            W.Q.append_data(r.Q, m_start - 1, m_end - 1, 0, r_last, row_complete);
                        }
                    }
                break;
            }

            offset += (sizes(j) - n_absorbing(abs_t));
        }
    }

    return W;

}
