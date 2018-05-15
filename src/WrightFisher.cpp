#include "WrightFisher.hpp"

double WrightFisher::psi_diploid(llong i, llong N, double s, double h, double u, double v) {

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

WrightFisher::Row WrightFisher::binom_row(llong i, llong Nx, llong Ny, double s, double h, double u, double v, double alpha) {

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
    r.Q.exp();
    // Re-weigh to sum to 1
    r.weight = r.Q.normalize();

    return r;

}

WrightFisher::Matrix WrightFisher::Equilibrium(llong N, double s, double h, double u, double v, double alpha, llong block_size) {
    llong N2 = 2 * N;
    WrightFisher::Matrix W(N2 + 1, N2 + 1, 0, WrightFisher::NEITHER);
    for(llong i = 0; i <= N2; i++) {
        WrightFisher::Row r = binom_row(i, N, N, s, h, u, v, alpha);
        // I - Q
        r.Q.negate();
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

WrightFisher::Matrix WrightFisher::Single(llong Nx, llong Ny, absorption_type abs_t, double s, double h, double u, double v, double alpha, llong block_size) {

   llong Nx2 = 2 * Nx; 
   llong Ny2 = 2 * Ny; 

    if(abs_t == NEITHER) {

        WrightFisher::Matrix W(Nx2 + 1, Ny2 + 1, 0, abs_t);
        for(llong i = 0; i <= Nx2; i++) {
            WrightFisher::Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);
            W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1);
        }
        return W;
    }

    else if (abs_t == EXTINCTION) {

        WrightFisher::Matrix W(Nx2, Ny2, 1, abs_t);
        for(llong i = 1; i <= Nx2; i++) {
            WrightFisher::Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);
            if (r.start == 0) {
                //                   m0       m1         r0 r1
                W.Q.append_data(r.Q, r.start, r.end - 1, 1, r.size - 1);
                W.R(i - 1, 0) = r.Q(0);
            } else {
                W.Q.append_data(r.Q, r.start - 1, r.end - 1, 0, r.size - 1);
            }
            
        }
        return W;
    }

    else if (abs_t == FIXATION) {

        WrightFisher::Matrix W(Nx2, Ny2, 1, abs_t);
        for(llong i = 0; i <= Nx2 - 1; i++) {
            WrightFisher::Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);
            if (r.end == Ny2) {
                //                   m0       m1         r0 r1
                W.Q.append_data(r.Q, r.start, r.end - 1, 0, r.size - 2);
                W.R(i, 0) = r.Q(r.size - 1);
            } else {
                W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1);
            }
            
        }
        return W;
    }

    else if (abs_t == BOTH) {

        WrightFisher::Matrix W(Nx2 - 1, Ny2 - 1, 2, abs_t);
        for(llong i = 1; i <= Nx2 - 1; i++) {
            WrightFisher::Row r = binom_row(i, Nx, Ny, s, h, u, v, alpha);

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
        return W;
    }

    else {
        throw runtime_error("WrightFisher::Single() unknown absorption_type");
    }


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
