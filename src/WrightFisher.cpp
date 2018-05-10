#include "WrightFisher.hpp"

double WrightFisher::psi_diploid(llong i, llong N, double s, double h, double u, double v)
{ //{{{
    llong j = (2 * N) - i;
    double w_11 = 1 + s;
    double w_12 = 1 + (s * h);
    double w_22 = 1;
    double a = w_11 * i * i;
    double b = w_12 * i * j;
    double c = w_22 * j * j;
    double w_bar = a + (2 * b) + c;
    return (((a + b) * (1 - u)) + ((b + c) * v)) / w_bar;
} //}}}

WrightFisher::Row WrightFisher::binom_row(llong i, llong Nx, llong Ny, double s, double h, double u, double v, double alpha)
{ //{{{
    llong Ny2 = 2 * Ny;
    double a = alpha / 2.0;

    // binomial sampling probability
    double p = psi_diploid(i, Nx, s, h, u, v);

    // start and end quantiles for covering 1 - alpha weight of the probability mass
    llong start = (llong)binom_tail_cover(a, Ny2, p, true);
    llong end = (llong)binom_tail_cover(a, Ny2, p, false);

    // make sure we didn't mess up
    assert(start < end);

    // Create row container
    llong row_size = end - start + 1;
    Row r(start, end, row_size);

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
} //}}}

WrightFisher::Matrix WrightFisher::Single(llong N, absorption_type mt, double s, double h, double u, double v, double alpha, llong block_size)
{ //{{{

    if(mt == NEITHER) {
        WrightFisher::Matrix W((2 * N) + 1, 0, mt);
        for(llong i = 0; i <= 2 * N; i++) {
            WrightFisher::Row r = binom_row(i, N, N, s, h, u, v, alpha);
            W.Q.append_data(r.Q, r.start, r.end, r.start, r.end);
        }
        return W;
    }
    else if (mt == EXTINCTION) {

        WrightFisher::Matrix W(2 * N, 1, mt);
        for(llong i = 1; i <= 2 * N; i++) {
            WrightFisher::Row r = binom_row(i, N, N, s, h, u, v, alpha);
            if (r.start == 0) {
                //                   m0       m1         r0 r1
                W.Q.append_data(r.Q, r.start, r.end - 1, 1, r.size - 1);
                W.R(i - 1, 0) = r.Q(0);
            } else {
                W.Q.append_data(r.Q, r.start, r.end, 0, r.size - 1);
            }
            
        }
        return W;
    }
    else if (mt == FIXATION) {

        WrightFisher::Matrix W(2 * N, 1, mt);
        return W;
    }
    else { // BOTH

        WrightFisher::Matrix W((2 * N) - 1, 2, mt);
        return W;
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

} //}}}