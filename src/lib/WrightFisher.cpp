#include "WrightFisher.hpp"

double WrightFisher::psi_diploid(const llong i, const llong N, const double s, const double h, const double u, const double v) 
{

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

WrightFisher::Row WrightFisher::binom_row(const llong size, const double p, const double alpha) 
{

    // start and end quantiles for covering 1 - alpha weight of the probability mass
    llong start = (llong)binom_tail_cover(alpha / 2, size, p, true);
    llong end = (llong)binom_tail_cover(alpha / 2, size, p, false);

    // make sure we didn't mess up
    assert(start < end);

    // Initialize row
    Row r(start, end);

    // Start iterative binomial calculation (WFES supplementary, eq 18,19)
    double d = ld_binom(start, size, p);
    double lc = log(p) - log(1 - p);
    r.Q(0) = d;

    // Iterative binomial (in log)
    for(llong j = start + 1; j <= end; j++) {
        d += log(size - j + 1) - log(j) + lc;
        r.Q(j - start) = d;
    }

    // Exponentiate
    vdExp(r.Q.size(), r.Q.data(), r.Q.data());
    // Re-weigh to sum to 1
    r.weight = r.Q.sum();
    r.Q /= r.weight;

    return r;

}

WrightFisher::Matrix WrightFisher::Equilibrium(
        const llong N, const double s, const double h, const double u, const double v, 
        const double alpha, const bool verbose, const llong block_size) 
{
    time_point t_start, t_end;
    if (verbose) t_start = std::chrono::system_clock::now();
    llong N2 = 2 * N;
    llong size = N2 + 1;
    WrightFisher::Matrix W(size, size, n_absorbing(WrightFisher::NON_ABSORBING));
    for(llong block_row = 0; block_row < size; block_row += block_size) {
        llong block_length = (block_row + block_size) < size ? block_size : size - block_row;
        std::deque<Row> buffer(block_length);

        #pragma omp parallel for
        for(llong b = 0; b < block_length; b++) {
            llong i = b + block_row;
            buffer[b] = binom_row(2 * N, psi_diploid(i, N, s, h, u, v), alpha);
            Row& r = buffer[b];
            // I - Q
            for(llong j = 0; j < r.Q.size(); j++) r.Q(j) = -r.Q(j);
            // r.Q = -r.Q;
            // diagonal is set in the sequential block - since it may require a structural change to the matrix
            //r.Q(i - r.start) += 1;
        }

        for(llong b = 0; b < block_length; b++) {
            Row& r = buffer[b];
            llong i = b + block_row;

            // diagnoal is left of chunk - insert new entry before chunk
            if (i < r.start) W.Q.append_value(1, i, i);
            // diagonal overlaps chunk - increment element
            if (i >= r.start && i <= r.end) r.Q(i-r.start) += 1;
            // update chunk if it contains last column
            if (r.end == N2) r.Q(r.size - 1) = 1;
            // append large chunk
            W.Q.append_chunk(r.Q, r.start, r.end, 0, r.size - 1);
            // diagonal is right of chunk - insert new entry after chunk
            if (i > r.end)   W.Q.append_value(1, i, i);
            // add a column of 1s on the end
            if (r.end != N2) W.Q.append_value(1, i, N2);

            W.Q.next_row();
        }
    }
    if (verbose) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Time to build matrix: " << dt.count() << " s" << std::endl;
    }
    return W;
}


WrightFisher::Matrix WrightFisher::Single(
        const llong Nx, const llong Ny, const absorption_type abs_t, const double s, const double h, const double u, const double v, 
        bool recurrent_mutation, const double alpha, const bool verbose, const llong block_size) 
{
    time_point t_start, t_end;
    if (verbose) t_start = std::chrono::system_clock::now();
    bool verify_diagonal = (Nx == Ny);
    llong Nx2 = 2 * Nx; 
    llong Ny2 = 2 * Ny; 
    llong size = Nx2 + 1;

    llong n_abs = n_absorbing(abs_t);

    Matrix W(Nx2 + 1 - n_abs, Ny2 + 1 - n_abs, n_abs);

    for(llong block_row = 0; block_row <= Nx2; block_row += block_size) {
        llong block_length = (block_row + block_size) < size ? block_size : size - block_row;
        std::deque<Row> buffer(block_length);

        #pragma omp parallel for
        for(llong b = 0; b < block_length; b++) {
            llong i = b + block_row;
            if(!recurrent_mutation && i != 0) {
                buffer[b] = binom_row(2 * Ny, psi_diploid(i, Nx, s, h, 0, 0), alpha);    
            } else {
                buffer[b] = binom_row(2 * Ny, psi_diploid(i, Nx, s, h, u, v), alpha);    
            }
        }

        for(llong b = 0; b < block_length; b++) {
            Row& r = buffer[b];
            llong i = b + block_row;
            llong r_last = r.size - 1;

            // diagonal is left of inserted chunk
            if (verify_diagonal && (i < r.start)) W.Q.append_value(0, i, i);

            switch(abs_t) {
                case NON_ABSORBING:
                    // Include full row
                    W.Q.append_chunk(r.Q, r.start, r.end, 0, r_last);
                break;

                case EXTINCTION_ONLY:
                    // Do not include 0th row and column
                    if (i == 0) continue;
                    else {
                        if (r.start == 0) {
                            W.Q.append_chunk(r.Q, r.start, r.end - 1, 1, r_last);
                            W.R(i - 1, 0) = r.Q(0);
                        } else {
                            W.Q.append_chunk(r.Q, r.start - 1, r.end - 1, 0, r_last);
                        }
                    }
                break;

                case FIXATION_ONLY:
                    // Do not include Nx2th row and column
                    if (i == Nx2) continue;
                    else {
                        if (r.end == Ny2) {
                            W.Q.append_chunk(r.Q, r.start, r.end - 1, 0, r_last - 1);
                            W.R(i, 0) = r.Q(r_last);
                        } else {
                            W.Q.append_chunk(r.Q, r.start, r.end, 0, r_last);
                        }
                    }
                break;

                case BOTH_ABSORBING:
                    // Do not include 0th and Nx2th row and column
                    if (i == 0 || i == Nx2) continue;
                    else {
                        if (r.start == 0 && r.end == Ny2) {
                            W.Q.append_chunk(r.Q, r.start, r.end - 2, 1, r_last - 1);
                            W.R(i - 1, 0) = r.Q(0);
                            W.R(i - 1, 1) = r.Q(r_last);
                        } else if (r.start == 0) {
                            W.Q.append_chunk(r.Q, r.start, r.end - 1, 1, r_last);
                            W.R(i - 1, 0) = r.Q(0);
                        } else if (r.end == Ny2) {
                            W.Q.append_chunk(r.Q, r.start - 1, r.end - 2, 0, r_last - 1);
                            W.R(i - 1, 1) = r.Q(r_last);
                        } else {
                            W.Q.append_chunk(r.Q, r.start - 1, r.end - 1, 0, r_last);
                        }
                    }
                break;
            }

            // diagonal on the right
            if (verify_diagonal && (i > r.end)) W.Q.append_value(0, i, i);
            W.Q.next_row();
        }
    }

    if (verbose) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Time to build matrix: " << dt.count() << " s" << std::endl;
    }
    return W;
}

std::deque<std::pair<llong, llong>> submatrix_indeces(const lvec& sizes) 
{
    llong i = 0;
    llong j = 0;

    llong size = sizes.sum();

    std::deque<std::pair<llong, llong>>idx (size);

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

WrightFisher::Matrix WrightFisher::Switching(
        const lvec& N, const absorption_type abs_t, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, 
        double alpha, const bool verbose, const llong block_size) 
{
    time_point t_start, t_end;
    if (verbose) t_start = std::chrono::system_clock::now();

    llong k = N.size();

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
    lvec sizes        = 2 * N + lvec::Ones(k);
    llong size        = sizes.sum();

    Matrix W(size - n_abs_total, size - n_abs_total, n_abs_total);
    std::deque<std::pair<llong, llong>> index = submatrix_indeces(sizes);

    for(llong block_row = 0; block_row <= size; block_row += block_size) {
        llong block_length = (block_row + block_size) < size ? block_size : size - block_row;
        std::deque<std::deque<Row>> buffer(block_length);
        for(llong b = 0; b < block_length; b++) buffer[b] = std::deque<Row>(k);

        #pragma omp parallel for
        for(llong b = 0; b < block_length; b++) {

            llong row = b + block_row;
            llong i   = index[row].first; // model index
            llong im  = index[row].second; // current index within model i

            for(llong j = 0; j < k; j++) {
                double p        = psi_diploid(im, N(i), s(j), h(j), u(j), v(j));
                buffer[b][j]    = binom_row(2 * N(j), p, alpha);
                buffer[b][j].Q *= switching(i, j);
            }
        }


        for(llong b = 0; b < block_length; b++) {
            llong row    = b + block_row;
            llong i      = index[row].first; // model index
            llong im     = index[row].second; // current index within model i
            llong offset = 0; // coordinate of the submodel start

            // BEGIN ITERATOR ROW
            for(llong j = 0; j < k; j++) {
                Row& r               = buffer[b][j];
                bool row_complete    = (j == (k - 1));
                llong m_start        = r.start + offset;
                llong m_end          = r.end + offset;
                llong r_last         = r.size - 1;
                bool verify_diagonal = (i == j);

                if (verify_diagonal && (im < r.start)) W.Q.append_value(0, im + offset, im + offset);

                switch(abs_t) {
                    case NON_ABSORBING:
                        W.Q.append_chunk(r.Q, m_start, m_end, 0, r_last);
                    break;

                    case EXTINCTION_ONLY:
                        if (im == 0) continue;
                        else {
                            if (r.start == 0) {
                                W.Q.append_chunk(r.Q, m_start, m_end - 1, 1, r_last);
                                W.R(row - (i + 1), j) = r.Q(0);
                            } else {
                                W.Q.append_chunk(r.Q, m_start - 1, m_end - 1, 0, r_last);
                            }
                        }
                    break;

                    case FIXATION_ONLY:
                        if (im == N(i) * 2) continue;
                        else {
                            if (r.end == N(j) * 2) {
                                W.Q.append_chunk(r.Q, m_start, m_end - 1, 0, r_last - 1);
                                W.R(row - i, j) = r.Q(r_last);
                            } else {
                                W.Q.append_chunk(r.Q, m_start, m_end, 0, r_last);
                            }
                        }
                    break;

                    case BOTH_ABSORBING:
                        if (im == 0 || im == N(i) * 2) continue;
                        else {
                            if (r.start == 0 && r.end == N(j) * 2) {
                                W.Q.append_chunk(r.Q, m_start, m_end - 2, 1, r_last - 1);
                                // TODO: why is this not `row` ?
                                W.R(row - i - i - 1, 2 * j) = r.Q(0);
                                W.R(row - i - i - 1, (2 * j) + 1) = r.Q(r_last);
                            } else if (r.start == 0) {
                                W.Q.append_chunk(r.Q, m_start, m_end - 1, 1, r_last);
                                W.R(row - i - i - 1, 2 * j) = r.Q(0);
                            } else if (r.end == N(j) * 2) {
                                W.Q.append_chunk(r.Q, m_start - 1, m_end - 2, 0, r_last - 1);
                                W.R(row - i - i - 1, (2 * j) + 1) = r.Q(r_last);
                            } else {
                                W.Q.append_chunk(r.Q, m_start - 1, m_end - 1, 0, r_last);
                            }
                        }
                    break;
                }
                if (verify_diagonal && (im > r.end)) W.Q.append_value(0, im + offset, im + offset);
                // Increment row offset
                offset += (sizes(j) - n_absorbing(abs_t));
                // This needs to be inside in case block_sizes are unbalanced
                if (row_complete) W.Q.next_row();
            } // END ROW
        }
    }
    
    if (verbose) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Time to build matrix: " << dt.count() << " s" << std::endl;
    }
    return W;
}

WrightFisher::Matrix WrightFisher::NonAbsorbingToFixationOnly(
        const llong N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, 
        const double alpha, const bool verbose, const llong block_size) 
{
    time_point t_start, t_end;
    if (verbose) t_start = std::chrono::system_clock::now();

    // TODO: proper error checking
    assert(s.size() == 2);
    // forward mutation rate should be above 0
    for(llong i = 0; i < 2; i++) assert(v(i) > 0);

    lvec sizes(2); sizes << (2 * N) + 1, 2 * N;
    llong size = sizes.sum();

    Matrix W(size, size, 1);
    std::deque<std::pair<llong, llong>> index = submatrix_indeces(sizes);
    
    for(llong block_row = 0; block_row < size; block_row += block_size) {
        llong block_length = (block_row + block_size) < size ? block_size : size - block_row;
        std::deque<Row> b_1(block_length);
        std::deque<Row> b_2(block_length);

        #pragma omp parallel for
        for(llong b = 0; b < block_length; b++) {
            llong row = block_row + b;
            llong i = index[row].first; // model index
            llong im = index[row].second; // current index within model i

            Row r_1 = binom_row(2 * N, psi_diploid(im, N, s(0), h(0), u(0), v(0)), alpha);
            r_1.Q *= switching(i, 0);
            b_1[b] = r_1;

            Row r_2 = binom_row(2 * N, psi_diploid(im, N, s(1), h(1), u(1), v(1)), alpha);
            r_2.Q *= switching(i, 1);
            b_2[b] = r_2;
        }

        for(llong b = 0; b < block_length; b++) {
            llong row = block_row + b;
            llong offset = (2 * N) + 1;

            W.Q.append_chunk(b_1[b].Q, b_1[b].start, b_1[b].end, 0, b_1[b].size - 1);

            if (b_2[b].end == N * 2) {
                W.Q.append_chunk(b_2[b].Q, b_2[b].start + offset, b_2[b].end + offset - 1, 0, b_2[b].size - 2);
                W.R(row, 0) = b_2[b].Q(b_2[b].size - 1);
            } else {
                W.Q.append_chunk(b_2[b].Q, b_2[b].start + offset, b_2[b].end + offset, 0, b_2[b].size - 1);
            }
            W.Q.next_row();
        }
    }

    if (verbose) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Time to build matrix: " << dt.count() << " s" << std::endl;
    }
    return W;
}

WrightFisher::Matrix WrightFisher::NonAbsorbingToBothAbsorbing(
        const llong N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, 
        const double alpha, const bool verbose, const llong block_size) 
{
    time_point t_start, t_end;
    if (verbose) t_start = std::chrono::system_clock::now();

    // TODO: proper error checking
    assert(s.size() == 2);
    // forward mutation rate should be above 0
    for(llong i = 0; i < 2; i++) assert(v(i) > 0);

    lvec sizes(2); sizes << (2 * N) + 1, (2 * N) - 1;
    llong size = sizes.sum();

    Matrix W(size, size, 2);
    std::deque<std::pair<llong, llong>> index = submatrix_indeces(sizes);
    
    for(llong block_row = 0; block_row < size; block_row += block_size) {
        llong block_length = (block_row + block_size) < size ? block_size : size - block_row;
        std::deque<Row> buffer_1(block_length);
        std::deque<Row> buffer_2(block_length);

        #pragma omp parallel for
        for(llong b = 0; b < block_length; b++) {
            llong row = block_row + b;
            llong i = index[row].first; // model index
            llong im = index[row].second + i; // current index within model i, correct for non-absorbing starting state - ugly

            Row r_1 = binom_row(2 * N, psi_diploid(im, N, s(0), h(0), u(0), v(0)), alpha);
            r_1.Q *= switching(i, 0);
            buffer_1[b] = r_1;

            Row r_2 = binom_row(2 * N, psi_diploid(im, N, s(1), h(1), u(1), v(1)), alpha);
            r_2.Q *= switching(i, 1);
            buffer_2[b] = r_2;
        }

        for(llong b = 0; b < block_length; b++) {
            llong row = block_row + b;
            llong offset = (2 * N) + 1;

            W.Q.append_chunk(buffer_1[b].Q, buffer_1[b].start, buffer_1[b].end, 0, buffer_1[b].size - 1);

            if (buffer_2[b].start == 0 && buffer_2[b].end == 2*N) {
                W.R(row, 0) = buffer_2[b].Q(0);
                W.R(row, 1) = buffer_2[b].Q(buffer_2[b].size - 1);
                W.Q.append_chunk(buffer_2[b].Q, buffer_2[b].start + offset,     buffer_2[b].end + offset - 2, 1, buffer_2[b].size - 2);
            } else if (buffer_2[b].start == 0) {
                W.Q.append_chunk(buffer_2[b].Q, buffer_2[b].start + offset,     buffer_2[b].end + offset - 1, 1, buffer_2[b].size - 1);
                W.R(row, 0) = buffer_2[b].Q(0);
            } else if (buffer_2[b].end == N * 2) {
                W.Q.append_chunk(buffer_2[b].Q, buffer_2[b].start + offset - 1, buffer_2[b].end + offset - 2, 0, buffer_2[b].size - 2);
                W.R(row, 1) = buffer_2[b].Q(buffer_2[b].size - 1);
            } else {
                W.Q.append_chunk(buffer_2[b].Q, buffer_2[b].start + offset - 1, buffer_2[b].end + offset - 1, 0, buffer_2[b].size - 1);
            }
            W.Q.next_row();
        }
    }

    if (verbose) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Time to build matrix: " << dt.count() << " s" << std::endl;
    }
    return W;
}