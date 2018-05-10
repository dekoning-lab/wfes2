#pragma once

#include <cmath>
#include <cfloat>

#define W_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define W_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define W_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p))  : W_D_Lval(p))

#define W_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) : W_D_Cval(p))

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define error_print(M, ...) fprintf(stderr, ANSI_COLOR_RED "[ERROR]" ANSI_COLOR_RESET " (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define DIST_TAIL_CORRECTION 15

double q_norm(double p, double mu, double sigma, bool lower_tail, bool log_p);
double ld_binom(double k, double n, double pr);
double d_binom(double k, double n, double pr);
double q_binom_approx(double p, double n, double pr, bool lower_tail);
double binom_tail_cover(double p, double n, double pr, bool lower_tail);
