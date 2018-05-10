#include "rdist.hpp"

double q_norm(double p, double mu, double sigma, bool lower_tail, bool log_p)
{
    double p_, q, r, val;

    if (std::isnan(p) || std::isnan(mu) || std::isnan(sigma)) return p + mu + sigma;
    if (log_p) {
        if(p > 0) return NAN;
        if(p == 0) return lower_tail ? INFINITY : -INFINITY;
        if(p == -INFINITY) return lower_tail ? -INFINITY : INFINITY;
    } else { /* !log_p */
        if(p < 0 || p > 1) return NAN;
        if(p == 0) return lower_tail ? -INFINITY : INFINITY;
        if(p == 1)	return lower_tail ? INFINITY : -INFINITY;
    }

    if(sigma  < 0)	return NAN;
    if(sigma == 0)	return mu;

    p_ = W_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;


/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
                q * (((((((r * 2509.0809287301226727 +
                           33430.575583588128105) * r + 67265.770927008700853) * r +
                         45921.953931549871457) * r + 13731.693765509461125) * r +
                       1971.5909503065514427) * r + 133.14166789178437745) * r +
                     3.387132872796366608)
                / (((((((r * 5226.495278852854561 +
                         28729.085735721942674) * r + 39307.89580009271061) * r +
                       21213.794301586595867) * r + 5394.1960214247511077) * r +
                     687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

        /* r = min(p, 1-p) < 0.075 */
        if (q > 0) r = W_DT_CIv(p);/* 1-p */
        else r = p_;/* = W_DT_Iv(p) ^=  p */

        r = sqrt(- ((log_p &&
                     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
                    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177) *
                       r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734)
                  / (((((((r *
                           1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                          r + .0151986665636164571966) * r +
                         .14810397642748007459) * r + .68976733498510000455) *
                       r + 1.6763848301838038494) * r +
                      2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                      r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
                  / (((((((r *
                           2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                          r + 1.8463183175100546818e-5) * r +
                         7.868691311456132591e-4) * r + .0148753612908506148525)
                       * r + .13692988092273580531) * r +
                      .59983220655588793769) * r + 1.);
        }

        if(q < 0.0) val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

double ld_binom(double k, double n, double pr)
{
    double log_choose = lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1));
    return log_choose + (k * log(pr)) + ((n - k) * log1p(-pr));
}

double d_binom(double k, double n, double pr)
{
    return exp(ld_binom(k, n, pr));
}

double q_binom_approx(double p, double n, double pr, bool lower_tail)
{
    double mu = n * pr;
    double sigma = sqrt(n * pr * (1.0 - pr));
    double q = q_norm(p, mu, sigma, lower_tail, false);
    if(q < 0) return 0;
    if(q > n) return n;
    else return round(q);
}

double binom_tail_cover(double p, double n, double pr, bool lower_tail)
{
    double q = q_binom_approx(p, n, pr, lower_tail) + (lower_tail ? -DIST_TAIL_CORRECTION : DIST_TAIL_CORRECTION);
    if(q < 0) return 0;
    if(q > n) return n;
    else return q;
}


