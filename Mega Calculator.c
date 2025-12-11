#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

/* ---------- small helpers ---------- */

void flush_stdin(void) {
    int c;
    while ((c = getchar()) != '\n' && c != EOF) {}
}

double read_double(const char *prompt) {
    double v;
    printf("%s", prompt);
    if (scanf("%lf", &v) != 1) { flush_stdin(); return 0.0; }
    return v;
}

long long read_ll(const char *prompt) {
    long long v;
    printf("%s", prompt);
    if (scanf("%lld", &v) != 1) { flush_stdin(); return 0; }
    return v;
}

/* ---------- basic utilities ---------- */

long long gcd_ll(long long a, long long b) {
    a = llabs(a); b = llabs(b);
    while (b) { long long t = a % b; a = b; b = t; }
    return a;
}

long long lcm_ll(long long a, long long b) {
    if (a == 0 || b == 0) return 0;
    return llabs(a / gcd_ll(a,b) * b);
}

long long factorial_ll(long long n) {
    if (n < 0) return -1;
    long long res = 1;
    for (long long i = 2; i <= n; ++i) res *= i;
    return res;
}

int is_prime_ll(long long n) {
    if (n <= 1) return 0;
    if (n <= 3) return 1;
    if (n % 2 == 0) return 0;
    for (long long i = 3; i * i <= n; i += 2)
        if (n % i == 0) return 0;
    return 1;
}

/* ---------- statistics: mean/median/mode/variance/stddev ---------- */

void statistics_menu(void) {
    int n;
    printf("Number of data points: ");
    if (scanf("%d", &n) != 1 || n <= 0) { flush_stdin(); printf("Invalid count.\n"); return; }
    double *a = malloc(sizeof(double) * n);
    for (int i = 0; i < n; ++i) { printf("x[%d]= ", i); scanf("%lf", &a[i]); }
    double sum = 0;
    for (int i = 0; i < n; ++i) sum += a[i];
    double mean = sum / n;

    /* median: sort copy */
    double *b = malloc(sizeof(double) * n);
    memcpy(b, a, sizeof(double) * n);
    for (int i = 0; i < n - 1; ++i)
        for (int j = i + 1; j < n; ++j)
            if (b[i] > b[j]) { double t = b[i]; b[i] = b[j]; b[j] = t; }

    double median = (n % 2) ? b[n/2] : (b[n/2 - 1] + b[n/2]) / 2.0;

    /* simple mode (one value) */
    double mode = b[0];
    int maxcount = 1;
    for (int i = 0; i < n; ) {
        int j = i + 1;
        while (j < n && fabs(b[j] - b[i]) < 1e-12) ++j;
        int cnt = j - i;
        if (cnt > maxcount) { maxcount = cnt; mode = b[i]; }
        i = j;
    }

    double varsum = 0;
    for (int i = 0; i < n; ++i) varsum += (a[i] - mean) * (a[i] - mean);
    double variance_pop = varsum / n;
    double variance_samp = (n > 1) ? varsum / (n - 1) : 0.0;

    printf("Mean = %.10g\nMedian = %.10g\nMode (one) = %.10g\nVariance(pop)= %.10g\nStdDev(pop)= %.10g\n",
           mean, median, mode, variance_pop, sqrt(variance_pop));

    free(a); free(b);
}

/* ---------- permutations/combinations (simple integer) ---------- */

long long nPr_ll(long long n, long long r) {
    if (r < 0 || n < 0 || r > n) return 0;
    long long res = 1;
    for (long long i = 0; i < r; ++i) res *= (n - i);
    return res;
}

long long nCr_ll(long long n, long long r) {
    if (r < 0 || n < 0 || r > n) return 0;
    if (r > n - r) r = n - r;
    long long res = 1;
    for (long long i = 1; i <= r; ++i) {
        res = res * (n - r + i) / i;
    }
    return res;
}

/* ---------- simple numeric and checks ---------- */

void palindrome_check(void) {
    long long v = read_ll("Enter integer: ");
    if (v < 0) { printf("Negative number: not considered palindrome here.\n"); return; }
    long long orig = v, rev = 0;
    while (v) { rev = rev * 10 + (v % 10); v /= 10; }
    printf("%lld is %s palindrome\n", orig, (rev == orig) ? "" : "not");
}

void prime_and_factors(void) {
    long long v = read_ll("Enter integer: ");
    printf("%lld is %sprime\n", v, is_prime_ll(v) ? "" : "not ");
    if (v == 0) { printf("0 has infinite factors.\n"); return; }
    long long m = llabs(v);
    printf("Prime factorization:\n");
    for (long long p = 2; p * p <= m; ++p) {
        int cnt = 0;
        while (m % p == 0) { m /= p; ++cnt; }
        if (cnt) printf("  %lld^%d\n", p, cnt);
    }
    if (m > 1) printf("  %lld^1\n", m);
}

/* ---------- simple algebra / equation ---------- */

void quadratic_solver(void) {
    double a = read_double("a= "); double b = read_double("b= "); double c = read_double("c= ");
    if (fabs(a) < 1e-15) {
        if (fabs(b) < 1e-15) { printf("No solution or infinite solutions.\n"); return; }
        printf("Linear solution: x = %.12g\n", -c / b); return;
    }
    double D = b*b - 4*a*c;
    if (D > 0) {
        printf("Two real roots: %.12g , %.12g\n", (-b + sqrt(D)) / (2*a), (-b - sqrt(D)) / (2*a));
    } else if (fabs(D) < 1e-15) {
        printf("One real root: %.12g\n", -b / (2*a));
    } else {
        double re = -b / (2*a), im = sqrt(-D) / (2*a);
        printf("Complex roots: %.12g + %.12gi , %.12g - %.12gi\n", re, im, re, im);
    }
}

/* ---------- simple geometry ---------- */

void polygon_area_coords(void) {
    int n;
    printf("Number of vertices (>=3): ");
    if (scanf("%d", &n) != 1 || n < 3) { flush_stdin(); printf("Invalid.\n"); return; }
    double *x = malloc(sizeof(double) * n);
    double *y = malloc(sizeof(double) * n);
    for (int i = 0; i < n; ++i) { printf("x y for vertex %d: ", i+1); scanf("%lf %lf", &x[i], &y[i]); }
    double area = 0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += x[i] * y[j] - x[j] * y[i];
    }
    area = fabs(area) * 0.5;
    printf("Polygon area = %.12g\n", area);
    free(x); free(y);
}

/* ---------- small matrix ops: add and multiply (simple) ---------- */

double *alloc_flat(int r, int c) {
    return malloc(sizeof(double) * r * c);
}

void matrix_ops_simple(void) {
    printf("1) Add 2) Multiply\nChoose: ");
    int op; if (scanf("%d", &op) != 1) { flush_stdin(); return; }
    int r1, c1, r2, c2;
    printf("A rows cols: "); scanf("%d %d", &r1, &c1);
    double *A = alloc_flat(r1, c1);
    printf("Enter A elements row-major:\n");
    for (int i = 0; i < r1*c1; ++i) scanf("%lf", &A[i]);

    printf("B rows cols: "); scanf("%d %d", &r2, &c2);
    double *B = alloc_flat(r2, c2);
    printf("Enter B elements row-major:\n");
    for (int i = 0; i < r2*c2; ++i) scanf("%lf", &B[i]);

    if (op == 1) {
        if (r1 != r2 || c1 != c2) { printf("Size mismatch.\n"); free(A); free(B); return; }
        printf("A + B =\n");
        for (int i = 0; i < r1; ++i) {
            for (int j = 0; j < c1; ++j) printf("%.6g ", A[i*c1 + j] + B[i*c1 + j]);
            printf("\n");
        }
    } else {
        if (c1 != r2) { printf("Size mismatch for multiply.\n"); free(A); free(B); return; }
        printf("A * B =\n");
        for (int i = 0; i < r1; ++i) {
            for (int j = 0; j < c2; ++j) {
                double s = 0;
                for (int k = 0; k < c1; ++k) s += A[i*c1 + k] * B[k*c2 + j];
                printf("%.6g ", s);
            }
            printf("\n");
        }
    }
    free(A); free(B);
}

/* ---------- currency converter (user supplies rate) ---------- */

void currency_converter_simple(void) {
    char from[16], to[16];
    printf("From currency code: "); scanf("%15s", from);
    printf("To currency code: "); scanf("%15s", to);
    double rate = read_double("Enter rate (1 from = ? to): ");
    double amt = read_double("Amount in source: ");
    printf("%.12g %s = %.12g %s (rate %.12g)\n", amt, from, amt * rate, to, rate);
}

/* ---------- menu and main ---------- */

void print_menu(void) {
    puts("\n--- Simple Calculator ---");
    puts("1) Basic stats (mean/median/mode/variance)");
    puts("2) Permutation/Combination");
    puts("3) Palindrome check");
    puts("4) Prime check & factorization");
    puts("5) Quadratic solver");
    puts("6) Polygon area (coords)");
    puts("7) Matrix add/multiply");
    puts("8) Currency converter (enter rate)");
    puts("9) GCD / LCM / factorial");
    puts("0) Exit");
    printf("Choose: ");
}

int main(void) {
    int running = 1;
    while (running) {
        print_menu();
        int opt;
        if (scanf("%d", &opt) != 1) { flush_stdin(); break; }
        switch (opt) {
            case 1: statistics_menu(); break;
            case 2: {
                long long n = read_ll("n= "); long long r = read_ll("r= ");
                printf("nPr = %lld\nnCr = %lld\n", nPr_ll(n,r), nCr_ll(n,r));
            } break;
            case 3: palindrome_check(); break;
            case 4: prime_and_factors(); break;
            case 5: quadratic_solver(); break;
            case 6: polygon_area_coords(); break;
            case 7: matrix_ops_simple(); break;
            case 8: currency_converter_simple(); break;
            case 9: {
                printf("1) GCD 2) LCM 3) Factorial\nChoose: ");
                int c; scanf("%d", &c);
                if (c == 1) { long long a = read_ll("a= "); long long b = read_ll("b= "); printf("gcd=%lld\n", gcd_ll(a,b)); }
                else if (c == 2) { long long a = read_ll("a= "); long long b = read_ll("b= "); printf("lcm=%lld\n", lcm_ll(a,b)); }
                else { long long n = read_ll("n= "); printf("%lld! = %lld\n", n, factorial_ll(n)); }
            } break;
            case 0: running = 0; break;
            default: printf("Unknown option.\n");
        }
    }
    printf("Bye.\n");
    return 0;
}
