#include <motion_primitive_library/primitive/math.h>
#include <motion_primitive_library/primitive/rpoly.h>

int factorial(int n) {
  int nf = 1;
  while(n > 0) {
    nf *= n;
    n --;
  }
  return nf;
}

decimal_t power(decimal_t t, int n) {
  decimal_t tn = 1;
  while(n > 0) {
    tn *= t;
    n --;
  }
  return tn;
  //return n <= 0 ? 1 : power(t, n-1);
}


/* **************************************************************** */
int quad(decimal_t b,
         decimal_t c,
         decimal_t d,
         std::array<decimal_t, 4>& out) {
  decimal_t p = c*c - 4*b*d;
  if(p < 0)
    return 0;
  else
  {
    out[0] = (-c-sqrt(p))/(2*b);
    out[1] = (-c+sqrt(p))/(2*b);
    return 2;
  }
}

/* **************************************************************** */
int cubic(decimal_t a,
          decimal_t b,
          decimal_t c,
          decimal_t d,
          std::array<decimal_t, 4>& out) {
  decimal_t a2 = b / a;
  decimal_t a1 = c / a;
  decimal_t a0 = d / a;
  //printf("a: %f, b: %f, c: %f, d: %f\n", a, b, c, d);

  decimal_t Q = (3*a1-a2*a2)/9;
  decimal_t R = (9*a1*a2-27*a0-2*a2*a2*a2)/54;
  decimal_t D = Q*Q*Q + R*R;
  //printf("R: %f, Q: %f, D: %f\n", R, Q, D);
  if(D > 0) {
    decimal_t S = std::cbrt(R+sqrt(D));
    decimal_t T = std::cbrt(R-sqrt(D));
    //printf("S: %f, T: %f\n", S, T);
    out[0] = -a2/3+(S+T);
    return 1;
  }
  else if(D == 0) {
    decimal_t S = std::cbrt(R);
    out[0] = -a2/3+S+S;
    out[1] = -a2/3-S;
    return 2;
  }
  else {
    decimal_t theta = acos(R/sqrt(-Q*Q*Q));
    out[0] = 2*sqrt(-Q)*cos(theta/3)-a2/3;
    out[1] = 2*sqrt(-Q)*cos((theta+2*M_PI)/3)-a2/3;
    out[2] = 2*sqrt(-Q)*cos((theta+4*M_PI)/3)-a2/3;
    return 3;
  }
}


/* **************************************************************** */
int quartic(decimal_t a,
            decimal_t b,
            decimal_t c,
            decimal_t d,
            decimal_t e,
            std::array<decimal_t, 4>& out) {
  decimal_t a3 = b / a;
  decimal_t a2 = c / a;
  decimal_t a1 = d / a;
  decimal_t a0 = e / a;

  std::array<decimal_t, 4> ys;
  cubic(1, -a2, a1*a3-4*a0, 4*a2*a0-a1*a1-a3*a3*a0, ys);
  decimal_t y1 = ys[0];
  //printf("y1: %f\n", y1);
  decimal_t r = a3*a3/4-a2+y1;
  //printf("r: %f\n", r);

  //printf("a = %f, b = %f, c = %f, d = %f, e = %f\n", a, b, c, d, e);
  if(r < 0)
    return 0;

  decimal_t R = sqrt(r);
  decimal_t D, E;
  if(R != 0) {
    D = sqrt(0.75*a3*a3-R*R-2*a2+0.25*(4*a3*a2-8*a1-a3*a3*a3)/R);
    E = sqrt(0.75*a3*a3-R*R-2*a2-0.25*(4*a3*a2-8*a1-a3*a3*a3)/R);
  }
  else {
    D = sqrt(0.75*a3*a3-2*a2+2*sqrt(y1*y1-4*a0));
    E = sqrt(0.75*a3*a3-2*a2-2*sqrt(y1*y1-4*a0));
  }

  int results = 0;
  if(!std::isnan(D)) {
    out[0] = -a3/4+R/2+D/2;
    out[1] = -a3/4+R/2-D/2;
    results += 2;
  }
  if(!std::isnan(E)) {
    out[results] = -a3/4-R/2+E/2;
    out[results + 1] = -a3/4-R/2-E/2;
    results += 2;
  }

  return results;
}

/* **************************************************************** */
int solve_fast(decimal_t a,
               decimal_t b,
               decimal_t c,
               decimal_t d,
               decimal_t e,
               std::array<decimal_t, 4>& out) {
  if(a != 0)
    return quartic(a, b, c, d, e, out);
  else if(b != 0)
    return cubic(b, c, d, e, out);
  else if(c != 0)
    return quad(c, d, e, out);
  else if(d != 0)
  {
    out[0] = -e/d;
    return 1;
  }
  else
    return 0;
}

std::vector<decimal_t> solve(decimal_t a,
                             decimal_t b,
                             decimal_t c,
                             decimal_t d,
                             decimal_t e) {
    std::array<decimal_t, 4> roots;
    int results = solve_fast(a, b, c, d, e, roots);

    std::vector<decimal_t> ts;
    for (int i = 0; i < results; i++) {
        ts.push_back(roots[i]);
    }

    return ts;
}

/* **************************************************************** */
int solve_fast(decimal_t a,
               decimal_t b,
               decimal_t c,
               decimal_t d,
               decimal_t e,
               decimal_t f,
               decimal_t g,
               std::array<decimal_t, 6>& out) {
  if(a == 0 && b == 0) {
    std::array<decimal_t, 4> out_small;
    int result = solve_fast(c, d, e, f, g, out_small);
    for (int i = 0; i < result; i++) {
      out[i] = out_small[i];
    }
    return result;
  }
  else {
    struct RPoly_State* rpoly_state = real_poly_alloc(6);

    double p[7];
    p[0] = a;
    p[1] = b;
    p[2] = c;
    p[3] = d;
    p[4] = e;
    p[5] = f;
    p[6] = g;

    double zeror[6];
    double zeroi[6];
    int root_n = real_poly_roots_compute(p, 6, rpoly_state, zeror, zeroi);
    int rpoly2_roots = 0;
    for (int i = 0; i < root_n; i++) {
        if (std::abs(zeroi[i]) < 0.0000001) {
            out[rpoly2_roots] = zeror[i];
            rpoly2_roots++;
        }
    }

    real_poly_release(rpoly_state);
    return rpoly2_roots;
  }
}

std::vector<decimal_t> solve(decimal_t a,
                             decimal_t b,
                             decimal_t c,
                             decimal_t d,
                             decimal_t e,
                             decimal_t f,
                             decimal_t g) {
  std::vector<decimal_t> ts;
  std::array<decimal_t, 6> ts_arr;
  int roots = solve_fast(a, b, c, d, e, f, g, ts_arr);
  for (int i = 0; i < roots; i++) {
      ts.push_back(ts_arr[i]);
  }
  return ts;
}
