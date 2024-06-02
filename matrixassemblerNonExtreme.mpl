# Matrices assembling for Lee-Wick black-holes' quasi-normal modes computation
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  s::numeric, # s : spin (0, 1, 2)
  L::numeric, # L : angular momentum, L >= s
  q::numeric, # q : Lee-Wick mass
  a::numeric, # a : real Lee-Wick mass
  p::string   # p : string containing the path where we save the assembled matrices
  )
  if L < s then
    error(1, "Angular momentum cannot be smaller than spin!");
  end if:
  local e::numeric, h::function, f::function, xh::numeric, at::numeric, fp::function, V::function, F::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Physical parameters:
  e := 1 - s^2: # self-explaining
  # Determination of the event horizon:
  h := x -> 1 - exp(-a*x)*((1 + 1/2*a*(q^2 + 1)*x)*cos(q*a*x) + 1/2*(1 - q^2 + a*(q^2 + 1)*x)*sin(q*a*x)/q):
  xh := fsolve(x - h(x) = 0, x = 2.0):
  # Other parameters appearing in the model:
  at := ((q*a*xh*(q^2 + 1) + 2*q)*cos(q*a*xh) + (a*xh*(q^2 + 1) + 1 - q^2)*sin(q*a*xh))/((q*a*xh*(q^2 + 1) + 2*q)*cos(q*a*xh) + (a*xh*(q^2 + 1)*(a*(q^2 + 1)*(xh - 1) + 1) + 1 - q^2)*sin(q*a*xh)):
  f := y -> 1/2*(2*xh - 1 + y)/xh + exp(-2*a*xh/(1 - y))*(1/2*((1 - y)/xh + a*(q^2 + 1))*cos(2*q*a*xh/(1 - y)) + 1/2*(1/2*(-q^2 + 1)*(1 - y)/xh + a*(q^2 + 1))*sin(2*q*a*xh/(1 - y))/q):
  fp := y -> 1/(2*xh) - 2*a*xh*exp(-2*a*xh/(1 - y))*((((1 - y)/xh + a*(q^2 + 1))*cos(2*q*a*xh/(1 - y)))/2 + ((((-q^2 + 1)*(1 - y))/(2*xh) + a*(q^2 + 1))*sin(2*q*a*xh/(1 - y)))/(2*q))/(1 - y)^2 + exp(-2*a*xh/(1 - y))*(-cos(2*q*a*xh/(1 - y))/(2*xh) - ((1 - y)/xh + a*(q^2 + 1))*q*a*xh*sin(2*q*a*xh/(1 - y))/(1 - y)^2 - ((-q^2 + 1)*sin(2*q*a*xh/(1 - y)))/(4*xh*q) + (((-q^2 + 1)*(1 - y))/(2*xh) + a*(q^2 + 1))*a*xh*cos(2*q*a*xh/(1 - y))/(1 - y)^2):
  V := y -> 1/16*f(y)*(1 - y)^2*(e*(1 - y)*fp(y) + L*(L + 1)):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -4*V(y)/((1 + y)*(1 - y)^2):
  L01 := y -> 1/4*f(y)/(1+y)*(-2*(1-y)*f(y) + (1-y)^2*fp(y)):
  L02 := y -> 1/4*(1 - y)^2*f(y)^2/(1 + y):
  L10 := y -> 1/2*f(y)*(fp(y)*(2*xh + 1 - y) - f(y) + xh*at*(1 - y)*((3 + y)*f(y) - (-y^2 + 1)*fp(y))/(1 + y)^2)/(1 + y):
  L11 := y -> f(y)^2*((2*xh + 1 - y)/(1 + y) - xh*at*(1 - y)^2/(1 + y)^2):
  L12 := y -> 0:
  L20 := y -> 4*xh^2/((1 + y)*(1 - y)^2) - ((1 + y)*(2*xh + 1 - y) - xh*at*(1 - y)^2)^2*f(y)^2/((1 + y)^3*(1 - y)^2):
  L21 := y -> 0:
  L22 := y -> 0:
  F := y -> simplify(add(a[j]*ChebyshevT(j, y), j=0..n-1)):
  M0 := Matrix(n):
  M1 := Matrix(n):
  M2 := Matrix(n):
  for i from 1 to n do
    xi := cos((2.0*i-1.0)*Pi/(2.0*n)); # Chebyshev roots collocation points
    expr0 := evalf(L00(xi)*F(xi) + L01(xi)*subs(x=xi, diff(F(x),x)) + L02(xi)*subs(x=xi, diff(F(x),x$2))):
    expr1 := evalf(L10(xi)*F(xi) + L11(xi)*subs(x=xi, diff(F(x),x)) + L12(xi)*subs(x=xi, diff(F(x),x$2))):
    expr2 := evalf(L20(xi)*F(xi) + L21(xi)*subs(x=xi, diff(F(x),x)) + L22(xi)*subs(x=xi, diff(F(x),x$2))):
    for j from 1 to n do
      M0[i,j] := coeff(expr0, a[j-1]):
      M1[i,j] := coeff(expr1, a[j-1]):
      M2[i,j] := coeff(expr2, a[j-1]):
    end do:
  end do:
  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(n, string);
  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
end proc: