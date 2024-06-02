# Matrices assembling for Lee-Wick black-holes' quasi-normal modes computation, the extreme case
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  s::numeric, # s : spin (0, 1, 2)
  L::numeric, # L : angular momentum, L >= s
  q::numeric, # q : Lee-Wick mass
  p::string   # p : string containing the path where we save the assembled matrices
  )
  if L < s then
    error(1, "Angular momentum cannot be smaller than spin!");
  end if:
  local e::numeric, h::function, f::function, xe::numeric, a::numeric, alpha::numeric, h1::numeric, fp::function, V::function, F::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string, eta::function, etap::function, etapp::function:
  with(LinearAlgebra):
  Digits := d:
  # Physical parameters:
  e := 1 - s^2: # self-explaining
  # Computation of the dimensionless mass and event horizon parameters:
  h := (xe, alpha) -> 1 - exp(-alpha*xe)*((1 + alpha/2*(1 + q^2)*xe)*cos(alpha*q*xe) + 1/(2*q)*(1 - q^2 + alpha*(1 + q^2)*xe)*sin(alpha*q*xe)):
  f := (xe, alpha) -> 1 - h(xe,alpha)/xe:
  sol := fsolve({f(xe,alpha) = 0, diff(f(xe,alpha),xe) = 0}, {xe = 0.5..0.8, alpha = 4.9..5.5}):
  assign(sol); # we assign the values of critical parameters
  printf("Found BH's event horizon xh       = %f\n", xe);
  printf("Found BH's mass generalized alpha = %f\n", alpha);
  # Other parameters appearing in the model:
  h1 := 1/2*(alpha^3*q^4*xe^3 - alpha^3*q^4*xe^2 + 2*alpha^3*q^2*xe^3 - 2*alpha^3*q^2*xe^2 + alpha^3*xe^3 + 2*alpha^2*q^2*xe^2 - alpha^3*xe^2 + 2*alpha^2*xe^2 - 2*alpha*q^2*xe + 2*alpha*xe - 2)/(alpha*q^2*xe + alpha*xe + 2):
  a := 1 + xe*(2*alpha^4*xe^3*(xe - 1)*(q^2 + 1)^2 - alpha^3*xe^2*(q^2 + 1)*(q^2 - 4*xe + 1) - 2*alpha^2*xe^2*(q^2 - 3) - 4*alpha*xe*q^2 - 6)/(6*h1^2*(alpha*xe*(q^2 + 1) + 2)):
  f := y -> (2*xe - 1 + y)/(2*xe) + exp(-2*alpha*xe/(1 - y))*(1/2*((1 - y)/xe + alpha*(q^2 + 1))*cos(2*q*alpha*xe/(1 - y)) + 1/(2*q)*(1/(2*xe)*(-q^2 + 1)*(1 - y) + alpha*(q^2 + 1))*sin(2*q*alpha*xe/(1 - y))):
  fp := y -> 1/(2*xe) - 2*alpha*xe*exp(-2*alpha*xe/(1 - y))*((((1 - y)/xe + alpha*(q^2 + 1))*cos(2*q*alpha*xe/(1 - y)))/2 + ((((-q^2 + 1)*(1 - y))/(2*xe) + alpha*(q^2 + 1))*sin(2*q*alpha*xe/(1 - y)))/(2*q))/(1 - y)^2 + exp(-2*alpha*xe/(1 - y))*(-cos(2*q*alpha*xe/(1 - y))/(2*xe) - ((1 - y)/xe + alpha*(q^2 + 1))*q*alpha*xe*sin(2*q*alpha*xe/(1 - y))/(1 - y)^2 - ((-q^2 + 1)*sin(2*q*alpha*xe/(1 - y)))/(4*q*xe) + (((-q^2 + 1)*(1 - y))/(2*xe) + alpha*(q^2 + 1))*alpha*xe*cos(2*q*alpha*xe/(1 - y))/(1 - y)^2):
  V := y -> (1 - y)^2/16*f(y)*(e*(1 - y)*fp(y) + L*(L + 1)):
  eta := y -> (1 + y)/(1 - y) + (1 - y)/(h1*(1 + y)):
  etap := y -> ((2*h1 - 2)*y^2 + (4*h1 + 4)*y + 2*h1 - 2)/((-1 + y)^2*h1*(1 + y)^2):
  etapp := y -> ((-4*h1 + 4)*y^3 + (-12*h1 - 12)*y^2 + (-12*h1 + 12)*y - 4*h1 - 4)/((-1 + y)^3*h1*(1 + y)^3):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -4*V(y)/(-y^2 + 1)^2:
  L01 := y -> 1/4*f(y)/(1+y)^2*(-2*(1-y)*f(y) + (1-y)^2*fp(y)):
  L02 := y -> (1 - y)^2/(4*(1 + y)^2)*f(y)^2:
  L10 := y -> xe/2*(1 - y)^2*f(y)*fp(y)*etap(y)/(1 + y)^2 + xe/2*(1 - y)*f(y)^2*((1 - y)*etapp(y) - 2*etap(y))/(1 + y)^2 + 1/2*(1 - y)*(2 - a*(1 - y))*f(y)*fp(y)/(1 + y)^3 - f(y)^2*(2 - 2*a*(1 - y) + 1/2*a*(1 - y)^2)/(1 + y)^4:
  L11 := y -> (1 - y)*f(y)^2*(xe*(-y^2 + 1)*etap(y) + 2 - a*(1 - y))/(1 + y)^3:
  L12 := y -> 0:
  L20 := y -> 4*xe^2/(-y^2 + 1)^2 - f(y)^2*(xe*(-y^2 + 1)*etap(y) + 2 - a*(1 - y))^2/(1 + y)^4:
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