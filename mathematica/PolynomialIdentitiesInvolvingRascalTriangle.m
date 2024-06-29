(* ::Package:: *)

BeginPackage["PolynomialIdentitiesInvolvingRascalTriangle`"]

RascalNumber::usage="Gives generalized rascal triangle via sum of product \\rascalNumber{n}{k}{i} = \\sum_{m=0}^{i} \\binom{n-k}{m} \\binom{k}{m}"
OneQBinomial::usage="Gives (1,q)-Binomial coefficient, see https://oeis.org/A096940"
ColumnIdentity1::usage="Validates column identity \\rascalNumber{n}{k}{i} = \\binom{n}{k}, \\quad 0 \\leq k \\leq i"
ColumnIdentity2::usage="Validates column identity \\rascalNumber{n}{i-j}{i}        &= \\binom{n}{i-j}"
RowIdentity1::usage="Validates column identity \\rascalNumber{2i+1-j}{k}{i}             &= \\binom{2i+1-j}{k}"

Begin["`Private`"]

Unprotect[Power];
Power[0|0., 0|0.] = 1;
Protect[Power];

RascalNumber[n_, k_, i_]:= Sum[Binomial[n-k, m]* Binomial[k, m], {m, 0, i}];
OneQBinomial[n_, m_, q_] := OneQBinomial[n, m, q] = Which[
  m > n, 0,
  n == 0 && m == 0, q,
  n >= 1 && m == 0, 1,
  True, OneQBinomial[n - 1, m, q] + OneQBinomial[n - 1, m - 1, q]
];
ColumnIdentity1[rowsNumber_, i_]:= Column[Table[RascalNumber[n, i-j, i]== Binomial[n, i-j], {n, 0, rowsNumber}, {j, 0, i}]];
ColumnIdentity2[rowsNumber_, i_]:= Column[Table[RascalNumber[n, n-i+j, i]== Binomial[n, n-i+j], {n, 0, rowsNumber}, {j, 0, i}]];
RowIdentity1[i_]:= Column[Table[RascalNumber[2i+1-j, k, i]== Binomial[2i+1-j, k], {j, 0, 2i+1}, {k, 0, 2i+1-j}]];

End[ ]

EndPackage[ ]







