(* ::Package:: *)

BeginPackage["PolynomialIdentitiesInvolvingRascalTriangle`"]

RascalNumber::usage="Gives generalized rascal triangle via sum of product \\rascalNumber{n}{k}{i} = \\sum_{m=0}^{i} \\binom{n-k}{m} \\binom{k}{m}"
OneQBinomial::usage="Gives (1,q)-Binomial coefficient, see https://oeis.org/A096940"
ColumnIdentity1::usage="Validates column identity \\rascalNumber{n}{k}{i} = \\binom{n}{k}, \\quad 0 \\leq k \\leq i"
ColumnIdentity2::usage="Validates column identity \\rascalNumber{n}{i-j}{i}        &= \\binom{n}{i-j}"
RowIdentity1::usage="Validates row identity \\rascalNumber{2i+1-j}{k}{i}             &= \\binom{2i+1-j}{k}"
RowIdentity2::usage="Validates row identity \\rascalNumber{t-j}{k}{t-i-1}            &= \\binom{t-j}{k}"
RowIdentity3::usage="Validates row identity \\rascalNumber{2i+1-j}{j}{i}         &= \\binom{2i+1-j}{j}"
RowIdentity4::usage="Validates row identity \\rascalNumber{t-j}{t-2j}{t-i-1}     &= \\binom{t-j}{t-2j}, \\quad t \\geq 2i+1"

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
RowIdentity2[t_, i_]:= Column[Table[RascalNumber[t-j, k, t-i-1]== Binomial[t-j, k], {j, 0, t}, {k, 0, t-j}]];
RowIdentity3[i_]:= Column[Table[RascalNumber[2i+1-j, j, i]== Binomial[2i+1-j, j], {j, 0, i}]];
RowIdentity4[t_, i_]:= Column[Table[RascalNumber[t-j, t-2j, t-i-1]== Binomial[t-j, t-2j], {j, 0, t-i-1}]];

End[ ]

EndPackage[ ]










