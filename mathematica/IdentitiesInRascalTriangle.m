(* ::Package:: *)

BeginPackage["PolynomialIdentitiesInvolvingRascalTriangle`"]

RascalNumber::usage="Gives generalized rascal triangle via sum of product \\rascalNumber{n}{k}{i} = \\sum_{m=0}^{i} \\binom{n-k}{m} \\binom{k}{m}"
OneQBinomial::usage="Gives (1,q)-Binomial coefficient, see https://oeis.org/A096940"
PQBinomial::usage="Gives (p,q)-Binomial coefficient, see https://oeis.org/A096940"
ColumnIdentity1::usage="Validates column identity \\rascalNumber{n}{k}{i} = \\binom{n}{k}, \\quad 0 \\leq k \\leq i"
ColumnIdentity2::usage="Validates column identity \\rascalNumber{n}{i-j}{i}        &= \\binom{n}{i-j}"
RowIdentity1::usage="Validates row identity \\rascalNumber{2i+1-j}{k}{i}             &= \\binom{2i+1-j}{k}"
PrintRowIdentity1::usage="Prints row identity \\rascalNumber{2i+1-j}{k}{i}             &= \\binom{2i+1-j}{k}"
RowIdentity2::usage="Validates row identity \\rascalNumber{t-j}{k}{t-i-1}            &= \\binom{t-j}{k}"
RowIdentity3::usage="Validates row identity \\rascalNumber{2i+1-j}{j}{i}         &= \\binom{2i+1-j}{j}"
RowIdentity4::usage="Validates row identity \\rascalNumber{t-j}{t-2j}{t-i-1}     &= \\binom{t-j}{t-2j}, \\quad t \\geq 2i+1"
RowColumnDifferenceIdentity1::usage="Validates row identity \\binom{2i+2+j}{i+1} - \\rascalNumber{2i+2+j}{i+1}{i} &= \\binom{i+1+j}{i+1} \\\\"
RowColumnDifferenceIdentity2::usage="Validates row identity \\binom{j+i}{j} - \\rascalNumber{j+i}{j}{i-1} &= \\binom{j}{j-i}"
OneQPascalIdentity1::usage="Validates row identity \\binom{2i+3+j}{i+2} - \\rascalNumber{2i+3+j}{i+2}{i} &= \\oneQBinomial{i+2+j}{i+2}{i+2}"
OneQPascalIdentity2::usage="Validates row identity \\binom{2t-1+j}{t} - \\rascalNumber{2t-1+j}{t}{t-2} &= \\oneQBinomial{t+j}{t}{t}"
BinomialCoefficientAsPolynomial::usage="Gives binomial coefficient in form of polynomial"
BinomialRascalDifferenceColumn::usage="Gives Binomial minus rascal column for i, k"
VandermondeIdentity::usage="Gives Vandermonde identity for Binom(a+b, r)"
DifferenceBinomialRascal::usage="Gives the difference binom(n,k)-rascal(n,k,i)"
RowSumsPowerOf2Identity::usage="Conjecture from iterated rascal triangles."
A::usage= "A[n, k] returns the real coefficient A of non-negative integers n, k such that n <= k. 
See https://kolosovpetro.github.io/pdf/AStudyOnDynamicEquations.pdf."
RascalTriangleDiagonal::usage= "Gives iterated rascal triangle R(i) n-th diagonal up to specific term."
ColumnDifferenceIdentity::usage= "Binom(n+2i, i)-R(n+2i, i) difference"
RascalRowSum::usage= "Gives n-th row sum of iterated rascal triangle"
RascalColumnSum::usage= "Gives sum of k-th column up to fixed value."
RascalInternalSum::usage= "Gives sum of rows k=2 to n of rascal number except 1."
RascalInternalSum2::usage= "Gives sum of rows k=2 to n of rascal number except 1."
RascalInternalSum3::usage= "Gives sum of rows k=2 to n of rascal number except 1."

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
PQBinomial[n_, m_, p_, q_] := OneQBinomial[n, m, p, q] = Which[
  m > n, 0,
  n == 0 && m == 0, q,
  n >= 1 && m == 0, p,
  True, PQBinomial[n - 1, m, p, q] + PQBinomial[n - 1, m - 1, p, q]
];
ColumnIdentity1[rowsNumber_, i_]:= Column[Table[RascalNumber[n, i-j, i]== Binomial[n, i-j], {n, 0, rowsNumber}, {j, 0, i}]];
ColumnIdentity2[rowsNumber_, i_]:= Column[Table[RascalNumber[n, n-i+j, i]== Binomial[n, n-i+j], {n, 0, rowsNumber}, {j, 0, i}]];
RowIdentity1[i_]:= Column[Table[RascalNumber[2i+1-j, k, i]== Binomial[2i+1-j, k], {j, 0, 2i+1}, {k, 0, 2i+1-j}]];
PrintRowIdentity1[i_]:= Column[Table[RascalNumber[2i+1-j, k, i], {j, 0, 2i+1}, {k, 0, 2i+1-j}]];
RowIdentity2[t_, i_]:= Column[Table[RascalNumber[t-j, k, t-i-1]== Binomial[t-j, k], {j, 0, t}, {k, 0, t-j}]];
RowIdentity3[i_]:= Column[Table[RascalNumber[2i+1-j, j, i]== Binomial[2i+1-j, j], {j, 0, i}]];
RowIdentity4[t_, i_]:= Column[Table[RascalNumber[t-j, t-2j, t-i-1]== Binomial[t-j, t-2j], {j, 0, t-i-1}]];
RowColumnDifferenceIdentity1[i_, rowNumber_]:= Column[Table[Binomial[n+2i, i] - RascalNumber[n+2i, i, i-1]== Binomial[n+i, i], {n, 0, rowNumber}]];
RowColumnDifferenceIdentity2[i_, rowNumber_]:= Column[Table[Binomial[j+i, j] - RascalNumber[j+i, j, i-1]== Binomial[j, j-i], {j, 0, rowNumber}]];
OneQPascalIdentity1[i_, rowNumber_]:= Column[Table[Binomial[2i+3+j, i+2] - RascalNumber[2i+3+j, i+2, i]== OneQBinomial[i+2+j, i+2, i+2], {j, 0, rowNumber}]];
OneQPascalIdentity2[i_, rowNumber_]:= Column[Table[Binomial[2i-1+j, i] - RascalNumber[2i-1+j, i, i-2]== OneQBinomial[i+j, i, i], {j, 0, rowNumber}]];
BinomialCoefficientAsPolynomial[n_, k_]:= Product[(n-i+1)/i, {i, 1, k}];
BinomialRascalDifferenceColumn[i_, k_, size_] := Table[Binomial[n,k]-RascalNumber[n,k, i],{n, 2i+1, 2i+1 + size}];
VandermondeIdentity[a_, b_, r_, upperLimit_] := Sum[Binomial[a, m] * Binomial[b, r-m], {m, 0, upperLimit}];
DifferenceBinomialRascal[n_, k_, i_]:=Sum[Binomial[n-k, m]*Binomial[k, k-m], {m, i+1, k}];
RowSumsPowerOf2Identity[i_]:= Sum[RascalNumber[4*i+3, d, i], {d, 0, 4i+3}];

A[n_, k_] := 0;
A[n_, k_] := (2k + 1) * Binomial[2k, k] * Sum[A[n, j] * Binomial[j, 2k + 1] * (-1)^(j - 1) / (j - k) * BernoulliB[2j - 2k], {j, 2k + 1, n}] /; 0 <= k < n;
A[n_, k_] := (2n + 1) * Binomial[2n, n] /; k == n;

RascalTriangleDiagonal[n_, i_, count_]:= Table[RascalNumber[n+d, d, i], {d, 0, count}];
ColumnDifferenceIdentity[i_, count_]:=Table[Binomial[n+2i,i]-RascalNumber[n+2i,i, i-1], {n, 0, count}];
RascalRowSum[n_, i_]:=Sum[RascalNumber[n, k, i], {k, 0, n}];
RascalColumnSum[k_, i_, max_]:= Sum[RascalNumber[n, k, i], {n, 0, max}];
RascalInternalSum[n_, i_]:=Sum[Sum[RascalNumber[k, t, i], {t, 1, k-1}], {k, 2, n}];
RascalInternalSum2[n_, i_]:=Sum[Sum[RascalNumber[k, t+1, i], {t, 0, k-2}], {k, 2, n}];
RascalInternalSum3[n_, i_]:=Sum[Sum[RascalNumber[k, k-2-t-1, i], {t, 0, k+2}], {k, 2, n}];

End[ ]

EndPackage[ ]






