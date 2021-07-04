# Zippel
Implementation of a Zippel's sparse polynomial interpolation algorithm [1].

The project involves study and analysis of sparse multivariate polynomial interpolation and GCD computation algorithms. The basic idea is to interpolate the GCD modulo several primes and recover the integer coefficients of the GCD using Chinese Remainder Algorithm. Several problems and issues arise while using the above approach which are analyzed and solved in the paper we studied. The polynomial to be interpolated is treated as a black box function. The first approach given by Zippel is a probabilistic interpolation algorithm which uses a combination of dense and sparse interpolation to construct the GCD of polynomials. The proof of correctness of this approach is given by Schwartz–Zippel lemma [1]. Ben-Or/Tiwari gives a deterministic interpolation algorithm based on Berlekamp-Massey algorithm from coding theory [2]. A combination of Newton interpolation and univariate version of Ben-Or/Tiwari given by Kaltofen et. al. can be used to get dense univariate result efficiently and it can be further utilized to interpolate multivariate polynomials using Zippel’s Approach [3]. A parallel algorithm which uses Kronecker substitution along with approaches of Zippel and Ben-Or/Tiwari gives extremely fast benchmark results [4].Performance of the GCD algorithms can be improved in two situations: 1) Modular GCD algorithm generates a monic gcd (as GCD is unique upto multiplication by a unit) which needs to be scaled appropriately to find the correct modular image which is an expensive task. 2) When GCD is non-monic, it is difficult to reconstruct GCD from modular images (Normalization Problem [5]).

[1] Zippel, Richard. ”Probabilistic algorithms for sparse polynomials.” Symbolic and algebraic computation (1979): 216-226.

[2] Ben-Or, Michael, and Prasoon Tiwari. ”A deterministic algorithm for sparse multivariate polynomial interpolation.” Proceedings of the twentieth annual ACM symposium on Theory of computing. ACM, 1988.

[3] Kaltofen, Erich, Wen-shin Lee, and Austin A. Lobo. ”Early termination in Ben-Or/Tiwari sparse interpolation and a hybrid of Zippel’s algorithm.” Proceedings of the 2000 international symposium on Symbolic and algebraic computation. ACM, 2000.

[4] Hu, Jiaxiong, and Michael Monagan. ”A Fast Parallel Sparse Polynomial GCD Algorithm.” Proceedings of the ACM on International Symposium on Symbolic and Algebraic Computation. ACM, 2016.

[5] de Kleine, Jennifer, Michael Monagan, and Allan Wittkopf. ”Algorithms for the non-monic case of the sparse modular GCD algorithm.” Proceedings of the 2005 international symposium on Symbolic and algebraic computation. ACM, 2005.
