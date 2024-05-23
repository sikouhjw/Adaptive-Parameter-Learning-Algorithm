function W = Noise_Complex_AWGN(M, K, nuw)
	W = sqrt(nuw / 2) * (randn(M, K) + 1j * randn(M, K));
end