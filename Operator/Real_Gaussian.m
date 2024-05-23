function ret = Real_Gaussian(x, m, v)
	ret = 1 ./ sqrt(2 * pi * v) .* exp(-(x - m).^2 ./ (2 * v));
	ret = real(ret);
end