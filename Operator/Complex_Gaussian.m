function ret = Complex_Gaussian(x, m, v)
  ret = 1 ./ (pi * v) .* exp(-abs(x - m).^2 ./ v);
  ret = real(ret);
end