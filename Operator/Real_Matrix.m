function [H] = Real_Matrix(H)
  [M, N] = size(H);
  if N == 1
    H = [real(H); imag(H)];
  else
    H = [real(H) -imag(H); imag(H) real(H)];
  end
end