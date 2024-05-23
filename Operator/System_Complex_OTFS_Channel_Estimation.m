function obj = System_Complex_OTFS_Channel_Estimation(Input)
  M = Input.M;
  N = Input.N;
  k_max = Input.k_max;
  l_max = Input.l_max;
  p     = Input.p;
  In  = Input.In;
  Out = Input.Out;

  x = In.Generation(M*N, 1);

  Q = (l_max + 1) .* (2 .* k_max + 1);
  l = [0 : l_max];
  l = repmat(l, 1, (2 .* k_max + 1))';

  k = [-k_max : k_max];
  k = repmat(k, (l_max + 1), 1);
  k = reshape(k, Q, 1);

  h_sparse = randperm(Q, p);
  h_sparse_zero = zeros(Q, 1);
  h_sparse_zero( h_sparse ) = 1;

  eta   = 1;
  eta   = eta / p;

  h = normrnd(0, sqrt(eta / 2), Q, 1) + j * normrnd(0, sqrt(eta / 2), Q, 1);
  h = h .* h_sparse_zero;

  D = zeros(M .* N, Q);

  Delta = exp( 1j .* 2 .* pi ./ (M .* N) ).^([1 : M*N]' - 1);

  for i = 1 : Q
    D(:, i) = circshift(eye(M .* N), l(i), 1) * diag( Delta.^(k(i)) ) * kron( dftmtx(N)' ./ sqrt(N), eye(M) ) * x;
  end

  z = D * h;

  obj.x = h;
  obj.H = D;
  obj.h = h;
  obj.eta = eta;
  obj.l = l;
  obj.k = k;
  obj.Q = Q;
  obj.z = z;
end