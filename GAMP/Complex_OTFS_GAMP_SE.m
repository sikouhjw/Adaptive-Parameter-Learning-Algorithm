function MSE_List = Complex_OTFS_GAMP_SE(Input, obj)
  addpath('../Operator');

  iterations = Input.iterations;
  damping = Input.damping;
  nuw = Input.nuw;
  p   = Input.p;
  Q   = obj.Q;
  In = Input.In;
  Out = Input.Out;
  
  x = obj.x;
  H = obj.H;
  y = obj.y;
  eta = obj.eta;
  [M, N] = size(H);

  rho = p / Q;

  Tz = 1;

  ALPHA = M / N;

  V_a = 1 - 1e-1;

  MSE_List = zeros(iterations, 1);

  for ii = 1 : iterations

    tilde_v_a = Out.MSE_SE(Tz, V_a, nuw);
    tilde_v_a = real(tilde_v_a);

    tau_a = (V_a - tilde_v_a) / (V_a^(2));

    tilde_v_i = 1 ./ (M * tau_a);

    hat_v = 1 ./ Q - rho.^2 .* integral(@(z) ...
      Real_Gaussian(z, 0, tilde_v_i + eta).^2 .* z.^2 ./ ( ...
        ( 1 + tilde_v_i ./ eta ).^2 .* ( ...
          ( 1 - rho ) .* Real_Gaussian(z, 0, tilde_v_i) + rho .* Real_Gaussian(z, 0, tilde_v_i + eta) ...
        ) ...
      ), -5, 5 ...
    );
    hat_v = real(hat_v);

    if ii == 1
      hat_v_old = hat_v;
    end

    [hat_v, hat_v_old] = Damping(hat_v, hat_v_old, damping);

    MSE = Q .* real(hat_v);
    MSE_List(ii, 1) = MSE;

    V_a = N .* hat_v;
  end
end