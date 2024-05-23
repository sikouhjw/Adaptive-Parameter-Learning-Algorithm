function [MSE_error, HIT_error, psi_error, term_scale_error] = Complex_OTFS_GAMP_EM(Input, obj)
  addpath('../Operator');

  iterations = Input.iterations;
  damping = Input.damping;
  nuw = Input.nuw;
  mode_size = Input.mode_size;
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

  hat_v_i = ones(N, 1);
  hat_m_i = zeros(N, 1);
  hat_s_a = zeros(M, 1);

  sqrH = real(abs(H).^(2));
  sqrHt = sqrH';
  Ht = H';

  Z_a = zeros(M, 1);
  V_a = ones(M, 1);
  V_a_old = ones(M, 1);
  Z_a_old = zeros(M, 1);
  hat_s_a_old = zeros(M, 1);
  hat_tau_a_old = zeros(M, 1);

  outFileName = sprintf(mfilename);

  MSE_error = zeros(iterations, 1);
  HIT_error = zeros(iterations, 1);
  MSE_old = 1e5;
  HIT_old = 0;
  MSE_fake_old = 1e5;
  hat_v_i_old = 1e5;

  term_scale = 0;
  term_scale_error = zeros(iterations, 1);

  psi_error = zeros(iterations, 1);
  psi = Input.psi;
  psi_old = psi;

  for ii = 1 : iterations

    [tilde_z_a, tilde_v_a] = Out.Estimation(y, psi, Z_a, V_a);

    hat_s_a = (tilde_z_a - Z_a) ./ V_a;
    hat_tau_a = (V_a - tilde_v_a) ./ abs(V_a).^(2);

    tilde_v_i = 1 ./ (sqrHt * hat_tau_a);
    tilde_m_i = hat_m_i + tilde_v_i .* (Ht * hat_s_a);

    hat_m_i = ( ...
      rho * Complex_Gaussian( 0, tilde_m_i, tilde_v_i + eta ) .* tilde_m_i .* eta ./ ( tilde_v_i + eta ) ...
    ) ./ ( ...
      rho * Complex_Gaussian( 0, tilde_m_i, tilde_v_i + eta ) + ( 1 - rho ) * Complex_Gaussian( 0, tilde_m_i, tilde_v_i ) ...
    );
    hat_v_i = ( ...
      rho * ( ...
        tilde_v_i .* eta ./ ( ...
          tilde_v_i + eta ...
        ) + abs( ...
          tilde_m_i .* eta ./ (...
            tilde_v_i + eta ...
          ) ...
        ).^2 ...
      ) .* Complex_Gaussian( 0, tilde_m_i, tilde_v_i + eta ) ...
    ) ./ ( ...
      rho * Complex_Gaussian( 0, tilde_m_i, tilde_v_i + eta ) + ( 1 - rho ) * Complex_Gaussian( 0, tilde_m_i, tilde_v_i ) ...
    ) - abs( hat_m_i ).^2;

    if ii == 1
      hat_m_i_old = hat_m_i;
      hat_v_i_old = hat_v_i;
    end

    [hat_m_i, hat_m_i_old] = Damping(hat_m_i, hat_m_i_old, damping);
    [hat_v_i, hat_v_i_old] = Damping(hat_v_i, hat_v_i_old, damping);

    hat_v_i_mean = mean( hat_v_i );
    hat_v_i_old_mean = mean( hat_v_i_old );

    MSE = norm(hat_m_i - x).^2 / norm(x).^2;
    MSE_fake = norm( y - Out.Quantization(H * hat_m_i) ).^2;

    [hat_m_i_number, hat_m_i_index]  = maxk(hat_m_i, p, 'ComparisonMethod','abs');
    HIT = length( intersect( hat_m_i_index, find(x) ) ) / p;

    if (sum(hat_v_i == 0))
      MSE_error(max(1, ii) : iterations) = MSE;
      HIT_error(max(1, ii) : iterations) = HIT;
      psi_error(max(1, ii) : iterations) = psi;
      break;
    end

    if (sum(isnan(hat_v_i)) > 0 || MSE > MSE_old )
      MSE_error(max(1, ii - 1) : iterations) = MSE_old;
      HIT_error(max(1, ii - 1) : iterations) = HIT_old;
      psi_error(max(1, ii - 1) : iterations) = psi_old;
      break;
    else
      MSE = norm(hat_m_i - x).^2 / norm(x).^2;
      HIT = length( intersect( hat_m_i_index, find(x) ) ) / p;
    end

    hat_v_i_old = hat_v_i;
    MSE_fake_old = MSE_fake;

    V_a = sqrH * hat_v_i;
    Z_a = H * hat_m_i - hat_s_a .* V_a;

    if Input.estimation_enable == 1
      if psi <= 1e-1
        if abs(MSE - MSE_old)^2 / MSE^2 <= 1e-3 || MSE > MSE_old
          [psi, term_scale] = Estimate_EM_ELBO(Input, obj, Z_a, V_a);
        end
      else
        if abs(MSE - MSE_old)^2 / MSE^2 <= 1e-3 || MSE > MSE_old
          [psi, term_scale] = Estimate_EM_Bethe(Input, obj, Z_a, V_a);
        end
      end
    elseif Input.estimation_enable == 2
      if abs(MSE - MSE_old)^2 / MSE^2 <= 1e-3 || MSE > MSE_old
        [psi, term_scale] = Estimate_EM_ELBO(Input, obj, Z_a, V_a);
      end
    elseif Input.estimation_enable == 3
      if abs(MSE - MSE_old)^2 / MSE^2 <= 1e-3 || MSE > MSE_old
        [psi, term_scale] = Estimate_EM_Bethe(Input, obj, Z_a, V_a);
      end
    else
      if abs(MSE - MSE_old)^2 / MSE^2 <= 1e-3 || MSE > MSE_old
        psi = Estimate_EM_AWGN(Input, obj, tilde_z_a, tilde_v_a);
        term_scale = 0;
      end
    end

    if (psi < 0) || (psi > 25)
      psi = psi_old;
    else
      psi_old = psi;
    end

    psi_error(ii, 1) = psi;
    term_scale_error(ii, 1) = term_scale;

    if (sum(isnan(hat_v_i)) == 0)
      MSE_old = MSE;
      HIT_old = HIT;
      MSE_error(ii, 1) = MSE;
      HIT_error(ii, 1) = HIT;
    end

  end
end