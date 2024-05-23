classdef Out_Complex_Quantization_Estimation_Rewrite

  properties
    ADC_enable;
    bit;
    quan_step;
    q;
    p;
    MAX_TERM;
  end

  methods
    function obj = Out_Complex_Quantization_Estimation_Rewrite(ADC_enable, bit)

      obj.MAX_TERM = 1e5;
      MAX_TERM = obj.MAX_TERM;
      obj.ADC_enable = ADC_enable;
      obj.bit = bit;
      quan_step = 1 / (2^bit - 1);

      q = zeros(2^bit + 1, 1);
      p = zeros(2^bit    , 1);
      q(1) = -MAX_TERM;
      q(2 : 1 : 2^bit) = (1 - 2^(bit - 1) + [0 : 1 : 2^bit-2]) * quan_step;
      q(2^bit + 1) = +MAX_TERM;
      p(1 : 1 : 2^bit) = (-2^(bit - 1) + 0.5 + [0 : 1 : 2^bit - 1]) * quan_step;

      obj.quan_step = quan_step;
      obj.q = q;
      obj.p = p;
    end

    function tilde_Y = Quantization(obj, Z)
      ADC_enable = obj.ADC_enable;
      bit = obj.bit;
      if 0 == ADC_enable
        tilde_Y = Z;
      elseif 1 == ADC_enable
        [M, K] = size(Z);
        Z = reshape(Z, M * K, 1);
        Y_real = real(Z);
        Y_imag = imag(Z);
        Y_R = zeros(M * K, 1);
        Y_I = zeros(M * K, 1);

        quan_step = obj.quan_step;
        q = obj.q;
        p = obj.p;

        for bIdx = 1 : 1 : 2^bit
          Y_R = Y_R + p(bIdx) * ( (q(bIdx) <= Y_real) & (Y_real < q(bIdx + 1)) );
        end

        for bIdx = 1 : 1 : 2^bit
          Y_I = Y_I + p(bIdx) * ( (q(bIdx) <= Y_imag) & (Y_imag < q(bIdx + 1)) );
        end

        tilde_Y = Y_R + 1j*Y_I ;
        tilde_Y = reshape(tilde_Y,M,K);
      else
        throw(MException('Foo:FatalError', 'Error in Out_Complex_Quantization_Estimation.Quantization(...)'));
      end
    end

    function [hatz, hatv] = Estimation(obj, init_Y, nuw, init_Z, init_V)
      ADC_enable = obj.ADC_enable;
      MAX_TERM = obj.MAX_TERM;
      bit = obj.bit;
      quan_step = obj.quan_step;
      q = obj.q;
      p = obj.p;

      if 0 == ADC_enable
        hatv = 1 ./ (1 ./ init_V + 1 / nuw);
        hatz = hatv .* (init_Z ./ init_V + init_Y / nuw);

      elseif 1 == ADC_enable
        [M, L] = size(init_Y);
        Y = reshape(init_Y, M * L, 1);
        Z = reshape(init_Z, M * L, 1);
        V = reshape(init_V, M * L, 1);
        
        Y = [real(Y); imag(Y)];
        Z = [real(Z); imag(Z)];
        V = [real(V); real(V)];

        y_up  = Y + quan_step / 2;
        y_low = Y - quan_step / 2;
        [pos1, ~] = find(Y == max(p));
        [pos2, ~] = find(Y == min(p));
        y_up(pos1)  = MAX_TERM;
        y_low(pos2) = -MAX_TERM;

        eta1 = (y_up - Z)   ./ sqrt((nuw + V) / 2);
        eta2 = (y_low - Z)  ./ sqrt((nuw + V) / 2);

        tem1 = normpdf(eta1) - normpdf(eta2);
        tem2 = normcdf(eta1) - normcdf(eta2) + eps;
        tem3 = eta1 .* normpdf(eta1) - eta2 .* normpdf(eta2);

        z_tem = Z - V ./ sqrt(2 * (nuw + V)) .* (tem1 ./ tem2);
        v_tem = V / 2 - V.^2 ./ (2 * (nuw + V)) .* (tem3 ./ tem2 + (tem1 ./ tem2).^2);

        hatz = z_tem(1 : M * L) + 1j * z_tem(M * L + 1 : end);
        hatv = max(v_tem(1 : M * L) + v_tem(M * L + 1 : end), 1e-10);
        hatz = reshape(hatz, M, L);
        hatv = reshape(hatv, M, L);
      else
        throw(MException('Foo:FatalError', 'Error in Out_Complex_Quantization_Estimation.Estimation(...)'));
      end
      hatv = real(hatv);
    end

    function [hatz] = Estimation_mean(obj, init_Y, nuw, init_Z, init_V)
      ADC_enable = obj.ADC_enable;
      MAX_TERM = obj.MAX_TERM;
      bit = obj.bit;
      quan_step = obj.quan_step;
      q = obj.q;
      p = obj.p;

      if 0 == ADC_enable
        hatv = 1 ./ (1 ./ init_V + 1 / nuw);
        hatz = hatv .* (init_Z ./ init_V + init_Y / nuw);

      elseif 1 == ADC_enable
        [M, L] = size(init_Y);
        Y = reshape(init_Y, M * L, 1);
        Z = reshape(init_Z, M * L, 1);
        V = reshape(init_V, M * L, 1);
        
        Y = [real(Y); imag(Y)];
        Z = [real(Z); imag(Z)];
        V = [real(V); real(V)];

        y_up  = Y + quan_step / 2;
        y_low = Y - quan_step / 2;
        [pos1, ~] = find(Y == max(p));
        [pos2, ~] = find(Y == min(p));
        y_up(pos1)  = MAX_TERM;
        y_low(pos2) = -MAX_TERM;

        eta1 = (y_up - Z)   ./ sqrt((nuw + V) / 2);
        eta2 = (y_low - Z)  ./ sqrt((nuw + V) / 2);

        tem1 = normpdf(eta1) - normpdf(eta2);
        tem2 = normcdf(eta1) - normcdf(eta2) + eps;
        tem3 = eta1 .* normpdf(eta1) - eta2 .* normpdf(eta2);

        z_tem = Z - V ./ sqrt(2 * (nuw + V)) .* (tem1 ./ tem2);

        hatz = z_tem(1 : M * L) + 1j * z_tem(M * L + 1 : end);
        hatz = reshape(hatz, M, L);
      else
        throw(MException('Foo:FatalError', 'Error in Out_Complex_Quantization_Estimation.Estimation(...)'));
      end
    end

    function [hatv] = Estimation_variance(obj, init_Y, nuw, init_Z, init_V)
      ADC_enable = obj.ADC_enable;
      MAX_TERM = obj.MAX_TERM;
      bit = obj.bit;
      quan_step = obj.quan_step;
      q = obj.q;
      p = obj.p;

      if 0 == ADC_enable
        hatv = 1 ./ (1 ./ init_V + 1 / nuw);
        hatz = hatv .* (init_Z ./ init_V + init_Y / nuw);

      elseif 1 == ADC_enable
        [M, L] = size(init_Y);
        Y = reshape(init_Y, M * L, 1);
        Z = reshape(init_Z, M * L, 1);
        V = reshape(init_V, M * L, 1);
        
        Y = [real(Y); imag(Y)];
        Z = [real(Z); imag(Z)];
        V = [real(V); real(V)];

        y_up  = Y + quan_step / 2;
        y_low = Y - quan_step / 2;
        [pos1, ~] = find(Y == max(p));
        [pos2, ~] = find(Y == min(p));
        y_up(pos1)  = MAX_TERM;
        y_low(pos2) = -MAX_TERM;

        eta1 = (y_up - Z)   ./ sqrt((nuw + V) / 2);
        eta2 = (y_low - Z)  ./ sqrt((nuw + V) / 2);

        tem1 = normpdf(eta1) - normpdf(eta2);
        tem2 = normcdf(eta1) - normcdf(eta2) + eps;
        tem3 = eta1 .* normpdf(eta1) - eta2 .* normpdf(eta2);

        v_tem = V / 2 - V.^2 ./ (2 * (nuw + V)) .* (tem3 ./ tem2 + (tem1 ./ tem2).^2);

        hatv = max(v_tem(1 : M * L) + v_tem(M * L + 1 : end), 1e-10);
        hatv = reshape(hatv, M, L);
      else
        throw(MException('Foo:FatalError', 'Error in Out_Complex_Quantization_Estimation.Estimation(...)'));
      end
      hatv = real(hatv);
    end

    function vz_sub = MSE_SE(obj, Tz, v1_plus, nuw)
      ADC_enable = obj.ADC_enable;
      MAX_TERM = obj.MAX_TERM;
      
      if 0 == ADC_enable
        vz_sub = (nuw * v1_plus) / (v1_plus + nuw);
      elseif 1 == ADC_enable
        quan_step = obj.quan_step;
        p = obj.p;
        Y_Up  = p + quan_step / 2;
        Y_Up(end) = Inf;
        Y_Low = p - quan_step / 2;
        Y_Low(1) = -Inf;
        alpha = zeros(length(p), 1);
        Gaussian = @(x, m, v) 1 ./ sqrt(2 * pi * v) .* exp(-(x - m).^(2) ./ (2 * v));
        for index = 1 : length(p)
          alpha(index) = integral(@(u) ...
            ( ...
              normpdf( ...
                ( ...
                  sqrt(2) * Y_Up(index) - sqrt(Tz - v1_plus) * u ...
                ) / sqrt(v1_plus + nuw) ...
                 ) - normpdf( ...
                   ( ...
                     sqrt(2) * Y_Low(index) - sqrt(Tz - v1_plus) * u ...
                   ) / sqrt(v1_plus + nuw) ...
                 ) ...
            ).^(2) ./ ( ...
              normcdf( ...
                ( ...
                  sqrt(2) * Y_Up(index) - sqrt(Tz - v1_plus) * u ...
                ) / sqrt(v1_plus + nuw) ...
              ) - normcdf( ...
                ( ...
                  sqrt(2) * Y_Low(index) - sqrt(Tz - v1_plus) * u ...
                ) / sqrt(v1_plus + nuw) ...
              ) + eps ...
            ) .* Gaussian(u, 0, 1), -Inf, Inf ...
          );
        end
        vz_sub = v1_plus - (v1_plus)^(2) / (v1_plus + nuw) * sum(alpha);
      end
    end

  end
end