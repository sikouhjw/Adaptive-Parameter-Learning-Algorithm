classdef In_Complex_Constellation_Estimation 

  properties
    mode_size;
    xo;
    normal_scal;
    Tx;
  end

  methods
    function obj = In_Complex_Constellation_Estimation(mode_size)
      obj.mode_size = mode_size;
      obj.xo = obj.Get_Constellation();
      obj.Tx = 1;
    end

    function X = Generation(obj, N, L)
      mode_size = obj.mode_size;
      informationBit = randi([0, 1], N * L * log2(mode_size), 1);
      informationSym = qammod(informationBit, mode_size, 'InputType', 'bit', 'UnitAveragePower', true);
      X = reshape(informationSym, N, L);
    end

    function [hat_x, hat_v] = Estimation(obj, x_mean, x_var, Init_flag, Init_x)
      flag = false;
      x = 0;
      switch nargin
        case 3
          flag = false;
        case 4
          flag = Init_flag;
          x = zeros(size(x_mean));
        case 5
          flag = Init_flag;
          x = Init_x;
        otherwise
          throw(MException('Foo:FatalError', 'Error in In_Complex_Constellation_Estimation.Estimation(...)'));
      end
      if ~flag
        Allsymbol = obj.xo;
        [N, L] = size(x_mean);
        if sum(size(x_var)) == 2
          x_var = x_var;
        else
          x_var = reshape(x_var, N * L, 1);
        end
        x_mean = reshape(x_mean, N * L, 1);
        
        log_posterior = bsxfun(@times, -1 ./ x_var, abs(bsxfun(@minus, x_mean, Allsymbol)).^(2));
        log_posterior = bsxfun(@minus, log_posterior, max(log_posterior, [], 2));
        posterior = exp(log_posterior); 
        posterior = bsxfun(@rdivide, posterior, sum(posterior, 2));
        
        hat_x = sum(bsxfun(@times, posterior, Allsymbol), 2);
        hat_v = sum(posterior .* abs(bsxfun(@minus, hat_x, Allsymbol)).^(2), 2);
        
        hat_x = reshape(hat_x, N, L);
        hat_v = reshape(hat_v, N, L);
      else
        hat_x = x;
        hat_v = ones(size(x));
      end
      hat_v = real(hat_v);
    end

    function [hat_x, hat_v] = Estimation_Real(obj, x_mean, x_var, Init_flag, Init_x)

      Allsymbol_Real = real(obj.xo);
      Allsymbol_Imag = imag(obj.xo);
      Allsymbol_Real = unique(Allsymbol_Real);
      Allsymbol_Imag = unique(Allsymbol_Imag);

      [N, L] = size(x_mean);
      if sum(size(x_var)) == 2
        x_var = x_var;
      else
        x_var = reshape(x_var, N * L, 1);
      end
      x_mean = reshape(x_mean, N * L, 1);
      
      log_posterior_Real = bsxfun(@times, -1 ./ ( 2 .* x_var(1 : N ./ 2) ), abs(bsxfun(@minus, x_mean(1 : N ./ 2), Allsymbol_Real)).^(2));
      log_posterior_Real = bsxfun(@minus, log_posterior_Real, max(log_posterior_Real, [], 2));
      posterior_Real = exp(log_posterior_Real); 
      posterior_Real = bsxfun(@rdivide, posterior_Real, sum(posterior_Real, 2));

      hat_x_Real = sum(bsxfun(@times, posterior_Real, Allsymbol_Real), 2);
      hat_v_Real = sum(posterior_Real .* abs(bsxfun(@minus, hat_x_Real, Allsymbol_Real)).^(2), 2);

      log_posterior_Imag = bsxfun(@times, -1 ./ ( 2 .* x_var(N ./ 2 + 1 : N) ), abs(bsxfun(@minus, x_mean(N ./ 2 + 1 : N), Allsymbol_Imag)).^(2));
      log_posterior_Imag = bsxfun(@minus, log_posterior_Imag, max(log_posterior_Imag, [], 2));
      posterior_Imag = exp(log_posterior_Imag); 
      posterior_Imag = bsxfun(@rdivide, posterior_Imag, sum(posterior_Imag, 2));

      hat_x_Imag = sum(bsxfun(@times, posterior_Imag, Allsymbol_Imag), 2);
      hat_v_Imag = sum(posterior_Imag .* abs(bsxfun(@minus, hat_x_Imag, Allsymbol_Imag)).^(2), 2);

      hat_x = [hat_x_Real; hat_x_Imag];
      hat_v = [hat_v_Real; hat_v_Imag];

      hat_x = reshape(hat_x, N, L);
      hat_v = reshape(hat_v, N, L);

      hat_v = real(hat_v);
    end

    function xo = Get_Constellation(obj)
      mode_size = obj.mode_size;
      xo = [qammod([0 : mode_size - 1]', mode_size, 'UnitAveragePower', true)]';
    end

    function hat_vx_plus = MSE_SE(obj, v0_sub)
      v0_sub_inv = 1 / v0_sub;
      hat_vx_plus = 1 - integral(@(u) ...
        tanh(v0_sub_inv + sqrt(v0_sub_inv) .* u) .* Real_Gaussian(u, 0, 1), -Inf, Inf ...
      );
    end
  end
end

