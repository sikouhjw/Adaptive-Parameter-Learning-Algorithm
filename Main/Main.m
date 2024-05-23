clc;
clear;

addpath('../Input');
addpath('../Output');
addpath('../Operator');
addpath('../GAMP');

%% Inputs Setting
M     = 32;
N     = 8;
k_max = 7;
l_max = 10;

damping               = 0.8;
iterations            = 30;
ADC_enable            = 1;
mode_size             = 4;% Constellation point modulation, e.g. QPSK = 4
estimation_enable     = 2;% Hyperparameter estimation algorithm, 1 = hybrid, 2 = proposed, 3 = SR-TSP15, 4 = VS-TSP13
New_data              = 1;% Generate new `obj' data

PP                    = 30;% path number
SNRR                  = 40;% SNR range
BITT                  = 3;% Number of bits
test_num              = 1e2;

GAMP_enable           = 1;
GAMP_SE_enable        = 1;
GAMP_EM_enable        = 1;

SNR_figure_enable     = 0;% The horizontal coordinate is SNR
One_figure_enable     = 1;
Averaging_SNR_enable  = 1;
MSE_psi_enable        = 1;% The noise estimation performance is represented by MSE

iter = 1 : iterations;

GAMP_p_MSE       = [];
GAMP_EM_p_MSE    = [];
GAMP_SE_p_MSE    = [];
GAMP_EM_p_psi    = [];

for pp = PP
  p = pp;

  for bitt = BITT

    GAMP_SNR_last_MSE           = [];
    GAMP_HIT_SNR_iter           = [];
    MSE_m_SNR_iter              = [];
    MSE_v_SNR_iter              = [];
    KL_SNR_iter                 = [];
    GAMP_SE_SNR_last_MSE        = [];
    GAMP_EM_SNR_last_MSE        = [];
    GAMP_EM_HIT_SNR_iter        = [];
    GAMP_EM_SNR_last_psi        = [];
    GAMP_EM_SNR_last_term_scale = [];
    GAMP_EM_SNR_last_MSE_psi    = [];

    bit         = bitt  ;

    %% Load Inputs
    Input.M     = M;
    Input.N     = N;
    Input.k_max = k_max;
    Input.l_max = l_max;
    Input.p     = p;

    Input.damping           = damping;
    Input.iterations        = iterations;
    Input.mode_size         = mode_size;
    Input.bit               = bit;
    Input.estimation_enable = estimation_enable;

    Input.In  = In_Complex_Constellation_Estimation(mode_size);
    Input.Out = Out_Complex_Quantization_Estimation_Rewrite(ADC_enable, bit);

    tic;
    if New_data
      Parfor_Progress(test_num * length(SNRR) * 2 + test_num);
      parfor kk = 1 : test_num
        obj{kk} = System_Complex_OTFS_Channel_Estimation(Input);
        Parfor_Progress;
      end
    else
      Parfor_Progress(test_num * length(SNRR) * 2);
      str = ["Main-",num2str(BITT),'bit.mat'];
      load(join(str,''));
    end
%     outFileName = sprintf(mfilename);
%     save([outFileName, '-', num2str(BITT), 'bit.mat'], 'obj');

    for snrr = SNRR
      snr         = snrr  ;
      Input.nuw   = 10.^(-snr / 10);

      parfor kk = 1 : test_num
        w{kk}     = Noise_Complex_AWGN(M.*N, 1, Input.nuw);
        y{kk}     = Input.Out.Quantization(obj{kk}.z + w{kk});
        obj{kk}.y = y{kk};
        Parfor_Progress;
      end

      quan_w    = Input.Out.Quantization(w{1});
      Input.psi = mean( abs( quan_w ).^2 );

      parfor kk = 1 : test_num
        if GAMP_enable
          [GAMP_MSE(:, kk), ~, MSE_m(:, kk), MSE_v(:, kk), KL(:, kk)] = Complex_OTFS_GAMP(Input, obj{kk});
        end

        if GAMP_SE_enable
          [GAMP_SE_MSE(:, kk)] = Complex_OTFS_GAMP_SE(Input, obj{kk});
        end

        if GAMP_EM_enable
          [GAMP_EM_MSE(:, kk), GAMP_EM_HIT(:, kk), GAMP_EM_psi(:, kk), GAMP_EM_term_scale(:, kk)] = Complex_OTFS_GAMP_EM(Input, obj{kk});
        end
        Parfor_Progress;
      end

      if GAMP_enable
        GAMP_SNR_original_MSE(:, :, snr + 11) = GAMP_MSE;
        MSE_m_SNR(:, :, snr + 11)             = MSE_m;
        MSE_v_SNR(:, :, snr + 11)             = MSE_v;
        KL_SNR(:, :, snr + 11)                = KL;
      end

      if GAMP_SE_enable
        GAMP_SE_SNR_original_MSE(:, :, snr + 11)  = GAMP_SE_MSE;
      end

      if GAMP_EM_enable
        GAMP_EM_SNR_original_MSE(:, :, snr + 11)        = GAMP_EM_MSE;
        GAMP_EM_SNR_original_psi(:, :, snr + 11)        = GAMP_EM_psi;
        GAMP_EM_SNR_original_term_scale(:, :, snr + 11) = GAMP_EM_term_scale;
      end
    end
    Parfor_Progress(0);
    toc;

    if GAMP_enable
      GAMP_SNR_mean_MSE = mean(GAMP_SNR_original_MSE, 2);
      MSE_m_iter        = mean(MSE_m_SNR, 2);
      MSE_v_iter        = mean(MSE_v_SNR, 2);
      KL_iter           = mean(KL_SNR, 2);
    end

    if GAMP_SE_enable
      GAMP_SE_SNR_mean_MSE  = mean(GAMP_SE_SNR_original_MSE, 2);
    end

    if GAMP_EM_enable
      GAMP_EM_SNR_mean_MSE        = mean(GAMP_EM_SNR_original_MSE, 2);
      GAMP_EM_SNR_mean_psi        = mean(GAMP_EM_SNR_original_psi, 2);
      GAMP_EM_SNR_mean_term_scale = mean(GAMP_EM_SNR_original_term_scale, 2);
      for snrr = SNRR
        GAMP_EM_SNR_mean_MSE_psi(:, :, snrr + 11) = abs(GAMP_EM_SNR_mean_psi(:, :, snrr + 11) - 10.^(-snrr / 10)).^2 / ( 10.^(-snrr / 10) )^2;
      end
    end

    if SNR_figure_enable
      if GAMP_enable
        for snrr = SNRR
          snr                   = snrr  ;
          GAMP_SNR_last_MSE     = [GAMP_SNR_last_MSE; GAMP_SNR_mean_MSE(iterations, 1, snr + 11)];
          MSE_m_SNR_iter        = [MSE_m_SNR_iter;    MSE_m_iter(iterations, 1, snr + 11)];
          MSE_v_SNR_iter        = [MSE_v_SNR_iter;    MSE_v_iter(iterations, 1, snr + 11)];
          KL_SNR_iter           = [KL_SNR_iter;       KL_iter(iterations, 1, snr + 11)];
        end
      end

      if GAMP_SE_enable
        for snrr = SNRR
          snr                   = snrr  ;
          GAMP_SE_SNR_last_MSE  = [GAMP_SE_SNR_last_MSE; GAMP_SE_SNR_mean_MSE(iterations, 1, snr + 11)];
        end
      end

      if GAMP_EM_enable
        for snrr = SNRR
          snr                         = snrr  ;
          GAMP_EM_SNR_last_MSE        = [GAMP_EM_SNR_last_MSE; GAMP_EM_SNR_mean_MSE(iterations, 1, snr + 11)];
          if Averaging_SNR_enable == 1
            GAMP_EM_SNR_last_psi      = [GAMP_EM_SNR_last_psi; GAMP_EM_SNR_mean_psi(iterations, 1, snr + 11)];
          else
            GAMP_EM_SNR_last_psi      = [GAMP_EM_SNR_last_psi; GAMP_EM_SNR_original_psi(iterations, :, snr + 11)];
          end
          GAMP_EM_SNR_last_term_scale = [GAMP_EM_SNR_last_term_scale; GAMP_EM_SNR_mean_term_scale(iterations, 1, snr + 11)];
          GAMP_EM_SNR_last_MSE_psi    = [GAMP_EM_SNR_last_MSE_psi;    GAMP_EM_SNR_mean_MSE_psi(iterations, 1, snr + 11)];
        end
      end
    end

    if SNR_figure_enable == 0
      if One_figure_enable
        figure(1);
      else
        figure;
      end
      for k = 1:3
        ax(k) = subplot(3,1,k);
      end
      subplot(ax(1));
      if GAMP_enable
        for snrr = SNRR
          snr         = snrr  ;
          plot(iter, 10*log10(GAMP_SNR_mean_MSE(:, 1, snr + 11)), '-o', 'DisplayName', 'oracle', 'Color', [0.4940 0.1840 0.5560]);
          hold on;
        end
      end
      if GAMP_SE_enable
        for snrr = SNRR
          snr         = snrr  ;
          plot(iter, 10*log10(GAMP_SE_SNR_mean_MSE(:, 1, snr + 11)), '*', 'DisplayName', 'SE', 'Color', [0.8500 0.3250 0.0980]);
          hold on;
        end
      end
      if GAMP_EM_enable
        for snrr = SNRR
          snr         = snrr  ;
          if estimation_enable <= 2
            plot(iter, 10*log10(GAMP_EM_SNR_mean_MSE(:, 1, snr + 11)), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
          elseif estimation_enable == 3
            plot(iter, 10*log10(GAMP_EM_SNR_mean_MSE(:, 1, snr + 11)), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
          else
            plot(iter, 10*log10(GAMP_EM_SNR_mean_MSE(:, 1, snr + 11)), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
          end
          hold on;
        end
      end
      xlabel('iter');
      hold on;
      ylabel('MSE-h (dB)');
      hold on;
      legend('Location', 'northeast');
      subplot(ax(2));
      if GAMP_EM_enable
        for snrr = SNRR
          snr         = snrr  ;
          if MSE_psi_enable == 0
            if Averaging_SNR_enable == 1
              if estimation_enable <= 2
                plot(iter, -10*log10(GAMP_EM_SNR_mean_psi(:, :, snr + 11)), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
              elseif estimation_enable == 3
                plot(iter, -10*log10(GAMP_EM_SNR_mean_psi(:, :, snr + 11)), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
              else
                plot(iter, -10*log10(GAMP_EM_SNR_mean_psi(:, :, snr + 11)), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
              end
            else
              if estimation_enable <= 2
                plot(iter, -10*log10(GAMP_EM_SNR_original_psi(:, :, snr + 11)), 'v', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
              elseif estimation_enable == 3
                plot(iter, -10*log10(GAMP_EM_SNR_original_psi(:, :, snr + 11)), 'p', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
              else
                plot(iter, -10*log10(GAMP_EM_SNR_original_psi(:, :, snr + 11)), 'h', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
              end
            end
            hold on;
            yline(snrr,'--','Color',[0.4940 0.1840 0.5560]);
            hold on;
            xlabel('iter');
            hold on;
            ylabel('estimated SNR (dB)');
            hold on;
          else
            if estimation_enable <= 2
              plot(iter, 10*log10(GAMP_EM_SNR_mean_MSE_psi(:, :, snr + 11)), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
            elseif estimation_enable == 3
              plot(iter, 10*log10(GAMP_EM_SNR_mean_MSE_psi(:, :, snr + 11)), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
            else
              plot(iter, 10*log10(GAMP_EM_SNR_mean_MSE_psi(:, :, snr + 11)), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
            end
            hold on;
            xlabel('iter');
            hold on;
            ylabel('MSE-Noise variance (dB)');
            hold on;
          end
          subplot(ax(3));
          plot(iter, 100*(GAMP_EM_SNR_mean_term_scale(:, 1, snr + 11)));
          hold on;
          ytickformat('percentage');
        end
      end
    end

    if SNR_figure_enable == 1
      if One_figure_enable
        figure(3);
      else
        figure;
      end
      for k = 1:3
        ax(k) = subplot(3,1,k);
      end
      subplot(ax(1));
      if GAMP_enable
        plot(SNRR, 10*log10(GAMP_SNR_last_MSE), '-o', 'DisplayName', 'oracle', 'Color', [0.4940 0.1840 0.5560]);
        hold on;
      end
      if GAMP_SE_enable
        plot(SNRR, 10*log10(GAMP_SE_SNR_last_MSE), '*', 'DisplayName', 'SE', 'Color', [0.8500 0.3250 0.0980]);
        hold on;
      end
      if GAMP_EM_enable
        if estimation_enable <= 2
          plot(SNRR, 10*log10(GAMP_EM_SNR_last_MSE), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
        elseif estimation_enable == 3
          plot(SNRR, 10*log10(GAMP_EM_SNR_last_MSE), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
        else
          plot(SNRR, 10*log10(GAMP_EM_SNR_last_MSE), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
        end
        hold on;
      end
      xlabel('SNR (dB)');
      hold on;
      ylabel('MSE-h (dB)');
      hold on;
      legend('Location', 'southwest');
      subplot(ax(2));
      if GAMP_EM_enable
        if MSE_psi_enable == 0
          if Averaging_SNR_enable == 1
            if estimation_enable <= 2
              plot(SNRR, -10*log10(GAMP_EM_SNR_last_psi), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
            elseif estimation_enable == 3
              plot(SNRR, -10*log10(GAMP_EM_SNR_last_psi), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
            else
              plot(SNRR, -10*log10(GAMP_EM_SNR_last_psi), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
            end
          else
            if estimation_enable <= 2
              plot(SNRR, -10*log10(GAMP_EM_SNR_last_psi), 'v', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
            elseif estimation_enable == 3
              plot(SNRR, -10*log10(GAMP_EM_SNR_last_psi), 'p', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
            else
              plot(SNRR, -10*log10(GAMP_EM_SNR_last_psi), 'h', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
            end
          end
          hold on;
          plot(SNRR, SNRR, '--','Color',[0.4940 0.1840 0.5560]);
          hold on;
          hold on;
          xlabel('SNR (dB)');
          hold on;
          ylabel('Estimated SNR (dB)');
          hold on;
        else
          if estimation_enable <= 2
            plot(SNRR, 10*log10(GAMP_EM_SNR_last_MSE_psi), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
          elseif estimation_enable == 3
            plot(SNRR, 10*log10(GAMP_EM_SNR_last_MSE_psi), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
          else
            plot(SNRR, 10*log10(GAMP_EM_SNR_last_MSE_psi), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
          end
          hold on;
          xlabel('SNR (dB)');
          hold on;
          ylabel('MSE-Noise variance (dB)');
          hold on;
        end
        subplot(ax(3));
        plot(SNRR, 100*(GAMP_EM_SNR_last_term_scale), 'Color', [0 0.4470 0.7410]);
        hold on;
        ytickformat('percentage')
        hold on;
        xlabel('SNR (dB)');
        hold on;
        ylabel('term 1 / term 2 + 3');
        hold on;
      end
    end
  end

  if GAMP_enable
    GAMP_p_MSE     = [GAMP_p_MSE, GAMP_SNR_mean_MSE(iterations, 1, SNRR + 11)];
  end
  if GAMP_EM_enable
    GAMP_EM_p_MSE  = [GAMP_EM_p_MSE, GAMP_EM_SNR_mean_MSE(iterations, 1, SNRR + 11)];
    GAMP_EM_p_psi  = [GAMP_EM_p_psi, GAMP_EM_SNR_mean_psi(iterations, 1, SNRR + 11)];
  end
  if GAMP_SE_enable
    GAMP_SE_p_MSE  = [GAMP_SE_p_MSE, GAMP_SE_SNR_mean_MSE(iterations, 1, SNRR + 11)];
  end
end

if SNR_figure_enable == 2
  GAMP_EM_p_MSE_psi = abs(GAMP_EM_p_psi - 10.^(-SNRR / 10)).^2 / ( 10.^(-SNRR / 10) )^2;
  figure(5);
  for k = 1:2
    ax(k) = subplot(2,1,k);
  end
  subplot(ax(1));
  if GAMP_enable
    plot(PP, 10*log10(GAMP_p_MSE), '-o', 'DisplayName', 'oracle', 'Color', [0.4940 0.1840 0.5560]);
    hold on;
  end
  if GAMP_SE_enable
    plot(PP, 10*log10(GAMP_SE_p_MSE), '*', 'DisplayName', 'SE', 'Color', [0.8500 0.3250 0.0980]);
    hold on;
  end
  if GAMP_EM_enable
    if estimation_enable <= 2
      plot(PP, 10*log10(GAMP_EM_p_MSE), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
    elseif estimation_enable == 3
      plot(PP, 10*log10(GAMP_EM_p_MSE), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
    else
      plot(PP, 10*log10(GAMP_EM_p_MSE), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
    end
    hold on;
  end
  xlabel('P');
  hold on;
  ylabel('MSE-h (dB)');
  hold on;
  legend('Location', 'southwest');
  subplot(ax(2));
  if GAMP_EM_enable
    if MSE_psi_enable == 0
      if estimation_enable <= 2
        plot(PP, -10*log10(GAMP_EM_p_psi), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
      elseif estimation_enable == 3
        plot(PP, -10*log10(GAMP_EM_p_psi), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
      else
        plot(PP, -10*log10(GAMP_EM_p_psi), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
      end
      hold on;
      xlabel('P');
      hold on;
      ylabel('Estimated SNR (dB)');
      hold on;
    else
      if estimation_enable <= 2
        plot(PP, 10*log10(GAMP_EM_p_MSE_psi), '-', 'DisplayName', 'ours', 'Color', [0 0.4470 0.7410]);
      elseif estimation_enable == 3
        plot(PP, 10*log10(GAMP_EM_p_MSE_psi), '-', 'DisplayName', 'SR-TSP15', 'Color', [0.9290 0.6940 0.1250]);
      else
        plot(PP, 10*log10(GAMP_EM_p_MSE_psi), '-', 'DisplayName', 'VS-TSP13', 'Color', [0.4660 0.6740 0.1880]);
      end
      hold on;
      xlabel('P');
      hold on;
      ylabel('MSE-Noise variance (dB)');
      hold on;
    end
  end

end