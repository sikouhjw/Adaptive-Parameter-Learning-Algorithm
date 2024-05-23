function psi = Estimate_AWGN_noise_EM(Input, obj, m_z, v_z)
  In = Input.In;
  Out = Input.Out;
  quan_step = Out.quan_step;
  y = obj.y;
  p = Out.p;
  MAX_TERM = Out.MAX_TERM;
  M = Input.M;

  psi = mean( abs( y - m_z ).^2 + v_z );

end