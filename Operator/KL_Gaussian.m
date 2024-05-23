function [KL, MSE_m, MSE_v] = KL_Gaussian(m_p, v_p, m_q, v_q)

  [N,~] = size(m_p);

  m_p_vector = Real_Matrix(m_p);
  v_p_vector = [v_p / 2; v_p / 2];
  m_q_vector = Real_Matrix(m_q);
  v_q_vector = [v_q / 2; v_q / 2];

  KL = 1 ./ 2 * ( ...
    log( prod( v_q_vector ./ v_p_vector ) ) ...
    - 2 .* N ...
    + ( m_p_vector - m_q_vector )' * diag( 1 ./ v_q_vector ) * ( m_p_vector - m_q_vector ) ...
    + sum( v_p_vector ./ v_q_vector ) ...
  );

  MSE_m = sum( abs( m_p - m_q ).^2 ) ./ sum( abs( m_p ).^2 );
  MSE_v = sum( abs( v_p - v_q ).^2 ) ./ sum( abs( v_p ).^2 );

end