function C = mimo_capacity(eqv_H, tx_cov, n_cov)
[r_eqv_H, c_eqv_H] = size(eqv_H); % size of eqv_H
% MIMO capacity for equivalent channel eqv_H, transmit covariance
% tx_cov and noise covariance n_cov
C = log(det(eye(r_eqv_H)+inv(n_cov)*eqv_H*tx_cov*eqv_H'))/log(2);
