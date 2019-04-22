function rlc(R, L, C, N)
  alpha = 1/(2*R*C)
  wo = 1/sqrt(L*C)

  A = [0, 1;
     -wo^2, -2*alpha]
  B = [0; wo^2]

  alpha^2/wo^2

  [V, D] = eig(A)

  % Como a tensão inicial no capacitor é v(0) = 10V
  % Temos que v(0) = L * iL'(0) 
  % e, portanto, iL'(0) = v(0)/L

  % Y = inv(V) * [iL(0), iL'(0)]
  i = 20*rand(N)-10;
  di = 300*rand(N)-150;
  
  ASv = inv(V)*[i(1); di(1)/L]
  Xt = V*diag(ASv)
  t = linspace(0,8*10^-6,wo/100);
  U = L*(Xt(2,1)*exp(D(1,1)*t) +Xt(2,2)*exp(D(2,2)*t));
  figure(1)
  plot(t, U);
  xlabel('t');
  ylabel('Tensão');
  %hold on;
  cor = Xt(1,1)*exp(D(1,1)*t) +Xt(1,2)*exp(D(2,2)*t);
  figure(2)
  plot(t, cor);
  xlabel('t');
  ylabel('Corrente');
  
  figure(3)
  for k = 1:1:N
      AS = inv(V)*[i(k); di(k)/L]

      Xt = V*diag(AS)

      %Para plotar
      t = linspace(0,8*10^-6,wo/100);
      X1 = Xt(1,1)*exp(D(1,1)*t) +Xt(1,2)*exp(D(2,2)*t);
      X2 = Xt(2,1)*exp(D(1,1)*t) +Xt(2,2)*exp(D(2,2)*t);  
      plot(X1,X2)
      xlabel('X1 - i(t)')
      ylabel("X2 - di(t)/dt")
      hold on
    endfor
endfunction