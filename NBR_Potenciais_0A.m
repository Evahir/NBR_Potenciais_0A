%% NBR 15751:2013 - Anexo B - Calculo das tensoes de passo e de toque
%% Metodologia simplificada para malhas reticuladas retangulares
%% ========================================================================
%% Aplicado as equacoes [24] a [30] do Anexo B da NBR 15751:2013 para
%% calcular a tensao de toque (Vt) e a tensao de passo (Vp) efetivas
%% produzidas pela malha, comparando com os limites permissiveis.
%% -----------------------------------------------------------------------
%% 1) Determinação de parâmetros iniciais
%% -----------------------------------------------------------------------
If = 12.5;              % Corrente de falta fase-terra, kA
t  = 0.3;                % Tempo de atuacao do sistema de protecao, s

Df = 1.085;              % Fator de decremento
Sf = 0.8;                % Fator de distribuicao
Cp = 1.1;                % Fator de projecao (crescimento do sistema)

Im_kA = Df*Sf*If*Cp;     % Corrente de malha simetrica eficaz [19], kA
Imcd  = Im_kA*1000;      % Convertida para A (curta duracao) - usada no Anexo B
fprintf('Corrente de malha (Im) = %.4f kA = %.1f A\n', Im_kA, Imcd);

Rch  = 1000;              % Resistencia do corpo humano, ohm
Ichcd = 0.116/sqrt(t);    % Corrente de choque de curta duracao, A [6]

rho_s = 3000;             % Resistividade da camada superficial (brita), ohm.m
rho_1 = 700;               % Resistividade da 1a camada de solo, ohm.m
rho_2 = 171.4597;          % Resistividade da 2a camada de solo, ohm.m

a2 = 1-(rho_1/rho_s);
a3 = 2*0.1+0.106;          % a = 0,106 m (constante da NBR, eq [14])
C  = 1 - 0.106*(a2/a3);    % Fator de reducao C

Valor_permissivel_passo = (Rch + 6*rho_s*C)*Ichcd;   %#ok<*NOPTS>  % [10]
Valor_permissivel_toque = (Rch + 1.5*rho_s*C)*Ichcd;              % [11]
fprintf('Tensao de passo permissivel = %.1f V\n', Valor_permissivel_passo);
fprintf('Tensao de toque permissivel = %.1f V\n', Valor_permissivel_toque);

%% -----------------------------------------------------------------------
%% 2) Resistividade do solo a usar no calculo dos potenciais (Anexo B)
%% -----------------------------------------------------------------------
% O Anexo B (texto introdutorio) recomenda:
%   a) usar rho_1 quando rho_2 < rho_1;
%   b) usar rho_a (aparente) quando rho_2 > rho_1.
% No nosso caso rho_2 (171,4597) < rho_1 (700) -> usa-se rho_1.
if rho_2 < rho_1
    rho_calc = rho_1;
    fprintf('rho_2 < rho_1 -> usando rho_1 = %.1f ohm.m no calculo de Vt/Vp\n', rho_calc);
else
    rho_calc = 329;  % rho_a (calculado no Anexo A) - ajustar se necessario
    fprintf('rho_2 > rho_1 -> usando rho_a = %.1f ohm.m no calculo de Vt/Vp\n', rho_calc);
end

%% -----------------------------------------------------------------------
%% 3) Geometria real da malha (extraida do DXF do AutoCAD)
%% -----------------------------------------------------------------------
Lt_horizontais = 11241.416;   % Comprimento dos condutores horizontais, m (do DXF)
n_hastes = 48;
comp_haste = 3.0;             % m

% B.1.1: Lt = comprimento total enterrado, incluindo hastes
Lt = Lt_horizontais + n_hastes*comp_haste;
fprintf('Lt (condutores + hastes) = %.3f m\n', Lt);

H = 0.6;                      % Profundidade da malha, m
Secao_utilizada_mm2 = 35;     % Secao nominal escolhida na 1a etapa, mm2
d = sqrt(4*(Secao_utilizada_mm2*1e-6)/pi);  % Diametro do condutor, m
fprintf('Diametro do condutor (secao %d mm2) = %.4f mm\n', Secao_utilizada_mm2, d*1000);

%% -----------------------------------------------------------------------
%% 4) Fator kii - hastes na periferia?
%% -----------------------------------------------------------------------
% kii = 1 para malhas COM hastes de aterramento locadas na periferia
% (nao exige que TODAS estejam la). A formula alternativa so vale para
% malhas sem hastes, ou com hastes apenas fora da periferia.
% Conferido nas coordenadas do DXF: das 48 hastes, 26 estao na periferia
% (e 22 sao internas, "ao longo de toda a malha") -> kii = 1 se aplica.
hastes_com_periferia = true;   % ajuste aqui se a geometria mudar

%% -----------------------------------------------------------------------
%% 5) Calculo de Km, Ki, Ks e Vt/Vp nas duas direcoes da malha
%%    (eq. [25],[26],[27],[30] - NBR 15751 Anexo B)
%% -----------------------------------------------------------------------
% Direcao 1: condutores horizontais, espacados 6,25 m (n = 25)
% Direcao 2: condutores verticais,   espacados 12,50 m (n = 25)
D = [6.25, 12.5];
n = [25, 25];

Vt = zeros(1,2);
Vp = zeros(1,2);
Km = zeros(1,2);
Ks = zeros(1,2);
Ki = zeros(1,2);
kii = zeros(1,2);

for k = 1:2
    Dk = D(k);
    nk = n(k);

    if hastes_com_periferia
        kii(k) = 1;
    else
        kii(k) = (2*nk)^(-2/nk);         % [26]
    end

    Ki(k) = 0.656 + 0.172*nk;             % [27]

    Km(k) = (1/(2*pi)) * ( ...
                log( (Dk^2)/(16*H*d) + ((Dk+2*H)^2)/(8*Dk*d) - H/(4*d) ) ...
                + (kii(k)/sqrt(1+H)) * log(8/(pi*(2*nk-1))) );   % [25]

    Ks(k) = (1/pi) * ( 1/(2*H) + 1/(Dk+H) + (1/Dk)*(1-0.5^(nk-2)) );  % [30]

    Vt(k) = (rho_calc * Imcd * Km(k) * Ki(k)) / Lt;   % [24]
    Vp(k) = (rho_calc * Imcd * Ks(k) * Ki(k)) / Lt;   % [29]

    fprintf('\n--- Direcao %d (D = %.2f m, n = %d) ---\n', k, Dk, nk);
    fprintf('  kii = %.4f | Ki = %.4f\n', kii(k), Ki(k));
    fprintf('  Km  = %.4f  ->  Vt = %.2f V\n', Km(k), Vt(k));
    fprintf('  Ks  = %.4f  ->  Vp = %.2f V\n', Ks(k), Vp(k));
end

%% -----------------------------------------------------------------------
%% 6) Piores casos (maiores valores encontrados nas duas direcoes)
%% -----------------------------------------------------------------------
[Vt_max, idx_t] = max(Vt);
[Vp_max, idx_p] = max(Vp);

fprintf('\n=========================================================\n');
fprintf('Tensao de toque calculada (pior caso, direcao %d) = %.2f V\n', idx_t, Vt_max);
fprintf('Tensao de toque PERMISSIVEL                       = %.2f V\n', Valor_permissivel_toque);
if Vt_max <= Valor_permissivel_toque
    fprintf('  -> OK: dentro do limite\n');
else
    fprintf('  -> ATENCAO: tensao de toque calculada EXCEDE o limite permissivel!\n');
end

fprintf('\nTensao de passo calculada (pior caso, direcao %d)  = %.2f V\n', idx_p, Vp_max);
fprintf('Tensao de passo PERMISSIVEL                        = %.2f V\n', Valor_permissivel_passo);
if Vp_max <= Valor_permissivel_passo
    fprintf('  -> OK: dentro do limite\n');
else
    fprintf('  -> ATENCAO: tensao de passo calculada EXCEDE o limite permissivel!\n');
end
fprintf('=========================================================\n');
