function [pa,un,fe] = flipStratification(pa,un,fe)

switch pa.strat_model{1}.name

    case 'double_pycnocline'
        % ---- Unpack stratification values
        s1=pa.s1;
        s2=pa.s2;
        Delta1=pa.Delta1;
        Delta2=pa.Delta2;
        H1=pa.H-pa.H2-pa.H3;

        % ---- Flip the stratification: the middle layer keeps its depth
        pa.s2=s1; pa.s1=s2;
        pa.Delta2=Delta1; pa.Delta1=Delta2;
        pa.H3=H1;

    case 'single_pycnocline'
        % ---- Flip the stratification
        pa.H2 = pa.H-pa.H2;

end

% ---- Unpack the streamfunction
psi = un.psi;
N = pa.N; M = pa.M;
Psi = reshape(psi,[],N);

% ---- Flip the streamfunction
Psi_new = zeros(M,N);
for j=2:N
    Psi_new(:,j) = Psi_new(:,j-1) + (Psi(:,N-j+2)-Psi(:,N-j+1));
end

un.psi = reshape(Psi_new,[],1);
un.Psi = Psi_new;

% ---- Change the location of the measure of amplitude
Amploc = fe.Amploc;
AmplocN = Amploc/pa.M;
fe.Amploc = round((1-AmplocN/N)*pa.N*pa.M);
