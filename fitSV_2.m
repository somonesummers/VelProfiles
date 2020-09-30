function [varargout] = fitSV_2(x, z, t, sMax, vMax, dr, m2_p)
    % [sv, m2] = fitVelocity(x, z, t, s, vMax, sMax, dr, m2_p) fits for
    % unknown slope-velocity from a prior slope-velocity estimate. 
    % [sv, m2, F] = fitVelocity(x, z, t, s, vMax, sMax, dr, m2_p) also returns
    % the fit values
    
    nOutputs = nargout;
    varargout = cell(1,nOutputs);
    sv = zeros(size(z)); % Vel out
    m2 = zeros(2,size(z,2)); % array of fits for each layer
    F = zeros(size(z)); % fit scores
    
    %Fit each profile line with a sine curve
    for i = size(x,2):-1:1
        yy = real(x(:,i));
        xx = t;
        % Define fit function
        fit = @(b,x)  1 .* (sin(2 * pi * xx * b(1) + 2 * pi * b(2))) + 0;    % Function to fit
        fcn = @(b) sum((fit(b,xx) - yy).^2);
        OPTIONS = optimset('Display','none');% Least-Squares cost function
        buff = 5; %buffer on bounded search
        % Seed 3 different starting phases to find global min
        [m_0  , fval_0  ] = fminsearchbnd(fcn, [m2_p(i);   0   ],[1/buff*m2_p(i) -2*pi],[buff*m2_p(i) 2*pi]);
        [m_pi , fval_pi ] = fminsearchbnd(fcn, [m2_p(i);   pi/2],[1/buff*m2_p(i) -2*pi],[buff*m2_p(i) 2*pi]);
        [m_npi, fval_npi] = fminsearchbnd(fcn, [m2_p(i);  -pi/2],[1/buff*m2_p(i) -2*pi],[buff*m2_p(i) 2*pi]);
        % Pick the best fit, disguard fit value f
        [f , I] = min([fval_0, fval_pi, fval_npi]);
        m_s = [m_0, m_pi, m_npi];
        m = m_s(:,I);
        m2(:,i) = m;
        F(i) = f;
        sv(i) = abs(m(1) * dr); %calc output sv profile
    end
    varargout{1} = sv;
    varargout{2} = m2;
    if(nOutputs ==3)
        varargout{3} = F;
    end
end