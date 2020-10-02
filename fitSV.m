function [varargout] = fitSV(varargin)
    % [sv, m2] = fitVelocity(x, z, t, vMax, sMax, dr) solves for velocity profile
    % of a region of known layer slopes (s(z)). velocity is assumed to be
    % roughly monotonically increasing with depth, with max velocity at the
    % bed. Estimated max velocity is an input vMax. Returns Vels and sin
    % fits. 
    % [sv, m2, F] = fitVelocity(x, z, t, vMax, sMax, dr) also returns the
    % score of the Fits
    % = fitVelocity(x, z, t, vMax, sMax, dr, v_zb) also takes a vertical
    % velocity profile into account.
    nOutputs = nargout;
    varargout = cell(1,nOutputs);
    nInputs = nargin;
    x = varargin{1};
    z = varargin{2};
    t = varargin{3};
    vMax = varargin{4};
    sMax = varargin{5};
    dr = varargin{6};
    if(nInputs == 7)
        v_zb = varargin{7};
    else
        v_zb = 0;
    end
    sv = zeros(size(z)); % Vel out
    m2 = zeros(2,size(z,2)); % array of fits for each layer
    m = abs([(vMax*sMax + v_zb)/dr 0]);% Lat fit value
    F = zeros(size(z)); % fit scores
    
    %Fit each profile line with a sine curve
    for i = size(x,2):-1:1
        yy = real(x(:,i));
        xx = t;
        yu = max(yy);
        yl = min(yy);
        yr = (yu-yl)/2;                               % Range of y
        yz = movmean(yy-yu+(yr),floor(length(yy)/20));
        zx = xx(yz .* circshift(yz,[1 0]) <= 0);     % Find zero-crossings
        % Estimate period to set initalization for fminsearch, manually
        % enlongate estimated period for small n cases
        if length(zx) > 3
            per = 2*mean(diff(zx));                     
        elseif length(zx) == 2
            per = 1*max(xx)-min(xx);
        elseif length(zx) == 1
            per = 2*max(xx)-min(xx);
        else
            per = 5*max(xx)-min(xx);
        end
        % Define fit function
        fit = @(b,x)  1 .* (sin(2 * pi * xx * b(1) + 2 * pi * b(2))) + 0;    % Function to fit
        fcn = @(b) sum((fit(b,xx) - yy).^2);
        OPTIONS = optimset('Display','none');% Least-Squares cost function
        buff = 2; %buffer on bounded search
        if i < size(x,2)
            bd = m(1) * 2; %bound search to look for lower velocities
        else
            bd = m(1); 
        end
        %limit range of these bounds manually
        if(abs(bd * buff *dr) > 1e-5)
            bd = 1e-5 / (buff * dr);
        elseif(abs(bd * buff *dr) < 1e-9)
            bd = 1e-9 / (buff * dr);
        end
        % Seed 3 different starting phases to find global min
        [m_0  , fval_0] = fminsearchbnd(fcn, [1/per;  0],[0 -2*pi],[abs(bd)*buff 2*pi]);
        [m_pi , fval_pi] = fminsearchbnd(fcn, [1/per;  pi/2],[0 -2*pi],[abs(bd)*buff 2*pi]);
        [m_npi, fval_npi] = fminsearchbnd(fcn, [1/per;  -pi/2],[0 -2*pi],[abs(bd)*buff 2*pi]);
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