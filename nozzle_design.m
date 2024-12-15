function [xw,yw,xcl,Mcl] = MinLenNozDes(yt,Me,gam,nlines,thi,oflag,pflag)
%
%
% Function which computes the contour of a minimum length diverging nozzle 
% that will expand the flow from sonic speed to the desired value of Me.
% The method is implemented for a two-dimensional nozzle, i.e., a nozzle
% that has uniform properties in the span-wise (in the direction normal to
% the plane of the contour) direction.
% 
% The user must input the value of gamma and the number of characteristic
% lines to be used. The user may also suppress the visual output by
% setting the output visulaizatoin flag oflag = 0 or activate it by setting
% oflag = 1 (default is 0).
%
% All angles in degrees
%
% Input parameters:
% yt...... half throat height; its dimensions sets the dimension (scale) of
%          the resulting nozzle
% Me...... exit mach nunber
% gam..... gamma
% nlines.. number of characteristic lines
% thi..... angle of nozzle throat, defines the angle of the first
%          characteristic line emerging from the throat; typically this is
%          a small value, such as 0.1deg
% oflag... flag to show the solution, set to 1 to show solution, default is
%          set to o (nothing is displayed)
% pflag... flag to print to file figures, if oflag set to 1 (both oflag and
%          pflag must be set to 1 to print figures to file
%
% Output variables:
% xw...... x-coordinate of nozzle wall 
% yw...... y-coordinate of nozzle wall
% xcl..... x-coordinate on centerline of nozzle
% Mcl..... Mach number of centerline of nozzle
%
% Other internal variables of interest:
% Prat.... Refers to the pressure along the nozzle, normalized by the
%          stagnation pressure at the throat
% Trat.... Refers to the temperature along the nozzle, normalized by the
%          stagnation temperature at the throat
%
% Example of use:
% Run default values: MinLenNozDes;
%
% Run by specifing input, but don't display results: MinLenNozDes(1, 3.3, 1.4, 11, 0.1);
%
% Run by specifing input, and display results: MinLenNozDes(1, 3.3, 1.4, 11, 0.1, 1);
%
%
%


%
%
% Credit: this function was originally developed by Logan White, adapted by
%         Mirko Gamba
% 


%%% Initialize input variables
if nargin == 0
    yt = 1; Me = 3.3; gam = 1.4; nlines = 11; thi = 0.1; oflag = 0; pflag = 0;
elseif nargin <= 5
    oflag = 0;
    pflag = 0;
elseif nargin <= 6
    pflag = 0;
end

%pflag = 0;


%%% Initialize solution

xt = 0;
[MACH, nue, MU] = flowprandtlmeyer(gam, Me, 'mach');
thmax = nue / 2; % deg
dth = (thmax - thi)/(nlines - 1); % thmax = 18.375;



%%%% Construct and solve through network of characteristics
% dth = 3;
th(1) = thi;
nu(1) = th(1);
[M(1), NU, mu(1)] = flowprandtlmeyer(gam, nu(1), 'nu'); 
Km(1) = th(1) + nu(1);
Kp(1) = th(1) - nu(1);
xmlines = NaN * ones(nlines+1); 
ymlines = NaN * ones(nlines+1); 
xmlines(1,:) = xt;
ymlines(1,:) = yt;
xplines = NaN * ones(nlines+1); 
yplines = NaN * ones(nlines+1); 
yplines(1,:) = 0;
xcl = zeros(1,nlines+1); 
Mcl = ones(1,nlines+1); 
Mcl(2) = M(1);


for i = 2:nlines
    th(i) = th(1) + (i - 1)*dth;
    nu(i) = th(i);
    [M(i), NU, mu(i)] = flowprandtlmeyer(gam, nu(i), 'nu'); 
    Km(i) = th(i) + nu(i);
    Kp(i) = th(i) - nu(i);
end
th(nlines+1) = th(nlines); 
nu(nlines+1) = nu(nlines); 
M(nlines+1)  = M(nlines); 
mu(nlines+1) = mu(nlines); 
Km(nlines+1) = Km(nlines); 
Kp(nlines+1) = Kp(nlines);
last = nlines + 1; 



for j = 1:(nlines-1)   
    th(last+1) = 0;
    Km(last+1) = Km(last-(nlines-j));
    nu(last+1) = Km(last+1)-th(last+1);
    Kp(last+1) = th(last+1)-nu(last+1);
    [M(last+1), NU, mu(last+1)] = flowprandtlmeyer(gam, nu(last+1), 'nu');
    Mcl(j+2) = M(last+1);
    
    for k = 1:(nlines-1-j)
        Km(last+1+k) = Km(last+1) + 2*k*dth;
        Kp(last+1+k) = Kp(last+1);
        th(last+1+k) = (Km(last+1+k)+Kp(last+1+k))/2;
        nu(last+1+k) = (Km(last+1+k)-Kp(last+1+k))/2;
        [M(last+1+k), NU, mu(last+1+k)] = flowprandtlmeyer(gam, nu(last+1+k), 'nu');
    end

    last = last + (nlines - j)+1;
    th(last) = th(last-1); nu(last) = nu(last-1); M(last) = M(last-1); mu(last) = mu(last-1); Km(last) = Km(last-1); Kp(last) = Kp(last-1);
end


out = [ [1:length(Km)]' Km' Kp' th' nu' M' mu'];



% Find locations of first Cplus points
y(1) = 0;
x(1) = yt / -tand(th(1)-mu(1));
xcl(2) = x(1);
xmlines(2,1) = x(1); 
ymlines(2,1) = y(1);

for i = 2:nlines
    mw = tand(th(i)-mu(i)); 
    bw = yt - mw*xt;
    mp = tand(th(i)+mu(i)); 
    bp = y(i-1) - mp*x(i-1);
    x(i) = (bw - bp)/(mp - mw); 
    y(i) = mw*x(i) + bw;
    xmlines(2,i) = x(i); 
    ymlines(2,i) = y(i);
end

last = nlines + 1;
mw = tand((th(last)+thmax)/2);
bw = yt - mw*xt;
    
mp = tand(th(last)+mu(last)); 
bp = y(last-1) - mp*x(last-1); 
x(last) = (bw - bp)/(mp - mw); 
y(last) = mw*x(last)+bw;
xplines(1:length(x),1) = x'; 
yplines(1:length(y),1) = y';
xw(1) = xt; 
yw(1) = yt;


for j = 1:(nlines-1) 
    y(last+1) = 0;
    x(last+1) = y(last+1-(nlines-j+1))/-tand(th(last+1-(nlines-j+1))-mu(last+1-(nlines-j+1)))+x(last+1-(nlines-j+1)); 
    xcl(j+2) = x(last+1);
    xplines(1,j+1) = x(last+1);
    xmlines(j+2,j+1) = x(last+1);
    ymlines(j+2,j+1) = y(last+1);
    
    for k = 1:(nlines-1-j)
        m1 = tand(th(last+k)+mu(last+k));
        b1 = y(last+k) - m1*x(last+k);
        m2 = tand(th(last+1+k-(nlines-j+1))-mu(last+1+k-(nlines-j+1))); 
        b2 = y(last+1+k-(nlines-j+1)) - m2*x(last+1+k-(nlines-j+1)); 
        x(last+1+k) = (b2 - b1)/(m1 - m2);
        y(last+1+k) = m1*x(last+1+k) + b1;
        xplines(1+k,j+1) = x(last+1+k); 
        yplines(1+k,j+1) = y(last+1+k);
        xmlines(j+2,j+k+1) = x(last+1+k); 
        ymlines(j+2,j+k+1) = y(last+1+k);
    end

    last = last + (nlines - j)+1;
    mw = tand((th(last)+th(last-(nlines-j+1)))/2);
    bw = y(last-(nlines-j+1)) - mw*x(last-(nlines-j+1)); 
    mp = tand(th(last-1)+mu(last-1));
    bp = y(last-1) - mp*x(last-1);
    x(last) = (bw - bp)/(mp - mw);
    y(last) = mw*x(last)+bw;
    xw(1+j) = x(last); 
    yw(1+j) = y(last);
    xplines(nlines+1-j,j+1) = x(last); 
    yplines(nlines+1-j,j+1) = y(last);
end



%%% Now plot data if required
if oflag == 1 
    close all
    
    disp(out)
    [X,Y] = meshgrid(0:0.5:x(end),0:0.5:y(end)); 
    Mq = griddata(x,y,M,X,Y);
    [xlen,ylen] = size(Mq);
    for i = 1:xlen
        numcheck = 0; 
        for j = 1:ylen
            if isnan(Mq(i,j)) == 0 
                numcheck = 1;
            end
            if (numcheck == 1)&&(isnan(Mq(i,j)) == 1)
                Mq(i,j) = Mq(i,j-1);
            end
        end
    end
    Trat = (1 + (gam - 1)/2*Mq.^2).^(-1); 
    Prat = Trat.^(gam/(gam-1));
    
    
   
    figure(19) % title('Nozzle Contour with Characteristics')
    for i = 1:nlines
        plot(xplines(:,i),  yplines(:,i),'k.-','LineWidth',1,'MarkerSize',12), hold on 
        plot(xplines(:,i), -yplines(:,i),'k.-','LineWidth',1,'MarkerSize',12) 
        plot(xmlines(:,i),  ymlines(:,i),'k.-','LineWidth',1,'MarkerSize',12) 
        plot(xmlines(:,i), -ymlines(:,i),'k.-','LineWidth',1,'MarkerSize',12)
    end
    plot(xw,  yw,'k-','LineWidth',3) 
    plot(xw, -yw,'k-','LineWidth',3)
    formatplot([])
    if pflag == 1; print(gcf, '-dpng', 'nozzle_minlength_characteristics', '-r300'); end
    
    figure(10) % title('Mach Number Distribution in Nozzle')
    pcolor(X, Y,Mq), hold on
    pcolor(X,-Y,Mq) 
    plot(xw,yw,'k-','LineWidth',2) 
    plot(xw,-yw,'k-','LineWidth',2) 
    shading interp
    formatplot('$M$')
    if pflag == 1; print(gcf, '-dpng', 'nozzle_minlength_mach', '-r300'); end

    figure(11)
    %title('Temperature Distribution in Nozzle') 
    pcolor(X, Y,Trat), hold on
    pcolor(X,-Y,Trat) 
    plot(xw, yw,'k-','LineWidth',2) 
    plot(xw,-yw,'k-','LineWidth',2)
    shading interp
    formatplot('$T / T_o$')
    if pflag == 1; print(gcf, '-dpng', 'nozzle_minlength_temperature', '-r300'); end
    
    figure(12) % title('Pressure Distribution in Nozzle')
    pcolor(X,Y,Prat)
    hold on
    pcolor(X,-Y,Prat)
    plot(xw,yw,'k-','LineWidth',2) 
    plot(xw,-yw,'k-','LineWidth',2)
    shading interp
    formatplot('$p / p_o$')
    if pflag == 1; print(gcf, '-dpng', 'nozzle_minlength_pressure', '-r300'); end
    
end

end



%%%%%%%%%%%%%% Internal functions

function [MACH, NU, MU] = flowprandtlmeyer(g, inp, type)

PM = @(M, g) ( sqrt((g+1)/(g-1)).*atan(sqrt((g-1)/(g+1).*(M.^2-1))) - atan(sqrt(M.^2-1)) ) / pi * 180 ;
T = @(M, g) (1 + 0.5*(g-1).*Me.^2)./(1 + 0.5*(g-1).*M.^2);
P = @(M, g) ((1 + 0.5*(g-1).*Me.^2)./(1 + 0.5*(g-1).*M.^2)).^(g./(g-1));
mach = @(M) asind(1./M);

switch type
    case 'nu' % invert function
        opt = optimset('MaxFunEvals',100000,'display','off','MaxIter',100000);
        [M, fval, exitflag, output, jacobian] = fsolve(@(x) PM(x, g)-inp, 2, opt);
    case 'mach' % given mach number
        M = inp;
end

MACH = M;
NU = PM(MACH, g);
MU = mach(MACH);

end



function formatplot(clab)
xlabel('Distance from throat, mm', 'interpreter', 'latex', 'fontsize',24,'fontname','times')
ylabel('Distance from axis, mm',   'interpreter', 'latex', 'fontsize',24,'fontname','times')

set(gca,'fontsize',24,'fontname','times')
set(gca,'linewidth',1.5,'box','off','ticklength',[.01 0])
set(gca,'tickdir','out')

axis equal

if ~isempty(clab)
    h = colorbar;
    ylabel(h, clab,'fontsize',24,'fontname','times','interpreter', 'latex')
end

set(gcf,'position',[100 100 900 600],'color',[1 1 1],'paperPositionMode','auto')
set(gca, 'position', [0.06,0.12,0.88, 0.85])
    
end
