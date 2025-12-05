% The matlab functions that describe the diffusion kinetics are separated
% into the cases for an impermeable and permeable film. In each of those
% cases, functions are provided for concentration gradients (as a function
% of y and t), integration concentrations (as a function of time) for a
% plasmonic interaction near the interface and over the entire film for the
% absorption and desorption processes.

%% PARAMETERS

% X1 (float): Analyte concentration
% D (float): Diffusion constant, in units of cm^2/s
% d (float): Thickness of membrane, in units of cm
% dy (float): Height of evanescent wave interaction, in units of cm
% t (array of float/int): Time, in units of s
% y (array of float/int): Height, in units of cm
% num_terms (int): Number of terms to sum over

X1=1;
%X2=0;
D=1e-6;  % or 1e-2
d=0.01 ; % or 1 
t=1:200; 
dy=d./10;

num_terms = 50;
%% Plotting example

figure
% ============================ Tiled layout figure ============================
tcl = tiledlayout(1,2);
tcl.TileSpacing="compact";
tcl.Padding="compact";
tcl.Subtitle.HorizontalAlignment = 'right'; 
    
nexttile(tcl)
    y=1:100; % better define y as a percentage of d: y=0:100 where 100 = 100% x d
    t=0.1;
    plot(y/100, Imp_Abs_conc(X1,D,d,t,d/100*y,num_terms), "Color", [1 0 0])
    hold on;
    for t=1:2:20
        plot(y/100, Imp_Abs_conc(X1,D,d,t,d/100*y,num_terms), "Color", [1 0 t/20])
    end
    for t=30:15:150
        plot(y/100, Imp_Abs_conc(X1,D,d,t,d/100*y,num_terms), "Color", [1-t/150 0 1])
    end
    hold off    
%    title('Impermeable film concentration gradients');
%    title('(A)')
    xlim([-.10 1.10])
    ylim ([-0.1 1.1])
   xlabel('Film thickness, y/d') 
   ylabel('Measurement, X_{1}')
    rectangle('Position',[-.10,0,.10,1],'FaceColor',[0 .5 .5],'LineWidth',1)
    rectangle('Position',[1.00,0,.10,1],'EdgeColor',[0 .5 .5],'LineWidth',1)
    line([.10 .10],[0 1],'Color','red','LineStyle','--');

nexttile(tcl)
    t=1:200;
    t = [0.01, t];
    plot((pi^2*D*t)/(4*d^2), Imp_Abs_ev(X1,D,d,t,dy,num_terms),lineWidth=1)
    hold on;
    plot((pi^2*D*t)/(4*d^2),Imp_Abs_film(X1,D,d,t,num_terms),lineWidth=1)
    hold off
%    title('Impermeable film uptake curves ');    
%    title('(B)')
    xlim([-.10 5.00])
    ylim ([-0.1 1.1])
   xlabel('Reduced Time, t/\tau = \pi^2Dt/4d^{2}')
   ylabel('Measurement, X_{1}')

%% *************** IMPERMEABLE FILM ***************
%% Concentration gradients for absorption and desorption
% equation (10): Concentration gradient across impermeable membrane during
% absorption
function f=Imp_Abs_conc(X1,D,d,t,y,num_terms)
    i=0;
    Terms(num_terms:length(y))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/(2*j+1)*cos((2*j+1)*pi*y/(2*d))*exp(-D.*(2*j+1)^2*pi^2*t/(4*d^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(1-4/pi*SumofTerms);
end

% equation (12): Concentration gradient across impermeable membrane during
% desorption
function f=Imp_Des_conc(X1,D,d,t,y, num_terms)
    i=0;
    Terms(num_terms:length(y))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/(2*j+1)*cos((2*j+1)*pi*y/(2*d))*exp(-D.*(2*j+1)^2*pi^2*t/(4*d^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(4/pi*SumofTerms);
end

%% Absorption and desorption processes through evanescent wave

% equation (16): Integrated concentration through evanescent wave during
% absorption as a function of time
function f=Imp_Abs_ev(X1,D,d,t,dy, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)= (-1)^j/(2*j+1)^2*sin((2*j+1)*pi*dy/(2*d))*exp(-D.*(2*j+1)^2*pi^2*t/(4*d.^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(1-8*d/(pi^2*dy)*SumofTerms);
end

% equation (18): Integrated concentration through evanescent wave during
% desorption as a function of time
function f=Imp_Des_ev(X1,D,d,t,dy, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)= (-1)^j/(2*j+1)^2*sin((2*j+1)*pi*dy/(2*d))*exp(-D.*(2*j+1)^2*pi^2*t/(4*d.^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(8*d/(pi^2*dy)*SumofTerms);
end

%% Absorption and desorption processes integrated over film

% equation (20): Integrated concentration in film during absorption as a
% function of time
function f=Imp_Abs_film(X1,D,d,t, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)= 1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/(4*d.^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(1-8/pi^2*SumofTerms);
end

% equation (22): Integrated concentration in film during desorption as a
% function of time
function f=Imp_Des_film(X1,D,d,t, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)= 1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/(4*d.^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(8/pi^2*SumofTerms);
end
%% *************** PERMEABLE FILM ***************
%% Concentration gradients for absorption and desorption
% equation (25): Concentration gradient across permeable membrane during
% absorption
function f=Per_Abs_conc(X1,D,d,t,y, num_terms)
    i=0;
    Terms(num_terms:length(y))=zeros;
        for j=1:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/j*sin(j*pi*y/d)*exp(-D.*j^2*pi^2*t/(d^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(y/100+2/pi*SumofTerms);
    %f=X1*(2/pi*SumofTerms);
end

% equation (27): Concentration gradient across permeable membrane during
% desorption assuming linear concentration gradient
function f=Per_Des_conc_grad(X1,D,d,t,y, num_terms)
    i=0;
    Terms(num_terms:length(y))=zeros;
        for j=1:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/(j)*sin(j*pi*y/d)*exp(-D.*j^2*pi^2*t/(d^2));
        end
    SumofTerms=sum(Terms);
    f=-X1*(2/pi*SumofTerms);
    %f=X1*(2/pi*SumofTerms);
end

% equation (28): Concentration gradient across permeable membrane during
% desorption assuming constant initial concentration
function f=Per_Des_conc_const(X1,D,d,t,y, num_terms)
    i=0;
    Terms(num_terms:length(y))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)*sin((2*j+1)*pi*y/d)*exp(-D.*(2*j+1)^2*pi^2*t/(d^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(4/pi*SumofTerms);
    %f=X1*(2/pi*SumofTerms);
end

%% Absorption and desorption processes through evanescent wave
% equation (39): Integrated concentration through evanescent wave during
% absorption as a function of time
function f=Per_Abs_ev(X1,D,d,t,dy, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=1:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/j^2*(cos(j*pi*dy/d)-1)*exp(-D.*j^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(dy/(d*2)-2*d/(dy*pi^2)*SumofTerms);
end

% equation (41): Integrated concentration through evanescent wave during
% desorption assuming constant initial concentration as a function of time
function f=Per_Des_ev_const(X1,D,d,t,dy, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*(1-cos((2*j+1)*pi*dy/d))*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(4*d/(dy*pi^2)*SumofTerms);
end

% equation (42): Integrated concentration through evanescent wave during
% desorption assuming linear concentration gradient as a function of time
function f=Per_Des_ev_grad(X1,D,d,t,dy, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=1:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/j^2*(cos(j*pi*dy/d)-1)*exp(-D.*j^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(2*d/(dy*pi^2)*SumofTerms);
end

%% Absorption and desorption processes integrated over film

% equation (44): Integrated concentration over film during absorption as a
% function of time
function f=Per_Abs_film(X1,D,d,t, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(1/2-4/pi^2*SumofTerms);
end

% equation (45): Integrated concentration over film during desorption
% assuming constant initial concentration as a function of time
function f=Per_Des_film_const(X1,D,d,t, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(8/(pi^2)*SumofTerms);
end

% equation (46): Integrated concentration over film during desorption
% assuming linear concentration gradient as a function of time
function f=Per_des_film_grad(X1,D,d,t, num_terms)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(4/(pi^2)*SumofTerms);
end

