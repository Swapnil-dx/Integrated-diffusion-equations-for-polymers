%% Supplementary material: Exact integrated equations to describe diffusion kinetics
% Swapnil Daxini, Jack A. Barnes, and Hans-Peter Loock

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

num_terms = 50;

%% *************** IMPERMEABLE FILM ***************
%% Concentration gradients for absorption and desorption
% equation (10): Concentration gradient across impermeable membrane during
% absorption
function f=Imp_Abs_conc(X1,D,d,t,y)
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
function f=Imp_Des_conc(X1,D,d,t,y)
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
function f=Imp_Abs_ev(X1,D,d,t,dy)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)= (-1)^j/(2*j+1)^2*sin((2*j+1)*pi*dy/(2*d))*exp(-D.*(2*j+1)^2*pi^2*t/(4*d.^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(1-8*d/(pi^2*dy)*SumofTerms);
end

% equation (17): Integrated concentration through evanescent wave during
% desorption as a function of time
function f=Imp_Des_ev(X1,D,d,t,dy)
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

% equation (18): Integrated concentration in film during absorption as a
% function of time
function f=Imp_Abs_film(X1,D,d,t)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)= 1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/(4*d.^2));
        end
    SumofTerms=sum(Terms);
    f=X1*(1-8/pi^2*SumofTerms);
end

% equation (19): Integrated concentration in film during desorption as a
% function of time
function f=Imp_Des_film(X1,D,d,t)
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
% equation (22): Concentration gradient across permeable membrane during
% absorption
function f=Per_Abs_conc(X1,D,d,t,y)
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

% equation (24): Concentration gradient across permeable membrane during
% desorption assuming linear concentration gradient
function f=Per_Des_conc_grad(X1,D,d,t,y)
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

% equation (25): Concentration gradient across permeable membrane during
% desorption assuming constant initial concentration
function f=Per_Des_conc_const(X1,D,d,t,y)
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
% equation (36): Integrated concentration through evanescent wave during
% absorption as a function of time
function f=Per_Abs_ev(X1,D,d,t,dy)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=1:num_terms
        i=i+1;
            Terms(i,:)=(-1)^j/j^2*(cos(j*pi*dy/d)-1)*exp(-D.*j^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(dy/(d*2)-2*d/(dy*pi^2)*SumofTerms);
end

% equation (38): Integrated concentration through evanescent wave during
% desorption assuming constant initial concentration as a function of time
function f=Per_Des_ev_const(X1,D,d,t,dy)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*(1-cos((2*j+1)*pi*dy/d))*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(4*d/(dy*pi^2)*SumofTerms);
end

% equation (39): Integrated concentration through evanescent wave during
% desorption assuming linear concentration gradient as a function of time
function f=Per_Des_ev_grad(X1,D,d,t,dy)
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

% equation (41): Integrated concentration over film during absorption as a
% function of time
function f=Per_Abs_film(X1,D,d,t)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(1/2-4/pi^2*SumofTerms);
end

% equation (42): Integrated concentration over film during desorption
% assuming constant initial concentration as a function of time
function f=Per_Des_film_const(X1,D,d,t)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(8/(pi^2)*SumofTerms);
end

% equation (43): Integrated concentration over film during desorption
% assuming linear concentration gradient as a function of time
function f=Per_des_film_grad(X1,D,d,t)
    i=0;
    Terms(num_terms:length(t))=zeros;
        for j=0:num_terms
        i=i+1;
            Terms(i,:)=1/(2*j+1)^2*exp(-D.*(2*j+1)^2*pi^2*t/d^2);
        end
    SumofTerms=sum(Terms);
    f=X1*(4/(pi^2)*SumofTerms);
end
