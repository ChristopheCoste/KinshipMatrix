%% inputs
%components of the projection matrix
F =[ 0.2300    0.5100    0.8400;  0   0   0;   0 0  0]
S = [ 0  0  0;  0.5700  0   0;   0    0.5400    0.4600]
%kinship distance investigated
gmax=3

% same-litter sisters
%Z= ?   %same-litter-sister matrix, enter it here (and uncomment) if there is more than one class of newborn and reproduction is neither Bernoulli nor Poisson
%particular cases:
Fstar=F; %if reproduction is Poisson
%Fstar=zeros(size(F)) ; %if reproduction is Bernoulli
%Fstar= ? % if reproduction is neither Poisson nor Bernoulli and there is sonly one newborn class

%% eigen-analysis
M=F+S %projection matrix
[w,lam]=eigs(M,1); %eigen analysis
lam % asymptotic growth rate
w=w/sum(w) %asymptotic relative abundances

%% same litter sisters 
Dw=diag(w);
Zdagger=Fstar*Dw*F'; % uncomment if in one of the particular cases
%Zdagger=lam*Z*Dw; % uncomment if not in one of the particular cases

%%  Generation structured matrices/tensors
s=size(F,1);

J1=diag(ones(1,gmax-1),-1);
J2=zeros(gmax);J2(2,2)=1;
Jy=kron(eye(gmax),ones(1,s));
J3=zeros(gmax,gmax);J3(1)=1;
J=kron(J3,eye(s));

bigF=kron(J1,F);
bigS=kron(eye(gmax),S);
bigM=bigF+bigS;

bigZdagger=kron(J2,Zdagger);

bigDw=kron(eye(gmax),Dw);

%% Direct Computation of kinship matrix 

vecKw=((lam*eye(s*s*gmax*gmax)- kron(bigM,bigM))^(-1)) * ((kron(bigS,bigF)+kron(bigF,bigS))*vec(J*bigDw) + vec(bigZdagger));
Kw=reshape(vecKw,[gmax*s,gmax*s]);
K=reshape(vecKw,[gmax*s,gmax*s])*(bigDw^(-1));
Ku=reshape(kron(Jy,Jy)*vecKw,[gmax,gmax]);

%% Displaying the kinship matrix and the unstructured kinship matrix 
format shortg
round(K,2)
round(Ku,2)

%% function
function [y] = vec(x)
    y=x(:);
end
