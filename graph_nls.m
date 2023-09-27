% Modeling of the stationary solution of the Nonlinear cubic 
% Schrodinger equation on the Y-junction graph. 
% Author: Adilbek Kairzhan


function nls_1
  
echo off, clc
clear all

global Nsave; global NSTART; global NADVANCE; global M; 
global DEALIAS;
Nsave    = 50;     % save solution after Nsave time steps
NSTART   = 0;
NADVANCE = 5000;
M = 60;            % Each edge of the graph is bounded by [0, M]
DEALIAS  = 1;      % should be 1 for dealising of the error 

Nmin = 9;          % Minimum number of collocation points
Nmax = 129;        % Maximum number of collocation points
dt = 1.0e-6;       % time step

N2 = Nmin;         
i = 1;

while (N2 <= Nmax)

X(1:N2) = cos( (0:N2-1) * pi / (N2-1)  ); % Gauss-Lobatto-Chebyshev points

% introducing the initial conditions
uIC = sqrt(2) .* sech( M .* X ./ 2 + M ./ 2 ); 
vIC = sqrt(2) .* sech( M .* X ./ 2 + M ./ 2 );
wIC = sqrt(2) .* sech( M .* X ./ 2 + M ./ 2 ); 

% finding the solutions and relative errors
[u, v, w, Err1, Err2, Err3] = Solve_NLS( uIC, vIC, wIC, N2, dt);
Err(i,1) = N2;
Err(i,2) = max([Err1,Err2,Err3]);

% plot the soluions and relative errors for Nmax collocation points
if (N2 == Nmax)
figure(1); surf(X,(1:NADVANCE/Nsave)*dt,abs(u));
title('Solution graph for u(x,t) as a function of x and time');
xlabel('x'); ylabel('time');
figure(2); surf(X,(1:NADVANCE/Nsave)*dt,abs(v));
title('Solution graph for v(x,t) as a function of x and time');
xlabel('x'); ylabel('time');
figure(3); surf(X,(1:NADVANCE/Nsave)*dt,abs(w)); 
title('Solution graph for w(x,t) as a function of x and time');
xlabel('x'); ylabel('time');
figure(4); plot(1:NADVANCE/Nsave,Err1(:));
title('L-inf of Relative error of u w.r.t. n');
xlabel('n time steps'); ylabel('RelErr');
figure(5); plot(1:NADVANCE/Nsave,Err2(:));
title('L-inf of Relative error of v w.r.t. n');
xlabel('n time steps'); ylabel('RelErr');
figure(6); plot(1:NADVANCE/Nsave,Err3(:));
title('L-inf of Relative error of w w.r.t. n');
xlabel('n time steps'); ylabel('RelErr');
end;
N2 = 2*N2-1;
i = i+1;
end;

%plot the values of relative errors for different numbers of collocation
%points
figure(7); loglog(Err(:,1),Err(:,2),'-o');
title('Max Relative error of solution w.r.t. number of collocation points');
xlabel('N points'); ylabel('RelErr');



% The function solves the NLS equation on the graph;
% Time-stepping is performed using first order explicit/implicit 
% Euler scheme
function [usave, vsave, wsave, Err1, Err2, Err3] = Solve_NLS( u0, v0, w0, N2, dt)
% INPUT
% u0,v0,w0  - initial conditions
% N2  - number of grid point (extended grid)
% dt  - time-step;
% OUTPUT
% u,v,w   - solutions saved every Nsave steps
% Err1, Err2, Err3 - appropriate relative errors
 
global DEALIAS; global Nsave; global NSTART; global NADVANCE; global M;

Err1 = 0;
Err2 = 0;
Err3 = 0;
% Determining the number of grid points
if ( DEALIAS )
  N = round( N2 / 2);
else
  N = N2;
end; 

coeff = M * M * j ./ dt;

% Computing the differentiation matrix
  D = cheb( N2-1 );
  D2 = 4 * D * D;
for i=2:N2-1;
    D2(i,i) = coeff + D2(i,i);
end;

% Constructing the 3Nx3N matrix for algebraic system
  D3(1:3*N2-3, 1:3*N2-3) = zeros(3*N2-3, 3*N2-3);
  %D3(1,1) = 1;
  D3(1:N2-2, 1:N2-1) = D2(2:N2-1,2:N2);
  D3(N2-1, N2-1) = 1; D3(N2-1, 2*N2-2) = -1;
  %D3(N2+1, N2+1) = 1;
  D3(N2:2*N2-3, N2:2*N2-2) = D2(2:N2-1, 2:N2);
  D3(2*N2-2, 2*N2-2) = 1; D3(2*N2-2, 3*N2-3) = -1;
  %D3(2*N2+1, 2*N2+1) = 1;
  D3(2*N2-1:3*N2-4, 2*N2-1:3*N2-3) = D2(2:N2-1, 2:N2);
  D3(3*N2-3, 1:N2-1) = D(N2, 2:N2);
  D3(3*N2-3, N2:2*N2-2) = D(N2, 2:N2);
  D3(3*N2-3, 2*N2-1:3*N2-3) = D(N2, 2:N2);
  
for i=NSTART:NADVANCE;
  
  if ( i == NSTART )
% at the initial time-step;   
    u(1:N2) = u0(1:N2);
    v(1:N2) = v0(1:N2);
    w(1:N2) = w0(1:N2);
  else    
% pseudospectral calculation of the nonlinear products;
    uuu(1:N2) =  u(1:N2) .* conj(u(1:N2)) .* u(1:N2);
    vvv(1:N2) =  v(1:N2) .* conj(v(1:N2)) .* v(1:N2);
    www(1:N2) =  w(1:N2) .* conj(w(1:N2)) .* w(1:N2);
    uuubar = Transform2Chebyshev( uuu, N, N2 );
    uuu(1:N2) = Transform2Complex(uuubar, N, N2);
    vvvbar = Transform2Chebyshev( vvv, N, N2 );
    vvv(1:N2) = Transform2Complex(vvvbar, N, N2);
    wwwbar = Transform2Chebyshev( www, N, N2 );
    www(1:N2) = Transform2Complex(wwwbar, N, N2);
% RHS vector of the algebraic system
    vector(1) = 0;
    vector(1:N2-2) = coeff .* u(2:N2-1) - M * M .* uuu(2:N2-1);
    vector(N2-1) = 0;
    vector(N2:2*N2-3) = coeff .* v(2:N2-1) - M * M .* vvv(2:N2-1);
    vector(2*N2-2) = 0;
    vector(2*N2-1:3*N2-4) = coeff .* w(2:N2-1) - M * M .* www(2:N2-1);
    vector(3*N2-3) = 0;

% Actual integration in time;
    vector(1:3*N2-3) = D3\transpose(vector);
% separation of answers on the edge solutions
    u(1) = 0;
    u(2:N2) = vector(1:N2-1);
    v(1) = 0;
    v(2:N2) = vector(N2:2*N2-2);
    w(1) = 0;
    w(2:N2) = vector(2*N2-1:3*N2-3);
  end;

% Saving data if necessary;
  if ( mod( i, Nsave ) )
    usave(floor(i/Nsave)+1, 1:N2 ) = u(1:N2);
    vsave(floor(i/Nsave)+1, 1:N2 ) = v(1:N2);
    wsave(floor(i/Nsave)+1, 1:N2 ) = w(1:N2);
% Relative errors   
    Err1(floor(i/Nsave)+1) = norm(abs(u(2:N2)/u0(2:N2)) - 1, inf);
    Err2(floor(i/Nsave)+1) = norm(abs(v(2:N2)/v0(2:N2)) - 1, inf);
    Err3(floor(i/Nsave)+1) = norm(abs(w(2:N2)/w0(2:N2)) - 1, inf);
    end;

end;
  
return;


%-------------------------------------------------------------
% The function performs the Chebyshev to Complex space transform;  
function [ur] = Transform2Complex (v, n, n2);
% INPUT
% v  - input vector (complex components);
% n  - length on the padded grid;
% n2 - length of the input vector (may be on padded grid);
% OUTPUT
% ur - the transformed complex vector;
utmp(1:n) = v(1:n);
x(1:n2) = cos( (0:n2-1) * pi ./ (n2-1)  );

for i=1:n2
T(1:n) = cos((0:n-1) * acos(x(i)));
ur(i) = utmp * transpose(T);
end

return;


%-------------------------------------------------------------
% The function carries out the complex to Chebyshev transform and
% removes the aliased elements if necessary;
function [Uch] = Transform2Chebyshev (v,n,n2);
% INPUT
% v - the input complex vector;
% n2 - its length (on the padded grid - N2);
% n  - size of the output vector (may be dealiased);
% OUTPUT
% Uch - the transformed vector of the Chebyshev coefficients;
% Checking of the number of grid points is consistent with FFT;
  if ( ~isposint(log2( 2 * (n2-1) )) )
    disp(' Wrong number of grid-points!')
    return;
  end;
  U1(1:n2) = v(1:n2)';
% Calculating the cosine transform
  U2 = fft([U1'; flipud(U1(2:n2-1)')]);
% Extracting the Chebyshev coefficients
  Uch(1:n) = U2(1:n).*[0.5; ones(n-2,1); 0.5] / (n2-1);
% Now padding zeros;
if ( n ~= n2 )
  Uch(n+1:n2) = 0.0;
end;  

return;


%-------------------------------------------------------------
% CHEB  compute D = differentiation matrix, x = Chebyshev grid
% Implemented by N.L. Trefethen

function [D,x] = cheb(N)
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries
  
  
%-------------------------------------------------------------
%ISPOSINT True for positive integer values. 
%       Mark Beale, 11-31-97
%       Copyright 1992-2001 The MathWorks, Inc.
% $Revision: 1.5 $
function f=isposint(v)
 
f = 1;
if ~isa(v,'double') | any(size(v) ~= [1 1]) | ...
  ~isreal(v) | v<0 | round(v) ~= v
  f = 0;
end