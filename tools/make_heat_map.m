function make_heat_map(A, B, T, Nx, Nu, R, M, loc, Tmax, myTitle, openLoop)
% makes heat maps for x and u according to the specified system
%   A,B      : system matrices
%   T        : finite impulse response horizon
%   Nx,Nu    : dimensions of x and u
%   R,M      : system responses of noise-to-state and noise-to-input
%   loc      : where the disturbance hits
%   Tmax     : number of steps to simulate for
%   myTitle  : overall title of the heat maps
%   openLoop : 1 if we want to simualte open loop, else 0; defaults to 0

% TODO: myTitle is actually only used if closed loop
% TODO: refactor so can be used for output feedback as well

if nargin == 10
    openLoop = 0;
end

w_d    = zeros(Nx,Tmax);
Tstart = T+1;
B1     = eye(Nx);

w_d(loc,Tstart) = 10;

x     = zeros(Nx,Tmax);
u     = zeros(Nu,Tmax);
x_ref = zeros(Nx,Tmax);
w_est = zeros(Nx,Tmax);

MATx = []; MATu = [];

for i=Tstart:1:Tmax-1

    w_est(:,i-1) = x(:,i) - x_ref(:,i);
    
    for jj=1:1:T
        if (openLoop==1)
           u(:,i) = zeros(Nu,1);
        else
           u(:,i) = u(:,i) + M{jj}*w_est(:,i-jj);
        end
    end
    
    x(:,i+1) = A*x(:,i) + B1*w_d(:,i)+ B*u(:,i);
    
    for jj=2:1:T
        if (openLoop==1)
           x_ref(:,i+1) = x_ref(:,i+1); 
        else
           x_ref(:,i+1) = x_ref(:,i+1) + R{jj}*w_est(:,i+1-jj);
        end
    end

    MATx = [MATx,log10(abs(x(:,i)))];
    MATu = [MATu,log10(abs(B*u(:,i)))];
    
end

if (openLoop == 0)
   figure; suptitle(myTitle); 
   subplot(1,2,1); heat_map_plotter(MATx, 'log10(|x|)', 'Time', 'Space');
   subplot(1,2,2); heat_map_plotter(MATu, 'log10(|u|)', 'Time', '');

else
   figure; heat_map_plotter(MATx, 'Open Loop: log10(|x|)', 'Time', 'Space');
end