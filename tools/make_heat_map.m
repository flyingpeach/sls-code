function make_heat_map(A,B,T,Nx,Nu,R,M,loc,Tmax,the_title,open_loop);

if nargin == 10
    open_loop = 0;
end

w_d = zeros(Nx,Tmax);
Tstart = T+1;
B1 = eye(Nx);
w_d(loc,Tstart) = 10;



x=zeros(Nx,Tmax);
u=zeros(Nu,Tmax);
x_ref=zeros(Nx,Tmax);
w_est=zeros(Nx,Tmax);


MATx = [];
MATu = [];
for i=Tstart:1:Tmax-1

    w_est(:,i-1) = x(:,i) - x_ref(:,i);
    
    for jj=1:1:T
        if (open_loop==1)
           u(:,i) = zeros(Nu,1);
        else
           u(:,i) = u(:,i) + M{jj}*w_est(:,i-jj);
        end
    end
    x(:,i+1) = A*x(:,i) + B1*w_d(:,i)+ B*u(:,i);
    for jj=2:1:T
        if (open_loop==1)
           x_ref(:,i+1) = x_ref(:,i+1); 
        else
           x_ref(:,i+1) = x_ref(:,i+1) + R{jj}*w_est(:,i+1-jj);
        end
    end
    

    MATx = [MATx,log10(abs(x(:,i)))];
    MATu = [MATu,log10(abs(B*u(:,i)))];
    
end

if (open_loop == 0)
   figure
   suptitle(the_title)
   subplot(1,2,1);
   handle = imagesc(((MATx)));
   colorbar
   colormap jet
   caxis([-4 0])
   title('log10(|x|)')
   xlabel('Time')
   ylabel('Space')
   set(gca,'FontSize',16,'fontWeight','bold')
   set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

   subplot(1,2,2);
   handle = imagesc(((MATu)));
   colorbar
   colormap jet
   caxis([-4 0])
   title('log10(|u|)')
   xlabel('Time')
   set(gca,'FontSize',16,'fontWeight','bold')
   set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
else
   figure
   title(the_title)
   handle = imagesc(((MATx)));
   colorbar
   colormap jet
   caxis([-4 0])
   title('Open Loop: log10(|x|)')
   xlabel('Time')

   set(gca,'FontSize',16,'fontWeight','bold')
   set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

end