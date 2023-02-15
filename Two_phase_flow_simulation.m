function two_phase_flow_simulator
global poi boi bwi krwi kroi swi ti plots;
global po_xt sw_xt Nt Nx dp max_diff xi A;
global max_iter lambda_w lambda_o info L q_r;
global q_c pvs dt swc sor phi q_w_0;
global bmin bmax dx loop_plot qrmin qrmax;
global pmin pmax;

po_xt(1:Nx,1)=poi;
sw_xt(1:Nx,1)=swi;
pp=zeros(Nx-1,1);
pn=zeros(Nx-1,1);
s_r=zeros(Nx,1);

df=zeros(Nx-1,Nx-1);

if(plots==1 && loop_plot~=0)
     figure;
end

for j=2:Nt
    f(1:Nx,1)=max_diff+1;
    krwi=krw(swi);
    kroi=kro(swi);
    poi_t=poi;
    it=1;
    while(max(abs(f))>max_diff && it<max_iter)
        boi=bo(poi);
        bwi=bw(poi);
        lambda_o=lambdao(poi);
        lambda_w=lambdaw(poi);
        f=f_calc(poi(1:Nx-1,1),poi_t(1:Nx-1,1));
        for m=1:Nx-1
            for i=1:Nx-1
                if(i==m)
                    pp(i,1)=poi(i,1)+dp;
                    pn(i,1)=poi(i,1)-dp;
                else
                    pp(i,1)=0;
                    pn(i,1)=0;
                end
            end
            fp=f_calc(pp,poi_t(1:Nx-1,1));
            fm=f_calc(pn,poi_t(1:Nx-1,1));
            df(1:Nx-1,m)=(fp-fm)/(2*dp);
        end
        dpt=-1*tdma(df,f);
        poi(1:Nx-1)=poi(1:Nx-1)+dpt;
        if(info)
            disp(['time ' num2str(ti(j)) ', iteration=' num2str(it) ', change in p ' num2str(max(abs(dpt)))]);
        end
        it=it+1;
    end
    if(it==max_iter)
        disp('max iterations exceeded.');
        keyboard;
    end
    
    swi_t=swi;
    swi=sw_calc(swi_t,poi,poi_t);
    
    qr=q_prod(poi,poi_t,swi,swi_t);
    q_r(j,1:2)=qr;
    q_c(j,1)=q_c(j-1,1)+dt*q_r(j,1);
    q_c(j,2)=q_c(j-1,2)+dt*q_r(j,2);
    
    s_r(1,1)=(1-swi(1,1)-sor)*dx/2*A*phi;
    for i=2:Nx-1
        s_r(i,1)=s_r(i-1,1)+(1-swi(i,1)-sor)*dx*A*phi;
    end
    s_r(Nx,1)=s_r(Nx-1,1)+(1-swi(Nx,1)-sor)*dx/2*A*phi;
    
    if(plots==1 && loop_plot~=0)
        if(j/loop_plot==round(j/loop_plot))
            subplot(2,2,1);plot(xi,poi);title(['Pressure @ t=' num2str(ti(j)/60,3) ' minutes']);xlabel('x (m)');ylabel('P_o (Pa)');ylim([pmin pmax]);
            subplot(2,2,2);plot(xi,bwi,xi,boi);title('B factor');xlabel('x (m)');ylabel('B');legend('water','oil','Location','SouthEast');ylim([bmin bmax]);
            subplot(2,2,3);plot(xi,swi,xi,1-swi);title('Saturation profile');xlabel('x (m)');ylabel('S');ylim([0 1]);
            subplot(2,2,4);plot(xi,s_r);title('Cummulative movable oil');xlabel('x (m)');ylabel('movable oil before x (m^3)');xlim([0 L]);ylim([0 phi*A*L*(1-sor-swc)]);
            drawnow;
        end
    end
    
    po_xt(1:Nx,j)=poi;
    sw_xt(1:Nx,j)=swi;
end

if(plots==1)
    figure;
    subplot(2,2,1:2);
    plot(ti/60,q_c(:,1),ti/60,q_c(:,2),ti/60,q_w_0*ti);
    title('Cummulative production profile');
    xlabel('time (min)');
    ylabel('Volume (m^3)');
    ylim([0 pvs*(1-sor-swc)*phi*A*L]);
    legend('water','oil','injected water','Location','EastOutside');
    subplot(2,2,3);
    plot(ti/60,po_xt(1,:),ti/60,po_xt(round(Nx/2),:),ti/60,po_xt(Nx,:));
    title('Pressure profile');
    xlabel('time (min)');
    ylabel('P (Pa)');
    legend('x=0','x=L/2','x=L','Location','NorthEast');
    ylim([pmin pmax]);
    subplot(2,2,4);
    plot(ti/60,q_r(:,1)*60,ti/60,q_r(:,2)*60);
    title('Production rate');
    xlabel('time (min)');
    ylabel('q (m^3/min)');
    legend('water','oil','Location','NorthEast');
    ylim([qrmin qrmax]);
end
    
function [ bo_t ] = bo(poi)
global c_o po_0 b_o_0
bo_t=b_o_0*exp(-1*c_o*(poi-po_0));

function [ bw_t ] = bw(poi)
global c_w po_0 b_w_0;
bw_t=b_w_0*exp(-1*c_w*(poi-po_0));

function [ kro ] = kro(swi)
global sor swc no kro0 Nx;
kro=zeros(Nx,1);
for i=1:Nx
    if(swi(i,1)==0)
        kro(i,1)=1;
    elseif(swi(i,1)>0 && swi(i,1)<swc)
        kro(i,1)=kro0;
    elseif(swi(i,1)>=1-sor)
        kro(i,1)=0;
    else
        kro(i,1)=kro0*(1-(swi(i,1)-swc)/(1-sor-swc))^no;
    end
end

function [ krw ] = krw(swi)
global sor swc nw krw0 Nx;
krw=zeros(Nx,1);
for i=1:Nx
    if(swi(i,1)<=swc)
        krw(i,1)=0;
    elseif(swi(i,1)==1)
        krw(i,1)=1;
    elseif(swi(i,1)>1-sor && swi(i,1)<1)
        krw(i,1)=krw0;
    else
        krw(i,1)=krw0*((swi(i,1)-swc)/(1-sor-swc))^nw;
    end
end

function [ lo ] = lambdao(p)
global kroi boi mu_o ki Nx;
lo=zeros(Nx+1,1);
lo(1,1)=(ki(1,1)*kroi(1,1))/(mu_o*boi(1,1));
for i=2:Nx
    if(p(i-1,1)>=p(i,1))
        lo(i,1)=2*ki(i-1,1)*ki(i,1)/(ki(i-1,1)+ki(i,1))*kroi(i-1,1)/(mu_o*boi(i-1,1));
    else
        lo(i,1)=2*ki(i-1,1)*ki(i,1)/(ki(i-1,1)+ki(i,1))*kroi(i,1)/(mu_o*boi(i,1));
    end
end
lo(Nx+1,1)=(ki(Nx,1)*kroi(Nx,1))/(mu_o*boi(Nx,1));

function [ lw ] = lambdaw(p)
global krwi bwi mu_w ki Nx;
lw=zeros(Nx+1,1);
lw(1,1)=(ki(1,1)*krwi(1,1))/(mu_w*bwi(1,1));
for i=2:Nx
    if(p(i-1,1)>=p(i,1))
        lw(i,1)=2*ki(i-1,1)*ki(i,1)/(ki(i-1,1)+ki(i,1))*krwi(i-1,1)/(mu_w*bwi(i-1,1));
    else
        lw(i,1)=2*ki(i-1,1)*ki(i,1)/(ki(i-1,1)+ki(i,1))*krwi(i,1)/(mu_w*bwi(i,1));
    end
end
lw(Nx+1,1)=(ki(Nx,1)*krwi(Nx,1))/(mu_w*bwi(Nx,1));

function [ f ] = f_calc(pn,po)
global boi bwi phi dx dt swi c_w A Nx po_0;
global lambda_w lambda_o c_o q_w_0 q_o_0;
D(1,1)=phi*((1-swi(1,1))*c_o+swi(1,1)*c_w)/dt;
T(1,1)=-2/dx^2*(bwi(1,1)*lambda_w(2,1)+boi(1,1)*lambda_o(2,1));
T(1,2)=2/dx^2*(bwi(1,1)*lambda_w(2,1)+boi(1,1)*lambda_o(2,1));
c(1,1)=-2/(A*dx)*(bwi(1,1)*q_w_0 + boi(1,1)*q_o_0);
for i=2:(Nx-1)
    D(i,i)=phi*((1-swi(i,1))*c_o+swi(i,1)*c_w)/dt;
    T(i,i-1)=(boi(i,1)*lambda_o(i,1)+bwi(i,1)*lambda_w(i,1))/dx^2;
    T(i,i)=-(boi(i,1)*(lambda_o(i,1)+lambda_o(i+1,1))+bwi(i,1)*(lambda_w(i,1)+lambda_w(i+1,1)))/dx^2;
    if(i==Nx-1)
        c(i,1)=-(boi(i,1)*lambda_o(i+1,1)+bwi(i,1)*lambda_w(i+1,1))/dx^2*po_0;
    else
        T(i,i+1)=(boi(i,1)*lambda_o(i+1,1)+bwi(i,1)*lambda_w(i+1,1))/dx^2;
        c(i,1)=0;
    end
end
f=T*pn-D*(pn-po)-c;

function [ swn ] = sw_calc(swo,pn,po)
global phi swc sor boi bwi c_w dx dt q_o_0 q_w_0 lambda_w Nx;
swn(1,1)=swc+(1-swc-sor)*(q_w_0*bwi(1,1))/(q_w_0*bwi(1,1)+q_o_0*boi(1,1));
for i=2:Nx-1
    swn(i,1)=1/(1+c_w*(pn(i,1)-po(i,1)))*(swo(i,1)+bwi(i,1)*dt/(phi*dx^2)*(lambda_w(i,1)*pn(i-1,1)-(lambda_w(i,1)+lambda_w(i+1,1))*pn(i,1)+lambda_w(i+1,1)*pn(i+1,1)));
end
swn(Nx,1)=1/(1+c_w*(pn(Nx,1)-po(Nx,1)))*(swo(Nx,1)+bwi(Nx,1)*dt/phi/dx^2*((3*lambda_w(Nx+1,1)-2*lambda_w(Nx,1))*pn(Nx,1)+(2*lambda_w(Nx,1)-4*lambda_w(Nx+1,1))*pn(Nx-1,1)+lambda_w(Nx+1,1)*pn(Nx-2,1)));

function [ qr ] = q_prod(pn,po,swn,swo)
global dx lambda_w lambda_o A Nx dt phi bwi boi;
qr(1)=-A*(lambda_w(Nx,1)*(pn(Nx,1)-pn(Nx-1,1))/dx + dx*phi/(2*bwi(Nx,1)*dt)*(swn(Nx,1)-swo(Nx,1)));
qr(2)=-A*(lambda_o(Nx,1)*(pn(Nx,1)-pn(Nx-1,1))/dx + dx*phi/(2*boi(Nx,1)*dt)*((1-swn(Nx,1))-(1-swo(Nx,1))));

function [x] = tdma(A,di)
N=size(A,1);
c=zeros(N,1);
d=zeros(N,1);
x=zeros(N,1);
if(size(di,1)==N)
    for i=1:N
        if(i==1)
            c(i)=A(i,i+1)/A(i,i);
            d(i)=di(i)/A(i,i);
        else
            if(i<N)
                c(i)=A(i,i+1)/(A(i,i)-c(i-1)*A(i,i-1));
            end
            d(i)=(di(i)-d(i-1)*A(i,i-1))/(A(i,i)-c(i-1)*A(i,i-1));
        end
    end
    for i=1:(N-1)
        i=N+1-i; %#ok<FXSET>
        if(i==N)
            x(i)=d(i);
        end
        x(i-1)=d(i-1)-c(i-1)*x(i);
    end
else
    disp('matrix dimensions must match');
    keyboard;
end