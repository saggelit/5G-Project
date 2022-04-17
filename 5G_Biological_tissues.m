                %  -------------------------------  %
                %   Computational E/M               %
                %   Angelitsi Sotiria, AEM:4366     %
                %   Dragatsikas Ioannis, AEM:4375   %
                %   Project - 5G                    %
                %  -------------------------------  %
            
clear; clc; clf;

j=sqrt(-1);
L=80;
N=100;     % Number of Cells
c=3e8;     % Speed of light
dx=10^-3;
dt=1.7e-12;
x=linspace(0,L,N);
updategamma=[1*10^-9 0.5*10^-9 0.1*10^-9 5.0*10^-9];  %bandwidth
    
% Free space
mi0=4.0*pi*1.0e-7;      % Permeability
epsilon0=1.0/(c*c*mi0); % Permittivity

%% Case ~0~
sigma=0.0*ones(1,N); %s --> 0.0
sigma(5)=0.540;sigma(6:40)=0.037;sigma(41:75)=0.747;sigma(76:N)=0.104;
epsilon_s=1.0*ones(1,N); %epsilon_s -->1.0
epsilon_s(5)=47.9;epsilon_s(6:40)=5.53;epsilon_s(41:75)=56.9;epsilon_s(76:N)=14.2;
epsilon_inf=1.0*ones(1,N); %epsilon_inf --> 1.0
epsilon_inf(5)=29.9;epsilon_inf(6:40)=4.00;epsilon_inf(41:75)=28.0;epsilon_inf(76:N)=7.36;
tau=0.0*ones(1,N); %tau --> 0
tau(5)=43.6*10^-12;tau(6:40)=23.6*10^-12;tau(41:75)=18.7*10^-12;tau(76:N)=34.1*10^-12;
        
mir=ones(N,1);
sim=zeros(N,1);

% Source
updatefreq=[2.0e+9 2.2e+9 2.4e+9 2.6e+9];   % Frequency

% Update Coefficients
for count1=1:8 % From 1-3 frequency changes
    if (count1<=4)
        freq=updatefreq(count1);
        gamma=1e-9;
    elseif (count1<=8)
        freq=2.0e+9;
        count2=count1-4;
        gamma=updategamma(count2);
    end
    
    epsilon_r=epsilon_inf-(epsilon_inf-epsilon_s)./(1+j*2*pi*freq.*tau)-j*sigma./(2*pi*freq*epsilon0);
    
    % Constants
    u1=(2*epsilon0.*epsilon_inf.*tau+2.*(epsilon0.*epsilon_s+sigma.*tau).*dt+sigma.*dt^2).^-1;
    u2=4*epsilon0.*epsilon_inf.*tau+2.*(epsilon0.*epsilon_s+sigma.*tau).*dt-sigma.*dt^2;
    u3=2*epsilon0.*epsilon_inf.*tau;
    u4=2.*(tau+dt);
    u5=2.*(2.*tau+dt);
    u6=2.*tau;
    
    
    Ez0=100.0;        % Initialization
    T=1/freq;
    nmax=17600;      % Total time for computing and animation
   
    
    %skin depth
    delta=c./(2*pi*freq*abs(imag(sqrt(epsilon_r))));
    
    % Initial Conditions
    ca=(1.0-(dt*sigma)./(2.0*epsilon0*real(epsilon_r)))./(1.0+(dt*sigma)./(2.0*epsilon0*real(epsilon_r)));
    cb=(dt/epsilon0./real(epsilon_r)/dx)./(1.0+(dt*sigma)./(2.0*epsilon0*real(epsilon_r)));
    da=(1.0-(dt*sim)./(2.0*mi0*mir))./(1.0+(dt*sim)./(2.0*mi0*mir));
    db=(dt/mi0./mir/dx)./(1.0+(dt*sim)./(2.0*mi0*mir));
    
    Ez=zeros(1,N);
    Eznew=zeros(1,N);
    Dznew=zeros(1,N);
    Ezold=Ez0*exp(-x./delta);
    Dzold=epsilon0.*real(epsilon_r).*Ezold;
    Hy=zeros(1,N);
    Dz=zeros(1,N);    
    emaxz(1:N)=0;
   
    % Time Stepping - Yee Algorithm - 1D FDTD Method

    for k=1:nmax
        Ez(3)=Ez(3)+Ez0*exp(-((dt*k-5*gamma)/gamma)^2)*cos(2*pi*freq*dt*k);
        
        Hy(1)=Hy(2);
        for i=2:N-1
            Hy(i)=da(i)*Hy(i)-db(i)*(Ez(i+1)-Ez(i));
        end
        
        Ez(N)=Ez(N-1);
        for i=2:N-1
             Ez(i)=ca(i)*Ez(i)-cb(i)*(Hy(i)-Hy(i-1));
        end
      
        for i=2:N-1
            Dz(i)= Ez(i)*epsilon0;
        end    
       
        for i=2:N-1
            Dznew(i)= Dz(i)+(dt/dx)*(Hy(i)-Hy(i-1));
        end 

        for i =2:N-1
            Eznew(i)=u1(i)*(u2(i)*Ez(i)-u3(i)*Ezold(i)+u4(i)*Dznew(i)-u5(i)*Dz(i)+u6(i)*Dzold(i));
        end    
        
        emaxz=max(Eznew,emaxz);
    end
    if (count1<=4) %gia allagh syxnothtas kai times gia s kai e apo case 1
        f1=figure(1);
        f1.Name=['Angelitsi-Dragatsikas - Plot:',num2str(1)];
        set(f1,'NumberTitle', 'off');
        plot(x,abs(emaxz));
        title({'The comparison of the electric field magnitude along the path by varying',...
            ' the carrier frequency. The pulse bandwidth is kept constant as 1.0 ns'});
        xlabel('Path (mm)')
        ylabel('Peak of |E_z| (V/m)')
        hold on
        legend({'2.0 GHz','2.2 GHz','2.4 GHz','2.6 GHz'},'location','east')
    else 
        f2=figure(2);
        f2.Name=['Angelitsi-Dragatsikas - Plot:',num2str(2)];
        set(f2,'NumberTitle', 'off');
        plot(x,abs(emaxz));
        title({'The comparison of the electric field magnitude along the path by varying'...
            'the pulse bandwidth. The carrier frequency is kept constant as 2.0 GHz'});
        xlabel('Path (mm)')
        ylabel('Peak of |E_z| (V/m)')
        hold on
        legend({'1.0 ns','0.5 ns','0.1 ns','5.0 ns'},'location','east')

    end
end

%% Case ~1~ & ~2~
for count=1:2
    if (count==1)
        sigma=0.0*ones(1,N); %s --> 0.0
        sigma(5)=0.540;sigma(6:15)=0.037;sigma(16:20)=0.747;sigma(21:52)=0.104;sigma(53:N)=0.747;
        epsilon_s=1.0*ones(1,N); %epsilon_s --> 1.0
        epsilon_s(5)=47.9;epsilon_s(6:15)=5.53;epsilon_s(16:20)=56.9;epsilon_s(21:52)=14.2;epsilon_s(53:N)=56.9;
        epsilon_inf=1.0*ones(1,N); %epsilon_inf --> 1.0
        epsilon_inf(5)=29.9;epsilon_inf(6:15)=4.00;epsilon_inf(16:20)=28.0;epsilon_inf(21:52)=7.36;epsilon_inf(53:N)=28.0;
        tau=0.0*ones(1,N); %tau --> 0.0
        tau(5)=43.6*10^-12;tau(6:15)=23.6*10^-12;tau(16:20)=18.7*10^-12;tau(21:52)=34.1*10^-12;tau(53:N)=18.7*10^-12;
    else
        sigma=0.0*ones(1,N); %s --> 0.0
        sigma(5)=0.540;sigma(6)=0.037;sigma(7:9)=0.747;sigma(10:42)=0.104;sigma(43:N)=0.747;
        epsilon_s=1.0*ones(1,N); %epsilon_s --> 1.0
        epsilon_s(5)=47.9;epsilon_s(6)=5.53;epsilon_s(7:9)=56.9;epsilon_s(10:42)=14.2;epsilon_s(43:N)=56.9;
        epsilon_inf=1.0*ones(1,N); %epsilon_inf --> 1.0
        epsilon_inf(5)=29.9;epsilon_inf(6)=4.00;epsilon_inf(7:9)=28.0;epsilon_inf(10:42)=7.36;epsilon_inf(43:N)=28.0;
        tau=0.0*ones(1,N); %tau --> 0.0
        tau(5)=43.6*10^-12;tau(6)=23.6*10^-12;tau(7:9)=18.7*10^-12;tau(10:42)=34.1*10^-12;tau(43:N)=18.7*10^-12;
    end
    
mir=ones(N,1);
sim=zeros(N,1);

% Update Coefficients
    freq=2.0e+9;
    gamma=1e-9;
    
    epsilon_r=epsilon_inf-(epsilon_inf-epsilon_s)./(1+j*2*pi*freq.*tau)-j*sigma./(2*pi*freq*epsilon0);
    
    % Constants
    u1=(2*epsilon0.*epsilon_inf.*tau+2.*(epsilon0.*epsilon_s+sigma.*tau).*dt+sigma.*dt^2).^-1;
    u2=4*epsilon0.*epsilon_inf.*tau + 2.*(epsilon0.*epsilon_s+sigma.*tau).*dt-sigma.*dt^2;
    u3=2*epsilon0.*epsilon_inf.*tau;
    u4=2.*(tau+dt);
    u5=2.*(2.*tau+dt);
    u6=2.*tau;
    
    Ez0=100.0;        % Initialization
    T=1/freq;
    nmax=17600;      % Total time for computing and animation
    
    %skin depth
    delta=c./(2*pi*freq*abs(imag(sqrt(epsilon_r))));
    
    % Initial Conditions
    ca=(1.0-(dt*sigma)./(2.0*epsilon0*real(epsilon_r)))./(1.0+(dt*sigma)./(2.0*epsilon0*real(epsilon_r)));
    cb=(dt/epsilon0./real(epsilon_r)/dx)./(1.0+(dt*sigma)./(2.0*epsilon0*real(epsilon_r)));
    da=(1.0-(dt*sim)./(2.0*mi0*mir))./(1.0+(dt*sim)./(2.0*mi0*mir));
    db=(dt/mi0./mir/dx)./(1.0+(dt*sim)./(2.0*mi0*mir));
    
    Ez=zeros(1,N);
    Eznew=zeros(1,N);
    Dznew=zeros(1,N);
    Ezold=Ez0*exp(-x./delta);
    Dzold=epsilon0.*real(epsilon_r).*Ezold;
    Hy=zeros(1,N);
    Dz=zeros(1,N);    
    emaxz(1:N)=0;
   
    % Time Stepping - Yee Algorithm - 1D FDTD Method

    for k=1:nmax
        Ez(3)=Ez(3)+Ez0*exp(-((dt*k-5*gamma)/gamma)^2)*cos(2*pi*freq*dt*k);
        
        Hy(1)=Hy(2);
        for i=2:N-1
            Hy(i)=da(i)*Hy(i)-db(i)*(Ez(i+1)-Ez(i));
        end
        
        Ez(N)=Ez(N-1);
        for i=2:N-1
             Ez(i)=ca(i)*Ez(i)-cb(i)*(Hy(i)-Hy(i-1));
        end
      
        for i=1:N-1
            Dz(i)= Ez(i)*epsilon0;
        end    
       
        for i=2:N-1
            Dznew(i)= Dz(i)+(dt/dx)*(Hy(i)-Hy(i-1));
        end 

        for i =2:N-1
            Eznew(i)=u1(i)*(u2(i)*Ez(i)-u3(i)*Ezold(i)+u4(i)*Dznew(i)-u5(i)*Dz(i)+u6(i)*Dzold(i));
        end
         
        emaxz=max(Eznew,emaxz);
    end

    f3=figure(3);
    f3.Name=['Angelitsi-Dragatsikas - Plot:',num2str(3)];
    set(f3,'NumberTitle', 'off');
    plot(x,abs(emaxz));
    title({'The variation of the magnitude of electric field along the path ',...
        'when varying thicknesses of Fat, Muscle, Bone'});
    xlabel('Path (mm)')
    ylabel('Peak of |E_z| (V/m)')
    hold on
    legend({'case 1','case 2'},'location','east')
end

%% Case ~3~
t1=10; t2=5; t3=32;
for t0=1:4
    sigma=0.0*ones(1,N); % s --> 0.0
    sigma(4+1:4+t0)=0.540;sigma(4+t0+1:4+t0+t1)=0.037;sigma(4+t0+t1+1:4+t0+t1+t2)=0.747;sigma(4+t0+t1+t2+1:4+t0+t1+t2+t3)=0.104;sigma(4+t0+t1+t2+t3+1:N)=0.747;
    epsilon_s=1.0*ones(1,N); %epsilon_s --> 1.0
    epsilon_s(4+1:4+t0)=47.9;epsilon_s(4+t0+1:4+t0+t1)=5.53;epsilon_s(4+t0+t1+1:4+t0+t1+t2)=56.9;epsilon_s(4+t0+t1+t2+1:4+t0+t1+t2+t3)=14.2;epsilon_s(4+t0+t1+t2+t3+1:N)=56.9;
    epsilon_inf=1.0*ones(1,N); %epsilon_inf -->1.0 
    epsilon_inf(4+1:4+t0)=29.9;epsilon_inf(4+t0+1:4+t0+t1)=4.00;epsilon_inf(4+t0+t1+1:4+t0+t1+t2)=28.0;epsilon_inf(4+t0+t1+t2+1:4+t0+t1+t2+t3)=7.36;epsilon_inf(4+t0+t1+t2+t3+1:N)=28.0;
    tau=0.0*ones(1,N); %tau --> 0.0
    tau(4+1:4+t0)=43.6*10^-12;tau(4+t0+1:4+t0+t1)=23.6*10^-12;tau(4+t0+t1+1:4+t0+t1+t2)=18.7*10^-12;tau(4+t0+t1+t2+1:4+t0+t1+t2+t3)=34.1*10^-12;tau(4+t0+t1+t2+t3+1:N)=18.7*10^-12;

    mir=ones(N,1);
    sim=zeros(N,1);

    % Update Coefficients
    freq=2.0e+9;
    gamma=1e-9;
    
    epsilon_r=epsilon_inf-(epsilon_inf-epsilon_s)./(1+j*2*pi*freq.*tau)-j*sigma./(2*pi*freq*epsilon0);
    
    % Constants
    u1=(2*epsilon0.*epsilon_inf.*tau+2.*(epsilon0.*epsilon_s+sigma.*tau).*dt+sigma.*dt^2).^-1;
    u2=4*epsilon0.*epsilon_inf.*tau + 2.*(epsilon0.*epsilon_s+sigma.*tau).*dt-sigma.*dt^2;
    u3=2*epsilon0.*epsilon_inf.*tau;
    u4=2.*(tau+dt);
    u5=2.*(2.*tau+dt);
    u6=2.*tau;
    
    Ez0=100.0;        % Initialization
    T=1/freq;
    nmax=17600;      % Total time for computing and animation
   
    %skin depth
    delta=c./(2*pi*freq*abs(imag(sqrt(epsilon_r))));
    
    % Initial Conditions
    ca=(1.0-(dt*sigma)./(2.0*epsilon0*real(epsilon_r)))./(1.0+(dt*sigma)./(2.0*epsilon0*real(epsilon_r)));
    cb=(dt/epsilon0./real(epsilon_r)/dx)./(1.0+(dt*sigma)./(2.0*epsilon0*real(epsilon_r)));
    da=(1.0-(dt*sim)./(2.0*mi0*mir))./(1.0+(dt*sim)./(2.0*mi0*mir));
    db=(dt/mi0./mir/dx)./(1.0+(dt*sim)./(2.0*mi0*mir));
    
    Ez=zeros(1,N);
    Eznew=zeros(1,N);
    Dznew=zeros(1,N);
    Ezold=Ez0*exp(-x./delta);
    Dzold=epsilon0.*real(epsilon_r).*Ezold;
    Hy=zeros(1,N);
    Dz=zeros(1,N);    
    emaxz(1:N)=0;
   
    % Time Stepping - Yee Algorithm - 1D FDTD Method

    for k=1:nmax
        Ez(3)=Ez(3)+Ez0*exp(-((dt*k-5*gamma)/gamma)^2)*cos(2*pi*freq*dt*k);
        
        Hy(1)=Hy(2);
        for i=2:N-1
            Hy(i)=da(i)*Hy(i)-db(i)*(Ez(i+1)-Ez(i));
        end
        
        Ez(N)=Ez(N-1);
        for i=2:N-1
             Ez(i)=ca(i)*Ez(i)-cb(i)*(Hy(i)-Hy(i-1));
        end
      
        for i=1:N-1
            Dz(i)= Ez(i)*epsilon0;
        end    
       
        for i=2:N-1
            Dznew(i)= Dz(i)+(dt/dx)*(Hy(i)-Hy(i-1));
        end 

        for i =2:N-1
            Eznew(i)=u1(i)*(u2(i)*Ez(i)-u3(i)*Ezold(i)+u4(i)*Dznew(i)-u5(i)*Dz(i)+u6(i)*Dzold(i));
        end
         
        emaxz=max(Eznew,emaxz);
    end

    f4=figure(4);
    f4.Name=['Angelitsi-Dragatsikas - Plot:',num2str(4)];
    set(f4,'NumberTitle', 'off');
    plot(x,abs(emaxz));
    title({'Comparison of the field magnitude when the thickness of the',...
        'first leyer (Skin) in case 1 was increased from 1mm to 4 mm'});
    xlabel('Path (mm)')
    ylabel('Peak of |E_z| (V/m)')
    hold on
    legend({'1 mm','2 mm','3 mm','4 mm'},'location','east')
end
