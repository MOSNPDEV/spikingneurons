
%%%%Time Constant of the gKd channel from Hemond et al 2008
clear
    celsius = 25;	
	gkdbar=.0 ;
    vhalfn=-33 ;  
    a0n=0.01  ;    
    zetan=3  ;  
    gmn=0.7 ;
	nmax=2  ;
	q10=1;
	sh = 0;
    
path='/home/saray/Desktop/taus/';

 ntau=[];   
 infn=[];
 nbet=[];
 nalp=[];
V=-100:1:50;

for i=1:length(V)
    
    v=V(i);

qt=q10^((celsius-24)/10);

alpn = exp(1e-3*zetan*(v-vhalfn-sh)*9.648e4/(8.315*(273.16+celsius)));

betn = exp(1.e-3*zetan*gmn*(v-vhalfn-sh)*9.648e4/(8.315*(273.16+celsius))); 

a = alpn;

nalp=[nalp,alpn];

nbet=[nbet,betn];

ninf=1/(1+a);

infn=[infn,ninf];

taun = betn/(qt*a0n*(1+a));

ntau=[ntau,taun];

end


figure
plot(V,ntau,'o')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Time Constant (ms)','FontSize',20)
title(' Tau-n gKd','FontSize',20)
saveas(gcf,[path,'kdtaun'],'jpg');


figure
plot(V,infn,'ro')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Inactivation variable n','FontSize',20)
title(' n-inf gKd','FontSize',20)
saveas(gcf,[path,'kdn'],'jpg');


figure
plot(V,nalp,'o')
hold on
plot(V,nbet,'ro')
legend('alpha','beta')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Alpha and Beta','FontSize',20)
title(' Tau-n gKd','FontSize',20)
saveas(gcf,[path,'kdtaun'],'jpg');



   






