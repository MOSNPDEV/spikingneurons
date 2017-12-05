%%%%Time Constant of the gKm channel from Hemond et al 2008
clear
    celsius = 25;	
	gbar=.0001; 
    vhalfl=-40;   	
	kl=-10;
    vhalft=-42 ;  	
     a0t=0.003 ;     	
     zetat=7  ;  
    gmt=.4  ; 	
	q10=5;
	b0=60;
	st=1;
	sh =24;
    
path='/home/saray/Desktop/taus/';

 mtau=[];   
 infm=[];
V=-100:1:50;

for i=1:length(V)
    
    v=V(i);

qt=q10^((celsius-35)/10);

  alpt = exp(0.0378*zetat*(v-vhalft-sh)) ;
  
    bett = exp(0.0378*zetat*gmt*(v-vhalft-sh)) ;
    

a = alpt;

minf = (1/(1 + exp((v-vhalfl-sh)/kl)));

infm=[infm,minf];

taum = b0 + bett/(a0t*(1+a));

mtau=[mtau,taum];

end


figure
plot(V,mtau,'o')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Time Constant (ms)','FontSize',20)
title(' Tau-m gKm','FontSize',20)
saveas(gcf,[path,'kmtaum'],'jpg');

figure
plot(V,infm,'ro')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Activation variable m','FontSize',20)
title(' m-inf gKm','FontSize',20)
saveas(gcf,[path,'km_m'],'jpg');


