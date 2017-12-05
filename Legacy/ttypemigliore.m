
%%%%Time Constant of the CaT channel from Hemond et al 2008
clear
    celsius = 25;	
	gcatbar=0.003 ;
	cai = 50.e-6 ;
	cao = 2 ;
	q10 = 5;
	mmin=0.2;
	hmin=10;
	a0h =0.015;
	zetah = 3.5;
	vhalfh = -75;
	gmh=0.6	;
	a0m =0.04;
	zetam = 2;
	vhalfm = -28;
	gmm=0.1	;
    
path='/home/saray/Desktop/taus/';
 mtau2=[];   
V=-100:1:50;
htau3=[];
%mtau1=[];
infm=[];
infh=[];

for i=1:length(V)
    
    v=V(i);

qt=q10^((celsius-25)/10);

alpmt = exp(0.0378*zetam*(v-vhalfm)) ;

betmt = exp(0.0378*zetam*gmm*(v-vhalfm));

mtau = betmt/(qt*a0m*(1+alpmt));

mtau2=[mtau2,mtau];
%mtau1=[mtau1,(1/(alpmt+betmt))];

    a = 0.2*(-1.0*v+19.26)/(exp((-1.0*v+19.26)/10.0)-1.0);
	b = 0.009*exp(-v/22.03);
	minf = a/(a+b);
    infm=[infm,minf];

end


figure
plot(V,mtau2,'o')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Time Constant (ms)','FontSize',20)
title(' Tau-m gCaT','FontSize',20)
saveas(gcf,[path,'ttypetaum'],'jpg');

figure
plot(V,infm,'ro')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Activation m','FontSize',20)
title('m-inf gCaT','FontSize',20)
saveas(gcf,[path,'ttypem'],'jpg');

% % figure
% % plot(V,mtau1,'o')
% % xlabel('Voltage (mV)')
% % ylabel('Time Constant (ms)')

for i=1:length(V)
    
    v=V(i);

alph = exp(0.0378*zetah*(v-vhalfh)); 
beth = exp(0.0378*zetah*gmh*(v-vhalfh)); 
htau = beth/(a0h*(1+alph));
htau3=[htau3,htau];

    a2 = 1.e-6*exp(-v/16.26);
	b2 = 1/(exp((-v+29.79)/10.)+1.);
	hinf = a2/(a2+b2);
    infh=[infh,hinf];
end


figure
plot(V,htau3,'o')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Time Constant h (ms)','FontSize',20)
title(' Tau-h gCaT','FontSize',20)
saveas(gcf,[path,'ttypetauh'],'jpg');

figure
plot(V,infh,'ro')
xlabel('Voltage (mV)','FontSize',20)
ylabel('Inactivation h','FontSize',20)
title('h-inf gCaT','FontSize',20)
saveas(gcf,[path,'ttypeh'],'jpg');



   






