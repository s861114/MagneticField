figname=sprintf('%s_%.3f(eV)_%d(T)_conductivity',deblank(wksp.Name(indc,:)),Ef,B);
figname2=sprintf('%s_%.3f(eV)_%d(T)_plasmon',deblank(wksp.Name(indc,:)),Ef,B);
figname=fullfile(cd,'data',figname);
figname2=fullfile(cd,'data',figname2);

cntrl_plot=plot(hwS,real(y));        
set(cntrl_plot,'Color','blue','LineWidth',2);
hold on;
cntrl_plot=plot(hwS,imag(y));        

set(cntrl_plot,'Color','red','LineWidth',2);
export_fig(figname,'-tif', '-cmyk', '-r200');

close;

%--------------------Plasmonic---------------------------
for hwn=1:size(hwS,2)
    hw=hwS(hwn);
    omg=hw/(wksp.hevbar);
    beta=omg/wksp.c*sqrt(1-(2/(wksp.imp*y(hwn)*unit))^2);
    propwav(hwn)=abs(real(beta))/abs(2*pi*imag(beta));
end

freqS=hwS/(wksp.hevbar*2*pi)/10^12;    
cntrl_plot=loglog(freqS,propwav);
set(cntrl_plot,'Color','black','LineWidth',2)   

axis([10^0 10^3 10^(-5) 10^2]);
export_fig(figname2,'-tif', '-cmyk', '-r200');    
clear figname;
clear figname2;
close;

