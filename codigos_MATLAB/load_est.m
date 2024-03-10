function [lat,lon,alt,name,time,velv,dirv,temp,hrel,humesp,pres,rglb] = load_est(est)
%% Funcion que carga los parámetros necesarios que utilizan todos los codigos
% También hace el control de calidad para evitar posteriores errores

lat = est.lat;
lon = est.lon;
alt = est.alt;
name = est.statname;
velv = est.vels010pro;
velv(est.vels010proQC == 4) = NaN;
dirv = est.dirv010pro;
dirv(est.dirv010proQC == 4) = NaN;
time = est.time;
temp = est.temp005;
temp(est.temp005QC == 4) = NaN;
hrel = est.hrel005;
hrel(est.hrel005QC == 4) = NaN;
pres = est.pres005;
pres(est.pres005QC == 4) = NaN;
rglb = est.rglb005;
rglb(est.rglb005QC == 4) = NaN;

t_0 = 273.15;
p0 = 1013;
temp_k = temp + t_0;
e_s = exp((17.67*(temp_k-t_0))./(temp_k - 29.65));
humesp = ((hrel.*e_s)./(26.3*p0))*10^3;

end