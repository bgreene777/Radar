%
% Program to read/plot sounding and calculate refractivity
%
clear all;
%
% load data and put into vectors
%-----------------------------------------------------------------------------
% PRES HGHT TEMP DWPT RELH MIXR DRCT SKNT THTA THTE THTV
% hPa m C C % g/kg deg knot K K K
%-----------------------------------------------------------------------------
%
sounding;
press=data(:,1);
alt=data(:,2)/1000;
temp=data(:,3);
td=data(:,4);
rh=data(:,5);
mixr=data(:,6);
direc=data(:,7);
speed=data(:,8);
%
% calculate mixing ratio and refractivity
% N: refractivity
% Ndry: refractivity w/o water vapor
% Nref: model refractivity (fig 2.7)
%
mixr=mixr/1000;
pw=press.*(mixr./(0.622+mixr));
N=(77.6./(temp+273)).*(press+4810*(pw./(temp+273)));
Ndry=(77.6./(temp+273)).*(press);
Nref=313*exp(-0.1439*alt);