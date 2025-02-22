%function [stat,x,y,Re,fltype,dstar,rlam,spanv,wAstat,sumw]=readxystats(filename)
function [stat,x,y,flow]=readxystats(filename)
%
% Read x-statistics file as defined in Simson.
% See wxys.f and boxxys.f for details.
%
% This version uses the vector form of fread
%
% The output "stat" is a structured object with
% all the statistics. (NOTE: function under construction).
%
%
rlam=0.0;
spanv=0.0;

%
% Open file
%
fid=fopen(filename,'r','ieee-le.l64');
eol=fread(fid,1,'int');

%
% If record is not 44 the file is either corrput or on big endian format
%
if eol ~= 44
  fclose(fid);
  disp(' ')
  disp(['Reading ' filename ' on big endian format'])
  fid=fopen(filename,'r','ieee-be.l64');
  eol=fread(fid,1,'int');
else
  disp(' ')
  disp(['Reading ' filename ' on little endian format'])
end
flow.Re=fread(fid,1,'float64');
flow.bau=fread(fid,1,'int');
flow.xl=fread(fid,1,'float64');
flow.zl=fread(fid,1,'float64');
flow.t=fread(fid,1,'float64');
flow.shift=fread(fid,1,'float64');eol=fread(fid,2,'int');
flow.A=fread(fid,1,'uint8=>char');
flow.mhd_n=fread(fid,1,'float64');
flow.b0=fread(fid,3,'float64');eol=fread(fid,2,'int');
flow.nx=fread(fid,1,'int');
flow.nyp=fread(fid,1,'int');
flow.nzc=fread(fid,1,'int');
flow.nfzsym=fread(fid,1,'int'); eol=fread(fid,2,'int');
flow.fltype=fread(fid,1,'int');
flow.dstar=fread(fid,1,'float64'); eol=fread(fid,2,'int');
if flow.fltype<0
     flow.rlam=fread(fid,1,'float64'); eol=fread(fid,2,'int');
elseif flow.fltype >=6
    flow.bstart=fread(fid,1,'float64');
    flow.blenght=fread(fid,1,'float64');
    flow.rlam=fread(fid,1,'float64');
    flow.spanv=fread(fid,1,'float64'); eol=fread(fid,2,'int');
end
flow.sumw=fread(fid,1,'float64');
flow.nxys=fread(fid,1,'int');
flow.nxysth=fread(fid,1,'int');
flow.scalar=fread(fid,1,'int');
%wAstat=logical(fread(fid,1,'int'));
wAstat = 0;
eol=fread(fid,1,'int'); %CHECKED

% List of velocity statistics.
Lvs={'u','v','w','u2','v2','w2',...
    'omX','omY','omZ','omX2','omY2','omZ2',...
    'uv','uw','vw',...
    'corrX_uu','corrX_vv','corrX_ww','corrY_uu','corrY_vv','corrY_ww',...
    'corrZ_uu','corrZ_vv','corrZ_ww',...
    're_eps11','re_eps22','re_eps33','re_eps12','re_eps13','re_eps23',...
    'p','p2','up','vp','wp','p_ux','p_vy','p_wz','p_uy','p_vz','p_wx','p_uz',...
    'umin','umax',...
    'Sij_Sij','S11','S12','S13','S22','S23','S33',...
    'p52addedinLES',...
    'tauij_Sij','tau_11','tau_12','tau_13','tau_22','tau_23','tau_33',...
    'nu_t_addedinLES',...
    'cond_sampling_forward_scatter','conditional_sampling_backscatter',...
    'u_skewness','v_skewness','w_skewness',...
    'u2_v','u2_w','v2_u','v2_w','w2_u','w2_v','u_v_w',...
    'p_vx','p_wy',...
    'u_flatness','v_flatness','w_flatness',...
    'p_skewness','p_flatness',...
    'tau11u_tau12v_tau13w','tau12u_tau22v_tau23w','tau13u_tau23v_tau33w',...
    'phi','phi2',...
    'j1','j2','j3','j1j1','j2j2','j3j3','j1_j2','j1_j3','j2_j3',...
    'phi_j1','phi_j2','phi_j3',...
    'umax','vmax','wmax','umin','vmin','wmin'};
% 
%last = 44; %nxys = last;
last=flow.nxys;
 if wAstat
     flow.nxys=flow.nxys+0;
 end
%rl=zeros(nyp*nx/2,1);iml=zeros(nyp*nx/2,1);
nxys=flow.nxys;
nx=flow.nx;
nyp=flow.nyp;
for i=1:1:nxys
    if i>last
        name = Lvs{end-6+mod(i,last)};
    else
        name = Lvs{i};
    end
    eol=fread(fid,1,'int');
    dummy=fread(fid,nx*nyp,'float64');
    vel0 = fliplr(reshape(dummy,nx,nyp)); 
    ix = [nx/2:nx, 1:nx/2-1]; vel0(:,:) = vel0(ix,:);
    stat.(name)=vel0';
    eol=fread(fid,1,'int');
end
fclose(fid);
flow.Re=flow.Re*flow.dstar;
x=(0:nx-1)/nx * flow.xl/flow.dstar;
y=(1-cos(linspace(0,pi,nyp))) * 1/2 * 2/flow.dstar;

% Correct value of umax, umin
flow.sumw
if nxys >= 43
    stat.umin = stat.umin*flow.sumw;
end
if nxys >= 44
    stat.umax = stat.umax*flow.sumw;
end

% Add root-mean-square velocity
if nxys >= 4
    stat.urms = sqrt(stat.u2-stat.u.^2);
end
if nxys >= 5
    stat.vrms = sqrt(stat.v2-stat.v.^2);
end
if nxys >= 6
    stat.wrms = sqrt(stat.w2-stat.w.^2);
end
 end
