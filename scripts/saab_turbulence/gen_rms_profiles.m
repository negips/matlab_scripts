% generate profiles

  addpath '/scratch/negi/git_repos/matlabscripts/scripts/';

  load 'rms.mat';

  nt = length(uu);

  
  col = ['kbrm'];
  mkr = ['--sd'];
  fs = 16;
  budgetfs = 12;
  
  destn='plots/';
  pcol=1;
  ifnorm = 0;
  ifrtf  = 0;
  ifprofile=0;
  ifplot=1;
  
  k=0.4;
  U0=1.0;
  chord=1;
  semichord=chord/2;
  omega=k*U0/semichord;
  Tosc=2*pi/omega;
  ptch_start=6.0;
  ptch_amp=1.0;
  ini_aoa=3.4;

  x0 = 0.8;


