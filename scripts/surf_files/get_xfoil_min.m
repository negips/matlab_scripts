% Select xfoil data

  clear
  clc
  close all

  xfiles{1}='saab750k_N7.9.t1.dat';     % file for Xfoil Data
  xfiles{2}='saab750k_N8.0.t1.dat';     % file for Xfoil Data
  xfiles{3}='saab750k_N8.3.t1.dat';     % file for Xfoil Data
  xfiles{4}='saab750k_N8.4.t1.dat';     % file for Xfoil Data
  xfiles{5}='saab750k_N8.5.t1.dat';     % file for Xfoil Data
  xfiles{6}='saab750k_N9.0.t1.dat';     % file for Xfoil Data


  tr24 = 0.766;
  tr34 = 0.659;
  tr44 = 0.436;

  tr_e = [tr24 tr34 tr44];

  nfiles=length(xfiles);

  for i=1:nfiles
    xfile=xfiles{i};
    re750 = importdata(xfile);
    alpha = re750.data(:,1);
    tr    = re750.data(:,7);

    trx     = interp1(alpha,tr,[2.4 3.4 4.4]);
    lsq(i)  = sum((trx-tr_e).^2);
  end        


