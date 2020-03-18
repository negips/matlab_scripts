function alag=getalpha(time,lag)

  kred         =  0.4;
  U0           =  1.0;
  chord        =  1.0;
  semichord    =  chord/2;
  omega        =  kred*U0/semichord;
  Tosc         =  2*pi/omega;
  ptch_amp     =  1.0;
  ptch_start   =  6.0;
  ini_aoa      =  3.4;
  ini_phase    = -pi/2;

  alag = ini_aoa + ptch_amp*sin(omega*(time-ptch_start)+ini_phase+lag);

return






