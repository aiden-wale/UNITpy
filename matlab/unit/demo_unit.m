clear;   % Make sure we start from scratch

global ch; ch=1;    % Selection character (number) from list
global dm; dm=1;    % Do we offer the keyboard at the end of each testvec routine (0=no)
global dsp; dsp=1;  % Display level for testvec routines (0=no display,
                    % 1=iteration display, 2=iteration display + echo of testvec routine to Matlab command window

while ch>0
 disp(' ')
 disp('-------------------------------------------------------------')
 disp('University of Newcastle Identification Toolbox (UNIT)')
 disp('-------------------------------------------------------------')
 disp( '  ');
 disp('           Available Demonstrations                          ');
 disp('   ');
 disp('Linear, discrete time, SISO Model Structure Estimation Demos')
 disp('-------------------------------------------------------------')
 disp('(1)   AR Model Structure Estimation')
 disp('(2)   ARMA Model Structure Estimation')
 disp('(3)   Non-parametric Estimation')
 disp('(4)   FIR Model Structure Estimation')
 disp('(5)   ARX Model Structure Estimation')
 disp('(6)   ARMAX Model Structure Estimation')
 disp('(7)   Output Error (OE) Model Structure Estimation')
 disp('(8)   Box-Jenkins Model Structure Estimation')
 disp('(9)   Recursive FIR model estimation')
 disp('(10)  Recursive ARX model estimation')
 disp('   ');
 disp('Linear, continuous time, Model Estimation Demos')
 disp('------------------------------------------------------------')
 disp('(11) Output Error Continous Time Model Structure Estimation')
 disp('(12) MIMO Continuous Time (CT) State Space Model Estimation ')
 disp(' ')
 disp('Nonlinear Model Structure Estimation Demos')
 disp('------------------------------------------------------------')
 disp('(13) Memoryless Nonlinearity')     
 disp('(14) Hammerstein/Output Error Model Estimation')
 disp('(15) Wiener/Output Error Model Estimation')
 disp('(16) Hammerstein-Wiener/Output Error Model Estimation')
 disp('(17) MISO Hammerstein/Output-Error Model Estimation')
 disp('(18) MISO Hammerstein--Wiener/Output-Error Model Estimation')
 disp('(19) MIMO Bilinear State Space Model Estimation') 
 disp('(20) SISO Nonlinear Growth Model Estimation')
 disp(' ')
 disp('Linear, discret time, MIMO/MISO Estimation from Time domain data demos')
 disp('------------------------------------------------------------')
 disp('(21) MISO Transfer Function Output-Error Model Estimation')
 disp('(22) Subspace Based State-Space Model Structure Estimation') 
 disp('(23) Multiple Output (Time Series) State Space Model Estimation') 
 disp('(24) MIMO State Space (SS) Model Estimation')
 disp('(25) MIMO grey-box parametrized SS Model Estimation ') 
 disp(' ')
 disp('Estimation from Frequency domain data demos')
 disp('------------------------------------------------------------')
 disp('(26) SISO ARX (linear least squares) model estimation')
 disp('(27) SISO Output Error (non-linear least squares) model estimation')
 disp('(28) MIMO State Space Model estimated via Subspace method')  
 disp('(29) MIMO State Space Model estimated via ML method') 
 disp(' ')
 disp('Filtering and Smoothing')
 disp('------------------------------------------------------------')
 disp('(30) Kalman Predictor/Filter/Smoother Signal Estimation')  
 disp('(31) Particle Filter Signal Estimation')   
 disp(' ')
 disp('Other demos')
 disp('------------------------------------------------------------')
 disp('(32) Compute Posterior Mean and Density of SISO OE Model')
 disp('(33) Test optimisation engine on Rosenbrock banana function')
 disp(' ')
 disp('Correct Installation Test')
 disp('------------------------------------------------------------')
 disp('(100) Run all the above simply to test for no crashes')
 disp(' ')
 disp('(0) Quit this demo programme')
 disp(' ')

 ch = input('Which demo would you like to run? :');

 if ch>0 
  % Translate from number selected to numerical identifier for a given
  % demo.  Looks messy, but this way if the list above is re-ordered or a
  % new demo is inserted, it only requires a change in the translation
  % array following, and not in all the code after that.
  
  global trans;
  trans(1)=29;
  trans(2)=5; 
  trans(3)=10; 
  trans(4)=1;
  trans(5)=2; 
  trans(6)=6; 
  trans(7)=3;
  trans(8)=4; 
  trans(9)=8; 
  trans(10)=9;
  trans(11)=33; 
  trans(12)=20; 
  trans(13)=14;
  trans(14)=11; 
  trans(15)=12; 
  trans(16)=13;
  trans(17)=16;
  trans(18)=30; 
  trans(19)=19;
  trans(20)=28; 
  trans(21)=15;
  trans(22)=7;
  trans(23)=17; 
  trans(24)=18; 
  trans(25)=32;
  trans(26)=21; 
  trans(27)=22; 
  trans(28)=31;
  trans(29)=23; 
  trans(30)=24; 
  trans(31)=34;  
  trans(32)=27;
  trans(33)=26; 
  
  if isempty(ch), ch = 0; end
  
  if ch>=100, dm=[]; dsp=0; end;  % Flag to routines to not pause
 
  warning off;
  if (ch==100)||(trans(ch)==1)
   disp('Running demo_fir');
   demo_fir;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==2)
   disp('Running demo_arx');
   demo_arx;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==3)
   disp('Running demo_oe');
   demo_oe;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==4)
   disp('Running demo_bj');
   demo_bj;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==5)
   disp('Running demo_arma');
   demo_arma;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==6)
   disp('Running demo_armax');
   demo_armax;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==7)
   disp('Running demo_sid');
   demo_sid;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==8)
   disp('Running demo_recfir');
   demo_recfir;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==9)
   disp('Running demo_recarx');
   demo_recarx;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==10)
   disp('Running demo_nonpar');
   demo_nonpar;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==11)
   disp('Running demo_hammer');
   demo_hammer;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==12)
   disp('Running demo_wiener');
   demo_wiener;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==13)
   disp('Running demo_hammwiener');
   demo_hammwiener;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==14)
   disp('Running demo_memoryless');
   demo_memoryless;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==15)
   disp('Running demo_miso');
   demo_miso;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==16)
   disp('Running demo_miso_hammer');
   demo_miso_hammer;
  end;
   
  global ch;
  if (ch==100)||(trans(ch)==17)
   disp('Running demo_mtseries');
   demo_mtseries;
  end;
  
  global ch;
  if (ch==100)||(trans(ch)==18)
   disp('Running demo_mimo');
   demo_mimo;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==19)
   disp('Running demo_bilin');
   demo_bilin;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==20)
   disp('Running demo_mimo_ct');
   demo_mimo_ct;
  end;
  
  global ch;
  if (ch==100)||(trans(ch)==21)
   disp('Running demo_farx');
   demo_farx;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==22)
   disp('Running demo_foe');
   demo_foe;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==23)
%   disp('Running demo_mimo_freq');
%   demo_mimo_freq;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==24)
   disp('Running demo_kf');
   demo_kf;
  end;
  
  global ch;
  if (ch==100)||(trans(ch)==26)
   disp('Running demo_rosenbrock');
   demo_rosenbrock;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==27)
   disp('Running demo_mcmc');
   demo_mcmc;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==28)
%   disp('Running demo_nlss');
%   demo_nlss;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==29)
   disp('Running demo_ar');
   demo_ar;
  end;
  
  global ch;
  if (ch==100)||(trans(ch)==30)
   disp('Running demo_miso_hammerwiener');
   demo_miso_hammer;
  end;
  
  global ch;
  if (ch==100)||(trans(ch)==31)
   disp('Running demo_fsid');
   demo_fsid;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==32)
   disp('Running demo_struct');
   demo_struct;
  end;
 
  global ch;
 
  if (ch==100)||(trans(ch)==33)
   disp('Running demo_tfcts');
   demo_struct;
  end;
 
  global ch;
  if (ch==100)||(trans(ch)==34)
   disp('Running demo_sir');
   demo_sir;
  end;
 
  global ch; global dm; global dsp; dsp=1;
 
  if (ch==100), dm=1; end
  
 end; % test on ch>0  
end;

clear global ch; clear global dm; clear global dsp; clear global trans;
