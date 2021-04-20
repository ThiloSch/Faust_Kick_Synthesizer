import("stdfaust.lib");
import("stdfaust.lib");
import("physmodels.lib");

// Trigger and global Gain
trigger  = button("gate");

// Amplitude Envelope
attack = 0.03;
decay = 1;
sustain = 0;
release  = hslider("/h:[1]envelope/Release",1, 0.05, 3, 0.001);

// Filter control
f_frequency = hslider("Cut Off", 900, 50, 4000, 1);
qfactor = hslider("Q", 1,0.1,20,0.1);

// Additional Controls
physmix = hslider("mix",-1,-1,1,0.01);                                                  // Classic/Physical Mix
sigType = hslider("Shape", 0, 0, 1.5, 0.01);                                            // Signal Type

// Pitch Controls
p_note  = hslider("/h:[0]Pitch/Pitch Note",250,40, 1000, 0.1);                          // Starting Pitch
frequency = hslider("/h:[0]Pitch/frequency", 50, 40, 100, 0.1);                         // Endfrequency
p_sustain = hslider("/h:[0]Pitch/Pitch Sustain", 0.02, 0.001, 0.05, 0.001);             // Pitch Sustain
p_release  = hslider("/h:[0]Pitch/Pitch Release",0.1, 0.001, 3, 0.001);                 // Pitch release


//OSCILLATOR SECTION//
//-----------------//
oscillator_sin(sig) = signal
with{
    swing = pitcher(p_note, frequency, p_sustain, p_release, trigger);
    coef = 1;
    signal = os.osc(swing) * coef ;
};

oscillator_tri(sig) = signal
with{
    swing = pitcher(p_note, frequency, p_sustain, p_release, trigger);
    coef = select2(sig<2,0,select2(sig<1,2-sig,sig));
    signal = os.triangle(swing)* coef  ;
};


oscillator_sqa(sig) = signal
with{
    swing = pitcher(p_note, frequency, p_sustain, p_release, trigger);
    coef = select2(sig>1,0,select2(sig<2,3-sig,-1+sig));
    signal =  os.square(swing)  * coef ;
};

// Classic Approach with Pitch Envelope
pitcher(start, end, sus, rel, trig) = freq
with {
//time raises/resets indefinitely when the trigger is pressed
  time = (raise*reset + upfront) ~ _
    with {
        upfront = trig > trig';
        reset = trig <= trig';
        raise(x) = (x + (x > 0));
    };

  susSamps = int(sus * 44100);

  target = select2( susSamps>time, end, start);

  relPhase = susSamps < time;

  pole = ba.tau2pole(rel/6.91);

  freq = target : select2(relPhase, start, si.smooth(pole));
};
//-----------------//




//PHYSICAL MODELLING SECTION//
//-----------------//
meshdecay = hslider("/h:[3]Physical Modelling/Damping",25,1,2000,0.01);
phystune = hslider("/h:[3]Physical Modelling/Pitch",1,0,1,0.01);
physfreq = frequency;  

striker = pm.impulseExcitation(trigger) : chain;

meshed = striker : drumModel(frequency,0,1,1,meshdecay) : *(en.ar(0.03,10,trigger)); //
composite = 8*meshed;
// Extracted Modes
drumModel(freq ,exPos,t60,t60DecayRatio,t60DecaySlope) = _ <: par(i,nModes,modeFilter(modesFreqs(i),modesT60s(i),modesGains(int(exPos),i))) :> /(nModes)
with{
nModes = 40;
nExPos = 5;
modesFreqRatios(n) = ba.take(n+1,(0.496487 ,0.498829 ,0.510539 ,1.000000 ,1.402810 ,1.407494 ,1.510539 ,1.512881 ,1.562061 ,1.805621 ,1.807963 ,1.814988 ,2.032787 ,2.072600 ,2.081967 ,2.086651 ,2.100703 ,2.271663 ,2.278689 ,2.510539 ,2.576112 ,2.580796 ,2.608899 ,2.648712 ,2.660422 ,2.702576 ,2.711944 ,2.714286 ,2.777518 ,2.789227 ,2.868852 ,2.873536 ,3.133489 ,3.222482 ,3.224824 ,3.302108 ,3.318501 ,3.416862 ,3.444965 ,3.449649));
modesFreqs(i) = freq/(2-phystune*1.5)*modesFreqRatios(i);
modesGains(p,n) = waveform{0.002800 ,0.002910 ,0.002828 ,1.000000 ,0.000329 ,0.261250 ,0.000724 ,0.004314 ,0.002937 ,0.009676 ,0.002475 ,0.003697 ,1.504010 ,0.030839 ,0.086771 ,0.000730 ,0.010672 ,0.000401 ,0.036587 ,0.002806 ,0.024838 ,0.000031 ,0.001831 ,0.000088 ,0.033158 ,0.004144 ,0.004803 ,0.002403 ,0.063744 ,0.000837 ,0.009838 ,0.001281 ,1.908791 ,0.005826 ,0.016356 ,0.000032 ,0.142409 ,0.004683 ,0.033018 ,0.001156},int(p*nModes+n) : rdtable : select2(modesFreqs(n)<(ma.SR/2-1),0);
modesT60s(i) = t60*pow(1-(modesFreqRatios(i)/3000)*t60DecayRatio,t60DecaySlope);
}; 
//-----------------//


process = ((oscillator_sin(sigType), oscillator_sqa(sigType), oscillator_tri(sigType)) :> _ * en.adsre(attack, decay, sustain, release, trigger))* (0.5*(1-physmix) : sqrt) ,  composite*(0.5*(1+physmix) : sqrt):> fi.resonlp(f_frequency,qfactor,1)  <:_,_;
