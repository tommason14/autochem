        program thermochemistry

        implicit none

!CC mass from IUPAC
      real*8:: mass(50)=(/1.00794,4.002602,6.941,9.012182,10.811,
     .  12.0107,14.0067,15.9994,18.9984032,20.1797,22.989770,
     . 24.3050,26.981538,28.0855,30.973761,32.065,35.453,39.948,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,72.64,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,118.71/)
!CC mass from Gaussian
!CC       real*8:: mass(18)=(/1.00783,0.0, 0.0, 0.0,0.0,12.0,14.00307,
!CC     . 15.9941,18.9984,0.0,0.0,0.0,0.0,0.0,0.0,31.97207,
!CC     . 34.96885,0.0/)
       real*4 hpitzer(24,20),spitzer(24,20),qpitzer(20),vpitzer(24)

      data hpitzer(1:24,1)/0.994,1.1822,1.3513,1.5011,1.6324,1.7460,
     & 1.9607,2.0934,2.1657,2.1971,2.2030,2.1944,2.1788,2.1607,2.1261,
     & 2.0984,2.0781,2.0634,2.0526,2.0382,2.0292,2.0229,2.0182,2.0147/
      data hpitzer(1:24,2)/0.994,1.142,1.300,1.437,1.556,1.660,1.856,
     & 1.971,2.031,2.049,2.043,2.024,1.998,1.971,1.918,1.875,1.840,
     & 1.811,1.787,1.749,1.717,1.690,1.666,1.646/
      data hpitzer(1:24,3)/0.994,1.106,1.249,1.374,1.482,1.576,1.753,
     & 1.854,1.900,1.909,1.893,1.864,1.829,1.794,1.727,1.670,1.623,
     & 1.583,1.548,1.492,1.441,1.401,1.363,1.329/
      data hpitzer(1:24,4)/0.994,1.074,1.200,1.311,1.411,1.495,1.654,
     & 1.742,1.779,1.777,1.753,1.715,1.673,1.631,1.552,1.484,1.427,
     & 1.379,1.335,1.264,1.202,1.150,1.102,1.061/
      data hpitzer(1:24,5)/0.994,1.050,1.151,1.251,1.340,1.418,1.561,
     & 1.636,1.662,1.651,1.621,1.577,1.529,1.481,1.392,1.315,1.251,
     & 1.196,1.147,1.067,0.997,0.937,0.886,0.841/
      data hpitzer(1:24,6)/0.994,1.032,1.106,1.190,1.272,1.344,
     . 1.472,1.536,1.550,1.535,1.497,1.448,1.394,1.344,1.247,1.164,
     . 1.095,1.035,0.982,0.896,0.823,0.760,0.707,0.660/
      data hpitzer(1:24,7)/0.994,1.022,1.073,1.138,1.211,1.275,1.385,
     . 1.440,1.448,1.426,
     . 1.382,1.329,1.273,1.218,1.115,1.029,0.955,0.892,0.838,0.745,
     . 0.672,0.613,0.561,0.515/
      data hpitzer(1:24,8)/0.994,1.015,1.051,1.099,1.157,1.211,1.306,
     . 1.350,1.351,1.321,
     . 1.275,1.221,1.162,1.104,0.999,0.908,0.833,0.768,0.715,0.624,
     . 0.551,0.493,0.443,0.399/
      data hpitzer(1:24,9)/0.994,1.008,1.036,1.072,1.114,1.155,1.230,
     . 1.265,1.260,1.224,
     . 1.176,1.121,1.061,1.002,0.893,0.802,0.725,0.661,0.608,0.519,
     . 0.450,0.394,0.347,0.307/
      data hpitzer(1:24,10)/0.994,1.004,1.025,1.049,1.077,1.106,1.164,
     . 1.190,1.179,1.140,
     . 1.088,1.030,0.968,0.909,0.799,0.708,0.631,0.569,0.515,0.431,
     . 0.365,0.314,0.271,0.236/
      data hpitzer(1:24,11)/0.994,1.000,1.015,1.030,1.048,1.065,1.103,
     . 1.120,1.104,1.060,
     . 1.006,0.947,0.884,0.824,0.714,0.624,0.549,0.488,0.437,0.356,
     . 0.295,0.249,0.211,0.181/
      data hpitzer(1:24,12)/0.994,0.996,1.006,1.014,1.026,1.038,1.059,
     . 1.057,1.032,0.988,
     . 0.933,0.872,0.810,0.750,0.644,0.554,0.480,0.421,0.370,0.296,
     . 0.240,0.198,0.164,0.138/
      data hpitzer(1:24,13)/0.994,0.994,0.999,1.004,1.009,1.014,1.019,
     . 1.005,0.972,0.924,
     . 0.868,0.806,0.744,0.685,0.580,0.491,0.420,0.363,0.314,0.244,
     . 0.195,0.157,0.128,0.105/
      data hpitzer(1:24,14)/0.994,0.994,0.994,0.995,0.996,0.996,0.987,
     .0.962,0.922,0.870,
     . 0.811,0.749,0.687,0.628,0.523,0.437,0.368,0.312,0.269,0.202,
     . 0.158,0.127,0.099,0.080/
      data hpitzer(1:24,15)/0.994,0.994,0.992,0.990,0.984,0.982,0.962,
     . 0.928,0.882,0.828,
     . 0.765,0.701,0.638,0.580,0.476,0.392,0.326,0.273,0.231,0.170,
     . 0.127,0.098,0.077,0.061/
      data hpitzer(1:24,16)/0.994,0.992,0.990,0.987,0.980,0.972,0.945,
     . 0.904,0.850,0.791,
     . 0.727,0.661,0.599,0.540,0.437,0.354,0.290,0.240,0.200,0.143,
     . 0.103,0.076,0.060,0.047/
      data hpitzer(1:24,17)/0.994,0.992,0.988,0.984,0.976,0.965,0.932,
     . 0.886,0.827,0.763,
     . 0.697,0.630,0.567,0.508,0.406,0.324,0.261,0.211,0.174,0.121,
     . 0.084,0.061,0.047,0.036/
      data hpitzer(1:24,18)/0.994,0.991,0.988,0.982,0.974,0.962,0.922, 
     . 0.873,0.811,0.744,
     . 0.676,0.609,0.545,0.485,0.383,0.302,0.239,0.191,0.154,0.104,
     . 0.072,0.051,0.036,0.028/
      data hpitzer(1:24,19)/0.994,0.990,0.986,0.980,0.972,0.960,0.916,
     . 0.864,0.801,0.732,
     . 0.663,0.595,0.531,0.470,0.368,0.286,0.223,0.176,0.140,0.091,
     . 0.062,0.044,0.029,0.022/
      data hpitzer(1:24,20)/0.994,0.989,0.985,0.979,0.971,0.959,0.915,
     . 0.860,0.796,0.728,
     . 0.659,0.590,0.526,0.465,0.361,0.279,0.215,0.168,0.132,0.084,
     . 0.056,0.038,0.026,0.018/

      data spitzer(1:24,1)/0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     . 0.00,0.00,0.00,0.00,0.00,0.00,
     . 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/
       data spitzer(1:24,2)/6.946,6.941,6.926,6.902,6.869,6.828,
     . 6.694,6.529,6.347,6.163,
     . 5.982,5.813,5.657,5.515,5.272,5.072,4.906,4.766,4.643,
     . 4.438,4.270,4.127,4.003,3.892/
       data spitzer(1:24,3)/5.569,5.565,5.551,5.526,5.492,5.452,
     . 5.319,5.154,4.975,4.792,
     . 4.612,4.443,4.289,4.148,3.907,3.709,3.545,3.406,3.285,
     . 3.084,2.919,2.781,2.659,2.552/
       data spitzer(1:24,4)/4.763,4.759,4.745,4.720,4.688,
     . 4.648,4.515,4.353,4.178,3.995,
     . 3.819,3.652,3.498,3.359,3.120,2.926,2.765,2.629,2.511,
     .2.316,2.156,2.023,1.908,1.807/
       data spitzer(1:24,5)/4.192,4.188,4.174,4.152,4.120,
     . 4.080,3.950,3.790,3.615,3.435,
     . 3.263,3.098,2.948,2.812,2.576,2.385,2.230,2.097,1.984,
     . 1.798,1.645,1.518,1.411,1.320/
       data spitzer(1:24,6)/3.748,3.743,3.730,3.709,3.679,
     . 3.638,3.512,3.355,3.180,3.008,
     . 2.838,2.678,2.528,2.396,2.166,1.983,1.830,1.703,1.593,
     . 1.417,1.275,1.157,1.058,0.975/
       data spitzer(1:24,7)/3.386,3.382,3.370,3.347,3.318,
     . 3.279,3.156,3.004,2.836,2.667,
     . 2.500,2.343,2.199,2.068,1.844,1.665,1.519,1.397,1.295,
     . 1.125,0.994,0.890,0.801,0.727/
       data spitzer(1:24,8)/3.079,3.076,3.065,3.043,3.013,
     . 2.974,2.854,2.709,2.548,2.380,
     . 2.180,2.069,1.926,1.798,1.585,1.411,1.272,1.156,1.060,
     . 0.904,0.783,0.688,0.609,0.542/
       data spitzer(1:24,9)/2.814,2.811,2.801,2.780,2.750,
     . 2.714,2.600,2.458,2.303,2.138,
     . 1.978,1.834,1.698,1.579,1.370,1.204,1.071,0.962,
     . 0.872,0.728,0.620,0.533,0.464,0.405/
       data spitzer(1:24,10)/2.580,2.578,2.568,2.547,2.519,
     . 2.485,2.376,2.241,2.091,1.933,
     . 1.782,1.643,1.511,1.392,1.192,1.033,0.906,0.804,
     . 0.719,0.588,0.492,0.414,0.353,0.303/
       data spitzer(1:24,11)/2.371,2.369,2.359,2.340,2.315,
     . 2.279,2.173,2.048,1.907,1.756,
     . 1.610,1.475,1.348,1.233,1.040,0.891,0.770,0.674,0
     . .596,0.476,0.388,0.322,0.270,0.228/
       data spitzer(1:24,12)/2.182,2.180,2.170,2.151,2.125,
     . 2.094,1.997,1.874,1.739,1.576,
     . 1.458,1.328,1.209,1.097,0.915,0.774,0.660,0.570,
     . 0.496,0.388,0.309,0.251,0.205,0.170/
       data spitzer(1:24,13)/2.009,2.003,1.996,1.980,1.957,
     . 1.928,1.833,1.718,1.589,1.456,
     . 1.323,1.199,1.086,0.982,0.808,0.672,0.566,0.483,
     . 0.414,0.315,0.247,0.196,0.158,0.129/
       data spitzer(1:24,14)/1.850,1.848,1.837,1.823,1.800,
     . 1.744,1.685,1.578,1.456,1.330,
     . 1.206,1.087,0.978,0.881,0.715,0.588,0.486,0.407,
     . 0.348,0.255,0.196,0.155,0.121,0.097/
       data spitzer(1:24,15)/1.703,1.701,1.691,1.677,1.654,
     . 1.629,1.552,1.450,1.335,1.217,
     . 1.100,0.988,0.884,0.794,0.637,0.516,0.422,0.350,
     . 0.293,0.213,0.157,0.119,0.093,0.073/
       data spitzer(1:24,16)/1.567,1.563,1.555,1.541,1.523,
     . 1.499,1.428,1.332,1.224,1.114,
     . 1.004,0.901,0.804,0.716,0.568,0.453,0.366,0.300,
     . 0.248,0.176,0.126,0.092,0.072,0.056/
       data spitzer(1:24,17)/1.438,1.433,1.428,1.415,1.399,
     . 1.377,1.310,1.224,1.126,1.021,
     . 0.919,0.821,0.730,0.648,0.509,0.401,0.320,0.258,
     . 0.211,0.146,0.100,0.075,0.056,0.042/
       data spitzer(1:24,18)/1.316,1.312,1.307,1.295,1.284,
     . 1.262,1.201,1.122,1.031,0.936,
     . 0.841,0.748,0.662,0.588,0.457,0.357,0.281,0.223,
     . 0.180,0.122,0.084,0.059,0.042,0.032/
       data spitzer(1:24,19)/1.203,1.196,1.193,1.184,1.171,
     . 1.153,1.094,1.024,0.942,0.855,
     . 0.769,0.683,0.607,0.535,0.412,0.319,0.248,0.195,
     . 0.154,0.101,0.069,0.048,0.034,0.024/
       data spitzer(1:24,20)/1.097,1.091,1.085,1.076,1.068,
     . 1.052,1.000,0.936,0.860,0.779,
     . 0.703,0.623,0.551,0.486,0.372,0.285,0.220,0.171,
     . 0.134,0.084,0.056,0.038,0.026,0.018/

       data vpitzer(1:24)/0.0,0.2,0.4,0.6,0.8,1.0,1.5,2.0,
     . 2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,
     . 8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0/
       data qpitzer(1:20)/0.0,0.05,0.10,0.15,0.20,0.25,
     . 0.30,0.35,0.40,0.45,0.50,0.55,0.60,
     . 0.65,0.70,0.75,0.80,0.85,0.90,0.95/

       real*8 ixx, izz, iyy, ixz, ixy, iyz
       real*8 coord(1000,3), x2(3), imatrix(3,3)
       real*8 imol(3), cosmol(3,3),mtot,xcm(3),coordt(1000,3)
       real*8 xcmrf(3,2), coordrf(1000,3,2)
       real*8 coordcm(1000,3), coordpra(1000,3,2)
       real*8 modvecr, vecr(3), rcm(3,2), nuni(3), veccm(3,2)
       real*8 modveccm(2),mrf(2)
       real(8) xdstar(1000,3,2), itoptot(2)
       real(8) tmp1, tmp2, itopstar(3,2), i23(1000)
       real*8 tmp3,tmp,imolkg(3),tk,a0
       real*8 eigtmp(3), eigvectmp(3,3),freqsp(1000)
       integer tmpminl, tmpmaxl, tmpmidl,geomcheck,geommax
       integer natomrf(1000,2),chargerf(1000,2)
       integer nrot,nrf(2),atomrf(1000,2),k,m,nsp,posfreq(1000)
       real*8 eigvalues(3), eigvectors(3,3)
       real*8 tcvibr, tcrot, tctotal, zpvetot, tcelec
       real*8 snusc(1000),theta(1000),qus(1000)
       real*8 snu(1000),quh(1000),hnusc(1000)
       real*8 hnu(1000),zpnusc(1000),quzp(1000)
       real*8 zpnu(1000)
       real(8) snuhr(1000),snufree(1000),hnuhr(1000)
       real(8) hnufree(1000),qunuhr(1000),qunufree(1000)
       real*8 qcorrfree,scorrfree,hcorrfree,qcorrhr,scorrhr
       real*8 hcorrhr,qnutotfree,qnutothr
       real*8 snutotfree,stotfree,htothr,htotfree,htot
       real*8 iredmom,stot,snutothr
       real(8) klevels,cm2kj,qrot,srot,qelec,selec
       real*8 strans,tctrans,qtrans,thetaxyz(3)
       real(8) qtothr,stothr,qtot,qnutot,snutot,pi
       real(8) qtotfree,hnutotfree,hnutothr
       real*8 tokg,tokg2
       integer i,nfreq,mult,sigma,j,j1,j2,i1
       real*8 freq(1000)
       real*8 na, kb,runi,hplanck,clight,patm,ppa
       real(kind=8) sfzpve,sfs,sffreq
       real*8 sftc
       real*8 v0(1000),invqunuhr(1000),valh,vals
       real*8 enthjk,enthj1k,enthjk1,enthj1k1,entrjk
       real*8 entrj1k,entrjk1,entrj1k1,uh,th
       character*80 afile,acommand,ifactor,alevel
       character*80 abasis,line,ired
       logical ioptlevel,ioptbasis
        character*80 filepart,a,b
        character*1 ipart
       integer natom, numatom(1000), charge(1000), chtmp
       integer parameter(1000), symint(1000),nim
       integer possp(1000)
       character*80 aawk,geom,imode,afiledat,ipitzer


C        print *,"Enter the frequency job file"
C        read *,afile
C        write(6,'("Is it (y) optimization followed by a frequency 
C     . calculation or (n) a frequency job?")')
C        read *,geom

C        write(afiledat,'(a,''.dat'')') trim(afile)

C        open(2,file=trim(afile),status='old')
C        open(10,file=trim(afiledat),status='unknown')

C        call system('rm -rf freq.out')
C        write(acommand,'(''awk -f ~/bin/freq.awk '',a)') trim(afile) 
C        call system(acommand)

        posfreq=0 
        possp=0
!CC parameters
        a0=0.5291772108
        na=6.0221367E+23
        kb=1.380658E-23
        runi=8.31441
        hplanck=6.6260755E-34
        clight=2.99792458E+10
        patm=1. 
        ppa=101325
        cm2kj=2625.5
        pi=3.1415926539
        tokg=1.6605402E-27
        tokg2=1.6605402E-47
        
        open(3,file='freq.out',status='unknown')
       
        i=1
  10    read(3,*,end=20) freq(i)
!CC        if(freq(i).eq.0.0) goto 20
        i=i+1 
        goto 10
  20    nfreq=i-1

C!CC get all necessary input information
!CC       print *,"Accepted levels of theory:
!CC     . HF, MP2-fu, MP2-fc, QCISD-fc, B3LYP,BLYP, PB86, B3-PW91"
       
!CC        print *,"Enter the level of theory"
!CC        read *,alevel
           alevel="M062X"
       
!CC        print *,"Possible basis sets: 
!CC     . 6-31G(d), 6-31G(d,p), 6-311G(d,p), 6-311G(df,p),6-31+G(d),
!CC     . 6-31G(2df,p), 6-311+G(3df,2p),vtz" 
!CC       print *,"Enter the basis set"
!CC        read *,abasis
        abasis="vdz"
        ioptbasis=.false.
        ioptlevel=.false.
 
        select case(alevel)
        case("HF")
        ioptlevel=.true. 
             select case(abasis)
             case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=0.9135
             sftc=0.8905
             sfs=0.8978
             sffreq=0.9061
             
             case("6-31G*")
             ioptbasis=.true. 
             sfzpve=0.9135
             sftc=0.8905
             sfs=0.8978
             sffreq=0.9061

             case("6-31+G(d)")
             ioptbasis=.true. 
             sfzpve=0.9153
             sftc=0.8945
             sfs=0.9027
             sffreq=0.9131

             case("6-31+G*")
             ioptbasis=.true. 
             sfzpve=0.9153
             sftc=0.8945
             sfs=0.9027
             sffreq=0.9131

             case("6-311G(d,p)")
             ioptbasis=.true. 
             sfzpve=0.9248
             sftc=0.8951
             sfs=0.9021
             sffreq=0.9110

             case("6-311G**)")
             ioptbasis=.true. 
             sfzpve=0.9248
             sftc=0.8951
             sfs=0.9021
             sffreq=0.9110
             

             case("6-311G(df,p)")
             ioptbasis=.true. 
             sfzpve=0.9247
             sftc=0.8908
             sfs=0.8981
             sffreq=0.9085

             case("6-31G(d,p)")
             sfzpve=0.9181
             sftc=0.8912
             sfs=0.8990
             sffreq=0.9089
             ioptbasis=.true. 

             case("6-31G**")
             sfzpve=0.9181
             sftc=0.8912
             sfs=0.8990
             sffreq=0.9089
             ioptbasis=.true. 

             case("3-21G(d)")
             ioptbasis=.true. 
             sfzpve=0.9207
             sftc=0.9444
             sfs=0.9666
             sffreq=1.0075
!CC end HF
             end select
             if(alevel.eq."HF") then
              if(.not.ioptbasis) then
             ioptbasis=.true. 
             sfzpve=0.8929
             sftc=0.8929
             sfs=0.8929
             sffreq=1.0
             endif
             endif
         case("B3LYP")
         ioptlevel=.true.
         select case(abasis)

         case("6-31G(2df,p)")              
             ioptbasis=.true. 
             sfzpve=0.9854
             sftc=0.9854
             sfs=0.9854
             sffreq=1.0

         case("6-311+G(3df,2p)")
             ioptbasis=.true. 
             sfzpve=0.985
             sftc=0.985
             sfs=0.985
             sffreq=1.0

         case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=0.9806
             sftc=0.9989
             sfs=1.0015
             sffreq=1.0013
           
         case("6-31G*")
             ioptbasis=.true. 
             sfzpve=0.9806
             sftc=0.9989
             sfs=1.0015
             sffreq=1.0013

         case("vtz")
             ioptbasis=.true. 
             sfzpve=0.985
             sftc=0.985
             sfs=0.985
             sffreq=1.0
!CC end B3LYP
         end select
             if(alevel.eq."B3LYP") then
              if(.not.ioptbasis) then
             ioptbasis=.true. 
             sfzpve=0.985
             sftc=0.985
             sfs=0.985
             sffreq=1.0
             endif
             endif

         case("BLYP")
             ioptlevel=.true. 
             select case(abasis)
             case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=1.0126
             sftc=1.0633
             sfs=1.067
             sffreq=1.062

             case("6-31G*")
             ioptbasis=.true. 
             sfzpve=1.0126
             sftc=1.0633
             sfs=1.067
             sffreq=1.062

            case("6-311G(df,p)")
             ioptbasis=.true. 
             sfzpve=1.0617
             sftc=1.0593
             sfs=1.0641
             sffreq=1.0667
!CC end BLYP                
            end select 
            
           case("QCISD-fc")
           ioptlevel=.true.
           if(abasis.eq."6-31G(d)".or.abasis.eq."6-31G*") then
             ioptbasis=.true. 
             sfzpve=0.9776
             sftc=1.0080
             sfs=1.0187
             sffreq=1.047
           endif
!CC end QCISD
           
           case("MP2-fu")
           ioptlevel=.true.
           select case(abasis)
           case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=0.9661
             sftc=1.0084
             sfs=1.0228
             sffreq=1.0214
           case("6-31G*")
             ioptbasis=.true. 
             sfzpve=0.9661
             sftc=1.0084
             sfs=1.0228
             sffreq=1.0214
!CC end MP2-FU
           end select

           case("MP2-fc")
           ioptlevel=.true.
           select case(abasis)
           case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=0.9670
             sftc=1.0211
             sfs=1.0444
             sffreq=1.0485
           case("6-31G*")
             ioptbasis=.true. 
             sfzpve=0.9670
             sftc=1.0211
             sfs=1.0444
             sffreq=1.0485
           case("6-31G(d,p)")             
             ioptbasis=.true. 
             sfzpve=0.9608
             sftc=1.0084
             sfs=1.0232
             sffreq=1.0229
           case("6-311G(d,p)")
             ioptbasis=.true. 
             sfzpve=0.9748
             sftc=1.0061
             sfs=1.0171
             sffreq=1.0127
!CC end "MP2-fc"
           end select  
           case("B3P86")
           ioptlevel=.true.
           select case(abasis)
           case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=0.9759
             sftc=0.9864
             sfs=0.9902
             sffreq=0.9923
           case("6-31G*")
             ioptbasis=.true. 
             sfzpve=0.9759
             sftc=0.9864
             sfs=0.9902
             sffreq=0.9923
!CC eend B3P86
          end select
           case("BP86")
           ioptlevel=.true.
           select case(abasis) 
           case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=1.0108
             sftc=1.0478
             sfs=1.0527
             sffreq=1.0512
           case("6-31G*")
             ioptbasis=.true. 
             sfzpve=1.0108
             sftc=1.0478
             sfs=1.0527
             sffreq=1.0512
!CC end BP86
           end select
           case("B3PW91")   
           ioptlevel=.true.
           select case(abasis)
           case("6-31G(d)")
             ioptbasis=.true. 
             sfzpve=0.9774
             sftc=0.9885
             sfs=0.9920
             sffreq=0.9930
           case("6-31G*")
             ioptbasis=.true. 
             sfzpve=0.9774
             sftc=0.9885
             sfs=0.9920
             sffreq=0.9930
!CC B3PW91
          end select
!CC end level
         end select
         if(.not.ioptlevel) then
             sfzpve=1.0
             sftc=1.0
             sfs=1.0
             sffreq=1.0
         endif
         if(.not.ioptbasis) then
             sfzpve=1.0
             sftc=1.0
             sfs=1.0
             sffreq=1.0
         endif


!CC         if(.not.ioptlevel) then
!CC             sfzpve=0.8929
!CC             sftc=0.8929
!CC             sfs=0.8929
!CC             sffreq=1.0
!CC         endif
!CC         if(.not.ioptbasis) then
!CC             sfzpve=0.8929
!CC             sftc=0.8929
!CC             sfs=0.8929
!CC             sffreq=1.0
!CC         endif
          

!CC         write(*,'(''For '',2x,a,''/'',a)') trim(alevel), trim(abasis)
!CC         write(*,'(''For your method/basis set selection'')')
         write(*,'(''Scaling factor for ZPVE 
     .  '',f6.4,2x,''Is it OK (y/n)?'')')
     .  sfzpve
         read *,ifactor
         if(ifactor.eq."n".or.ifactor.eq."N") then
            print *,"Enter your own value"
            read *,sfzpve
         endif  
         ifactor=""
         write(*,'(''Scaling factor for Temp correction
     .  '',f6.4,2x,''Is it OK (y/n)?'')') sftc
         read *,ifactor
         if(ifactor.eq."n".or.ifactor.eq."N") then
            print *,"Enter your own value"
            read *,sftc
         endif  

         ifactor=""
         write(*,'(''Scaling factor for Entropy 
     . '',f6.4,2x,''Is it OK (y/n)?'' )') sfs
         read *,ifactor
         if(ifactor.eq."n".or.ifactor.eq."N") then
            print *,"Enter your own value"
            read *,sfs
         endif  

!CC         ifactor=""
!CC         write(*,'(''Scaling factor for Vibrational frequencies '',f6.4,2x,''Is it OK (y/n)?'')') sftc
!CC         read *,ifactor
!CC         if(ifactor.eq."no".or.ifactor.eq."NO".or.ifactor.eq."n".or.ifactor.eq."N") then
!CC            print *,"Enter your own value"
!CC            read *,sffreq
!CC         endif  

         print *,"Enter multiplicity"
         read *,mult
!CC         print *,"Enter rotational symmetry number of the molecule"
!CC         read *,sigma
         sigma=1
         print *,"Enter temperature"
         read *,tk

         
!CC get coordinates and molecular mass
C         close(2) 

C         open(2,file=trim(afile),status='old')

C         open(4,file='geom.input',status='unknown')
!CC         write(*,*) geom

C        if(geom.eq."n".or.geom.eq."N") then

C 100    read(2,'(a)',end=110) line
C        if(line(26:46)=="Standard orientation:") then
C          do i=1,4
C             read(2,'(a)') line
C          enddo
C  120      read(2,'(a)') line
C          if(line(1:10).eq." ---------") goto 110
C             write(4,'(a)') line
C             goto 120
C        else 
C        goto 100
C        endif
C        endif
C  110   continue
C
C        if(trim(geom).eq."y".or.trim(geom).eq."Y") then
C           geomcheck=0
C 130      read(2,'(a)',end=140) line
C          if(line(1:23)==" Optimization completed") then
C 131          read(2,'(a)') line
C              if(line(26:46)=="Standard orientation:".or.
C     . line(26:46)=="Z-Matrix orientation:") then
C!CC                 write(*,*) "I am in"
C                 geomcheck=geomcheck+1
C                if(geomcheck.ne.2) goto 131 
C                 do i=1,4
C                   read(2,'(a)') line
C                 enddo
C  132            read(2,'(a)') line
C                 if(line(1:10).eq." ---------") goto 140
C                 write(4,'(a)') line
C                 goto 132
C                 else
C                 goto 131
C               endif
C            else
C            goto 130
C            endif
C         endif
C  140   close(4)
        

!CC calculate moment of inertia 

       call mominertia(natom,coord,charge,xcm,imol,cosmol,
     . imolkg,a0,numatom,mtot)
!CC       write(*,*) 'mtot',mtot
       mtot=mtot*tokg 
!CC check which modes should be treated specially
       j=0
       nim=0
       do i=1,nfreq
          if(freq(i).le.0.) then
             nim=nim+1
          endif
       enddo

         do i=nim+1,nfreq 
          if(freq(i).le.300.) then
       goto 900 
       write(*,'(''Do you want to treat 
     . this vibr. mode '',f8.4,'' 
     . specially (y/n) ? '')') freq(i)
             read *,imode
             if(imode.eq."Y".or.imode.eq."y") then
                j=j+1
                freqsp(j)=freq(i)
                posfreq(i)=j
                possp(j)=i
                i23(j)=0.0
        call redmoment(natom,i23(j),coord,charge,xcm,a0,
     . imol,cosmol,freqsp(j),j,numatom)
        print *,"Enter the internal symmetry number of this rotation"
        read *,symint(j)
        iredmom=i23(j)*tokg2
        qunufree(j)=1/(symint(j)*hplanck)*
     . sqrt(8.*pi**3*i23(j)*tokg2*kb*tk)
        snufree(j)=runi*(0.5*(log(8.*pi**3*i23(j)*1.6605402*1.380658*tk)
     . -161.18095651)-log(symint(j)*hplanck)+0.5)
                hnufree(j)=0.5*runi*tk/1000
         write(*,'(''Do you want use Piters tables (y/n)?'')')
                read *,ipitzer
                if(ipitzer.eq."n".or.ipitzer.eq."N") then
        write(*,'(''Read partition function from file?'')')
                  read *,ipart
                  if(ipart.eq.'y'.or.ipart.eq.'Y') then
        write(*,'(''Enter the file name for this vibr mode'')')
                    read *,filepart
                    open(100,file=filepart,status='old')
                    read(100,*)
                    read(100,*)
                    read(100,*)
                    read(100,*)
                    read(100,*) a,b,qunuhr(j)
                    read(100,*) a,b,hnuhr(j)
                    read(100,*) a,b,snuhr(j)
                    close(100)
                  else
                  print *,"Enter the partition function Qt"
                  read *,qunuhr(j)
                  print *,"Enter the enthalpy H"
                  read *,hnuhr(j)
                  print *,"Enter the entropy S"
                  read *,snuhr(j)
                  endif
                endif
                if(ipitzer.eq."y".or.ipitzer.eq."Y") then
                  print *,"Enter your rotational barrier"
                  read *,v0(j)
                  v0(j)=v0(j)*cm2kj/(runi*tk)*1000
      qunuhr(j)=2.7935/symint(j)*sqrt(i23(j)*tk*1.66035*0.01)
                  invqunuhr(j)=1/qunuhr(j)
                j1=0.0
                j2=0.0
                do i1=1,20   
                if(invqunuhr(j).ge.qpitzer(i1)) j1=i1
                enddo
                do i1=1,24
                  if(v0(j).ge.vpitzer(i1)) j2=i1
                enddo
      th=(invqunuhr(j)-qpitzer(j1))/(qpitzer(j1+1)-qpitzer(j1))
                uh=(v0(j)-vpitzer(j2))/(vpitzer(j2+1)-vpitzer(j2))
                enthjk=hpitzer(j2,j1)
                enthj1k=hpitzer(j2,j1+1)
                enthjk1=hpitzer(j2+1,j1+1)
                enthj1k1=hpitzer(j2+1,j1)
                entrjk=spitzer(j2,j1)
                entrj1k=spitzer(j2,j1+1)
                entrjk1=spitzer(j2+1,j1+1)
                entrj1k1=spitzer(j2+1,j1)

                valh=(1-th)*(1-uh)*enthjk+th*(1-uh)*enthj1k+
     . th*uh*enthjk1+(1-th)*uh*enthj1k1
                hnuhr(j)=tk*4.184*valh/1000
                vals=(1-th)*(1-uh)*entrjk+th*(1-uh)*entrj1k+
     . th*uh*entrjk1+(1-th)*uh*entrj1k1
                snuhr(j)=4.184*vals
                endif 
             endif
           endif
 900      enddo

         nsp=j


!CC calculate vibrational contribution     

       klevels=0.69503877
      
       do i=nim+1,nfreq
          snusc(i)=freq(i)*clight*sfs
          theta(i)=snusc(i)*hplanck/kb
          tmp=exp(-theta(i)/tk)
          qus(i)=1/(1-tmp)
          snu(i)=runi*((theta(i)/tk)/(exp(theta(i)/tk)-1)-log(1-tmp))
          hnusc(i)=freq(i)*sftc*clight
          quh(i)=hnusc(i)*hplanck/kb
          hnu(i)=runi*(quh(i)/(exp(quh(i)/tk)-1))/1000 
          zpnusc(i)=freq(i)*clight*sfzpve
          quzp(i)=zpnusc(i)*hplanck/kb
          zpnu(i)=runi*quzp(i)*0.5/1000
!CC          write(*,*) "zpnu",zpnu(i)
       enddo

!CC temperature correction

          tcvibr=sum(hnu(nim+1:nfreq))     
          tctrans=2.5*runi*tk/1000
!CC check for non-linearity
!CC          write(*,*) 'natom',natom 
          select case(natom)
          case(1:2)
          tcrot=runi*tk/1000
          case(3:)
          tcrot=1.5*runi*tk/1000
          end select
          tcelec=0.0
!CC          tctotal=(tcvibr+tctrans+tcrot+tcelec)/cm2kj
          tctotal=tcvibr+tctrans+tcrot+tcelec 
!CC ZPVE          
!CC          zpvetot=sum(zpnu(nim+1:nfreq))/cm2kj
          zpvetot=sum(zpnu(nim+1:nfreq))
!CC          write(*,*) 'zpve',zpvetot

!CC Entropy
          qtrans=(2*pi*mtot*kb*tk/hplanck**2)**1.5*kb*tk/ppa
          strans=runi*(log(qtrans)+2.5)

          qelec=mult
          selec=runi*log(qelec)

          do i=1,3
             thetaxyz(i)=hplanck**2/(8*pi**2*kb*imolkg(i))
          enddo

          qrot=tk**1.5/(thetaxyz(1)*thetaxyz(2)*thetaxyz(3))**0.5*pi**0.5/sigma
          srot=runi*(log(qrot)+1.5)

          qnutot=product(qus(nim+1:nfreq))
          snutot=sum(snu(nim+1:nfreq))

          qtot=qtrans*qelec*qrot*qnutot
          stot=strans+selec+srot+snutot
!CC calculating total S for free rotor and Hindered rotor
          if(nsp.ne.0) then
          snutotfree=0. 
          snutothr=0.
          do i=nim+1,nfreq
             if(posfreq(i).ne.0) then
             snutotfree=snutotfree+snufree(posfreq(i))
             snutothr=snutothr+snuhr(posfreq(i))
             else
             snutotfree=snutotfree+snu(i)
             snutothr=snutothr+snu(i)
             endif
          enddo
          endif

!CC Output 
           write(10,*)
           write(10,'(''Thermochemistry was done at '',
     . 1x,f10.2,1x,''(k) and 1(Atm)'')') tk
           write(10,*)
           write(10,'(''Total vibr frequencies        
     .                 ='',i3)') nfreq
          write(10,'(''Imaginary frequencies           
     .               ='',i3)') nim
           write(10,'(''Low freqs (<300 cm^-1) with 
     . special treatment  ='',i3)') nsp

           write(10,*)
           write(10,'(''Scaling factors'')')   
           write(10,'(''ZPVE ='',1x,f8.4)') sfzpve
           write(10,'(''TC   ='',1x,f8.4)') sftc
           write(10,'(''Svib ='',1x,f8.4)') sfs
           write(10,*)
           
             
      write(10,'(''   Frequency     ZPVE(vi)      H(vi)      S(vi)'')')
           do i=nim+1,nfreq          
              if(posfreq(i).ne.0) then
                 write(10,'(f10.4,2x,3(f15.6,2x),2x,''FR'')') 
     . freq(i),zpnu(i),hnufree(posfreq(i)),snufree(posfreq(i)) 
                 write(10,'(f10.4,2x,3(f15.6,2x),2x,''HR'')')
     . freq(i),zpnu(i),hnuhr(posfreq(i)),snuhr(posfreq(i)) 
                 write(10,'(f10.4,2x,3(f15.6,2x),2x,''HO'')')
     .  freq(i),zpnu(i),hnu(i),snu(i) 
                 write(10,*)
              else
                 write(10,'(f10.4,2x,3(f15.6,2x),2x,''HO'')') 
     . freq(i),zpnu(i),hnu(i),snu(i)
              endif
           enddo
!CC HR correction
      if(nsp.ne.0) then
      qcorrfree=1.
      scorrfree=0.
      hcorrfree=0. 
      qcorrhr=1.
      scorrhr=0.
      hcorrhr=0. 
       do i=1,nsp
        qcorrhr=qcorrhr*qunuhr(i)/qus(possp(i))  
        scorrhr=scorrhr+(snuhr(i)-snu(possp(i)))
        hcorrhr=hcorrhr+(hnuhr(i)-hnu(possp(i))) 
        qcorrfree=qcorrfree*qunufree(i)/qus(possp(i))  
        scorrfree=scorrfree+(snufree(i)-snu(possp(i)))
        hcorrfree=hcorrfree+(hnufree(i)-hnu(possp(i))) 
!CC        zpvetot=zpvetot+zpnu(i)
      enddo
!CC      hcorrhr=hcorrhr/cm2kj
!CC      hcorrfree=hcorrfree/cm2kj
      hcorrhr=hcorrhr
      hcorrfree=hcorrfree
      qtothr=qtot*qcorrhr
      qtotfree=qtot*qcorrfree
      stothr=stot+scorrhr 
      htothr=tctotal+hcorrhr 
      stotfree=stot+scorrfree 
      htotfree=tctotal+hcorrfree 
      write(10,*)
      write(10,'(''ZPVE                        ='',f20.5,'' kJ/mol'')') 
     . zpvetot
      write(10,*)
      write(10,*)
      write(10,'(''Thermal correction in kJ/mol'')') 
      write(10,*)
      write(10,'(''TC Free rotor               ='',f20.5)') htotfree
      write(10,'(''TC Hindered rotor           ='',f20.5)') htothr
      write(10,'(''TC harmonic oscillator      ='',f20.5)') tctotal
      write(10,*)
      write(10,'(5x,''Entropies in J/(mol K)'')')
      write(10,'(''S electronic                ='',f20.5)') selec
      write(10,'(''S translational             ='',f20.5)') strans     
      write(10,'(''S rotational                ='',f20.5)') srot
      write(10,'(''S vibrational FR            ='',f20.5)') snutotfree
      write(10,'(''S vibrational HR            ='',f20.5)') snutothr
      write(10,'(''S vibrational HO            ='',f20.5)') snutot
      write(10,*)
      write(10,'(''Stotal Free rotor           ='',f20.5)') stotfree        
      write(10,'(''Stotal Hindered rotor       ='',f20.5)') stothr        
      write(10,'(''Stotal Harmonic oscillator  ='',f20.5)') stot        
      write(10,*)
      write(10,'(''TC-TdeltaS in kJ/mol        ='',f20.5)') 
     . tctotal-tk*stot/1000.0 
      write(10,*)   

      else

      write(10,*)
      write(10,'(''ZPVE                        ='',f20.5,'' kJ/mol'')') 
     . zpvetot
      write(10,*)
      write(10,'(''Thermal correction in kJ/mol'')') 
      write(10,*)
      write(10,'(''TC harmonic oscillator      ='',f20.5)') tctotal 
      write(10,*)
      write(10,'(5x,''Entropies in J/(mol K)'')')
      write(10,'(''S electronic                ='',f20.5)') selec   
      write(10,'(''S translational             ='',f20.5)') strans   
      write(10,'(''S rotational                ='',f20.5)') srot
      write(10,'(''S vibrational HO            ='',f20.5)') snutot   
      write(10,*)
      write(10,'(''Stotal Harmonic oscillator  ='',f20.5)') stot
      write(10,*)   
      write(10,'(''TC-TdeltaS in kJ/mol        ='',f20.5)') 
     . tctotal-tk*stot/1000.0 
      write(10,*)   
      endif

!CC print out infor on low modes
      do i=1,nsp
        write(10,'(''Low mode'',1x,i3,1x,''
     . ('',f8.4,2x,''1/cm )'')') i,freq(possp(i))
        write(10,'(''Reduced moment of inertia    ='',
     . f15.4,'' (amu*A**2)'',f20.4,'' (amu*bohr**2)'')') 
     . i23(i),i23(i)/(a0*a0)
        if(ipitzer.eq."n".or.ipitzer.eq."N") then  
        write(10,'(''Partition function (1/Qf)    
     . ='',f15.4)') 1/qunuhr(i)
        write(10,'(''Internal symmetry number     
     . ='',i10)') symint(i)
        else
        write(10,'(''Partition function (1/Qf)    
     . ='',f15.4)') 1/qunuhr(i)
        write(10,'(''Internal symmetry number     
     . ='',i10)') symint(i)
        write(10,'(''Rotational potential (V/RT)  
     . ='',f15.4)') v0(i)
           
        endif
      enddo  
       end program

       subroutine mominertia(natom,coord,charge,xcm,imol,
     . cosmol,imolkg,a0,numatom,mtot)

       implicit none

!CC mass from IUPAC
      real*8:: mass(50)=(/1.00794,4.002602,6.941,9.012182,10.811,
     .  12.0107,14.0067,15.9994,18.9984032,20.1797,22.989770,
     . 24.3050,26.981538,28.0855,30.973761,32.065,35.453,39.948,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,72.64,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,118.71/)
!CC mass from Gaussian
!CC       real*8:: mass(18)=(/1.00783,0.0, 0.0, 0.0,0.0,12.0,14.00307,
!CC     . 15.9941,18.9984,0.0,0.0,0.0,0.0,0.0,0.0,31.97207,
!CC     . 34.96885,0.0/)

       real*8 ixx, izz, iyy, ixz, ixy, iyz
       real*8 coord(1000,3), x2(3), imatrix(3,3)
       real*8 imol(3), cosmol(3,3),mtot,xcm(3),coordt(1000,3)
       real*8 eigtmp(3), eigvectmp(3,3),imolkg(3),a0
       integer tmpminl, tmpmaxl, tmpmidl,nrot,i,j
       integer natom, numatom(1000), charge(1000)
       integer chtmp, parameter(1000) 
       character*2 label(1000)

       open(4,file='geom.input',status='old') 
       open(9,file='moments',status='unknown')

       mtot=0.0
       i=1
  300  read(4,*,end=310) label(i),charge(i),
       coord(i,1),coord(i,2),coord(i,3)
       mtot=mtot+mass(charge(i))
       numatom(i)=i
       i=i+1
       goto 300
  310  natom=i-1
       ixx=0.0
       iyy=0.0
       izz=0.0
       ixy=0.0
       ixz=0.0
       iyz=0.0
       coord=coord/a0
       xcm=0.0
       do i=1,natom
          do j=1,3
          xcm(j)=xcm(j)+mass(charge(i))*coord(i,j)
          enddo
       enddo
       xcm=xcm/mtot
       write(9,'(''Center of mass'',5x,3(f20.6,3x))') xcm*a0
       do i=1,natom
          do j=1,3
          coordt(i,j)=coord(i,j)-xcm(j)
          enddo
       enddo
       do i=1,natom
          chtmp=charge(i)
          do j=1,3
             x2(j)=coordt(i,j)*coordt(i,j)
          enddo
          if(chtmp.eq.9999) then
          ixx=ixx+9999.*(x2(2)+x2(3))
          iyy=iyy+9999.*(x2(1)+x2(3))
          izz=izz+9999.*(x2(1)+x2(2))
          ixy=ixy-9999.*coordt(i,1)*coordt(i,2)
          ixz=ixz-9999.*coordt(i,1)*coordt(i,3)
          iyz=iyz-9999.*coordt(i,2)*coordt(i,3)
          else
          ixx=ixx+mass(chtmp)*(x2(2)+x2(3))
          iyy=iyy+mass(chtmp)*(x2(1)+x2(3))
          izz=izz+mass(chtmp)*(x2(1)+x2(2))
          ixy=ixy-mass(chtmp)*coordt(i,1)*coordt(i,2)
          ixz=ixz-mass(chtmp)*coordt(i,1)*coordt(i,3)
          iyz=iyz-mass(chtmp)*coordt(i,2)*coordt(i,3)
          endif
       enddo

        imatrix(1,1)=ixx
        imatrix(2,2)=iyy
        imatrix(3,3)=izz
        imatrix(1,2)=ixy
        imatrix(1,3)=ixz
        imatrix(2,3)=iyz
        imatrix(2,1)=imatrix(1,2)
        imatrix(3,1)=imatrix(1,3)
        imatrix(3,2)=imatrix(2,3)

        call jacobi(imatrix,3,3,imol,cosmol,nrot)

        tmpminl=minloc(imol,1)
        tmpmaxl=maxloc(imol,1)
        eigtmp(1)=imol(tmpminl)
        eigtmp(3)=imol(tmpmaxl)

        do i=1,3
           if(i.ne.tmpminl.and.i.ne.tmpmaxl) then
              eigtmp(2)=imol(i)
              tmpmidl=i
           endif
        enddo
        do j=1,3
           eigvectmp(1,j)=cosmol(tmpminl,j)
           eigvectmp(2,j)=cosmol(tmpmidl,j)
           eigvectmp(3,j)=cosmol(tmpmaxl,j)
        enddo
        imol=eigtmp
        cosmol=eigvectmp

        write(9,'(''Moments of inertia in amu*bohr^2:
     . '',5x,3(f18.5,2x))') imol(1:3)
        write(9,'(''Direction cosines'')')
        do i=1,3
           write(9,'(3(f10.6,3x))') cosmol(i,1:3)
        enddo

        imolkg=imol*1.66054202E-27*(5.29177249E-11)**2

        xcm=xcm*a0
        coord=coord*a0
        imol=imol*a0*a0
        write(9,'(''Moments of inertia in amu*A^2
     . '',2x,3(f20.5))') imol(1:3)
 
        end subroutine


       subroutine redmoment(natom,i23,coord,charge,xcm,
     . a0,imol,cosmol,freq,nsp,numatom)
             
       implicit none
       real*8 ixx, izz, iyy, ixz, ixy, iyz
       real*8 coord(1000,3), x2(3), imatrix(3,3)
       real*8 imol(3), cosmol(3,3),mtot,xcm(3),coordt(1000,3)
       real*8 xcmrf(3,2), coordrf(1000,3,2),freq
       real*8 coordcm(1000,3), coordpra(1000,3,2)
       real*8 modvecr, vecr(3), rcm(3,2), nuni(3)
       real*8 veccm(3,2), modveccm(2),mrf(2)
       real*8 xdstar(1000,3,2), itoptot(2), tmp1 
       real*8 tmp2, itopstar(3,2), i13(2), i23
       real*8 i33(2), tmp3,tmp,a0
       real*8 eigtmp(3), eigvectmp(3,3)
       integer tmpminl, tmpmaxl, tmpmidl,i,j,nsp
       integer natomrf(1000,2),chargerf(1000,2)
       integer nrot,nrf(2),atomrf(1000,2),k,m
       real*8 eigvalues(3), eigvectors(3,3)
       integer natom, numatom(1000), charge(1000)
       integer chtmp, parameter(1000)
       character*80 ired,afile
!CC mass from IUPAC
      real*8:: mass(50)=(/1.00794,4.002602,6.941,9.012182,10.811,
     .  12.0107,14.0067,15.9994,18.9984032,20.1797,22.989770,
     . 24.3050,26.981538,28.0855,30.973761,32.065,35.453,39.948,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,72.64,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,118.71/)
!CC mass from Gaussian
!CC       real*8:: mass(18)=(/1.00783,0.0, 0.0, 0.0,0.0,12.0,14.00307,
!CC     . 15.9941,18.9984,0.0,0.0,0.0,0.0,0.0,0.0,31.97207,
!CC     . 34.96885,0.0/)


       print *,"Do you want to (y) read data for the rotating fragment from file or .
     . (n) enter them manually?"
       read *,ired

       if(ired.eq."y".or.ired.eq."Y") then
       if(nsp.lt.10) then
         write(afile,'(''rotfrag.'',i1)') nsp
       else
         write(afile,'(''rotfrag.'',i2)') nsp
       endif
         
       open(8,file=trim(afile),status='unknown')
       read(8,*) nrf(1)
       read(8,*) atomrf(1:nrf(1),1)
       else
       print *,"Enter the number of atoms
     .  in the first rotating fragment"
       read *,nrf(1)
       print *,"Enter the atom numbers"
       read *,atomrf(1:nrf(1),1)
       endif

       nrf(2)=natom-nrf(1)
       if(nrf(2).lt.1) stop 'wrong number of the rotating fragments'

       call rotfrag(natom,coord,nrf,atomrf,coordrf,numatom,
     . xcmrf,natomrf,charge,xcm,mrf, chargerf)
!CC calculate Rcm

       do i=1,2
          do j=1,nrf(i)
             do k=1,3
                coordrf(j,k,i)=coordrf(j,k,i)-xcm(k)
          enddo
       enddo
       enddo
              coordpra=0.0
       do i=1,2
          do j=1,nrf(i)
             do k=1,3
                tmp=0
                do m=1,3
                   tmp=tmp+coordrf(j,m,i)*cosmol(m,3-k+1)
                enddo
              coordpra(j,k,i)=tmp
             enddo
           enddo
       enddo

      rcm=0.0
!CC      write(*,*) xcmrf
      do m=1,2
         do j=1,3
            do k=1,3
               rcm(j,m)=rcm(j,m)+xcmrf(k,m)*cosmol(k,3-j+1)
            enddo
         enddo
      enddo
       modvecr=0.0
       do j=1,3
          vecr(j)=rcm(j,2)-rcm(j,1)
          modvecr=modvecr+vecr(j)*vecr(j)
      enddo
      modvecr=sqrt(modvecr)
      nuni=vecr/modvecr
      modveccm=0.0
      do i=1,2
          veccm(1,i)=nuni(2)*rcm(3,i)-nuni(3)*rcm(2,i)
          veccm(2,i)=nuni(1)*rcm(3,i)-nuni(3)*rcm(1,i)
          veccm(3,i)=nuni(1)*rcm(2,i)-nuni(2)*rcm(1,i)
          do j=1,3
             modveccm(i)=modveccm(i)+veccm(j,i)*veccm(j,i)
          enddo
          modveccm(i)=sqrt(modveccm(i))
      enddo
!CC  check out for veccm being almost zero!
!CC  calculate I(1,3) and I(3,3)
!CC for both fragments
              itoptot=0.0
       do m=1,2
        do i=1,nrf(m)
           do j=1,3
           xdstar(i,j,m)=coordpra(i,j,m)-rcm(j,m)
           enddo
        enddo
!CC      do i=1,nrf(m)
!CC         write(*,*) coordpra(i,1:3,m)
!CC         write(*,*) xdstar(i,1:3,m)
!CC      enddo
!CC calculate I_top
        do i=1,nrf(m)
           tmp1=0.0
           tmp2=0.0
           do j=1,3
              tmp1=tmp1+xdstar(i,j,m)*xdstar(i,j,m)
              tmp2=tmp2+xdstar(i,j,m)*nuni(j)
           enddo
           if(chargerf(i,m).eq.9999) then
           itoptot(m)=itoptot(m)+9999.*(tmp1-tmp2*tmp2)
           else
           itoptot(m)=itoptot(m)+mass(chargerf(i,m))*(tmp1-tmp2*tmp2)
           endif
         enddo
!CC        write(*,*) 'itop',m,itoptot(m)
         tmp1=0.0
         do i=1,3
            tmp1=+veccm(i,m)*veccm(i,m)
         enddo
         i13(m)=itoptot(m)-tmp1*mrf(m)
         tmp1=0.0
         do i=1,3
          itopstar(i,m)=i13(m)*nuni(i)*nuni(i)/imol(i)
          tmp1=tmp1+itopstar(i,m)
         enddo
         i33(m)=i13(m)*(1-tmp1)
!CC next fragment
         enddo
         i23=1/(1/i13(1)+1/i13(2))
         write(9,*)
         write(9,'(''Reduced moments (in amu A^2) 
     . for the vibrational mode'',2x,f20.2)') freq
         write(9,*) 
         write(9,'(''I(2,3)'',2x,f20.6,2x,''
     . in A^2'')') i23
!CC         i23=i23/(a0*a0)
         write(9,'(''I(2,3)'',2x,f20.6,2x,''
     . in bohr^2'')') i23/(a0*a0)
         write(9,*)
         write(9,'(''I(1,3) for the first rot. 
     . fragment'',2x,f20.6)') i13(1)
         write(9,'(''I(3,3) for the first rot. 
     . fragment'',2x,f20.6)') i33(1)
         write(9,'(''I(1,3) for the second rot.
     . fragment'',2x,f20.6)') i13(2)
         write(9,'(''I(3,3) for the second rot.
     . fragment'',2x,f20.6)') i33(2)

       close(8)

       end subroutine

         subroutine rotfrag(natom,coord,nrf,atomrf,coordrf,
     . numatom,xcmrf,natomrf,charge,xcm,mrf,chargerf)

      implicit none

      integer natom, nrf(2), numatom(1000),atomrf(1000,2)
      integer i,j,j1,j2,k,m
      integer temp(1000),charge(1000),chargerf(1000,2)
      real*8 coord(1000,3),coordrf(1000,3,2)
      real*8 xcmrf(3,2),mrf(2),xcm(3)
      real*8:: mass(50)=(/1.00794,4.002602,6.941,9.012182,10.811,
     .  12.0107,14.0067,15.9994,18.9984032,20.1797,22.989770,
     . 24.3050,26.981538,28.0855,30.973761,32.065,35.453,39.948,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,72.64,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     . 0.0,0.0,0.0,0.0,0.0,0.0,0.0,118.71/)
      integer natomrf(1000,2)

      temp=numatom

      do m=1,2
      j1=0
      j2=0
      if(m.eq.1) then
      do i=1,nrf(m)
         do j=1,natom
            if(temp(j).eq.atomrf(i,m)) then
               j1=j1+1
               do k=1,3
                  coordrf(j1,k,m)=coord(j,k)
                  chargerf(j1,m)=charge(j)
                  natomrf(j1,m)=temp(j)
               enddo
               temp(j)=0
            endif
         enddo
       enddo
      else
      do i=1,nrf(m)
         do j=1,natom
            if(temp(j).ne.0) then
               j2=j2+1
               do k=1,3
                  coordrf(j2,k,m)=coord(j,k)
                  chargerf(j2,m)=charge(j)
                  natomrf(j2,m)=temp(j)
               enddo
            endif
         enddo
       enddo
      endif
      enddo

      xcmrf=0.0
      mrf=0.0
      do m=1,2
         do i=1,nrf(m)
            if(chargerf(i,m).eq.9999) then
            mrf(m)=mrf(m)+9999.
            else
            mrf(m)=mrf(m)+mass(chargerf(i,m))
            endif
             do j=1,3
                if(chargerf(i,m).eq.9999) then
                xcmrf(j,m)=xcmrf(j,m)+9999.*coordrf(i,j,m)
                else
                xcmrf(j,m)=xcmrf(j,m)+mass(chargerf(i,m))*coordrf(i,j,m)
                endif
         enddo
      enddo
            do i=1,3
         xcmrf(i,m)=xcmrf(i,m)/mrf(m)-xcm(i)
      enddo
      enddo
!CC      write(*,*) 'center of mass'
!CC      do i=1,2
!CC         write(*,*) mrf(i),xcmrf(1:3,i)
!CC      enddo

      end subroutine   

       SUBROUTINE jacobi(a,n,np,d,v,nrot)
       INTEGER n,np,nrot,NMAX
       REAL*8 a(np,np),d(np),v(np,np)
       real*8 l1,l2,l3,l4
       PARAMETER (NMAX=500)
!CC    Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of 
!CC    size n by n, stored in a physical np by np array. On output, elements of a above the 
!CC    diagonal are destroyed. d returns the eigenvalues of a in its first n elements. v is a 
!CC    matrix with the same logical and physical dimensions as a, whose columns contain, 
!CC    on output, the normalized eigenvectors of a. nrot returns the number of Jacobi 
!CC    rotations that were required. 

       INTEGER i,ip,iq,j 
       REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
       do  ip=1,n
!CC    Initialize to the identity matrix.
           do iq=1,n
              v(ip,iq)=0.
           enddo
           v(ip,ip)=1.
       enddo

       do ip=1,n
             b(ip)=a(ip,ip)
!CC Initialize band d to the diagonal of a.
             d(ip)=b(ip)
             z(ip)=0.
!CC This vector will accumulate terms of the form tapq as in equation(11.1.14)
       enddo

       nrot=0
       do  i=1,50
             sm=0.
             do  ip=1,n-1
!CC Sum off-diagonal elements.
                   do  iq=ip+1,n
                         sm=sm+abs(a(ip,iq))
                   enddo 
             enddo 
             if(sm.eq.0.) return
!CC The normal return, which relies on quadratic convergence to machine underflow.
             if(i.lt.4) then 
                tresh=0.2*sm/n**2 
!CC ...on the first three sweeps.
             else
                tresh=0.
!CC ...thereafter.
             endif

             do ip=1,n-1
                do iq=ip+1,n
                      g=100.*abs(a(ip,iq))
!CC After four sweeps, skip the rotation if the off-diagonal element is small.
      l1=abs(d(ip))+g
      l2=abs(d(ip))
      l3=abs(d(iq))+g
      l4=abs(d(iq))
      if((i.gt.4).and.(l1.eq.l2).and.(l3.eq.l4)) then
             a(ip,iq)=0.0
           else if(abs(a(ip,iq)).gt.tresh)then
                              h=d(iq)-d(ip)
                     if(abs(h)+g.eq.abs(h))then
                                 t=a(ip,iq)/h
                              else
                                 theta=0.5*h/a(ip,iq)
                     t=1./(abs(theta)+sqrt(1.+theta**2))
                                if(theta.lt.0.)t=-t
                              endif
                           c=1./sqrt(1+t**2)
                           s=t*c
                           tau=s/(1.+c)
                           h=t*a(ip,iq)
                           z(ip)=z(ip)-h
                           z(iq)=z(iq)+h
                           d(ip)=d(ip)-h
                           d(iq)=d(iq)+h
                           a(ip,iq)=0.
                          do j=1,ip-1
!CC  Case of rotations 1<=j<p.
                                g=a(j,ip)
                                h=a(j,iq)
                                a(j,ip)=g-s*(h+g*tau)
                               a(j,iq)=h+s*(g-h*tau)
                           enddo
                           do j=ip+1,iq-1
!CC Caseof rotations p<j<q.
                                 g=a(ip,j)
                                 h=a(j,iq)
                                 a(ip,j)=g-s*(h+g*tau)
                                 a(j,iq)=h+s*(g-h*tau)
                            enddo
                            do j=iq+1,n
!CC Case of rotations q<j<=n
                                  g=a(ip,j)
                                  h=a(iq,j)
                                  a(ip,j)=g-s*(h+g*tau)
                                  a(iq,j)=h+s*(g-h*tau)
                             enddo
                             do j=1,n
                                   g=v(j,ip)
                                   h=v(j,iq)
                                   v(j,ip)=g-s*(h+g*tau)
                                   v(j,iq)=h+s*(g-h*tau)
                             enddo
                             nrot=nrot+1
                        endif
                     enddo
                  enddo
                  do ip=1,n
                        b(ip)=b(ip)+z(ip)
                        d(ip)=b(ip)
                        z(ip)=0.
!CC Updated with the sum of tapq and reinitialize z.
                 enddo 
              enddo
              return
              END 
          
