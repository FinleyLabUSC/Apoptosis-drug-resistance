function dsdt = coreFile_DISC(t,y,params)

DISC    =   y(1);
C8	    =   y(2); 
C8a	    =	y(3); 
C3	    =	y(4); 
C3a	    =	y(5); 
XIAP	=	y(6); 
C3aXIAP	=	y(7); 
BAR	    =	y(8); 
C8aBAR  =	y(9); 
Bid     =   y(10);
tBid    =   y(11);
tBid_Bax =  y(12);
Cytcm    =  y(13);
Cytc     =  y(14);
SMACm    =  y(15);
SMAC     =  y(16);
SMACXIAP =  y(17);
Apop     =  y(18);
C9       =  y(19);
C9a      =  y(20);
C9aXIAP  =  y(21);

%reaction velocities
k1  = params(1,1);    
k2  = params(2,1);    
k3  = params(3,1);    
k4  = params(4,1); 
k4r = params(5,1); 
k5  = params(6,1);    
k6  = params(7,1);  
k6r = params(8,1);
k7  = params(9,1);  
k8  = params(10,1);  
k9  = params(11,1);  
k10 = params(12,1);  
k10r = params(13,1); 
k11  = params(14,1);  
k12  = params(15,1);  
k13  = params(16,1);  
k14  = params(17,1);  
k14r = params(18,1); 
sC8  = params(19,1);  
sC3  = params(20,1); 
sXIAP = params(21,1);  
sBAR  = params(22,1);  
sBid  = params(23,1);  
sCytcm  = params(24,1);  
sSMACm  = params(25,1); 
sC9     = params(26,1); 
dDISC   = params(27,1); 
dC8     = params(28,1); 
dC8a    = params(29,1);
dC3     = params(30,1);
dC3a    = params(31,1);
dXIAP   = params(32,1);
dC3aXIAP = params(33,1);
dBAR    = params(34,1);
dC8aBAR = params(35,1);
dBid    = params(36,1);
dtBid   = params(37,1);
dtBid_Bax = params(38,1);
dCytcm  = params(39,1);
dCytc   = params(40,1);
dSMACm  = params(41,1);
dSMAC   = params(42,1);
dSMACXIAP = params(43,1);
dApop   = params(44,1);
dC9     = params(45,1);
dC9a    = params(46,1);
dC9aXIAP = params(47,1); 
fDISC    = params(48,1);
ftBidBax = params(49,1);
fApop    = params(50,1);
Toco     = params(51,1);
Narc     = params(52,1);
Cele     = params(53,1);
beta     = params(54,1);



R1	=	k1*DISC*C8 ;%V1
R2	=	k2*C3a*C8 ;%V2
R3	=	Toco+k3*C8a*C3;%V3  Toco effect 1 of 2
R4	=	k4*C3a*XIAP - k4r*C3aXIAP; %V4 
R5	=	k5*C3a*XIAP; %V5	
R6	=	k6*C8a*BAR - (1+Narc)*k6r*C8aBAR; %V6 Narc effect 
R7	=	k7*C8a*Bid; %V7
R8	=	k8*tBid_Bax*Cytc; %V8
R9	=	k9*tBid_Bax*SMAC;	%V9
R10	=	k10*SMAC*XIAP - k10r*SMACXIAP;	%V10
R11	=	k11*Apop*C9;%V11
R12	=	k12*C3a*C9;%V12
R13	=	k13*C9a*C3; %V13
R14	=	k14*C9a*XIAP - k14r*C9aXIAP;%V14 


dsdt(1,1)	= dDISC*fDISC - DISC; %DISC 
dsdt(2,1)	= -R1-R2+sC8-dC8*C8; %C8
dsdt(3,1)	= R1+R2-R6-dC8a*C8a; %C8a
dsdt(4,1)	= -R3-R13+sC3-dC3*C3; %C3
dsdt(5,1)	= R3 -R4 +R13+ Toco*exp(-t/8.64e4) - (beta*dC3a)*C3a; %C3a Toco effect 2 of 2
dsdt(6,1)	= Cele*(-R4 - R5 - R10 - R14 - dXIAP*XIAP) + sXIAP; %XIAP Cele effect
dsdt(7,1)	= R4 - dC3aXIAP*C3aXIAP; %C3aXIAP
dsdt(8,1)	=  -R6 + sBAR - dBAR*BAR; %BAR
dsdt(9,1)	= R6 - dC8aBAR*C8aBAR; %C8aBAR
dsdt(10,1)	= -R7 + sBid - dBid*Bid; %Bid
dsdt(11,1)	= R7 - dtBid*tBid; %tBid
dsdt(12,1)  = dtBid_Bax*ftBidBax - tBid_Bax ; %tBid-Bax*
dsdt(13,1)  = 0; %Cytcm
dsdt(14,1)  = R8 - dCytc*Cytc; %Cytc
dsdt(15,1)  = -R9 + sSMACm - dSMACm*SMACm; %SMACm
dsdt(16,1)  = R9 - R10 - dSMAC*SMAC; %SMAC
dsdt(17,1)  = R10 - dSMACXIAP*SMACXIAP;%SMACXIAP
dsdt(18,1)  = dApop*fApop - Apop;%Apop*
dsdt(19,1)  = -R11 - R12 + sC9 - dC9*C9; %C9 
dsdt(20,1)  = R11 + R12 - R14 - dC9a*C9a;%C9a
dsdt(21,1)  = R14 - dC9aXIAP*C9aXIAP;%C9aXIAP




return
