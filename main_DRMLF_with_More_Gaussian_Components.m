 clear all; close all;

L = 5000;
p = 5;
type=1;
lamda_New_ATC=0.986;
lamda_New_CTA=0.986;
lamda_New_ATC_QPD=0.985;
lamda_New_CTA_QPD=0.985;
lamda_New_ATC_6PD=0.986;
lamda_New_CTA_6PD=0.985;
lamda_New_ATC_8PD=0.987;
lamda_New_CTA_8PD=0.985;

orders=2;
orders1=4;
orders2=6;
orders3=8;
num=20;    %number of node
%% Initial
load('C.mat');


UU=[];
for i=1:num
    UU=[UU randn(p,L)];
end
wo1 = randn(p,1);
wo = [ kron(wo1, ones(1,L)) ];

  %% Mean by EM
    MIU_EM(1,:)=[0.1,0];
    MIU_EM(2,:)=[0.2,0.3];
    MIU_EM(3,:)=[-0.1,0];
    MIU_EM(4,:)=[0,0.2];

    MIU_EM2(1,:)=[-0.23,-0.08,0.11,0.15];
    MIU_EM2(2,:)=[-0.1,-0.1,0.1,0.07];
    MIU_EM2(3,:)=[-0.12,-0.1,0.1,0.06];
    MIU_EM2(4,:)=[-0.1,-0.08,0.08,0.1];

    MIU_EM3(1,:)=[-0.15,-0.11,0.04,0.13,0.15,0.03];
    MIU_EM3(2,:)=[-0.1,-0.09,0.03,0.12,0.14,0.03];
    MIU_EM3(3,:)=[-0.09,-0.11,0.05,0.12,0.1,0.01];
    MIU_EM3(4,:)=[-0.09,-0.12,0.02,0.12,0.12,0];

    MIU_EM4(1,:)=[-0.21,-0.12,-0.1,-0.03,0.01,0.11,0.2,0.25];
    MIU_EM4(2,:)=[-0.19,-0.1,-0.08,-0.02,0.03,0.12,0.16,0.22];
    MIU_EM4(3,:)=[-0.2,-0.09,-0.11,0,0.02,0.14,0.18,0.25];
    MIU_EM4(4,:)=[-0.18,-0.12,-0.12,0,0.02,0.12,0.2,0.22];


    MIU_EM(5,:)=[0.2,0.3];
    MIU_EM(6,:)=[-0.1,0.2];
    MIU_EM(7,:)=[0.2,0.2];
    MIU_EM(8,:)=[-0.3,0.2];

    MIU_EM2(5,:)=[-0.1,-0.06,0.03,0.05];
    MIU_EM2(6,:)=[-0.06,-0.08,0.02,0];
    MIU_EM2(7,:)=[-0.11,-0.05,0.02,0.01];
    MIU_EM2(8,:)=[-0.05,-0.04,0.02,0.03];

    MIU_EM3(5,:)=[-0.08,-0.01,0.03,0.03,0.07,0.10];
    MIU_EM3(6,:)=[-0.07,-0.02,0.02,0.03,0.06,0.14];
    MIU_EM3(7,:)=[-0.05,0,0.01,0.03,0.1,0.11];
    MIU_EM3(8,:)=[-0.07,0,0.01,0.03,0.02,0.11];

    MIU_EM4(5,:)=[-0.07,-0.04,-0.01,0,0.03,0.05,0.11,0.13];
    MIU_EM4(6,:)=[-0.05,-0.02,-0.02,0.03,0.05,0.02,0.15,0.16];
    MIU_EM4(7,:)=[-0.05,-0.03,-0.03,0.09,0.03,0.05,0.1,0.14];
    MIU_EM4(8,:)=[-0.04,-0.08,-0.01,0.1,0.06,0.04,0.1,0.12];

    MIU_EM(9,:)=[-2.7,3.1];
    MIU_EM(10,:)=[-1.8,3.1];
    MIU_EM(11,:)=[-1.2,2.3];
    MIU_EM(12,:)=[-1.6,2.7];

    MIU_EM2(9,:)=[-7.4 -3.1 2.3 7.2];
    MIU_EM2(10,:)=[-4.1  -1.5   1.5   5];
    MIU_EM2(11,:)=[-2.5   0.2    3.1  3.9];
    MIU_EM2(12,:)=[-4.0   -2.2   1.2   4.8];

    MIU_EM3(9,:)=[-9.3  -7.1  -2.6  0.2    3.8    7.6];
    MIU_EM3(10,:)=[-3.1   -4.5   -1.5    2.2   6.2  8.3];
    MIU_EM3(11,:)=[-2.5  -1.2 0.3   1.2 3.1  4.2];
    MIU_EM3(12,:)=[ -3.2   -3.8   -1   1.1  3.8  8.2];


    MIU_EM4(9,:)=[-8.5  -6.1  -3.7   0.3   4.4   8.2   10.5   12.1];
    MIU_EM4(10,:)=[-4.1 -2.8 -1.5  0.2 1.5  3.3  4.8 9.1];
    MIU_EM4(11,:)=[-4.2  -3.6 -1.6    -0.1    1    2    3.7  5.2];
    MIU_EM4(12,:)=[-4.4   -3.7   -2.1   -0.7   2.7    4.3   6.5   7.6];

    MIU_EM(13,:)=[-2.4,2.7];
    MIU_EM(14,:)=[-5.2,5];
    MIU_EM(15,:)=[-7.2,7.6];
    MIU_EM(16,:)=[-10.2,9.9];

    MIU_EM2(13,:)=[-4.4   -2.6    2.8    4.5];
    MIU_EM2(14,:)=[-7.6   -3.1    3.2    8.0 ];
    MIU_EM2(15,:)=[-12.3   -5.5    5.2   13.2];
    MIU_EM2(16,:)=[-17.2   -6.5    7.2  16.4];

    MIU_EM3(13,:)=[-4.5  -3.1 -2.1    1.9  2.7  4.5];
    MIU_EM3(14,:)=[-7.8  -5.1 -3.2    3.1  5.2  8.3 ];
    MIU_EM3(15,:)=[-13.4  -9.2 -4.8    5.4 8.6  128];
    MIU_EM3(16,:)=[-17.5  -12.1 -7.2    7.1 12.2 17.3];

    MIU_EM4(13,:)=[-4.3   -3.1  -2.1  -0.5 0.7   1.6   2.9  4.5];
    MIU_EM4(14,:)=[-8.2  -4.8 -3 -1.5 1.5  3  5  8 ];
    MIU_EM4(15,:)=[-13.2  -8.8 -5  -2 2  5 9  13];
    MIU_EM4(16,:)=[-17.1  -12.2 -7.1 -2 2   6.9 12.5 17.8];

    MIU_EM(17,:)=[-4.8,5.3];
    MIU_EM(18,:)=[-4.5,5.3];
    MIU_EM(19,:)=[-8.3,8.1];
    MIU_EM(20,:)=[-7.6,8.1];

    MIU_EM2(17,:)=[-5.5,-6.1,4.6,5.3];
    MIU_EM2(18,:)=[-5.1,-6.1,4.8,5.1];
    MIU_EM2(19,:)=[-7.8,-7.9,8.2,8.1];
    MIU_EM2(20,:)=[-7.6,-7.9,8.2,8.3];

    MIU_EM3(17,:)=[-5.3,-6.4,-5.2,5.7,4.4,6.1];
    MIU_EM3(18,:)=[-5,-6.2,-5.2,5.5,4.7,6];
    MIU_EM3(19,:)=[-7.6,-7.7,-7.9,8.3,8.2,8.3,];
    MIU_EM3(20,:)=[-7.8,-7.9,-7.7,8.1,8.3,8.2];

    MIU_EM4(17,:)=[-5.1,-4.8,-4.7,-6.2,6.3,5.1,4.8,5];
    MIU_EM4(18,:)=[-5,-4.8,-4.8,-6,6.3,5.3,4.6,5];
    MIU_EM4(19,:)=[-8.2,-8,-8,-8,7.8,8.4,7.7,7.9];
    MIU_EM4(20,:)=[-8.3,-8.1,-8,-8,8,8.1,8,8.2];
    %% Sigma by EM

    SIGMA_EM(1,:)=[0.6,4.6];
    SIGMA_EM(2,:)=[0.7,8.2];
    SIGMA_EM(3,:)=[1.2,11.3];
    SIGMA_EM(4,:)=[0.7,15.4];

    SIGMA_EM2(1,:)=[0.5,0.6,8.9,11.3];
    SIGMA_EM2(2,:)=[0.5,11.8,13.9,20.6];
    SIGMA_EM2(3,:)=[0.8,10.3,15.3,20.7];
    SIGMA_EM2(4,:)=[0.6,12.3,14.8,21.6];

    SIGMA_EM3(1,:)=[0.54,0.58,0.45,0.5,9.32,10.82];
    SIGMA_EM3(2,:)=[0.3,1,5,9,15,20];
    SIGMA_EM3(3,:)=[0.6,1,5,8,17,21];
    SIGMA_EM3(4,:)=[0.8,1.5,6,11,14,22];

    SIGMA_EM4(1,:)=[0.47,0.62,0.69,0.49,9.64,11.2,10.4,11.5];
    SIGMA_EM4(2,:)=[0.5,0.7,1,7,11,15,20,25];
    SIGMA_EM4(3,:)=[0.8,0.9,1,6,10,13,20,24];
    SIGMA_EM4(4,:)=[0.6,0.8,1,7,11,15,22,24];

    SIGMA_EM(5,:)=[1.2,1.1];
    SIGMA_EM(6,:)=[1.8,2.1];
    SIGMA_EM(7,:)=[3,2.9];
    SIGMA_EM(8,:)=[3.7,3.9];

    SIGMA_EM2(5,:)=[0.9,1,1.2,1];
    SIGMA_EM2(6,:)=[2.1,2.3,2.1,2.1];
    SIGMA_EM2(7,:)=[3.2,3,3.1,2.8];
    SIGMA_EM2(8,:)=[4.1,3.8,4.1,4];


    SIGMA_EM3(5,:)=[0.8,1.1,0.9,1,0.9,1];
    SIGMA_EM3(6,:)=[2.1,1.9,1.8,2,2,2.2];
    SIGMA_EM3(7,:)=[2.8,3,2.9,3,3.1,2.6];
    SIGMA_EM3(8,:)=[4.2,3.9,4,4.1,4.1,4.1];

    SIGMA_EM4(5,:)=[1.1,1,0.9,0.9,1.1,0.9,1,1];
    SIGMA_EM4(6,:)=[1.8,2.3,1.7,2,2,1.8,2,2];
    SIGMA_EM4(7,:)=[3.2,2.7,3,3.2,3,3.1,3,3.1];
    SIGMA_EM4(8,:)=[4.1,3.9,3.8,4.1,3.9,4.1,4.1,4.1];

    SIGMA_EM(9,:)=[3.1,5.3];
    SIGMA_EM(10,:)=[2.2,4.3];
    SIGMA_EM(11,:)=[1.9,3];
    SIGMA_EM(12,:)=[1.9,3.8];

    SIGMA_EM2(9,:)=[1.4   2.5  3.2  5];
    SIGMA_EM2(10,:)=[1.1    1.4   2.6    3.8];
    SIGMA_EM2(11,:)=[1.6    1.7    2.2    3.2];
    SIGMA_EM2(12,:)=[1.1   1.6   2.4    4.2];

    SIGMA_EM3(9,:)=[0.6 1.1 1.9 2.4 3.6 4.8];
    SIGMA_EM3(10,:)=[1.3    0.6   1.5    2.1  1.8   4.5];
    SIGMA_EM3(11,:)=[1.6   1.1 1.6  2  1.8    2.8];
    SIGMA_EM3(12,:)=[1.1   0.7   1.1  1.9   3.1 5.2];

    SIGMA_EM4(9,:)=[0.8    1.2    1.5    2.1    3.2    2.5    0.5  2.9];
    SIGMA_EM4(10,:)=[1.3   1.4 0.5  1.5 1.5    2.2  2.1  4.5];
    SIGMA_EM4(11,:)=[1.5  1.2  1.5 1.3  1.5 2.1   2.3 3.1];
    SIGMA_EM4(12,:)=[1.2 0.6   1.1   1.1   1.5  0.5  3.3  4.2];

    SIGMA_EM(13,:)=[1.5,1.7];
    SIGMA_EM(14,:)=[2.5,3];
    SIGMA_EM(15,:)=[4.1,4.7];
    SIGMA_EM(16,:)=[5.5,6.1];

    SIGMA_EM2(13,:)=[0.8   1.8   1.5   0.9];
    SIGMA_EM2(14,:)=[1.3   2.6   2.3   1.2];
    SIGMA_EM2(15,:)=[1.5   4   4   1.5];
    SIGMA_EM2(16,:)=[2.1   5.2    4.6    2.2];

    SIGMA_EM3(13,:)=[0.8   1.2 1.1   1.9 2.1  0.5];
    SIGMA_EM3(14,:)=[0.5   1 2.5   2.6 1.1  0.7];
    SIGMA_EM3(15,:)=[1.6  1.8  4.1   4.3 2.1  1.6];
    SIGMA_EM3(16,:)=[2.2  3.2 4.6    5.2   3.3 2.1];

    SIGMA_EM4(13,:)=[0.4 0.8 1.1 1.7   1.6 1 0.7 0.5];
    SIGMA_EM4(14,:)=[1.2 1.5 2  2.5   2.5 2 1.5 0.9 ];
    SIGMA_EM4(15,:)=[1.5 2.2 3.1 3.7   4.2  3 2 1.5 ];
    SIGMA_EM4(16,:)=[2.2 3.7 4.1 5.1   4.8  4 3.2  1.9];



    SIGMA_EM(17,:)=[2.2,1.8];
    SIGMA_EM(18,:)=[2,1.7];
    SIGMA_EM(19,:)=[2.9,3];
    SIGMA_EM(20,:)=[3.2,2.7];

    SIGMA_EM2(17,:)=[1.8,2.2,2.4,1.6];
    SIGMA_EM2(18,:)=[1.9,2.2,2.2,1.7];
    SIGMA_EM2(19,:)=[3.1,2.8,2.9,3.1];
    SIGMA_EM2(20,:)=[3.1,3,3.2,2.8];

    SIGMA_EM3(17,:)=[2.8,1.6,2.4,1.9,2.3,3];
    SIGMA_EM3(18,:)=[2.5,1.9,2.4,1.9,2.1,2.5];
    SIGMA_EM3(19,:)=[3.1,3.4,2.9,3,3,2.6];
    SIGMA_EM3(20,:)=[3.2,2.9,3,2.8,3,3.1];

    SIGMA_EM4(17,:)=[2.1,2.3,1.4,1.8,2.4,3.1,2.2,2.4];
    SIGMA_EM4(18,:)=[2.3,1.9,1.5,1.9,2.1,3,2.2,2.5];
    SIGMA_EM4(19,:)=[3.2,3,3.3,3,3,2.7,2.8,3];
    SIGMA_EM4(20,:)=[3.1,2.8,3,2.9,3,2.9,3.2,3.1];
    %% Weight by EM
    WEIGHT_EM(1,:)=[0.92,0.08];
    WEIGHT_EM(2,:)=[0.89,0.11];
    WEIGHT_EM(3,:)=[0.84,0.16];
    WEIGHT_EM(4,:)=[0.8,0.2];
    
    WEIGHT_EM2(1,:)=[0.42,0.45, 0.07,0.06];
    WEIGHT_EM2(2,:)=[0.86,0.05, 0.05,0.04];
    WEIGHT_EM2(3,:)=[0.83,0.07, 0.05,0.05];
    WEIGHT_EM2(4,:)=[0.87,0.03, 0.06,0.04];

    WEIGHT_EM3(1,:)=[0.41,0.24,0.19,0.07,0.06,0.03];
    WEIGHT_EM3(2,:)=[0.65,0.12,0.08,0.05,0.05,0.05];
    WEIGHT_EM3(3,:)=[0.7,0.1,0.06,0.05,0.04,0.05];
    WEIGHT_EM3(4,:)=[0.6,0.2,0.08,0.02,0.03,0.02];

    WEIGHT_EM4(1,:)=[0.35,0.21,0.2,0.15,0.04,0.02,0.02,0.01];
    WEIGHT_EM4(2,:)=[0.6,0.1,0.1,0.04,0.04,0.05,0.03,0.04];
    WEIGHT_EM4(3,:)=[0.45,0.15,0.15,0.09,0.04,0.04,0.04,0.04];
    WEIGHT_EM4(4,:)=[0.5,0.2,0.1,0.04,0.04,0.04,0.05,0.03];



    WEIGHT_EM(5,:)=[0.52,0.48];
    WEIGHT_EM(6,:)=[0.46,0.54];
    WEIGHT_EM(7,:)=[0.5,0.5];
    WEIGHT_EM(8,:)=[0.51,0.49];
    WEIGHT_EM(9,:)=[0.5,0.5];
    WEIGHT_EM(10,:)=[0.5,0.5];
    WEIGHT_EM(11,:)=[0.53,0.47];
    WEIGHT_EM(12,:)=[0.6,0.4];
    WEIGHT_EM(13,:)=[0.52,0.48];
    WEIGHT_EM(14,:)=[0.5,0.5];
    WEIGHT_EM(15,:)=[0.53,0.47];
    WEIGHT_EM(16,:)=[0.46,0.54];
    WEIGHT_EM(17,:)=[0.51,0.49];
    WEIGHT_EM(18,:)=[0.53,0.47];
    WEIGHT_EM(19,:)=[0.5,0.5];
    WEIGHT_EM(20,:)=[0.45,0.55];
    
    WEIGHT_EM2(5,:)=[0.24,0.25,0.26,0.25];
    WEIGHT_EM2(6,:)=[0.24,0.28,0.22,0.26];
    WEIGHT_EM2(7,:)=[0.22,0.28,0.26,0.24];
    WEIGHT_EM2(8,:)=[0.25,0.25,0.25,0.25];
    WEIGHT_EM2(9,:)=[0.15   0.35   0.4   0.1];
    WEIGHT_EM2(10,:)=[0.15   0.35    0.35   0.1];
    WEIGHT_EM2(11,:)=[0.26    0.49    0.18    0.07];
    WEIGHT_EM2(12,:)=[0.1   0.3    0.4    0.2];
    WEIGHT_EM2(13,:)=[0.14    0.35   0.34   0.17];
    WEIGHT_EM2(14,:)=[0.15   0.35    0.35    0.15];
    WEIGHT_EM2(15,:)=[0.1   0.42  0.38   0.1];
    WEIGHT_EM2(16,:)=[0.12   0.41  0.38   0.09];
    WEIGHT_EM2(17,:)=[0.35,0.15,0.27,0.23];
    WEIGHT_EM2(18,:)=[0.3,0.2,0.25,0.25];
    WEIGHT_EM2(19,:)=[0.24,0.26,0.27,0.23];
    WEIGHT_EM2(20,:)=[0.22,0.3,0.2,0.28];

    WEIGHT_EM3(5,:)=[0.16,0.17,0.16,0.15,0.16,0.16];
    WEIGHT_EM3(6,:)=[0.17,0.15,0.16,0.15,0.17,0.16];
    WEIGHT_EM3(7,:)=[0.2,0.16,0.17,0.13,0.15,0.16];
    WEIGHT_EM3(8,:)=[0.16,0.18,0.16,0.14,0.16,0.16];
    WEIGHT_EM3(9,:)=[0.05    0.05    0.25    0.35    0.25    0.05];
    WEIGHT_EM3(10,:)=[0.21    0.05    0.34    0.28  0.07    0.05];
    WEIGHT_EM3(11,:)=[0.22   0.08 0.45  0.1  0.1    0.05];
    WEIGHT_EM3(12,:)=[0.2   0.05   0.25    0.3    0.15  0.05];
    WEIGHT_EM3(13,:)=[0.06   0.16    0.25   0.31    0.14    0.08];
    WEIGHT_EM3(14,:)=[0.1   0.15    0.25   0.25    0.15    0.1];
    WEIGHT_EM3(15,:)=[0.05   0.17    0.28   0.3    0.15    0.05];
    WEIGHT_EM3(16,:)=[0.05   0.15    0.31   0.29    0.15    0.05];
    WEIGHT_EM3(17,:)=[0.18,0.2,0.16,0.12,0.16,0.18];
    WEIGHT_EM3(18,:)=[0.17,0.2,0.18,0.12,0.17,0.16];
    WEIGHT_EM3(19,:)=[0.16,0.21,0.17,0.12,0.16,0.18];
    WEIGHT_EM3(20,:)=[0.17,0.2,0.16,0.12,0.16,0.19];


    WEIGHT_EM4(5,:)=[0.13,0.12,0.12,0.13,0.12,0.12,0.14,0.12];
    WEIGHT_EM4(6,:)=[0.11,0.12,0.12,0.15,0.12,0.1,0.14,0.14];
    WEIGHT_EM4(7,:)=[0.13,0.11,0.13,0.13,0.13,0.12,0.13,0.12];
    WEIGHT_EM4(8,:)=[0.12,0.13,0.12,0.13,0.12,0.13,0.13,0.12];
    WEIGHT_EM4(9,:)=[0.05   0.1    0.2   0.3    0.3    0.02   0.02   0.01];
    WEIGHT_EM4(10,:)=[0.14  0.06  0.07    0.28 0.03    0.25 0.07  0.05];
    WEIGHT_EM4(11,:)=[0.05    0.1    0.3   0.15    0.15    0.1    0.1    0.05];
    WEIGHT_EM4(12,:)=[0.06    0.12    0.27    0.14    0.26    0.05    0.05   0.05];
    WEIGHT_EM4(13,:)=[0.06 0.09   0.15    0.25   0.21    0.19    0.11 0.04];
    WEIGHT_EM4(14,:)=[0.08 0.1   0.12    0.2   0.2    0.15    0.1 0.05];
    WEIGHT_EM4(15,:)=[0.05 0.09   0.15    0.2   0.2    0.15    0.1 0.06];
    WEIGHT_EM4(16,:)=[0.03 0.11   0.16    0.21   0.22    0.13    0.09 0.05];
    WEIGHT_EM4(17,:)=[0.21,0.12,0.09,0.12,0.11,0.17,0.1,0.08];
    WEIGHT_EM4(18,:)=[0.2,0.13,0.08,0.13,0.12,0.16,0.09,0.09];
    WEIGHT_EM4(19,:)=[0.21,0.13,0.08,0.11,0.11,0.18,0.1,0.08];
    WEIGHT_EM4(20,:)=[0.2,0.15,0.09,0.1,0.1,0.16,0.12,0.08];


    MIU_Rand=MIU_EM;
    SIGMA_Rand=SIGMA_EM;
    WEIGHT_Rand=WEIGHT_EM;

    MIU_Rand2=MIU_EM2;
    SIGMA_Rand2=SIGMA_EM2;
    WEIGHT_Rand2=WEIGHT_EM2;

    MIU_Rand3=MIU_EM3;
    SIGMA_Rand3=SIGMA_EM3;
    WEIGHT_Rand3=WEIGHT_EM3;

    MIU_Rand4=MIU_EM4;
    SIGMA_Rand4=SIGMA_EM4;
    WEIGHT_Rand4=WEIGHT_EM4;

tic
P=50;  %% Independent Simulation
for mm = 1 : P

    un1=zeros(p,num);
    dd=zeros(num,L);
    Fai=zeros(p,num);
    Fai_QPD=zeros(p,num);
    Fai_cta_QPD=Fai;
    Fai_cta=Fai;
    FaiR=zeros(p,num);
    FaiR_cta=Fai;
    Fai_6PD=Fai;
    Fai_8PD=Fai;


       %% Noise
    %% Impulse Noise
    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV1(n)=normrnd(0,0.5);
        else
            VV1(n)=normrnd(0,5);
        end
    end

    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV2(n)=normrnd(0,0.5);
        else
            VV2(n)=normrnd(0,10);
        end
    end
 
    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV3(n)=normrnd(0,1);
        else
            VV3(n)=normrnd(0,10);
        end
    end

    VV4=alpha_noise(1.97,0,2,0,1,L);

    v(1,:)=VV1;
    v(2,:)=VV2;
    v(3,:)=VV3;
    v(4,:)=VV4;

    %% Gaussian Noise
    v(5,:)=1*randn(1,L);
    v(6,:)=2*randn(1,L);
    v(7,:)=3*randn(1,L);
    v(8,:)=4*randn(1,L);
    %% Skewed Noise
    v(9,:)=raylrnd(8,1,L);
    v(10,:)=gamrnd(3,2,1,L);
    v(11,:)=poissrnd(6,1,L);
    v(12,:)=chi2rnd(6,1,L);
    %% Uniform Noise
    v(13,:)=10*rand(1,L)-5;
    v(14,:)=20*rand(1,L)-10;
    v(15,:)=30*rand(1,L)-15;
    v(16,:)=40*rand(1,L)-20;
    %% Multi-peak Noise
    for n=1:L
        Epsilon=rand;
        if Epsilon<0.5
            VV17(n)=normrnd(-5,2);
        else
            VV17(n)=normrnd(5,2);
        end
    end

    mu=-5;
    sigma=4;
    b=sigma/sqrt(2);
    mu1=5;
    sigma1=4;
    b1=sigma1/sqrt(2);
    weight_laplace=[0.5 0.5];

    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_laplace(1)
            a=rand(1,1)-0.5;
            VV18(n)=mu-b*sign(a).*log(1-2*abs(a));
        else
            a=rand(1,1)-0.5;
            VV18(n)=mu1-b1*sign(a).*log(1-2*abs(a));
        end
    end

    for n=1:L
        Epsilon=rand;
        if Epsilon<0.5
            VV19(n)=normrnd(-8,3);
        else
            VV19(n)=normrnd(8,3);
        end
    end

    mu=-8;
    sigma=9; 
    b=sigma/sqrt(2);
    mu1=8;
    sigma1=9;
    b1=sigma1/sqrt(2);
    weight_laplace=[0.5 0.5];

    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_laplace(1)
            a=rand(1,1)-0.5;
            VV20(n)=mu-b*sign(a).*log(1-2*abs(a));
        else
            a=rand(1,1)-0.5;
            VV20(n)=mu1-b1*sign(a).*log(1-2*abs(a));
        end
    end

    v(17,:)=VV17;
    v(18,:)=VV18;
    v(19,:)=VV19;
    v(20,:)=VV20;



    for i=1:num
        vv(i,:)=v(i,:)-mean(v(i,:));
        Sig(i)=var (vv(i,1:L));
    end


    for i=1:num
        for ii = 1 : L
            dd(i,ii) = wo(:,ii)' * UU(:,(i-1)*L+ii) + vv(i,ii);
        end
    end
  
    C=Weight(num,C,type,Sig,orders);

    %% Iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ATC
    w_ATC=randn(p,num);
    w_ATC_RLS=w_ATC;
    w_ATC_RLS_BPD=w_ATC;
    w_ATC_RLS_QPD=w_ATC;
    w_ATC_RLS_6PD=w_ATC;
    w_ATC_RLS_8PD=w_ATC;
    for cnt=1:num
        Pn_RLS_BPD(:,:,cnt) = eye(p);
    end
    Pn_RLS_QPD =Pn_RLS_BPD;
    Pn_RLS_6PD =Pn_RLS_BPD;
    Pn_RLS_8PD =Pn_RLS_BPD;
    %% CTA
    w_CTA=w_ATC;
    w_CTA_RLS=w_ATC;
    w_CTA_RLS_BPD=w_ATC;
    w_CTA_RLS_QPD=w_ATC;
    w_CTA_RLS_6PD=w_ATC;
    w_CTA_RLS_8PD=w_ATC;
    w_M_RLS=w_ATC;
    w_M_CTA_RLS=w_ATC;
    Pn_RLS_BPD_cta=Pn_RLS_BPD;
    Pn_RLS_QPD_cta=Pn_RLS_BPD;
    Pn_RLS_6PD_cta=Pn_RLS_BPD;
    Pn_RLS_8PD_cta=Pn_RLS_BPD;
    %% ATC
    for ii = 1 : L
        AddR=zeros(p,num);
        AddR_QPD=zeros(p,num);
        AddR_6PD=zeros(p,num);
        AddR_8PD=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            Err_ATC_RLS_BPD(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_BPD(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_BPD(:,cnt));
            Err_ATC_RLS_QPD(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_QPD(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_QPD(:,cnt));
            Err_ATC_RLS_6PD(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_6PD(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_6PD(:,cnt));
            Err_ATC_RLS_8PD(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_8PD(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_8PD(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);%
            %% ATC-RLS-BPD
            ek_Rand1(cnt)=dn(cnt)-w_ATC_RLS_BPD(:,cnt)'*un1(:,cnt);
            for k=1:orders
                P_Rand(ii,k+orders*(cnt-1) )=exp(-1*(ek_Rand1(cnt)-MIU_Rand(cnt,k))^2/(2*SIGMA_Rand(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand(cnt,k));
            end
            for k=1:orders
                V_Rand(ii,k+orders*(cnt-1) )=WEIGHT_Rand(cnt,k)*P_Rand(ii,k+orders*(cnt-1))/(WEIGHT_Rand(cnt,:)*P_Rand(ii,orders*(cnt-1)+1:orders*(cnt-1)+orders)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders
                R1=R1+V_Rand(ii,k+orders*(cnt-1))/SIGMA_Rand(cnt,k)^2;
                R2=R2+V_Rand(ii,k+orders*(cnt-1))*MIU_Rand(cnt,k)/SIGMA_Rand(cnt,k)^2;
            end
            kn(:,:,cnt) = Pn_RLS_BPD(:,:,cnt) * un1(:,cnt) / ( lamda_New_ATC+ R1*un1(:,cnt)' * Pn_RLS_BPD(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_BPD(:,:,cnt) = 1/lamda_New_ATC * ( Pn_RLS_BPD(:,:,cnt) - R1*kn(:,:,cnt) * un1(:,cnt)' * Pn_RLS_BPD(:,:,cnt));
            FaiR(:,cnt) = w_ATC_RLS_BPD(:,cnt) +kn(:,:,cnt) * (R1*ek_Rand1(cnt)-R2);
            %% ATC-RLS-QPD
            ek_Rand2(cnt)=dn(cnt)-w_ATC_RLS_QPD(:,cnt)'*un1(:,cnt);

            for k=1:orders1
                P_Rand2(ii,k+orders1*(cnt-1) )=exp(-1*(ek_Rand2(cnt)-MIU_Rand2(cnt,k))^2/(2*SIGMA_Rand2(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand2(cnt,k));
            end
            for k=1:orders1
                V_Rand2(ii,k+orders1*(cnt-1) )=WEIGHT_Rand2(cnt,k)*P_Rand2(ii,k+orders1*(cnt-1))/(WEIGHT_Rand2(cnt,:)*P_Rand2(ii,orders1*(cnt-1)+1:orders1*(cnt-1)+orders1)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders1
                R1=R1+V_Rand2(ii,k+orders1*(cnt-1))/SIGMA_Rand2(cnt,k)^2;
                R2=R2+V_Rand2(ii,k+orders1*(cnt-1))*MIU_Rand2(cnt,k)/SIGMA_Rand2(cnt,k)^2;
            end
            kn2(:,:,cnt) = Pn_RLS_QPD(:,:,cnt) * un1(:,cnt) / (lamda_New_ATC_QPD+ R1*un1(:,cnt)' * Pn_RLS_QPD(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_QPD(:,:,cnt) = 1/lamda_New_ATC_QPD * ( Pn_RLS_QPD(:,:,cnt) - R1*kn2(:,:,cnt) * un1(:,cnt)' * Pn_RLS_QPD(:,:,cnt));
            Fai_QPD(:,cnt) = w_ATC_RLS_QPD(:,cnt) +kn2(:,:,cnt) * (R1*ek_Rand2(cnt)-R2);
            %% ATC-RLS-6PD
            ek_Rand3(cnt)=dn(cnt)-w_ATC_RLS_6PD(:,cnt)'*un1(:,cnt);

            for k=1:orders2
                P_Rand3(ii,k+orders2*(cnt-1) )=exp(-1*(ek_Rand3(cnt)-MIU_Rand3(cnt,k))^2/(2*SIGMA_Rand3(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand3(cnt,k));
            end
            for k=1:orders2
                V_Rand3(ii,k+orders2*(cnt-1) )=WEIGHT_Rand3(cnt,k)*P_Rand3(ii,k+orders2*(cnt-1))/(WEIGHT_Rand3(cnt,:)*P_Rand3(ii,orders2*(cnt-1)+1:orders2*(cnt-1)+orders2)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders2
                R1=R1+V_Rand3(ii,k+orders2*(cnt-1))/SIGMA_Rand3(cnt,k)^2;
                R2=R2+V_Rand3(ii,k+orders2*(cnt-1))*MIU_Rand3(cnt,k)/SIGMA_Rand3(cnt,k)^2;
            end
            kn3(:,:,cnt) = Pn_RLS_6PD(:,:,cnt) * un1(:,cnt) / (lamda_New_ATC_6PD+ R1*un1(:,cnt)' * Pn_RLS_6PD(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_6PD(:,:,cnt) = 1/lamda_New_ATC_6PD * ( Pn_RLS_6PD(:,:,cnt) - R1*kn3(:,:,cnt) * un1(:,cnt)' * Pn_RLS_6PD(:,:,cnt));
            Fai_6PD(:,cnt) = w_ATC_RLS_6PD(:,cnt) +kn3(:,:,cnt) * (R1*ek_Rand3(cnt)-R2);
            %% ATC-RLS-8PD
            ek_Rand4(cnt)=dn(cnt)-w_ATC_RLS_8PD(:,cnt)'*un1(:,cnt);

            for k=1:orders3
                P_Rand4(ii,k+orders3*(cnt-1) )=exp(-1*(ek_Rand4(cnt)-MIU_Rand4(cnt,k))^2/(2*SIGMA_Rand4(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand4(cnt,k));
            end
            for k=1:orders3
                V_Rand4(ii,k+orders3*(cnt-1) )=WEIGHT_Rand4(cnt,k)*P_Rand4(ii,k+orders3*(cnt-1))/(WEIGHT_Rand4(cnt,:)*P_Rand4(ii,orders3*(cnt-1)+1:orders3*(cnt-1)+orders3)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders3
                R1=R1+V_Rand4(ii,k+orders3*(cnt-1))/SIGMA_Rand4(cnt,k)^2;
                R2=R2+V_Rand4(ii,k+orders3*(cnt-1))*MIU_Rand4(cnt,k)/SIGMA_Rand4(cnt,k)^2;
            end
            kn4(:,:,cnt) = Pn_RLS_8PD(:,:,cnt) * un1(:,cnt) / (lamda_New_ATC_8PD+ R1*un1(:,cnt)' * Pn_RLS_8PD(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_8PD(:,:,cnt) = 1/lamda_New_ATC_8PD * ( Pn_RLS_8PD(:,:,cnt) - R1*kn4(:,:,cnt) * un1(:,cnt)' * Pn_RLS_8PD(:,:,cnt));
            Fai_8PD(:,cnt) = w_ATC_RLS_8PD(:,cnt) +kn4(:,:,cnt) * (R1*ek_Rand4(cnt)-R2);
        end
        for cnt=1:num
            for jj=1:num
                AddR(:,cnt)=C(cnt,jj)*FaiR(:,jj)+AddR(:,cnt);
                AddR_QPD(:,cnt)=C(cnt,jj)* Fai_QPD(:,jj)+AddR_QPD(:,cnt);
                AddR_6PD(:,cnt)=C(cnt,jj)* Fai_6PD(:,jj)+AddR_6PD(:,cnt);
                AddR_8PD(:,cnt)=C(cnt,jj)* Fai_8PD(:,jj)+AddR_8PD(:,cnt);
            end
            w_ATC_RLS_BPD(:,cnt) = AddR(:,cnt);
            w_ATC_RLS_QPD(:,cnt) = AddR_QPD(:,cnt);
            w_ATC_RLS_6PD(:,cnt) = AddR_6PD(:,cnt);
            w_ATC_RLS_8PD(:,cnt) = AddR_8PD(:,cnt);
        end
        %% CTA
        AddR_cta=zeros(p,num);
        AddR_QPD_cta=zeros(p,num);
        AddR_6PD_cta=zeros(p,num);
        AddR_8PD_cta=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            for jj=1:num
                AddR_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS_BPD(:,jj)+AddR_cta(:,cnt);
                AddR_QPD_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS_QPD(:,jj)+AddR_QPD_cta(:,cnt);
                AddR_6PD_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS_6PD(:,jj)+AddR_6PD_cta(:,cnt);
                AddR_8PD_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS_8PD(:,jj)+AddR_8PD_cta(:,cnt);
            end
            FaiR_cta(:,cnt) = AddR_cta(:,cnt);
            Fai_cta_QPD(:,cnt) = AddR_QPD_cta(:,cnt);
            Fai_cta_6PD(:,cnt) = AddR_6PD_cta(:,cnt);
            Fai_cta_8PD(:,cnt) = AddR_8PD_cta(:,cnt);
        end
        for cnt=1:num
            Err_CTA_RLS_BPD(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_BPD(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_BPD(:,cnt));
            Err_CTA_RLS_QPD(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_QPD(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_QPD(:,cnt));
            Err_CTA_RLS_6PD(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_6PD(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_6PD(:,cnt));
            Err_CTA_RLS_8PD(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_8PD(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_8PD(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);%
            %% CTA-RLS-BPD
            ek_Rand1_cta(cnt)=dn(cnt)-FaiR_cta(:,cnt)'*un1(:,cnt);
            for k=1:orders
                P_Rand_cta(ii,k+orders*(cnt-1) )=exp(-1*(ek_Rand1_cta(cnt)-MIU_Rand(cnt,k))^2/(2*SIGMA_Rand(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand(cnt,k));
            end
            for k=1:orders
                V_Rand_cta(ii,k+orders*(cnt-1) )=WEIGHT_Rand(cnt,k)*P_Rand_cta(ii,k+orders*(cnt-1))/(WEIGHT_Rand(cnt,:)*P_Rand_cta(ii,orders*(cnt-1)+1:orders*(cnt-1)+orders)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders
                R1=R1+V_Rand_cta(ii,k+orders*(cnt-1))/SIGMA_Rand(cnt,k)^2;
                R2=R2+V_Rand_cta(ii,k+orders*(cnt-1))*MIU_Rand(cnt,k)/SIGMA_Rand(cnt,k)^2;
            end
            kn_cta(:,:,cnt) = Pn_RLS_BPD_cta(:,:,cnt) * un1(:,cnt) / ( lamda_New_CTA+ R1*un1(:,cnt)' * Pn_RLS_BPD_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_BPD_cta(:,:,cnt) = 1/lamda_New_CTA * ( Pn_RLS_BPD_cta(:,:,cnt) - R1*kn_cta(:,:,cnt) * un1(:,cnt)' * Pn_RLS_BPD_cta(:,:,cnt));
            w_CTA_RLS_BPD(:,cnt) = FaiR_cta(:,cnt) +kn_cta(:,:,cnt) * (R1*ek_Rand1_cta(cnt)-R2);
            %% CTA-RLS-QPD
            ek_Rand2_cta(cnt)=dn(cnt)-Fai_cta_QPD(:,cnt)'*un1(:,cnt);
            for k=1:orders1
                P_Rand_cta2(ii,k+orders1*(cnt-1) )=exp(-1*(ek_Rand2_cta(cnt)-MIU_Rand2(cnt,k))^2/(2*SIGMA_Rand2(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand2(cnt,k));
            end
            for k=1:orders1
                V_Rand_cta2(ii,k+orders1*(cnt-1) )=WEIGHT_Rand2(cnt,k)*P_Rand_cta2(ii,k+orders1*(cnt-1))/(WEIGHT_Rand2(cnt,:)*P_Rand_cta2(ii,orders1*(cnt-1)+1:orders1*(cnt-1)+orders1)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders1
                R1=R1+V_Rand_cta2(ii,k+orders1*(cnt-1))/SIGMA_Rand2(cnt,k)^2;
                R2=R2+V_Rand_cta2(ii,k+orders1*(cnt-1))*MIU_Rand2(cnt,k)/SIGMA_Rand2(cnt,k)^2;
            end
            kn_cta2(:,:,cnt) = Pn_RLS_QPD_cta(:,:,cnt) * un1(:,cnt) / (lamda_New_CTA_QPD+ R1*un1(:,cnt)' * Pn_RLS_QPD_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_QPD_cta(:,:,cnt) = 1/lamda_New_CTA_QPD * ( Pn_RLS_QPD_cta(:,:,cnt) - R1*kn_cta2(:,:,cnt) * un1(:,cnt)' * Pn_RLS_QPD_cta(:,:,cnt));
            w_CTA_RLS_QPD(:,cnt) =    Fai_cta_QPD(:,cnt) +kn_cta2(:,:,cnt) * (R1*ek_Rand2_cta(cnt)-R2);
            %% CTA-RLS-6PD
            ek_Rand3_cta(cnt)=dn(cnt)-Fai_cta_6PD(:,cnt)'*un1(:,cnt);
            for k=1:orders2
                P_Rand_cta3(ii,k+orders2*(cnt-1) )=exp(-1*(ek_Rand3_cta(cnt)-MIU_Rand3(cnt,k))^2/(2*SIGMA_Rand3(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand3(cnt,k));
            end
            for k=1:orders2
                V_Rand_cta3(ii,k+orders2*(cnt-1) )=WEIGHT_Rand3(cnt,k)*P_Rand_cta3(ii,k+orders2*(cnt-1))/(WEIGHT_Rand3(cnt,:)*P_Rand_cta3(ii,orders2*(cnt-1)+1:orders2*(cnt-1)+orders2)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders2
                R1=R1+V_Rand_cta3(ii,k+orders2*(cnt-1))/SIGMA_Rand3(cnt,k)^2;
                R2=R2+V_Rand_cta3(ii,k+orders2*(cnt-1))*MIU_Rand3(cnt,k)/SIGMA_Rand3(cnt,k)^2;
            end
            kn_cta3(:,:,cnt) = Pn_RLS_6PD_cta(:,:,cnt) * un1(:,cnt) / (lamda_New_CTA_6PD+ R1*un1(:,cnt)' * Pn_RLS_6PD_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_6PD_cta(:,:,cnt) = 1/lamda_New_CTA_6PD * ( Pn_RLS_6PD_cta(:,:,cnt) - R1*kn_cta3(:,:,cnt) * un1(:,cnt)' * Pn_RLS_6PD_cta(:,:,cnt));
            w_CTA_RLS_6PD(:,cnt) =    Fai_cta_6PD(:,cnt) +kn_cta3(:,:,cnt) * (R1*ek_Rand3_cta(cnt)-R2);
            %% CTA-RLS-8PD
            ek_Rand4_cta(cnt)=dn(cnt)-Fai_cta_8PD(:,cnt)'*un1(:,cnt);
            for k=1:orders3
                P_Rand_cta4(ii,k+orders3*(cnt-1) )=exp(-1*(ek_Rand4_cta(cnt)-MIU_Rand4(cnt,k))^2/(2*SIGMA_Rand4(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand4(cnt,k));
            end
            for k=1:orders3
                V_Rand_cta4(ii,k+orders3*(cnt-1) )=WEIGHT_Rand4(cnt,k)*P_Rand_cta4(ii,k+orders3*(cnt-1))/(WEIGHT_Rand4(cnt,:)*P_Rand_cta4(ii,orders3*(cnt-1)+1:orders3*(cnt-1)+orders3)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders3
                R1=R1+V_Rand_cta4(ii,k+orders3*(cnt-1))/SIGMA_Rand4(cnt,k)^2;
                R2=R2+V_Rand_cta4(ii,k+orders3*(cnt-1))*MIU_Rand4(cnt,k)/SIGMA_Rand4(cnt,k)^2;
            end
            kn_cta4(:,:,cnt) = Pn_RLS_8PD_cta(:,:,cnt) * un1(:,cnt) / (lamda_New_CTA_8PD+ R1*un1(:,cnt)' * Pn_RLS_8PD_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_8PD_cta(:,:,cnt) = 1/lamda_New_CTA_8PD * ( Pn_RLS_8PD_cta(:,:,cnt) - R1*kn_cta4(:,:,cnt) * un1(:,cnt)' * Pn_RLS_8PD_cta(:,:,cnt));
            w_CTA_RLS_8PD(:,cnt) =    Fai_cta_8PD(:,cnt) +kn_cta4(:,:,cnt) * (R1*ek_Rand4_cta(cnt)-R2);
        end
    end

    disp(mm);
end

Err_ATC_RLS_BPD1=mean(Err_ATC_RLS_BPD,3);
Err_ATC_RLS_QPD1=mean(Err_ATC_RLS_QPD,3);
Err_ATC_RLS_6PD1=mean(Err_ATC_RLS_6PD,3);
Err_ATC_RLS_8PD1=mean(Err_ATC_RLS_8PD,3);

Err_CTA_RLS_BPD1=mean(Err_CTA_RLS_BPD,3);
Err_CTA_RLS_QPD1=mean(Err_CTA_RLS_QPD,3);
Err_CTA_RLS_6PD1=mean(Err_CTA_RLS_6PD,3);
Err_CTA_RLS_8PD1=mean(Err_CTA_RLS_8PD,3);

for cnt=1:num
    Err_TH_RLS_New(cnt) = mean(Err_ATC_RLS_BPD(:,L,cnt));
    Err_TH_CTA_RLS_BPD(cnt) = mean(Err_CTA_RLS_BPD(:,L,cnt));
    Err_TH_RLS_New_QPD(cnt) = mean(Err_ATC_RLS_QPD(:,L,cnt));
    Err_TH_CTA_RLS_QPD(cnt) = mean(Err_CTA_RLS_QPD(:,L,cnt));
    Err_TH_RLS_New_6PD(cnt) = mean(Err_ATC_RLS_6PD(:,L,cnt));
    Err_TH_CTA_RLS_6PD(cnt) = mean(Err_CTA_RLS_6PD(:,L,cnt));
    Err_TH_RLS_New_8PD(cnt) = mean(Err_ATC_RLS_8PD(:,L,cnt));
    Err_TH_CTA_RLS_8PD(cnt) = mean(Err_CTA_RLS_8PD(:,L,cnt));
end


toc
figure(1),hold on;
plot(10*log10(mean(Err_CTA_RLS_BPD1)),'g');
plot(10*log10(mean(Err_ATC_RLS_BPD1)),'c');
plot(10*log10(mean(Err_CTA_RLS_QPD1)),'y');
plot(10*log10(mean(Err_ATC_RLS_QPD1)),'k');
plot(10*log10(mean(Err_CTA_RLS_6PD1)),'b');
plot(10*log10(mean(Err_ATC_RLS_6PD1)),'r');
plot(10*log10(mean(Err_CTA_RLS_8PD1)),'m');
plot(10*log10(mean(Err_ATC_RLS_8PD1)));
legend('CTA-RMLF(K=2)','ATC-RMLF(K=2)','CTA-RMLF(K=4)','ATC-RMLF(K=4)','CTA-RMLF(K=6)','ATC-RMLF(K=6)','CTA-RMLF(K=8)','ATC-RMLF(K=8)');
xlabel('Iterations');ylabel('MSD');
ylim([-20,10]);
box on;grid on;
hold off

figure(2),hold on
plot(10*log10((Err_TH_CTA_RLS_BPD)),'-gp','LineWidth',1.5);
plot(10*log10((Err_TH_RLS_New)),'-co','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RLS_QPD)),'-yh','LineWidth',1.5);
plot(10*log10((Err_TH_RLS_New_QPD)),'-kv','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RLS_6PD)),'-b>','LineWidth',1.5);
plot(10*log10((Err_TH_RLS_New_6PD)),'-r<','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RLS_8PD)),'-m*','LineWidth',1.5);
plot(10*log10((Err_TH_RLS_New_8PD)),'-d','LineWidth',1.5);
legend('CTA-RMLF(K=2)','ATC-RMLF(K=2)','CTA-RMLF(K=4)','ATC-RMLF(K=4)','CTA-RMLF(K=6)','ATC-RMLF(K=6)','CTA-RMLF(K=8)','ATC-RMLF(K=8)');
xlabel('Node number,k');ylabel('Steady State MSD (dB)');
ylim([-20,-6]);
xlim([1,20]);
box on;grid on;
hold off

