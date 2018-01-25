%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes an image and obtains the mean prediction error (detail coefficient energy) 
% when the proposed Update/Prediction solution is used in a multiresolution
% framework (4 levels of the transform)
%     
% References: "Optimized Update/Prediction Assignment for
% Lifting Transforms on Graphs", Eduardo Martinez-Enriquez, Jesus Cid-Sueiro, 
% Fernando Diaz-de-Maria, and Antonio Ortega
%
% Author:
%  - Eduardo Martinez-Enriquez <emenriquez@tsc.uc3m.es>
% 
%     Copyright (C)  2017 Eduardo Martínez-Enríquez.
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image=imread('lena_md.png','png');
% image=imread('cameraman_sm.png','png');
% image=imread('barbara.png','png');

  if length(size(image))==3
        image= rgb2gray(image);
  end

% VARIABLES INITIALIZATION  
reconst_5_perc=zeros(size(image));
reconst_10_perc=zeros(size(image));
reconst_20_perc=zeros(size(image));
reconst_30_perc=zeros(size(image));
reconst_40_perc=zeros(size(image));
reconst_50_perc=zeros(size(image));
reconst_5_percWMC=zeros(size(image));
reconst_10_percWMC=zeros(size(image));
reconst_20_percWMC=zeros(size(image));
reconst_30_percWMC=zeros(size(image));
reconst_40_percWMC=zeros(size(image));
reconst_50_percWMC=zeros(size(image));
reconst_5_percRAN=zeros(size(image));
reconst_10_percRAN=zeros(size(image));
reconst_20_percRAN=zeros(size(image));
reconst_30_percRAN=zeros(size(image));
reconst_40_percRAN=zeros(size(image));
reconst_50_percRAN=zeros(size(image));

figure
imshow(image,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu_ep=0;
nu_et=0;
var_et=10;
var_ep=var(double(image(:)));
c=mean(double(image(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORK IN BLOCKS TO REDUCE THE COMPUTATIONAL COMPLEXITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % size_block=25; % 9*9 blocks of 25x25 in Lena (225x225); 5 blocks 45
% Configuration used in the paper:
size_block=32; % 4*4 blocks of 32x32 in Lena and cameraman (128x128); 16*16 blocks of 32x32 in Barbara (512x512);
% size_block=16;

% MAXIMUM PERCENTAGE OF |U| NEIGHBORS ALLOWED FOR EVERY LEVEL.
restriction=0.5;
restriction_G2=0.5;
restriction_G3=0.5;
restriction_G4=0.5;

N=size_block^2;
num_blocks=length(image)/size_block;
num_it=num_blocks^2;
it=1;

prop_ev_acum_MC_matrix=zeros(num_it,N);
prop_ev_acum_MAM_matrix=zeros(num_it,N);
prop_ev_acum_RAN_matrix=zeros(num_it,N);
D0_acum_MC_matrix=zeros(num_it,N);
D0_acum_MAM_matrix=zeros(num_it,N);
D0_acum_RAN_matrix=zeros(num_it,N);
Energia_final_bloques_J1=zeros(num_it,1);
Energia_final_bloques_J2=zeros(num_it,1);
Energia_final_bloques_J3=zeros(num_it,1);
Energia_final_bloques_J4=zeros(num_it,1);


% % % % % LOOP OVER BLOCKS

for ver_block=1:num_blocks
    for hor_block=1:num_blocks


    partial_image=image(size_block*(ver_block-1)+1:size_block*(ver_block), size_block*(hor_block-1)+1:size_block*(hor_block));
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE THE WEIGHTED GRAPH 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[W,x]=generate_signal_W_imagen3(size_block, partial_image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE Q MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Q] = Generate_Q_matrix(W);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GREEDY ALGORITHM FOR THE PROPOSED (MAM) UPA SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
label=ones(N,1);

tic

[energia_acum D label_MAM U_set P_set D0_acum  prop_ev_acum] = greedy_MAM_image_low(Q, var_ep,var_et,nu_et,label,W,c,x,size_block,ver_block,hor_block,[],[],[],[],[],[], restriction,1);
energia_final=energia_acum/D
Energia_final_bloques_J1(it)=energia_final;

a=toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level of the transform j=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A_G2 A_odds A_evens]=construct_Gj_mat(W, label_MAM);

a=sum(A_G2,1);
nodes_save=find(a~=0);
matrix_aux=A_G2(nodes_save,nodes_save);
x_G2=x(nodes_save);
label_MAM_G2=ones(length(nodes_save),1);

[Q_G2] = Generate_Q_matrix(matrix_aux);
[energia_acum_G2 D_G2 label_MAM_G2 U_set_G2 P_set_G2 D0_acum_G2  prop_ev_acum_G2] = greedy_MAM_image_low(Q_G2, var_ep,var_et,nu_et,label_MAM_G2,matrix_aux,c,x_G2,size_block,ver_block,hor_block,[],[],[],[],[],[], restriction_G2,0);

P_set_G2=nodes_save(P_set_G2);
U_set_G2=nodes_save(U_set_G2);

P_final=union(P_set,P_set_G2);
energia_final=(energia_acum+energia_acum_G2)/(D+D_G2)
Energia_final_bloques_J2(it)=energia_final;

label_MAM_G2=ones(size(A_G2,1),1);
label_MAM_G2(U_set_G2)=-1;
label_MAM_G2(P_final)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level of the transform j=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A_G3 A_odds_G2 A_evens_G2]=construct_Gj_mat(A_G2, label_MAM_G2);
a=sum(A_G3,1);
nodes_save=find(a~=0);
matrix_aux=A_G3(nodes_save,nodes_save);
x_G3=x(nodes_save);
label_MAM_G3=ones(length(nodes_save),1);

[Q_G3] = Generate_Q_matrix(matrix_aux);
[energia_acum_G3 D_G3 label_MAM_G3 U_set_G3 P_set_G3 D0_acum_G3  prop_ev_acum_G3] = greedy_MAM_image_low(Q_G3, var_ep,var_et,nu_et,label_MAM_G3,matrix_aux,c,x_G3,size_block,ver_block,hor_block,[],[],[],[],[],[], restriction_G3,0);

P_set_G3=nodes_save(P_set_G3);
U_set_G3=nodes_save(U_set_G3);

P_final=union(P_final,P_set_G3);
energia_final=(energia_acum+energia_acum_G2+energia_acum_G3)/(D+D_G2+D_G3)
Energia_final_bloques_J3(it)=energia_final;

label_MAM_G3=ones(size(A_G3,1),1);
label_MAM_G3(U_set_G3)=-1;
label_MAM_G3(P_final)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level of the transform j=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A_G4 A_odds_G3 A_evens_G3]=construct_Gj_mat(A_G3, label_MAM_G3);
a=sum(A_G4,1);
nodes_save=find(a~=0);
matrix_aux=A_G4(nodes_save,nodes_save);
x_G4=x(nodes_save);
label_MAM_G4=ones(length(nodes_save),1);

[Q_G4] = Generate_Q_matrix(matrix_aux);
[energia_acum_G4 D_G4 label_MAM_G4 U_set_G4 P_set_G4 D0_acum_G4  prop_ev_acum_G4] = greedy_MAM_image_low(Q_G4, var_ep,var_et,nu_et,label_MAM_G4,matrix_aux,c,x_G4,size_block,ver_block,hor_block,[],[],[],[],[],[], restriction_G4,0);

P_set_G4=nodes_save(P_set_G4);
U_set_G4=nodes_save(U_set_G4);

P_final=union(P_final,P_set_G4);
energia_final=(energia_acum+energia_acum_G2+energia_acum_G3+energia_acum_G4)/(D+D_G2+D_G3+D_G4)


label_MAM_G4=ones(size(A_G4,1),1);
label_MAM_G4(U_set_G4)=-1;
label_MAM_G4(P_final)=1;


a=min(prop_ev_acum);



figure, plot(prop_ev_acum,sqrt(D0_acum));
hold on
plot(prop_ev_acum_G2,sqrt(D0_acum_G2),'r');
plot(prop_ev_acum_G3,sqrt(D0_acum_G3),'g');
plot(prop_ev_acum_G4,sqrt(D0_acum_G4),'k');


% D0_acum_MAM_matrix(it,1:length(prop_ev_acum))=D0_acum_G2;
Energia_final_bloques_J4(it)=energia_final;


it=it+1;
close all
    end
end

sqrt(mean(Energia_final_bloques_J1))
sqrt(mean(Energia_final_bloques_J2))
sqrt(mean(Energia_final_bloques_J3))
sqrt(mean(Energia_final_bloques_J4))






