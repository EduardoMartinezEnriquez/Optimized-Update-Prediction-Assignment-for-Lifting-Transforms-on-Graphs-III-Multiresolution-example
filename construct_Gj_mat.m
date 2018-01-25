function [A_Gj A_odds A_evens]=construct_Gj_mat(A, label)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM FOR THE CONSTRUCTION OF THE ADJACENCY MATRIX AT THE (j)th LEVEL OF 
% THE TRANSFORM FROM THE ADJACENCY MATRIX AT THE (j-1)th LEVEL
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


NoOfSensors=length(label);
% load ej_G2

% Find the Evens of G2

% Initialize the Adj Matrix in G3

A_Gj=speye(NoOfSensors);
A_Gj=A_Gj*0;

Evens=find(label==-1);
Odds=find(label==1);

A_evens=A;
A_evens(Odds,:)=0;
A_evens(:,Evens)=0;

% A_odds=A;
% A_odds(Evens,:)=0;
% A_odds(:,Odds)=0;

A_odds=A_evens';

A_Gj_nonorm=A_evens*A_odds;


% Evens neighbors of Even nodes

A_evens_evens=A;
A_evens_evens(Odds,:)=0;
A_evens_evens(:,Odds)=0;

A_Gj=A_evens_evens+A_Gj_nonorm;


% Normalization
 NZ=A_evens~=0;
 NZ=double(NZ);
 
 NZevev=A_evens_evens~=0;
 NZevev=double(NZevev);
 
 C=NZ*NZ';
 
 % Include the nodes directly connected 
 C=C+NZevev;
% Avoid singularities
% D=C==0;
nz=find(C~=0);
% G=C+D;


A_Gj(nz)=A_Gj(nz)./C(nz);





%Remove autoloops (conections of nodes with themselves) 
d=diag(A_Gj);
A_Gj=A_Gj-diag(d);


% for c=1:NoOfSensors
%     A_G2(c,c)=0;
% end

    

