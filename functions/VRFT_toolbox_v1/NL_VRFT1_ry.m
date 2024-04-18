%
%     theta = NL_VRFT1_ry(u,y,Mr,B,k,preFilt)
% 
% Design a 1 degree of freedom (1 d.o.f.) nonlinear controller so that the
% closed-loop output (see figure below) in response to the reference signal 
% r(t) = Mr(z)^-1 y(t)
% is as close as possible to y(t). Refer to paper [1] for the iterative use 
% of this procedure so as to design a closed-loop that resembles Mr(z) when
% the reference is a desired signal r°(t).
% Reference: 
% [1] M.C. Campi and S.M. Savaresi (2006).
%     Direct nonlinear control design: the Virtual Refrence Feedback Tuning (VRFT) approach.
%     IEEE Trans. on Automatic Control, Vol.51, n.1, pp.14-27.
% 
% Inputs (compulsory)--------------------------------------------------------------
% u:  column vector (Nx1) that contains the INPUT data collected from the plant.
% y:  column vector (Nx1) that contains the OUTPUT data collected from the plant.
%     If y is a Nx2 matrix, the two columns contain the output data collected
%     in 2 different experiments (both experiments are made with the same input u(t); 
%     the two noise realizations must be uncorrelated).
% Mr: tf-object that represents the discrete transfer function of the reference 
%     model behaviour. The reference model Mr(z) describes the desired closed-loop behaviour 
%     from the reference r(t) to the output y(t).
% B:  cell array (K x 1) of strings. Each B{i}, i=1,...,K, is a string 
%     that describes a controller block.
%
%     VRFT constructs a controller whose output u(t) is the linear combination with the coefficients
%     given by theta of the outputs of the blocks.
%     Each controller block is a nonlinear system with input sequence "e" 
%     and output sequence "b", in the following relation:
%     b(t)= function(e(t),...,e(t-Te),b(t-1),...,b(t-Tu)).
%     B{i} contains the symbolic description of the right-hand side of 
%     the above system equation. 
%     
%     EXAMPLE:
%     By setting
%     B{1}='e(t)'
%     the first block is 
%     b(t)= e(t)
%    
%     By setting 
%     B{2}='0.1*b(t-1)+0.9*e(t-1)^1/3'
%
%     the second block is 
%     b(t)= 0.1*b(t-1)+0.9*e(t-1)^1/3
%
%     In this case, 
%     B={'e(t)','0.1*b(t-1)+0.9*e(t-1)^1/3'}
%
%     Please note that the initial conditions of each block are always set to zero.
%   
%
% Inputs (optional)---------------------------------------------------------------
% k:  this parameter must be used only if the measured output y is noisy, and a 
%     single experiment is available; in this case this parameter sets the order 
%     of an ARX(k,k) model used to make an approximate model of the plant. Otherwise 
%     this parameter must be empty: [].
% preFilt: if this parameter is set to 'n', the optimal filter is disabled, and the 
%     filter F(z) is set to 1. If this parameter is empty ([]) the function uses 
%     the filter 1-Mr(z)P(z), where P(z) is a rough linear estimate of the plant.
%          
% Outputs -------------------------------------------------------------------------
% theta:  weighting coefficients to combine the building blocks of the controller into
% the designed controller.     
%      
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                   ____________             ____________         d(t)
%                  |            |           |            |         |
% r(t)       e(t)  |      C     |   u(t)    |     P      |         |  y(t)
% -------->O------>|            |---------->|            |-------->O------->           
%          |-      |____________|           |____________|             |
%          |                                                           |
%          |___________________________________________________________|
%          
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


function [theta,ref] = NL_VRFT1_ry(u,y,Mr,B,k,preFilt)

% Passo 1: controllo della correttezza dei parametri inseriti

errore = 0; %inizializzazioni variabli

if( nargin ~= 6 )   %controllo numero argomenti passati alla funzione
    errore = 1;
    dispMsg(errore)
    return
end

[nRu nCu] = size(u);
[nRy nCy] = size(y);
[nRB nCB] = size(B);
[nRk nCk] = size(k);

if( isempty(u) )    %controllo u
    errore = 2;
    dispMsg(errore)
    return
elseif( nRu == 1 & nCu >  nRu ) %correggo se u viene inserito come vettore riga
    u = u';
    [nRu nCu] = size(u);
elseif( nCu ~= 1)
    errore = 3;
    dispMsg(errore)
    return
end

if( isempty(y) )    %controllo y
    errore = 4;
    dispMsg(errore)
    return
elseif( nRy ~= nRu )    %controllo consistenza con u
    errore = 5;
    dispMsg(errore)
    return
elseif( nCy > 2)
    errore = 6;
    dispMsg(errore)
    return
elseif( nRy > nCy & nCy == 2)   %ricavo y1 e y2
    y1 = y(:,1);
    y2 = y(:,2);
elseif(nRy > nCy & nCy == 1)    %no noise / noise IV sim. experim.
    y1 = y;
    y2 = [];
end

if( isempty(Mr) )   %controllo Mr
    errore = 7;
    dispMsg(errore)
    return
elseif( isstatic(Mr) && Mr.Numerator{1}==1 )
    errore = 24;
    dispMsg(errore)
    return
elseif( Mr.variable ~= 'z^-1')  %Mr puo' essere inserito in qualunque forma (z,z^-1,q)
    Mr.variable = 'z^-1';
end

if( isempty(B) )    %controllo B
    errore = 8;
    dispMsg(errore)
    return
elseif(nRB < nCB)
    B=B';
    %Further check-ups to be done!!!
end

%MODIFICATO PER IL RITARDO 15/10/2004
if( ~isempty(k) )   %controllo k
    if(nRk > nCk)   %sistemo eventuali vettori colonna
        k = k';
    end
    if(nCk > 3)
        errore = 11;
        dispMsg(errore)
        return
    elseif(nCk == 1)
        k = [k k 0];  %la funzione VRFT_engine si aspetta un vettore di due elementi
    end
end

if( isempty(preFilt) )
    Tsampling = Mr.Ts;
    
    P_approx=identifica(u(:,1),y(:,1),Tsampling,[4,4,0]);    
    L = minreal( (1 - Mr)*P_approx);  %calcolo il filtro (diretto) di default
else
    try
        if( preFilt == 'n')
        L = tf([1],[1],Mr.Ts,'variable','z^-1');
        end
    catch   %FUNZIONALITA' NON DICHIARATA: permette di filtrare i dati con il filtro proposto
        L = preFilt;
        if(L.Ts ~= Mr.Ts)
            errore = 12;
            dispMsg(errore)
            return
        else
            L.variable = 'z^-1';
        end
    end
end


% Passo 2: chiamo il VRFT_engine se non ci sono errori
M = minreal(Mr/(1 -Mr));

[ThetaR, VirtualRef, An, Fn] = NL_vrft_engine(u,y1,y2,M,[],B,[],L,[],k,'n');
theta = ThetaR;
ref=VirtualRef;

