%
%              C = VRFT1_dy_du(u,y,Md,Mu,B,Wd,Wu,lambda,k,preFilt)
%     [C, theta] = VRFT1_dy_du(u,y,Md,Mu,B,Wd,Wu,lambda,k,preFilt)
% 
% Design a 1 d.o.f. linear controller so as to match the d(t) to y(t) output  
% sensitivity with the model reference Md and the d(t) to u(t) input sensitivity
% with the model reference Mu (see figure below).
%
% Reference: 
% [1] M.C. Campi, A. Lecchini and S.M. Savaresi. 
%     An application of the virtual reference feedback tuning method to a benchmark problem. 
%     European Journal of Control, Special Issue on  "Design and Optimisation of Restricted Complexity Controllers", 1:66-76, 2003.
%
% 
% Inputs (compulsory)--------------------------------------------------------------
% u:  column vector (Nx1) that contains the INPUT data collected from the plant.
% y:  matrix (Nx2) that contain the OUTPUT data collected from the plant:
%     y(:,1) contains the output of the plant when it is fed by u
%     y(:,2) contains the output of the plant when it is fed by y(:,1).
%     If y is a Nx4 matrix, the columns 1 and 2 contain the OUTPUT data collected in
%     the first experiment and the columns 3 and 4 contain the OUTPUT data collected
%     in the second experiment (the noise realizations in the two experiments must be
%     uncorrelated). 
% Md: tf-object that represents the discrete transfer function of the reference model.
%     The reference model Md(z) describes the desired closed loop behaviour from the 
%     signal d(t) to the output y(t) (output sensitivity).
% Mu: tf-object that represents the discrete trasfer function Mu(z) that describes 
%     the desired dynamical relationship from d(t) to u(t) (input sensitivity).
% B:  column vector of tf-objects. The linear controller has the following structure: 
%     C(z,theta)= B'*theta, where B is a column vector of transfer functions, and 
%     theta is the vector of parameters.
% 
% Inputs (optional)----------------------------------------------------------------
% Wd: tf-object of the weighting function Wd(z). If this parameter is empty [],
%     the function automatically sets Wd(z) = 1. 
% Wu: tf-object of the weighting function Wu(z). If this parameter is empty [], 
%     the function automatically sets Wu(z) = 1.   
% lambda: this parameter is used to balance the emphasis on Md
%     (lambda close to 0) or on Mu (lambda close to 1).  Note that 0<=lambda<=1. If this 
%     parameter is empty [], the function automatically sets lambda = 0.5.
% k:  this parameter must be used only if the measured y is noisy, and a single experiment
%     is avaible; in this case this parameter sets the order of an ARX(k,k) model used to 
%     make an approximate model of the plant. Otherwise this parameter must be empty: [].
% preFilt: if this parameter is set to 'n', the optimal filter is disabled, and the filter 
%     L(z) is set to 1. If this parameter is empty ([]) the function uses the optimal VRFT filter. 
%     
% Outputs -------------------------------------------------------------------------
% C:  tf-object, which represents the transfer function of the designed controller.
% theta (optional): the vector such that C = theta'*B     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                   ____________             ____________         d(t)
%                  |            |           |            |         |
% r(t)       e(t)  | C(z,theta) |   u(t)    |    P(z)    |         |  y(t)
% -------->O------>|            |---------->|            |-------->O------->           
%          |-      |____________|           |____________|             |
%          |                                                           |
%          |___________________________________________________________|
% 
%          
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


function [C, theta_vector] = VRFT1_dy_du(u,y,Md,Mu,B,Wd,Wu,lambda,k,preFilt)

% Passo 1: controllo della correttezza dei parametri inseriti

errore = 0; %inizializzazioni variabli

if( nargin ~= 10 )   %controllo numero argomenti passati alla funzione
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
elseif( nCy > 4)
    errore = 18;
    dispMsg(errore)
    return
elseif( nRy > nCy & nCy == 2)   %ricavo y1 e y2 #no noise / noise IV sim. experim.
    y1 = y(:,1);
    y2 = y(:,2);
    y3 = [];
    y4 = [];
elseif( nRy > nCy & nCy == 4)   %ricavo y1 - y2 - y3 - y4 #noise IV rep. experim.
    y1 = y(:,1);
    y2 = y(:,2);
    y3 = y(:,3);
    y4 = y(:,4);    
elseif( nCy == 1 | nCy == 3)
    errore = 19;
    dispMsg(errore)
    return
end

if( isempty(Md) | isempty(Mu) )   %controllo Md e Mu
    errore = 7;
    dispMsg(errore)
    return
elseif(Md.Ts ~= Mu.Ts)
    errore = 21;
    dispMsg(errore)
    return
elseif( isstatic(Md) && Md.Numerator{1}==1 )
    errore = 24;
    dispMsg(errore)
    return
elseif( Md.variable ~= 'z^-1' | Mu.variable ~= 'z^-1')  %Md e Mu puo' essere inserito in qualunque forma (z,z^-1,q)
    Md.variable = 'z^-1';
    Mu.variable = 'z^-1';
end

if( isempty(B) )    %controllo B
    errore = 8;
    dispMsg(errore)
    return
elseif(nRB >= nCB)
    for h = 1:nRB   %per ogni elemento di B
        fdt = B(h,1);
        if(fdt.variable ~= 'z^-1')  %sistemo la variabile
            fdt.variable = 'z^-1';
            B(h,1) = fdt;
        end
        if(fdt.Ts ~= Md.Ts) %controllo il tempo di campionamento
            errore = 9;
            dispMsg(errore)
        return
        end
    end
end

if( isempty(Wd) )    %controllo Wd
    Tsampling = Md.Ts;
    Wd = tf([1],[1],Tsampling,'variable','z^-1');
elseif( ~isempty(Wd) )
    if( Wd.Ts ~= Md.Ts)
        errore = 10;
        dispMsg(errore)
        return
    end
    Wd.variable = 'z^-1';
end

if( isempty(Wu) )    %controllo Wu
    Tsampling = Md.Ts;
    Wu = tf([1],[1],Tsampling,'variable','z^-1');
elseif( ~isempty(Wu) )
    if( Wu.Ts ~= Md.Ts)
        errore = 10;
        dispMsg(errore)
        return
    end
    Wu.variable = 'z^-1';
end
  
if( isempty(lambda) )   %controllo lambda
    lambda = 0.5;
elseif(lambda < 0 | lambda > 1)
    errore = 21;
    dispMsg(errore)
    return
end

if( ~isempty(k) )   %controllo k
    if(nRk > nCk)   %sistemo eventuali vettori colonna
        k = k';
    end
    if(nCk > 2)
        errore = 11;
        dispMsg(errore)
        return
    elseif(nCk == 1)
        k = [k k];  %la funzione VRFT_engine si aspetta un vettore di due elementi
    end
end

if( isempty(preFilt) )  %controllo filtro
    Tsampling = Md.Ts;
    U = stima_U(u,Tsampling);
    L = minreal( (1 - Md)*Md*Wd*(U\1) );  %calcolo il filtro ottimo
    Lu = [];    %la funzione InSens calcola il filtro ottimo
else
    try
        if( preFilt == 'n')
        L = tf([1],[1],Md.Ts,'variable','z^-1');
        Lu = L; 
        end
    catch   %FUNZIONALITA' NON DICHIARATA: permette di filtrare i dati con i filtri proposti --> preFilt = [L ; Lu]
        L = preFilt(1,1);
        Lu = preFilt(2,1);
        if(L.Ts ~= Md.Ts | Lu.Ts ~= Md.Ts)
            errore = 12;
            dispMsg(errore)
            return
        else
            L.variable = 'z^-1';
            Lu.variable = 'z^-1';
        end
    end
end    


% Passo 2: chiamo il VRFT_engine e InSens se non ci sono errori

[Cr, Cy, An_d, Fn_d] = vrft_engine(u,y1,y3,[],Md,[],B,[],L,k,'n');  %calcolo Cy

Mu = Mu * (-1);    % ### da valutare ###
[Cu, An_u, Fn_u] = InSens(u,y1,y3,y2,y4,Mu,B,Wu,Lu,k);  %calcolo Cu

%combinazione convessa dei parametri  
An = ((1 - lambda)*An_d) + (lambda*An_u);
Fn = ((1 - lambda)*Fn_d) + (lambda*Fn_u);
            
theta_vector = An\Fn; 
%C = minreal(theta_vector' * B);
if (length(theta_vector)==0 | length(theta_vector)<=2) %note that when length(theta_vector)==0 (theta_vector' * B) is the empty transfer function
C = minreal(theta_vector' * B);
else
    C=theta_vector(1)*B(1,:);
  for j=2:length(theta_vector)
    C=minreal(C+theta_vector(j)*B(j,:));
  end
end


