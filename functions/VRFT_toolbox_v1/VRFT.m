%
% VRFT Toolbox - Version 1.1 - 02-Aug-2019
% ---------------------------------------------------------------------------------
% 
% <a href="matlab:eval('help VRFT1_ry')":>VRFT1_ry</a>:    designs a 1 degree of freedom (1 d.o.f.) linear controller so as to
%             match the r(t) to y(t) closed-loop transfer function with the model 
%             reference Mr (see Fig.1). 
%             
% <a href="matlab:eval('help VRFT1_dy')":>VRFT1_dy</a>:    designs a 1 d.o.f. linear controller so as to match the d(t) to y(t) 
%             output sensitivity with the model reference Md (see Fig.1).
%             
% <a href="matlab:eval('help VRFT1_ry_ru')":>VRFT1_ry_ru</a>: designs a 1 d.o.f. linear controller so as to match the r(t) to y(t) 
%             closed-loop transfer function with the model reference Mr and the
%             r(t) to u(t) input sensitivity with the model reference Mu (see Fig.1).
%              
% <a href="matlab:eval('help VRFT1_dy_du')":>VRFT1_dy_du</a>: designs a 1 d.o.f. linear controller so as to match the d(t) to y(t) 
%             closed-loop transfer function with the model reference Md and the
%             d(t) to u(t) input sensitivity with the model reference Mu (see Fig.1).
%             
% <a href="matlab:eval('help VRFT2_ry_dy')":>VRFT2_ry_dy</a>: designs a 2 d.o.f. linear controller so as to match the r(t) to y(t)
%             transfer function and the d(t) to y(t) transfer function (see Fig.2).
%          
%
% <a href="matlab:eval('help NL_VRFT1_ry')":>NL_VRFT1_ry</a>: design a nonlinear controller to shape the closed-loop response 
% to the virtual reference.
%
%
% Type "help Name_Function" for more details.
%             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% - Figure 1: 1 DEGREE OF FREEDOM SCHEME
% 
%                   ____________             ____________         d(t)
%                  |            |           |            |         |
% r(t)       e(t)  | C(z,theta) |   u(t)    |    P(z)    |         |  y(t)
% -------->O------>|            |---------->|            |-------->O------->           
%          |-      |____________|           |____________|             |
%          |                                                           |
%          |___________________________________________________________|
%          
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% - Figure 2: 2 DEGREES OF FREEDOM SCHEME
% 
%          ______________                    ______________         d(t)
%         |              |                  |              |         |
% r(t)    | Cr(z,theta_r)|           u(t)   |     P(z)     |         |   y(t)
% ------->|              |------->O-------->|              |-------->O------->           
%         |______________|        |-        |______________|             |
%                                 |                                      |
%                                 |          ______________              |
%                                 |         |              |             |
%                                 |_________| Cy(z,theta_y)|_____________| 
%                                           |              |
%                                           |______________|
% 
% -----------------------------------------------------------------------------------
% 
% VRFT_GUI starts the Graphical User Interface for the VRFT toolbox.
% 
% -----------------------------------------------------------------------------------
%
% Toolbox reference:
%     - A. Carè, F. Torricelli, M.C. Campi, S. Savaresi. 
%       <a href="matlab:web('http://marco-campi.unibs.it/pdf-pszip/VRFT-toolbox.pdf')">A Toolbox for Virtual Reference Feedback Tuning (VRFT).</a>
%       In Proc. of the European Control Conference 2019, Napoli, Italy, 4252-4257, 2019.     
%
% Main references:
%     - Campi M.C., A. Lecchini, S.M. Savaresi (2002). 
%       Virtual Reference Feedback Tuning: a Direct Method for the Design of Feedback Controllers. 
%       Automatica, Vol.38, n.8, pp.1337-1346.
%
%     - S. Formentin, M.C. Campi, A. Carè, S.M. Savaresi (2019)
%       Deterministic continuous-time Virtual Reference Feedback Tuning (VRFT) with application to PID design.
%       Systems & Control Letters, vol. 127, pp. 25-34.
%     - M.C. Campi and S.M. Savaresi (2006).
%       Direct nonlinear control design: the Virtual Refrence Feedback Tuning (VRFT) approach.
%       IEEE Trans. on Automatic Control, Vol.51, n.1, pp.14-27.
%     - Campi M.C., Lecchini A., Savaresi S.M. (2003). 
%       An application of the Virtual Reference Feedback Tuning (VRFT) method to a benchmark active suspension system. 
%       European Journal of Control, Vol.9, pp.66-76. 
%     - Lecchini A., M.C. Campi, S.M. Savaresi (2002). 
%       Virtual reference feedback tuning for two degree of freedom controllers. 
%       International Journal on Adaptive Control and Signal Processing, vol.16, n.5, pp.355-371. 
%     - A. Lecchini, M.C. Campi and S.M Savaresi (2001). 
%       Sensitivity shaping via Virtual Reference Feedback Tuning. 
%       In Proc. 40th Conf. on Decision and Control, Orlando, pages 750-755.
%     - Guardabassi G.O., Savaresi S.M. (2000). 
%       Virtual Reference Direct Design method: an off-line approach to data-based control system design. 
%       IEEE Transactions on Automatic Control, Vol.45, n.5, pp.954-959. 
% -----------------------------------------------------------------------------------
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                               LICENSE AGREEMENT 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We (the licensee) understand that VRFT toolbox is supplied "as is", without expressed or implied warranty.
% We agree on the following:
% 
% - The licensers do not have any obligation to provide any maintenance or consulting help with respect 
%   to "VRFT toolbox".
% - The licensers neither have any responsibility for the correctness of systems designed using 
%   "VRFT toolbox", nor for the correctness of "VRFT toolbox" itself. 
% - We will never distribute or modify any part of the "VRFT toolbox" code without a written permission 
%   from Prof. Marco Campi (University of Brescia) or Prof. Sergio Savaresi (University of Milano). 
% - We will only use "VRFT toolbox" for non-profit research purposes. This implies that neither 
%   "VRFT toolbox" nor any part of its code should be used or modified for any commercial software product.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

