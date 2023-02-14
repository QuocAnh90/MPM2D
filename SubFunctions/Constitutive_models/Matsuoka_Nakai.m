function Particle=Matsuoka_Nakai(SolidModel,Particle,Time)

%% Update particle's velocity strain Rate
% Formula dESP = (L_sp + L_sp')/2*dt;
Particle.strainInc(:,1)    = Particle.Gradvelocity(:,1) * Time.timestep;                                   % particle strain increment xx
Particle.strainInc(:,2)    = Particle.Gradvelocity(:,4) * Time.timestep;                                   % particle strain increment yy
Particle.strainInc(:,3)    = (Particle.Gradvelocity(:,2)+Particle.Gradvelocity(:,2))/2 * Time.timestep;    % particle strain increment xy

%% Parameter of the model
E           = SolidModel.Young_modul;                           % Young's modoulus
nu          = SolidModel.nu;                                    % Poisson's ratio

% N           = SolidModel.N;                                     % state variable N
M           = SolidModel.M;                                     % state variable M
lambda_c    = SolidModel.lambda_c;                              % Dilatancy parameter
% Szz         = SolidModel.stressZZ;                              % Stress ZZ


Particle.stressZZ       = zeros(Particle.Count,1);                     % Stress out of plane
Particle.N              = zeros(Particle.Count,1);                     % Dilatancy state variables

% Elastic matrix 2x2
C = zeros(4,4);                                                                       
	C(1,1) = 1.0 - nu ; C(1,2) = nu		  ; C(1,3) = nu;
	C(2,1) = nu			; C(2,2) = 1.0 - nu ; C(2,3) = nu;
	C(3,1) = nu			; C(3,2) = nu		  ; C(3,3) = 1.0 - nu;
    C(4,4) = (1-2*nu)/2;
	C = E/((1+nu)*(1-2*nu)) * C; 

% Elastic matrix 3x3
D = E/(1+nu)/(1-2*nu)*[ 1-nu nu nu 0 0 0 ...
                      ; nu 1-nu nu 0 0 0 ...
                      ; nu nu 1-nu 0 0 0 ...
                      ; 0 0 0 (1-2*nu)/2 0 0 ...
                      ; 0 0 0 0 (1-2*nu)/2 0 ...
                      ; 0 0 0 0 0 (1-2*nu)/2];
                    
%% Loop all particles
for p = 1:Particle.Count
%     p
%     if p ==1
%         Particle.Gradvelocity(p,:)
%     end
% Plastic multiplier 
lambda = 0;

% State variables
N           = Particle.N(p,1);            % state variable N
Szz         = Particle.stressZZ(p,1);     % Stress ZZ


%% -----------------------------------------Trial step------------------------------------------%
de_sp    = zeros(4,1);      % strain increment vector 1x4
de_sp(1) = -Particle.strainInc(p,1);
de_sp(2) = -Particle.strainInc(p,2);
de_sp(3) = 0;
de_sp(4) = -Particle.strainInc(p,3);
Particle.stress(:,p) = Particle.stress(:,p)*-1;

deps = [de_sp;0;0];         % strain increment vector 1x4

dSigma = C*de_sp;           % stress increment
Sigma  = [Particle.stress(1,p) ; Particle.stress(2,p) ; Szz ; Particle.stress(3,p)] + dSigma;

%% --- Value of Yield function f = X - (M+N)
Sigma_Tensor3x3 = [Sigma(1) Sigma(4) 0 ; Sigma(4) Sigma(2) 0 ; 0 0 Sigma(3)];
  
I1 = trace(Sigma_Tensor3x3);
I2 = 0.5*(trace(Sigma_Tensor3x3)^2 - trace(Sigma_Tensor3x3^2));
I3 = det (Sigma_Tensor3x3);
X = sqrt (I1*I2/I3 - 9);
f = X - (M+N);                % Matsuoka Nakai Yield function

% Cut off condition
if I3 <= 0
   f = 0;
   Sigma(1:4) = 0;
end
        
if I1*I2/I3 <= 9
   f = 0;
   Sigma(1:4) = 0;
end
        
%----------------------------------------------------------

  if (f > 0)      
    %% elaso_plastic condition
      i = 0;
        while abs(f) > 0.00000000001
%             i = i + 1
            
%             if i > 100
%                 de_sp
%             end
            
            S = Sigma_Tensor3x3; %Just to simplify notation
            % Vector B
            Sigma_Inverse = inv(S);
            Sigma_Inverse_Dev = Sigma_Inverse - 1/3*eye(3,3)*trace(Sigma_Inverse);
            Sigma_Inverse_Dev_Vector = [Sigma_Inverse_Dev(1,1); Sigma_Inverse_Dev(2,2); Sigma_Inverse_Dev(3,3); Sigma_Inverse_Dev(1,2); Sigma_Inverse_Dev(1,3); Sigma_Inverse_Dev(2,3)];
            N_vector = [N/3;N/3;N/3;0;0;0];
            dQdSigma = -(N_vector + I1 / X * Sigma_Inverse_Dev_Vector);
            
            % A
            A = lambda_c * (M + N) * N;
            
            % Derivative
            dYdI1       = I2/(I3*2*X);
            dYdI2       = I1/(I3*2*X);
            dYdI3       = -I1*I2/(I3^2*2*X);
            
            dI1dSigma   = eye(3,3);
            
            dI2dSigma   = [ S(2,2)+S(3,3)   -2*S(1,2)       -2*S(1,3) ...
                          ; -2*S(1,2)       S(1,1)+S(3,3)   -2*S(2,3) ...
                          ; -2*S(1,3)       -2*S(2,3)       S(2,2)+S(1,1)];
                      
            dI3dSigma   = [ S(2,2)*S(3,3)-S(2,3)^2          2*(S(2,3)*S(1,3)-S(3,3)*S(1,2)) 2*(S(2,3)*S(1,2)-S(2,2)*S(1,3)) ...
                          ; 2*(S(2,3)*S(1,3)-S(3,3)*S(1,2)) S(1,1)*S(3,3)-S(1,3)^2          2*(S(1,3)*S(1,2)-S(1,1)*S(2,3)) ...
                          ; 2*(S(2,3)*S(1,2)-S(2,2)*S(1,3)) 2*(S(1,3)*S(1,2)-S(1,1)*S(2,3)) S(2,2)*S(1,1)-S(1,3)^2];
            
            % Tensor 3x3
            dYdSigma = dYdI1 * dI1dSigma + dYdI2 * dI2dSigma + dYdI3 * dI3dSigma;
            
            % Vector
            dYdSigma = [dYdSigma(1,1) dYdSigma(2,2) dYdSigma(3,3) dYdSigma(1,2) dYdSigma(1,3) dYdSigma(2,3)];
            
            dlambda = f / (A + dYdSigma * D * dQdSigma);
            
            lambda = lambda + dlambda;
            
            % Flow rule           
            depsp = lambda * dQdSigma;
            
            % Update stress
            depse = deps - depsp;
            dSigma = C*depse(1:4);           % stress increment
            Sigma  = [Particle.stress(1,p) ; Particle.stress(2,p) ; Szz ; Particle.stress(3,p)] + dSigma;

            % Update state variables
            Epsilon_v_plastic = depsp(1) + depsp(2) + depsp(3);
            N =  Particle.N(p,1) + lambda_c*(M+N)*Epsilon_v_plastic;
              
%             if sig(1) <= 10 
%                 Sigma(1:4) = 0;
%             end
%             
%             if sig(2) <= 10 
%                 Sigma(1:4) = 0;
%             end
%             
%             if sig(3) <= 10 
%                 Sigma(1:4) = 0;
%             end
            
            % Recompute yield function
            % Stress invariants
            Sigma_Tensor3x3 = [Sigma(1) Sigma(4) 0 ; Sigma(4) Sigma(2) 0 ; 0 0 Sigma(3)];
  
            I1 = trace(Sigma_Tensor3x3);
            I2 = 0.5*(trace(Sigma_Tensor3x3)^2 - trace(Sigma_Tensor3x3^2));
            I3 = det (Sigma_Tensor3x3);
            X = sqrt (I1*I2/I3 - 9);
            f = X - (M+N);                % Matsuoka Nakai Yield function
            
            if I3 <= 0
                f = 0;
                Sigma(1:4) = 0;
            end
            
            if I1*I2/I3 <= 9
                f = 0;
                Sigma(1:4) = 0;
            end
            
            % Liquefaction cut-off on p
            mean = Sigma(1)+Sigma(2)+Sigma(3);         
            if mean <= 1 
                Sigma(1:4) = 0;
                f = 0;
            end
        end      
  end
  
  %% Out of loop update
  Particle.stress(1,p) = -Sigma(1);
  Particle.stress(2,p) = -Sigma(2);
  Particle.stress(3,p) = -Sigma(4);
  
  % State variables
  Particle.N(p,1) = N;
  Particle.stressZZ(p,1) = -Sigma(3);
  
end

  

 
 
 
 

 

        
        
        
 