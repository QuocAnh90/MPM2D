function Particle=Mohr_Coulomb(SolidModel,Particle,Time)

%% Update particle's velocity strain Rate
% Formula dESP = (L_sp + L_sp')/2*dt;
Particle.strainRate(:,1)    = Particle.Gradvelocity(:,1) * Time.timestep;                                   % particle strain rate xx
Particle.strainRate(:,2)    = Particle.Gradvelocity(:,4) * Time.timestep;                                   % particle strain rate yy
Particle.strainRate(:,3)    = (Particle.Gradvelocity(:,2)+Particle.Gradvelocity(:,2))/2 * Time.timestep;    % particle strain rate xy

%% Parameter of the model
E       = SolidModel.Young_modul;                           % Young modoulus
nu      = SolidModel.nu;                                    % Poissons ratio
phi     = SolidModel.Friction*pi/180;                       % Friction angle in radian
psi     = SolidModel.Dilation*pi/180;                       % Dilation angle in radian
c       = SolidModel.cohesion;

k       = (1.0+sin(phi))/(1.0-sin(phi));                    % Friction parameter in principal stress state 
m       = (1.0+sin(psi))/(1.0-sin(psi));                    % Dilation parameter in principal stress state 
comp    = 2*c*sqrt(k);                                      % Uniaxial compressive strength
PlasPar = [k , comp,  m];                                   % Plastic Parameter array

for p = 1:Particle.Count

%% -----------------------------------------Trial step------------------------------------------%
% Mohr_Coulomb is based on plane strain formulation
de_sp    = zeros(4,1);     % change of strain
de_sp(1) = Particle.strainInc(p,1);
de_sp(2) = Particle.strainInc(p,2);
de_sp(3) = 0.0;
de_sp(4) = Particle.strainInc(p,3);
SigP_up  = zeros(3,1);

% Elastic matrix
C = zeros(4,4);                                                                       
	C(1,1) = 1.0 - nu ; C(1,2) = nu		  ; C(1,3) = nu;
	C(2,1) = nu			; C(2,2) = 1.0 - nu ; C(2,3) = nu;
	C(3,1) = nu			; C(3,2) = nu		  ; C(3,3) = 1.0 - nu;
    C(4,4) = (1-2*nu);
	C = E/((1+nu)*(1-2*nu)) * C; 
dSigma = C*de_sp;    % change of stress
Sigma  = [Particle.stress(1,p) ; Particle.stress(2,p) ; nu * (Particle.stress(1,p) + Particle.stress(2,p)) ; Particle.stress(3,p)] + dSigma;

%% ----------------------------------Principal Stress & Angle-----------------------------------%
  sig_av = 0.5*(Sigma(1) + Sigma(2));
  sig_hj = sqrt((0.5*(Sigma(1) - Sigma(2)))^2 + Sigma(4)^2);
  
  SigP=zeros(3,1); 
  SigP(1) = sig_av + sig_hj;
  SigP(2) = sig_av - sig_hj;
  SigP(3) = Sigma(3); % The out-of-plane stress

  if (Sigma(3) > SigP(1)) 
		SigP(3) = SigP(2)  ;  SigP(2) = SigP(1)  ; SigP(1) = Sigma(3);
	elseif (Sigma(3) > SigP(2) )
		SigP(3) = SigP(2)  ;  SigP(2) = Sigma(3); 
  end
  %-- Sorting done --
  %-- Principal angle ------------------------
        if (Sigma(1) > Sigma(2)  &&  Sigma(4) >= 0.0) 
            psi = 0.5*atan(2*Sigma(4)/(Sigma(1)-Sigma(2)));
        elseif (Sigma(1) < Sigma(2)  &&   Sigma(4) >= 0.0) 
            psi = 0.5*(pi - atan(2*Sigma(4)/(Sigma(2)-Sigma(1))));
        elseif (Sigma(1) < Sigma(2)  &&   Sigma(4) < 0.0) 
            psi = 0.5*(atan(-2*Sigma(4)/(Sigma(2)-Sigma(1))) + pi);
        elseif (Sigma(1) > Sigma(2)  &&   Sigma(4) < 0.0) 
            psi = 0.5*(2*pi - atan(-2*Sigma(4)/(Sigma(1)-Sigma(2))));
        elseif (Sigma(1) == Sigma(2) &&  Sigma(4) > 0.0) 
            psi = 0.25*pi;
        elseif (Sigma(1) == Sigma(2) &&  Sigma(4) < 0.0) 
            psi = 0.75*pi;
        elseif (Sigma(1) == Sigma(2) &&  Sigma(4) == 0.0) 
            psi = 0.0;
        end
%--Principal angle done  ------------------------

%% --- Value of Yield function f = k*sigP1 - sigP3 - comp ---  
  f = PlasPar(1)*SigP(1) - SigP(3) - PlasPar(2); % Mohr-Coulomb yield criterion
%----------------------------------------------------------
  
  if (f <= 0)      %elastic condition
      Sigma_up = Sigma;
        Particle.stress(1,p) = Sigma_up(1);
        Particle.stress(2,p) = Sigma_up(2);
        Particle.stress(3,p) = Sigma_up(4);
      
  else             %elaaso_plastic condition
      
     %---- Position of out-of plane or in principal stress vector SigP ----	
	   if (Sigma(3) == SigP(1)) 
			ouplP = 1 ; %s1 = 2 ; s2 = 3;
	   elseif (Sigma(3) == SigP(2)) 
			ouplP = 2 ; %s1 = 1 ; s2 = 3;
       elseif (Sigma(3) == SigP(3)) 
			ouplP = 3 ; %s1 = 1 ; s2 = 2;
       end 
	 %-----Positioning done------------------------------------------------          
     %---- Stress return --------------------------------------------------	     
        %--- Preliminary parameters ----------------
        apex = [comp/(k-1.0) comp/(k-1.0) comp/(k-1.0)]'; % Stress coordinate of the criterions apex
        D = C(1:3,1:3);% Relates to normal stresses
        
        den = k*(D(1,1)*m-D(1,3)) - D(3,1)*m + D(3,3); % denominator a'*D*b
        Rp=zeros(3,1);
        Rp(1) = (D(1,1)*m-D(1,3))/den; % Rp is the scaled direction of
        Rp(2) = (D(2,1)*m-D(2,3))/den; % the update direction,
        Rp(3) = (D(3,1)*m-D(3,3))/den; % Rp = D*b/(a'*D*b) a, b gradient of
							  % yield surface and plastic potential, respectively.
        %--Preliminary parameters calculated---------
        %--Boundary planes --------------------------        
        SigPapex = SigP - apex; % Vector from predictor stress to the apex

        % Boundary plane between regions I and II
        NI_II=zeros(3,1);
        NI_II(1) = Rp(2)*k - Rp(3);   % NI_II = cross(Rp,R1)
        NI_II(2) = Rp(3)   - Rp(1)*k; % R1 = [1 1 k], direction
        NI_II(3) = Rp(1)   - Rp(2);   % vector of line 1
	
        pI_II = NI_II(1)*SigPapex(1) + NI_II(2)*SigPapex(2) + NI_II(3)*SigPapex(3);

        % Boundary plane between regions I and III
        NI_III=zeros(3,1);
        NI_III(1) = Rp(2)*k - Rp(3)*k; % NI_III = cross(Rp,R2)
        NI_III(2) = Rp(3)   - Rp(1)*k; % R2 = [1 k k], direction
        NI_III(3) = Rp(1)*k - Rp(2);   % vector of line 2
	
        pI_III = NI_III(1)*SigPapex(1) + NI_III(2)*SigPapex(2) + NI_III(3)*SigPapex(3);
        %--Boundary planes calculated ----------------
        
        %--- t-paramters for region determination --
	    % secondary surface in region II a = [0 k -1], b = [0 m -1]
        den = k*(D(2,2)*m-D(2,3)) - D(3,2)*m + D(3,3); % denominator a'*D*b
        Rp2=zeros(3,1);
        Rp2(1) = (D(1,2)*m-D(1,3))/den; % Rp is the scaled direction of
        Rp2(2) = (D(2,2)*m-D(2,3))/den; % the update direction,
        Rp2(3) = (D(3,2)*m-D(3,3))/den; % Rp = D*b/(a'*D*b) a, b gradient of
							   % yield surface and plastic potential, respectively.	
        N2=zeros(3,1);
        N2(1) = Rp(2)*Rp2(3) - Rp(3)*Rp2(2); % N2 = cross(Rp,Rp2)
        N2(2) = Rp(3)*Rp2(1) - Rp(1)*Rp2(3);
        N2(3) = Rp(1)*Rp2(2) - Rp(2)*Rp2(1);
	
        den_t = N2(1) + N2(2) + k*N2(3); % N2'*R1, R1 = [1 1 k] Direction of line 1
        % t-parameter of line 1
        t1 = (N2(1)*SigPapex(1) + N2(2)*SigPapex(2) + N2(3)*SigPapex(3))/den_t;

        % secondary surface in region III a = [k -1 0], b = [m -1 0]
        den = k*(D(1,1)*m-D(1,2)) - D(2,1)*m + D(2,2); % denominator a'*D*b
        Rp3=zeros(3,1);
        Rp3(1) = (D(1,1)*m-D(1,2))/den; % Rp is the scaled direction of
        Rp3(2) = (D(2,1)*m-D(2,2))/den; % the update direction,
        Rp3(3) = (D(3,1)*m-D(3,2))/den; % Rp = D*b/(a'*D*b) a, b gradient of
							   % yield surface and plastic potential, respectively.
	    N3=zeros(3,1);
        N3(1) = Rp(2)*Rp3(3) - Rp(3)*Rp3(2); % N3 = cross(Rp,Rp3)
        N3(2) = Rp(3)*Rp3(1) - Rp(1)*Rp3(3);
        N3(3) = Rp(1)*Rp3(2) - Rp(2)*Rp3(1);
	
        den_t = N3(1) + k*N3(2) + k*N3(3); % N3'*R2, R2 = [1 k k] Direction of line 2
        % t-parameter of line 2
        t2 = (N3(1)*(SigP(1)-apex(1)) + N3(2)*(SigP(2)-apex(2)) + N3(3)*(SigP(3)-apex(3)))/den_t;
        %---t-paramters for region determination calculated -------
        
        %--- Region determination and stress update ------- 
        if (t1 > 0.0 && t2 > 0.0)
% 	    region = 4;
	    SigP_up = apex;
        elseif (pI_II < 0.0) 
% 	    region = 2;
	    SigP_up(1) = t1   + apex(1); % SigP_up = t1*R1 + apex
	    SigP_up(2) = t1   + apex(2);  % R1 = [1 1 k], direction
	    SigP_up(3) = t1*k + apex(3);  % vector of line 1
        elseif (pI_III <= 0.0) 
% 	    region = 1;
	    SigP_up(1) = SigP(1) - f*Rp(1); % SigP_up = SigP - SiPla
	    SigP_up(2) = SigP(2) - f*Rp(2); % SiPla is the plastic corrector
	    SigP_up(3) = SigP(3) - f*Rp(3); % given by f*Rp
        else
% 	    region = 3;
	    SigP_up(1) = t2   + apex(1) ; %SigP_up = t2*R2 + apex
	    SigP_up(2) = t2*k + apex(2) ; % R2 = [1 k k], direction
	    SigP_up(3) = t2*k + apex(3) ; % vector of line 2
        end                     % Regions and update
        %---Region determination and stress update completed-       
        %--- Tranformation matrix A ------------------
        sin_psi  = sin(psi)  ; cos_psi = cos(psi);
	    sin_psi2 = sin_psi^2 ; cos_psi2 = cos_psi^2;
	    sin_cos_psi = sin_psi*cos_psi; sin_2psi=sin(2*psi);
        
        if (ouplP == 2) 
		A=[cos_psi2 0 sin_psi2 sin_cos_psi;0 1 0 0;sin_psi2 0 cos_psi2 -sin_cos_psi;-sin_2psi 0 sin_2psi (cos_psi2-sin_psi2)];
	    elseif (ouplP == 1) 
		A=[1 0 0 0;0 cos_psi2 sin_psi2 sin_cos_psi;0 sin_psi2 cos_psi2 -sin_cos_psi;0 -sin_2psi sin_2psi (cos_psi2-sin_psi2)];
	    elseif (ouplP == 3) 
		A=[cos_psi2 sin_psi2 0 sin_cos_psi;sin_psi2 cos_psi2 0 -sin_cos_psi;0 0 1 0;-sin_2psi sin_2psi 0 (cos_psi2-sin_psi2)];
        end 

        %--- Tranformation matrix A calculated ----------
        SigP_upextend=[SigP_up(1);SigP_up(2);SigP_up(3);0];
        Sigma_up=(A')*SigP_upextend;  
               
        if (ouplP == 2) 
		    Particle.stress(1,p) = Sigma_up(1);
            Particle.stress(2,p) = Sigma_up(3);
            Particle.stress(3,p) = Sigma_up(4);
	    elseif (ouplP == 1) 
		    Particle.stress(1,p) = Sigma_up(2);
            Particle.stress(2,p) = Sigma_up(3);
            Particle.stress(3,p) = Sigma_up(4);
	    elseif (ouplP == 3) 
		    Particle.stress(1,p) = Sigma_up(1);
            Particle.stress(2,p) = Sigma_up(2);
            Particle.stress(3,p) = Sigma_up(4);
        end              
  end
end

  

 
 
 
 

 

        
        
        
 