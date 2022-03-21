function [Particle]=Neo_Hookean_elastic(SolidModel,Particle)

% This is the hyper-elastic Neo-Hookean model
% Parameter of the model
Lambda = SolidModel.Young_modul*SolidModel.nu/(1+SolidModel.nu)/(1-2*SolidModel.nu);
Mu     = SolidModel.Young_modul/2/(1+SolidModel.nu);

% Stress tensor = Lambda*log(J)/J*eye(2,2) + Mu/J*(Defgrad tensor*Defgrad tensor'-eye(2,2));
Particle.stress(1,:) = Lambda*log(Particle.Jacobian)    +  Mu./Particle.Jacobian.*...
    (Particle.defgrad (:,1).*Particle.defgrad (:,1) + Particle.defgrad (:,2).*Particle.defgrad (:,2) - 1);


Particle.stress(2,:) = Lambda*log(Particle.Jacobian)    +  Mu./Particle.Jacobian.*...
    (Particle.defgrad (:,2).*Particle.defgrad (:,3) + Particle.defgrad (:,4).*Particle.defgrad (:,4) - 1); 

Particle.stress(3,:) =              0                   +  Mu./Particle.Jacobian.*...
    (Particle.defgrad (:,1).*Particle.defgrad (:,2) + Particle.defgrad (:,3).*Particle.defgrad (:,4) - 0); 