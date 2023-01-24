classdef Mantle_Unit_Properties
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        rho                % density
        alpha              % Thermal expansion
        Cp                 % Heat Capacity
        Bd                 % Compliance Diffusion
        Bn                 % Compliance Dislocation
        Ed                 % Diffusion Energy activaction
        En                 % Dislocation Energy activaction
        n                  % Stress exponent 
        R                  % Gas constant
        g                  % gravity acceleration
    end

     methods
         function obj = Mantle_Unit_Properties(rho,alpha,Cp,Bd,Bn,Ed,En,n)
             obj.rho = rho; 
             obj.alpha = alpha;
             obj.Cp = Cp;
             obj.Bd = Bd; 
             obj.Bn = Bn; 
             obj.Ed = Ed; 
             obj.En = En; 
             obj.n  = n ; 
             obj.R  = 8.314;
             obj.g = 9.81; 
         end

         function [Cd,Cn,phi] = Compute_Cd_Cn(obj,Pr,Vn,Vd,Tp)
             %UNTITLED2 Construct an instance of this class
             %   Detailed explanation goes here
            phi = obj.alpha./(obj.rho*obj.Cp);
            Cn = (obj.En*phi-Vn*(1-phi*Pr))/(obj.R*Tp);
            Cd = (obj.Ed*phi-Vd*(1-phi*Pr))/(obj.R*Tp);
         end
 
     end
end