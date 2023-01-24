classdef Problem_type
    %PROBLEM_TYPE Class that describe if the problem is linear or non
    %linear upper mantle. 
   enumeration
      Linear, NonLinear
   end
    
    methods
        function tf = islinear(obj,~)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = Problem_type.Linear == obj;        
        end
    end
end
