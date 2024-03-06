classdef Problem_type
    %PROBLEM_TYPE Class that describe if the problem is linear or non
    %linear upper mantle. 
   properties
        iteration 
        cut_off  
        Linear
    end
   
      
   
    
    methods
        function tf = islinear(obj,~)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tf = obj.Linear;        
        end
    end
end
