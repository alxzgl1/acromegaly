% The mode to be used in this library
classdef Mode < uint32
    enumeration
        % Matrix without null-space correction
        GMLVQ (0)
        
        % Matrix with null-space correction
        GMLVQNS (1)
        
        % Diagonal matrix (discouraged)
        GRLVQ (2)
        
        % GLVQ with Euclidian distance
        GLVQ (3)
    end
end

