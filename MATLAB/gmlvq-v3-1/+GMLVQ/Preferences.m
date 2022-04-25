classdef Preferences
    %PREFERENCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ConfMatStyle = 'default'
    end
    
    properties (Constant)
        AvailableStyles = {
            'default'
            'plotconfusion'
            'confusionchart'
            'custom'
        };
    end
    
    methods
        function obj = Preferences()
            %PREFERENCES Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = set.ConfMatStyle(obj, style)
            if ~any(strcmp(obj.AvailableStyles,style))
                fprintf(2,'The style %s is not available, choose one of the following:\n', style);
                for i = 1:length(obj.AvailableStyles)
                    fprintf('%s ', obj.AvailableStyles{i});
                end
                fprintf('\n');                
            else
                obj.ConfMatStyle = style;                
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

