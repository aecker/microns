%{
vangogh.RFMethod (lookup) #   methods for computing receptive fields from noise
rf_method :tinyint
-----
method :varchar(20)    # the name of the method
%}

classdef RFMethod < dj.Relvar
    methods
        function fill(self)
            self.inserti({
                1   'STA'
                })
        end
    end
end