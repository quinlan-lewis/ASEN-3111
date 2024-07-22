function [c] = ChordArrayFunc(Span,N)
%This function takes in the span of the wing of the plane and creates a
%hald ellipse function for the chord length of the wings
    x = linspace(-Span, Span, N);
    
    b = .23;
    a = 1.53;
    
    c = (b/a)*sqrt(a^2 - x.^2);
end

