function [M] = Fenetre(Mi,borneinf,bornesup, varargin)
%varargin : numéro(s) (ordonnés) du ou des muscles exclu(s)
    if nargin==3 % A tester avec length(varargin)==0
        M=Mi(:,borneinf:bornesup);
    else
        for i=1:length(varargin)
            Mi(varargin{i}-(i-1),:)=[];
        end    
        M=Mi(:,borneinf:bornesup);
    end    
end


