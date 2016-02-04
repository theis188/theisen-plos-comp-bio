function [ FuncString ] = getMatlabFunction( functionDefinition )
% Converting the SBML function into a matlab function

    MathString = functionDefinition.math; %The function as defined in SBML
    MathString = MathString(1:end-1); % Remove the last parenthesis which closes lambda
    %We need to figure out when the last input variable is defined
    %Whenever non-comma symbols (operators or parenthesis) start being used then the expression is
    %being written. ( must be used in an expression before a , is used on
    %it so this is safe.
    EarlyExp = regexp(MathString, '\(|\*|\+|/|-|\^');
    if length(EarlyExp)==1,
        EarlyExp = length(MathString);
    else
        EarlyExp = EarlyExp(2);
    end
    CommaPos = strfind(MathString, ',');
    LastVarComma = CommaPos(max(find(CommaPos<EarlyExp)));
    MathString(LastVarComma) = ')';
    MathString = strrep(MathString, 'lambda', '@');
    FuncString = strcat(functionDefinition.id, '=', MathString, ';');
    
end

