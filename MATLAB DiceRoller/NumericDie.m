classdef NumericDie < Die
    % Represents a die with numeric sides.

    properties (Dependent)
        ExpectedValue (1,1) double {mustBePositive}
    end

    methods
        function obj = NumericDie(sides, color, weight)
            arguments
                sides (:,1) double {mustBePositive}
                color (:,3) double {mustBeNonnegative} = [0 0 0] % RGB triplets.
                weight (:,1) double {mustBeNonnegative} = [1]
            end

            obj@Die(sides, color, weight);
        end

        function value = get.ExpectedValue(obj)
            value = sum(obj.Sides .* obj.Weight) / length(obj.Sides);
        end

        function [result, color] = roll(obj, x, history)
            % Rolls the die x times and returns the result(s). 
            % Calling with no arguments rolls the die once.
            arguments (Input)
                obj
                x (1,1) uint32 {mustBePositive} = 1
                history (1,1) logical = true
            end
            arguments (Output)
                result (:,1) double
                color (:,3) double
            end

            indices = probability.randi_weight(x, obj.Weight);
            result = obj.Sides(indices);
            color = obj.Color(indices,:);
            if history
                obj.add_to_history(result);
            end
        end

        function [result, color] = rollmod(obj, x, mod)
            % Rolls the die x times and returns the result(s). 
            % A modifier can be applied to each result.
            arguments (Input)
                obj
                x (1,1) uint32 {mustBePositive}
                mod (1,1) double
            end
            arguments (Output)
                result (:,1) double
                color (:,3) double
            end

            [result, color] = obj.roll(x);
            result = result + mod;
        end
    end
end