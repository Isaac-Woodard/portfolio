classdef CategoricDie < Die
    % Represents a die with categorical sides.

    methods
        function obj = CategoricDie(sides, color, weight)
            arguments
                sides (:,1) string
                color (:,3) double {mustBeNonnegative} = [0 0 0] % RGB triplets.
                weight (:,1) double {mustBeNonnegative} = [1]
            end

            sides = categorical(sides);
            obj@Die(sides, color, weight);
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
                result (:,1) string
                color (:,3) double
            end

            indices = probability.randi_weight(x, obj.Weight);
            result = obj.Sides(indices);
            color = obj.Color(indices,:);
            if history
                obj.add_to_history(result);
            end
        end
    end
end