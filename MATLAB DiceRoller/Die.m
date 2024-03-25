classdef (Abstract) Die < handle
    % Abstract base class for a die. 
    % 
    % Note: Inheriting Handle causes instances to be passed 
    % by reference, rather than by value. In many languages,
    % pass by reference is the default.

    properties (SetAccess=protected)
        Sides (:,1)
        Color (:,3) % RGB triplets.
        Weight (:,1)
        History (:,1)
    end

    methods (Abstract)
        [result, color] = roll(obj, x, history)
            % Rolls the die x times and returns the result(s).
            % Calling with no arguments rolls the die once.
            %
            % The following arguments blocks are recommended:
            %
            % arguments (Input)
            %     obj
            %     x (1,1) uint32 {mustBePositive} = 1
            %     history (1,1) bool = true
            % end
            % 
            % arguments (Output)
            %     result (:,1)
            %     color (:,3) double
            % end
    end

    methods (Sealed)
        function clear_history(obj)
            obj.History = [];
        end
    end

    methods
        function obj = Die(sides, color, weight)
            arguments
                sides (:,1)
                color (:,3) double {mustBeNonnegative} % RGB triplets.
                weight (:,1) double {mustBeNonnegative}
            end

            a = length(sides);
            if size(color, 1) == 1
                color = repmat(color, [a 1]);
            end
            if length(weight) == 1
                weight = ones([a 1]);
            end
            if any([size(color, 1) ~= a length(weight) ~= a])
                error("The number of colors and weights should be 1 or equal to the number of sides.");
            end
            if any([any(color > 1) any(color < 0)])
                error("Color triplets must have values between 0 and 1.")
            end

            weight = weight ./ sum(weight); % Normalize.
            
            obj.Color = color;
            obj.Weight = weight;
            obj.Sides = sides;
        end
    end

    methods (Access=protected)
        function add_to_history(obj, result)
            % Adds result to history.
            arguments (Input)
                obj
                result (:,1)
            end

            if isempty(obj.History)
                obj.History = result;
            else
                obj.History = [obj.History; result];
            end
        end
    end
end