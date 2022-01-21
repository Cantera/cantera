classdef ReactorSurface < handle

    properties
        surf_id
        area
        reactor
    end

    methods
        %% ReactorSurface class constructor

        function s = ReactorSurface(kleft, reactor, area)
            % A surface on which heterogeneous reactions take place. The
            % mechanism object (typically an instance of class 'Interface')
            % mustb be constructed so that it is properly linked to the
            % object representing the fluid in the reactor. The surface
            % temperature on each side is taken to be equal to the
            % temprature of the reactor.
            %
            % Note: all of the arguments are optional and can be activated
            % after initial construction by using the various methods of
            % the 'ReactorSurface' class.
            %
            % parameter kleft:
            %    Surface reaction mechanisms for the left-facing surface.
            %    This must bean instance of class 'Kinetics', or of a class
            %    derived from 'Kinetics', such as 'Interface'.
            % parameter reactor:
            %    Instance of class 'Reactor' to be used as the adjacent
            %    bulk phase.
            % parameter area:
            %    The area of the surface in m^2. Defaults to 1.0 m^2 if not
            %    specified.
            % return:
            %    Instance of class 'ReactorSurface'.

            checklib;

            s.surf_id = calllib(ct, 'reactorsurface_new', 0);
            s.reactor = -1;
%             if r.id < 0
%                 error(geterr);
%             end

            if nargin >= 1
                s.setKinetics(kleft);
            end

            if nargin >= 2
                if isa(reactor, 'Reactor')
                    s.install(reactor);
                else
                    warning('Reactor was not installed due to incorrect type');
                end
            end

            if nargin >= 3
                if isnumeric(area)
                    s.area = area;
                else
                    warning('Area was not a number and was not set');
                end
            end

        end

        %% Utility methods

        function clear(s)
            % Clear the ReactorSurface object from the memory.
            checklib;
            calllib(ct, 'reactorsurface_del', s.surf_id);
        end

        function install(s, r)
            % Install a ReactorSurface in a Reactor.
            checklib;
            s.reactor = r;
            calllib(ct, 'reactorsurface_install', s.surf_id, r.id);
        end

        %% ReactorSurface get methods

        function a = get.area(s)
            % Get the areaof the reactor surface in m^2.
            checklib;
            a = calllib(ct, 'reactorsurface_area', s.surf_id);
        end

        %% ReactorSurface set methods

        function set.area(s, a)
            % Set the area of a reactor surface
            checklib;
            calllib(ct, 'reactorsurface_setArea', s.surf_id, a);
        end

        function setKinetics(s, kin)
            % Setthe surface reaction mechanism on a reactor surface.
            % parameter kin:
            %    Instance of class 'Kinetics' (or another object derived
            %    from kin) to be used as the kinetic mechanism for this
            %    surface. Typically an instance of class 'Interface'.

            checklib;

            ikin = 0;
            if isa(kin, 'Kinetics')
                ikin = kin.kin_id;
            end

            calllib(ct, 'reactorsurface_setkinetics', s.surf_id, ikin);
        end

    end
end
