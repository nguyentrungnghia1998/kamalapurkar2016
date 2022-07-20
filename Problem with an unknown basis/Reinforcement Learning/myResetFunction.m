function [InitialObservation,LoggedSignal] = myResetFunction()

LoggedSignal.State = [1;1;0;0;0];
InitialObservation = [1;1;0;0];
end