function [EdgesOfCube] = BuildCube

r = 1000;
Vertex1 = [0 0 0 1]';
Vertex8 = [r 0 r 1]';
Vertex6 = [0 r r 1]';
Vertex7 = [r r r 1]';
Vertex3 = [r r 0 1]';
Vertex4 = [r 0 0 1]';
Vertex2 = [0 r 0 1]';
Vertex5 = [0 0 r 1]';


%C = [0 0 0; r 0 r; 0 r r; r r r; r r 0; r 0 0; 0 r 0; 0 0 r];

EdgesOfCube = [Vertex1 Vertex2 Vertex2 Vertex3 Vertex3 Vertex4 Vertex4 ...
    Vertex1 Vertex6 Vertex7 Vertex7 Vertex8 Vertex8 Vertex5 ...
    Vertex5 Vertex6 Vertex1 Vertex5 Vertex2 Vertex6 Vertex3 ...
    Vertex7 Vertex4 Vertex8];

    


   
