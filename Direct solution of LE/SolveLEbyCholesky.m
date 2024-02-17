%need matrix A and b beforehand
    n = size(A,1);
    A=Cholesky(A);
    y=SolveLowTriangle(A,b);
    x=SolveUpTriangle(A',y)
